#!/usr/bin/env python
import argparse
import sys
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

def parse_paf_line(line):
    """
    Parse one line of a PAF file into a dictionary.
    """
    fields = line.strip().split("\t")
    if len(fields) < 12:
        return None
    paf = {
        "query_name": fields[0],
        "query_len": int(fields[1]),
        "query_start": int(fields[2]),
        "query_end": int(fields[3]),
        "strand": fields[4],
        "ref_name": fields[5],
        "ref_len": int(fields[6]),
        "ref_start": int(fields[7]),
        "ref_end": int(fields[8]),
        "matches": int(fields[9]),
        "aln_len": int(fields[10]),
        "mapq": int(fields[11]),
    }
    # Optional tags
    for f in fields[12:]:
        try:
            tag, typ, val = f.split(":", 2)
            paf[tag] = val
        except ValueError:
            print(f"Warning: could not parse optional field: {f} for {line.strip()}", file=sys.stderr)
            continue
    return paf

def read_fasta_lengths(fasta_path):
    """
    Read sequence lengths from a FASTA file.
    Returns a dict {seq_name: length}.
    """
    lengths = {}
    with open(fasta_path) as fh:
        name = None
        seq_len = 0
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]  # take first word as ID
                seq_len = 0
            else:
                seq_len += len(line)
        if name is not None:
            lengths[name] = seq_len
    return lengths


def filtering_identity_counts(paf_entry, use_cigar=False, gap_threshold=10):
    """Return identity counts and whether they came from an explicit CIGAR.

    Without ``use_cigar``, use standard PAF identity. With ``use_cigar``,
    require an explicit =/X CIGAR and omit insertions and deletions at least
    ``gap_threshold`` bases long from the denominator.
    """
    if gap_threshold < 1:
        raise ValueError("gap_threshold must be at least 1")

    if not use_cigar:
        return paf_entry["matches"], paf_entry["aln_len"], False

    cigar = paf_entry.get("cg")
    if not cigar:
        raise ValueError(
            "PAF record has no cg:Z: CIGAR tag, but --use-cigar was requested; "
            "regenerate the PAF with minimap2 -c --eqx"
        )

    tokens = re.findall(r'(\d+)([MIDNSHP=XB])', cigar)
    if not tokens or ''.join(length + op for length, op in tokens) != cigar:
        raise ValueError(f"Could not parse CIGAR string: {cigar}")

    matches = 0
    denominator = 0
    for length_string, operation in tokens:
        length = int(length_string)
        if operation == "=":
            matches += length
            denominator += length
        elif operation == "X":
            denominator += length
        elif operation in ("I", "D", "N"):
            if length < gap_threshold:
                denominator += length
        elif operation == "M":
            raise ValueError(
                "CIGAR contains M operations, but --use-cigar requires "
                "explicit =/X operations; regenerate the PAF with --eqx"
            )
        elif operation == "B":
            raise ValueError(f"Unsupported CIGAR B operation: {cigar}")
        # Clipping and padding do not contribute.

    return matches, denominator, True


def union_length(intervals):
    """Return the number of bases covered by the union of half-open intervals."""
    if not intervals:
        return 0

    covered = 0
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            covered += current_end - current_start
            current_start, current_end = start, end
    return covered + current_end - current_start


def edge_nodes(contig_name):
    """Return both endpoint node names from node1.id1_node2.id2."""
    edge_parts = contig_name.split('_')
    if len(edge_parts) != 2:
        raise ValueError(
            f"Expected contig name node1.id1_node2.id2, got: {contig_name}"
        )

    nodes = set()
    for edge_part in edge_parts:
        if "." not in edge_part:
            raise ValueError(
                f"Expected node.id component in contig name, got: {edge_part}"
            )
        node, edge_id = edge_part.rsplit(".", 1)
        if not node or not edge_id:
            raise ValueError(
                f"Expected node.id component in contig name, got: {edge_part}"
            )
        nodes.add(node)
    return nodes


def filter_alignments_with_identity(
    paf_file_path,
    cigar_paf_file_path,
    output_path,
    coverage_threshold=85.0,
    shared_node_coverage_threshold=50.0,
    identity_threshold=70.0,
    gap_threshold=10,
    min_cigar_coverage=20.0,
):
    """Combine broad-span and CIGAR alignments per unordered sequence pair.

    The ordinary PAF supplies union coverage. The CIGAR PAF supplies explicit
    no-long-gap identity and the amount of base-level alignment support.
    """
    def read_pair_alignments(path, require_cigar):
        pairs = {}
        with open(path) as fh:
            for line in fh:
                paf = parse_paf_line(line)
                if paf is None or paf["query_name"] == paf["ref_name"]:
                    continue

                query_name = paf["query_name"]
                ref_name = paf["ref_name"]
                pair = tuple(sorted((query_name, ref_name)))
                if pair not in pairs:
                    pairs[pair] = {
                        'lengths': {
                            pair[0]: lengths[pair[0]],
                            pair[1]: lengths[pair[1]],
                        },
                        'intervals': {
                            pair[0]: [],
                            pair[1]: [],
                        },
                        'paf_matches': 0,
                        'paf_alignment_length': 0,
                        'identity_matches': 0,
                        'identity_denominator': 0,
                        'record_count': 0,
                    }

                aggregate = pairs[pair]
                aggregate["intervals"][query_name].append(
                    (paf["query_start"], paf["query_end"])
                )
                aggregate["intervals"][ref_name].append(
                    (paf["ref_start"], paf["ref_end"])
                )
                aggregate["paf_matches"] += paf["matches"]
                aggregate["paf_alignment_length"] += paf["aln_len"]
                if require_cigar:
                    identity_matches, identity_denominator, _ = (
                        filtering_identity_counts(
                            paf, use_cigar=True, gap_threshold=gap_threshold
                        )
                    )
                    aggregate["identity_matches"] += identity_matches
                    aggregate["identity_denominator"] += identity_denominator
                aggregate["record_count"] += 1
        return pairs

    normal_pairs = read_pair_alignments(paf_file_path, require_cigar=False)
    cigar_pairs = read_pair_alignments(cigar_paf_file_path, require_cigar=True)

    # output summary
    print("Candidate\tCandidate Len\tCovered Bases\tCoverage (%)\tOther\tOther Len\tIdentity\tCIGAR Covered Bases\tCIGAR Coverage (%)\tCIGAR Identity\tIdentity No-long-gap\tNormal Records\tCIGAR Records\tShared Node\tStatus")
    removed_edges = set()
    with open(output_path, 'w') if output_path else sys.stdout as out_edges:
        for pair, aggregate in normal_pairs.items():
            name1, name2 = pair
            length1 = aggregate["lengths"][name1]
            length2 = aggregate["lengths"][name2]
            if length1 < length2:
                candidate_name, other_name = name1, name2
            elif length2 < length1:
                candidate_name, other_name = name2, name1
            else:
                candidate_name, other_name = name1, name2

            candidate_length = aggregate["lengths"][candidate_name]
            other_length = aggregate["lengths"][other_name]
            covered_bases = union_length(
                aggregate["intervals"][candidate_name]
            )
            coverage = (
                covered_bases / candidate_length * 100
                if candidate_length else 0.0
            )
            paf_alignment_length = aggregate["paf_alignment_length"]
            identity = (
                aggregate["paf_matches"] / paf_alignment_length * 100
                if paf_alignment_length else 0.0
            )

            cigar_aggregate = cigar_pairs.get(pair)
            if cigar_aggregate:
                cigar_covered_bases = union_length(
                    cigar_aggregate["intervals"][candidate_name]
                )
                cigar_coverage = (
                    cigar_covered_bases / candidate_length * 100
                    if candidate_length else 0.0
                )
                cigar_paf_alignment_length = (
                    cigar_aggregate["paf_alignment_length"]
                )
                cigar_identity = (
                    cigar_aggregate["paf_matches"]
                    / cigar_paf_alignment_length * 100
                    if cigar_paf_alignment_length else 0.0
                )
                identity_denominator = cigar_aggregate["identity_denominator"]
                identity_nogap = (
                    cigar_aggregate["identity_matches"]
                    / identity_denominator * 100
                    if identity_denominator else 0.0
                )
                cigar_record_count = cigar_aggregate["record_count"]
            else:
                cigar_covered_bases = 0
                cigar_coverage = 0.0
                cigar_identity = 0.0
                identity_nogap = 0.0
                cigar_record_count = 0

            candidate_id = candidate_name.split('_', 1)[0]
            other_id = other_name.split('_', 1)[0]
            share_node = bool(
                edge_nodes(candidate_name) & edge_nodes(other_name)
            )
            required_coverage = (
                shared_node_coverage_threshold
                if share_node
                else coverage_threshold
            )
            required_cigar_bases = (
                candidate_length * min_cigar_coverage / 100
            )
            is_cognate = (
                candidate_length < other_length
                and coverage >= required_coverage
                and identity_nogap >= identity_threshold
                and cigar_covered_bases >= required_cigar_bases
            )
            status = "Cognate" if is_cognate else "None"

            if is_cognate:
                removed_edges.update(candidate_name.split('_', 1))

            print(
                f'{candidate_id}\t{candidate_length}\t'
                f'{covered_bases}\t{coverage:.2f}\t'
                f'{other_id}\t{other_length}\t'
                f'{identity:.2f}\t'
                f'{cigar_covered_bases}\t{cigar_coverage:.2f}\t'
                f'{cigar_identity:.2f}\t'
                f'{identity_nogap:.2f}\t'
                f'{aggregate["record_count"]}\t{cigar_record_count}\t'
                f'{"Yes" if share_node else "No"}\t{status}'
            )

        for edge in sorted(removed_edges):
            out_edges.write(edge + '\n')

    return normal_pairs, cigar_pairs

def generate_paf(
    query_fasta,
    ref_fasta,
    output_paf,
    threads=8,
    minimap2_preset="asm20",
    use_P=False,
    use_X=False,
    use_cigar=False,
):
    """
    Generate paf file from query vs reference using minimap2.
    """
    extra_opts = []
    if use_P and not use_X:
        extra_opts = ["-p", "0.1"]
    elif use_X and not use_P:
        extra_opts = ["-X"]
    elif use_P and use_X:
        extra_opts = ["-p", "0.1", "-X"]

    cmd = [
        "minimap2", "-t", str(threads), "-x", minimap2_preset,
    ]
    if use_cigar:
        cmd.extend(["-c", "--eqx"])
    cmd.extend(extra_opts)
    cmd.extend([ref_fasta, query_fasta])

    output_paf_tmp = f"{output_paf}.tmp"
    with open(output_paf_tmp, "w") as out:
        try:
            subprocess.run(cmd, stdout=out, check=True)
        except subprocess.CalledProcessError as error:
            raise RuntimeError(
                f"PAF generation failed for {output_paf}"
            ) from error
    os.replace(output_paf_tmp, output_paf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True, help="Query FASTA file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--paf", required=True, help="Ordinary minimap2 PAF")
    parser.add_argument(
        "--cigar-paf",
        help="CIGAR minimap2 PAF (default: <paf>.cigar)"
    )
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--output", help="Output file for summary")
    parser.add_argument(
        "--coverage-threshold", type=float, default=85.0,
        help="Minimum aggregate query coverage percentage (default: 85)"
    )
    parser.add_argument(
        "--shared-node-coverage-threshold",
        "--same-start-coverage-threshold",
        dest="shared_node_coverage_threshold",
        type=float,
        default=20.0,
        help=(
            "Minimum aggregate query coverage percentage when query and "
            "reference share either endpoint node (default: 50)"
        )
    )
    parser.add_argument(
        "--identity-threshold", type=float, default=70.0,
        help="Minimum aggregate filtering identity percentage (default: 70)"
    )
    parser.add_argument(
        "--gap-threshold", type=int, default=10,
        help=(
            "With CIGAR, exclude insertions/deletions of at least this length "
            "from the identity denominator (default: 10)"
        )
    )
    parser.add_argument(
        "--min-cigar-coverage", type=float, default=20.0,
        help="Minimum CIGAR-supported candidate percentage (default: 20)"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-P", action="store_true", help="Use minimap2 with -p 0.1")
    group.add_argument("-X", action="store_true", help="Use minimap2 with -X")
    args = parser.parse_args()
    cigar_paf = args.cigar_paf or args.paf + ".cigar"

    missing_pafs = []
    if not os.path.exists(args.paf):
        missing_pafs.append((args.paf, False))
    if not os.path.exists(cigar_paf):
        missing_pafs.append((cigar_paf, True))

    if len(missing_pafs) == 2:
        normal_threads = max(1, args.threads // 2)
        cigar_threads = max(1, args.threads - normal_threads)
        job_threads = [normal_threads, cigar_threads]
    else:
        job_threads = [args.threads] * len(missing_pafs)

    if missing_pafs:
        with ThreadPoolExecutor(max_workers=len(missing_pafs)) as executor:
            futures = []
            for (paf_path, use_cigar), threads in zip(
                missing_pafs, job_threads
            ):
                futures.append(executor.submit(
                    generate_paf,
                    args.query,
                    args.ref,
                    paf_path,
                    threads=threads,
                    use_P=args.P,
                    use_X=args.X,
                    use_cigar=use_cigar,
                ))
            for future in futures:
                future.result()

    global lengths
    lengths = read_fasta_lengths(args.ref)
    lengths.update(read_fasta_lengths(args.query))

    filter_alignments_with_identity(
        args.paf,
        cigar_paf,
        args.output,
        coverage_threshold=args.coverage_threshold,
        shared_node_coverage_threshold=args.shared_node_coverage_threshold,
        identity_threshold=args.identity_threshold,
        gap_threshold=args.gap_threshold,
        min_cigar_coverage=args.min_cigar_coverage,
    )
