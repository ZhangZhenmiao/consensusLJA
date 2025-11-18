#!/usr/bin/env python
import argparse
import sys


def calculate_identity_paf(n_match, aln_len):
    if aln_len <= 0:
        return 0.0
    return n_match / aln_len


def filter_alignments_with_identity(paf_file_path, threshold=0.0):
    high_identity = {}

    with open(paf_file_path) as paf:
        for line in paf:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            qname = fields[0]
            try:
                qstart = int(fields[2])
                qend = int(fields[3])
                strand = fields[4]
                tname = fields[5]
                tlen = int(fields[6])
                tstart = int(fields[7])
                tend = int(fields[8])
                n_match = int(fields[9])
                aln_len = int(fields[10])
            except:
                continue

            identity = calculate_identity_paf(n_match, aln_len)
            if identity < threshold:
                continue

            is_forward = (strand == "+")
            if is_forward:
                ref_id = tname
                ref_start = tstart
                ref_end = tend
            else:
                ref_id = "-" + tname
                ref_start = tlen - tend
                ref_end = tlen - tstart

            high_identity.setdefault(qname, []).append({
                "ref_id": ref_id,
                "ref_start": ref_start,
                "ref_end": ref_end,
                "ref_len": tlen,
                "identity": identity,
                "span": abs(ref_end - ref_start),
            })

    return high_identity


def read_fasta_lengths(fasta_path):
    lengths = {}
    name = None
    size = 0

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = size
                name = line[1:].split()[0]
                size = 0
            else:
                size += len(line)

        if name is not None:
            lengths[name] = size

    return lengths


def compute_covered_length(intervals):
    """compute union of reference intervals"""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    cur_s, cur_e = intervals[0]
    total = 0
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    total += cur_e - cur_s
    return total, cur_s, cur_e


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Base-level coverage per contig (keep output format unchanged).")
    parser.add_argument("paf_file")
    parser.add_argument("fasta_file")
    parser.add_argument("-t", "--threshold", type=float, default=0.0)
    parser.add_argument("-o", "--output")
    return parser.parse_args()


def main():
    args = parse_arguments()

    contig_lengths = read_fasta_lengths(args.fasta_file)
    aligns = filter_alignments_with_identity(args.paf_file, args.threshold)

    output_lines = []

    for qname, alns in aligns.items():
        if qname not in contig_lengths:
            continue
        qlen = contig_lengths[qname]

        # group intervals by ref ignoring "-"
        ref_groups = {}
        for a in alns:
            ref_base = a["ref_id"]
            ref_groups.setdefault(ref_base, []).append(a)

        # compute base-level coverage
        best_ref = None
        best_cov = -1
        best_info = None

        for ref_base, group in ref_groups.items():
            intervals = [(a["ref_start"], a["ref_end"]) for a in group]
            covered, start, end = compute_covered_length(intervals)

            if covered > best_cov:
                best_cov = covered
                best_ref = ref_base
                ref_len = group[0]["ref_len"]

                # weighted identity
                total_w = sum(a["span"] for a in group)
                if total_w > 0:
                    w_identity = sum(a["identity"] * a["span"] for a in group) / total_w
                else:
                    w_identity = 0.0

                best_info = (start, end, ref_len, w_identity)

        if best_info is None:
            continue

        start, end, ref_len, identity = best_info

        # KEEP EXACT OUTPUT FORMAT
        output_lines.append(
            f"{qname}\t{best_ref}\t{qlen}\t{start}\t{end}\t{ref_len}\t{identity:.2f}"
        )

    if args.output:
        with open(args.output, "w") as out:
            out.write("\n".join(output_lines))
    else:
        print("\n".join(output_lines))


if __name__ == "__main__":
    main()
