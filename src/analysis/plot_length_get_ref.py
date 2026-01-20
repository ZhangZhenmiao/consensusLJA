#!/usr/bin/env python
import argparse
import sys
import re
from collections import defaultdict


def merge_two_alignments(aln1, aln2):
    """
    Merge aln2 into aln1, updating coordinates and computing weighted identity.
    """

    # Determine reference interval of each
    ref_start1 = aln1["ref_start"]
    ref_end1 = aln1["ref_end"]
    ref_start2 = aln2["ref_start"]
    ref_end2 = aln2["ref_end"]

    # Compute overlap
    overlap_start = max(ref_start1, ref_start2)
    overlap_end = min(ref_end1, ref_end2)
    overlap_len = max(0.0, overlap_end - overlap_start)

    # Unique lengths
    unique1_len = max(0.0, (ref_end1 - ref_start1) - overlap_len)
    unique2_len = max(0.0, (ref_end2 - ref_start2) - overlap_len)

    # For weighting, approximate via reference length:
    weight_A = unique1_len
    weight_B = overlap_len
    weight_C = unique2_len

    # Degenerate case (complete overlap)
    if weight_A == 0 and weight_B == 0 and weight_C == 0:
        weight_B = 1.0

    total_weight = weight_A + weight_B + weight_C

    # Compute weighted identity and identity_nogap
    identity = (
        aln1["identity"] * weight_A +
        ((aln1["identity"] + aln2["identity"]) / 2.0) * weight_B +
        aln2["identity"] * weight_C
    ) / total_weight

    identity_nogap = (
        aln1["identity_nogap"] * weight_A +
        ((aln1["identity_nogap"] + aln2["identity_nogap"]) / 2.0) * weight_B +
        aln2["identity_nogap"] * weight_C
    ) / total_weight

    # Update ref and query boundaries
    aln1["ref_start"] = min(ref_start1, ref_start2)
    aln1["ref_end"] = max(ref_end1, ref_end2)
    aln1["ref_start_cood"] = min(aln1["ref_start_cood"], aln2["ref_start_cood"])
    aln1["ref_end_cood"] = max(aln1["ref_end_cood"], aln2["ref_end_cood"])
    aln1["query_start"] = min(aln1["query_start"], aln2["query_start"])
    aln1["query_end"] = max(aln1["query_end"], aln2["query_end"])
    aln1["query_start_cood"] = min(aln1["query_start_cood"], aln2["query_start_cood"])
    aln1["query_end_cood"] = max(aln1["query_end_cood"], aln2["query_end_cood"])

    aln1["length_query"] += aln2["length_query"]
    aln1["length_ref"] += aln2["length_ref"]

    aln1["identity"] = identity
    aln1["identity_nogap"] = identity_nogap


def ref_overlap(aln1, aln2):
    """
    Return overlap proportion of shorter alignment over the longer one
    based on reference percentage coordinates.
    """
    start1, end1 = aln1["ref_start"], aln1["ref_end"]
    start2, end2 = aln2["ref_start"], aln2["ref_end"]

    inter_start = max(start1, start2)
    inter_end = min(end1, end2)
    overlap_len = max(0.0, inter_end - inter_start)

    len1 = end1 - start1
    len2 = end2 - start2
    shorter = min(len1, len2)

    return overlap_len / shorter if shorter > 0 else 0.0


def query_overlap(a, b):
    # Compute overlap fraction over the shorter query span (percent scale)
    start = max(a["query_start"], b["query_start"])
    end = min(a["query_end"], b["query_end"])
    overlap = max(0.0, end - start)
    len_a = a["query_end"] - a["query_start"]
    len_b = b["query_end"] - b["query_start"]
    return overlap / min(len_a, len_b) if min(len_a, len_b) > 0 else 0.0


def merge_alignments(alignments):
    alignments = sorted(
        alignments,
        key=lambda x: (x["ref_start"], -x["identity_nogap"], -x["identity"])
    )
    merged = []
    i = 0
    while i < len(alignments):
        current = alignments[i]
        j = i + 1
        while j < len(alignments):
            nxt = alignments[j]

            # Check reference overlap
            ro = ref_overlap(current, nxt)

            if ro >= 0.9:
                # Choose longer ref span
                span1 = current["ref_end"] - current["ref_start"]
                span2 = nxt["ref_end"] - nxt["ref_start"]
                if span1 >= span2:
                    j += 1
                    continue
                else:
                    current = nxt
                    j += 1
                    continue

            # Else: check if gap is acceptable
            gap = nxt["ref_start"] - current["ref_end"]
            if gap <= 5.0:  # still allow small gap or overlap
                merge_two_alignments(current, nxt)
                j += 1
                continue
            break  # too far to merge

        merged.append(current)
        i = j
    return merged


def dedup_by_query_span(alignments):
    alignments = sorted(
        alignments,
        key=lambda x: (-(x["query_end"] - x["query_start"]),
                       -x["identity_nogap"], -x["identity"])
    )
    used = [False] * len(alignments)
    deduped = []

    for i, aln_i in enumerate(alignments):
        if used[i]:
            continue
        deduped.append(aln_i)
        for j in range(i + 1, len(alignments)):
            if used[j]:
                continue
            aln_j = alignments[j]
            if query_overlap(aln_i, aln_j) >= 0.9:
                used[j] = True

    deduped = [
        aln for aln in deduped
        if (aln['query_end'] - aln['query_start'] > 10.0) or
           (aln['length_query'] > 1_000_000)
    ]
    return deduped


def filter_alignments_with_identity(paf_file_path, threshold=0.0):
    """
    Read PAF and return dict: query_name -> list of merged high-identity alignments.
    Identity is computed as n_match / aln_len from columns 10 and 11 (0-100%).
    """
    high_identity_alignments = {}

    # same id_map as in your original code
    id_map = {
        "chr1_mat_hsa1": "1M",
        "chr2_mat_hsa3": "2M",
        "chr3_mat_hsa4": "3M",
        "chr4_mat_hsa5": "4M",
        "chr5_mat_hsa6": "5M",
        "chr6_mat_hsa7": "6M",
        "chr7_mat_hsa8": "7M",
        "chr8_mat_hsa10": "8M",
        "chr9_mat_hsa11": "9M",
        "chr10_mat_hsa12": "10M",
        "chr11_mat_hsa9": "11M",
        "chr12_mat_hsa2a": "12M",
        "chr13_mat_hsa2b": "13M",
        "chr14_mat_hsa13": "14M",
        "chr14_mat_hsa13_random_utig4-822": "14M",
        "chr14_mat_hsa13_random_utig4-823": "14M",
        "chr14_mat_hsa13_random_utig4-824": "14M",
        "chr14_mat_hsa13_random_utig4-825": "14M",
        "chr14_mat_hsa13_random_utig4-826": "14M",
        "chr14_mat_hsa13_random_utig4-827": "14M",
        "chr14_mat_hsa13_random_utig4-2241": "14M",
        "chr14_mat_hsa13_random_utig4-2242": "14M",
        "chr15_mat_hsa14": "15M",
        "chr16_mat_hsa15": "16M",
        "chr17_mat_hsa18": "17M",
        "chr18_mat_hsa16": "18M",
        "chr19_mat_hsa17": "19M",
        "chr20_mat_hsa19": "20M",
        "chr21_mat_hsa20": "21M",
        "chr22_mat_hsa21": "22M",
        "chr23_mat_hsa22": "23M",
        "chrX_mat_hsaX": "X",
        "chr1_pat_hsa1": "1P",
        "chr2_pat_hsa3": "2P",
        "chr3_pat_hsa4": "3P",
        "chr4_pat_hsa5": "4P",
        "chr5_pat_hsa6": "5P",
        "chr6_pat_hsa7": "6P",
        "chr7_pat_hsa8": "7P",
        "chr8_pat_hsa10": "8P",
        "chr9_pat_hsa11": "9P",
        "chr10_pat_hsa12": "10P",
        "chr11_pat_hsa9": "11P",
        "chr12_pat_hsa2a": "12P",
        "chr13_pat_hsa2b": "13P",
        "chr14_pat_hsa13": "14P",
        "chr15_pat_hsa14": "15P",
        "chr16_pat_hsa15": "16P",
        "chr17_pat_hsa18": "17P",
        "chr18_pat_hsa16": "18P",
        "chr19_pat_hsa17": "19P",
        "chr20_pat_hsa19": "20P",
        "chr21_pat_hsa20": "21P",
        "chr22_pat_hsa21": "22P",
        "chr23_pat_hsa22": "23P",
        "chr1522_pat_hsa1421_random_utig4-95": "1522P",
        "chr1522_pat_hsa1421_random_utig4-96": "1522P",
        "chr1522_pat_hsa1421_random_utig4-97": "1522P",
        "chr1522_pat_hsa1421_random_utig4-99": "1522P",
        "chr1522_pat_hsa1421_random_utig4-100": "1522P",
        "chr1522_pat_hsa1421_random_utig4-145": "1522P",
        "chr1522_pat_hsa1421_random_utig4-147": "1522P",
        "chr1522_pat_hsa1421_random_utig4-325": "1522P",
        "chr1522_pat_hsa1421_random_utig4-327": "1522P",
        "chr1522_pat_hsa1421_random_utig4-839": "1522P",
        "chr1522_pat_hsa1421_random_utig4-996": "1522P",
        "chr1522_pat_hsa1421_random_utig4-997": "1522P",
        "chr1522_pat_hsa1421_random_utig4-1011": "1522P",
        "chr1522_pat_hsa1421_random_utig4-2055": "1522P",
        "chrY_pat_hsaY": "Y",
        "chr1_mat": "1M",
        "chr2_mat": "2M",
        "chr3_mat": "3M",
        "chr4_mat": "4M",
        "chr5_mat": "5M",
        "chr6_mat": "6M",
        "chr7_mat": "7M",
        "chr8_mat": "8M",
        "chr9_mat": "9M",
        "chr10_mat": "10M",
        "chr11_mat": "11M",
        "chr12_mat": "12M",
        "chr13_mat": "13M",
        "chr14_mat": "14M",
        "chrX_mat": "X",
        "chr1_pat": "1P",
        "chr2_pat": "2P",
        "chr3_pat": "3P",
        "chr4_pat": "4P",
        "chr5_pat": "5P",
        "chr6_pat": "6P",
        "chr7_pat": "7P",
        "chr8_pat": "8P",
        "chr9_pat": "9P",
        "chr10_pat": "10P",
        "chr11_pat": "11P",
        "chr12_pat": "12P",
        "chr13_pat": "13P",
        "chr14_pat": "14P",
        "chrY_pat": "Y",
        "chrM": "mtDNA",
        "CM075330.1": "Chr1",
        "CM075331.1": "Chr2",
        "CM075332.1": "Chr3",
        "CM075333.1": "Chr4",
        "CM075334.1": "Chr5",
        "CM075335.1": "Chr6",
        "CM075336.1": "Chr7",
        "CM075337.1": "Chr8",
        "CM075338.1": "Chr9",
        "CM075339.1": "Chr10",
        "CM075340.1": "Chr11",
        "CM075341.1": "Chr12",
        "CM075342.1": "Chr13",
        "CM075343.1": "Chr14",
        "NC_091245.1": "Chr1",
        "NC_091246.1": "Chr2",
        "NC_091247.1": "Chr3",
        "NC_091248.1": "Chr4",
        "NC_091249.1": "Chr5",
        "NC_091250.1": "Chr6",
        "NC_091251.1": "Chr7",
        "NC_091252.1": "Chr8",
        "NC_091253.1": "Chr9",
        "NC_091254.1": "Chr10",
        "NC_091255.1": "Chr11",
        "NC_091256.1": "Chr12",
        "NC_091257.1": "Chr13",
        "NC_091258.1": "Chr14",
        "NC_091259.1": "Chr15",
        "NC_091260.1": "Chr16",
        "NC_091261.1": "Chr17",
        "NC_091262.1": "Chr18",
        "NC_091263.1": "Chr19",
        "NC_091264.1": "Chr20",
        "NC_091265.1": "Chr21",
        "NC_091266.1": "Chr22",
        "NC_091267.1": "Chr23",
        "NC_091268.1": "Chr24",
        "NC_091269.1": "Chr25",
        "NC_091270.1": "Chr26",
        "NC_091271.1": "Y",
        "NC_091727.1": "X"
    }

    with open(paf_file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 12:
                continue

            qname = fields[0]
            try:
                qlen = int(fields[1])
                qstart = int(fields[2])
                qend = int(fields[3])
                strand = fields[4]
                tname = fields[5]
                tlen = int(fields[6])
                tstart = int(fields[7])
                tend = int(fields[8])
                n_match = int(fields[9])
                aln_len = int(fields[10])
            except ValueError:
                continue

            if aln_len <= 0:
                continue

            # Identity in percent
            identity = (n_match / aln_len) * 100.0
            identity_nogap = identity  # no gap info in plain PAF

            if identity_nogap < threshold * 100.0:
                continue

            ref_id = tname
            if tname in id_map:
                ref_id = id_map[tname]

            is_forward = (strand == '+')
            is_reverse = not is_forward

            # Percent coordinates on reference
            if is_forward:
                ref_start_pct = 100.0 * (tstart / tlen)
                ref_end_pct = 100.0 * (tend / tlen)
                ref_start_cood = float(tstart)
                ref_end_cood = float(tend)
            else:
                ref_start_pct = 100.0 * ((tlen - tend) / tlen)
                ref_end_pct = 100.0 * ((tlen - tstart) / tlen)
                ref_start_cood = float(tlen - tend)
                ref_end_cood = float(tlen - tstart)

            # Percent coordinates on query
            if is_forward:
                query_start_pct = 100.0 * (qstart / qlen)
                query_end_pct = 100.0 * (qend / qlen)
                query_start_cood = float(qstart)
                query_end_cood = float(qend)
            else:
                query_start_pct = 100.0 * ((qlen - qend) / qlen)
                query_end_pct = 100.0 * ((qlen - qstart) / qlen)
                query_start_cood = float(qlen - qend)
                query_end_cood = float(qlen - qstart)

            length_query = float(qend - qstart)
            length_ref = float(tend - tstart)

            # First "orientation" (like your first dict)
            qname_primary = qname.split('_')[0]
            if qname_primary not in high_identity_alignments:
                high_identity_alignments[qname_primary] = []

            high_identity_alignments[qname_primary].append({
                'identity': identity,
                'identity_nogap': identity_nogap,
                'ref_id': ref_id if is_forward else '-' + ref_id,
                'ref_start': ref_start_pct,
                'ref_end': ref_end_pct,
                'ref_start_cood': ref_start_cood,
                'ref_end_cood': ref_end_cood,
                'query_start': query_start_pct,
                'query_end': query_end_pct,
                'query_start_cood': query_start_cood,
                'query_end_cood': query_end_cood,
                'reverse': not is_forward,
                'length_entire_query': float(qlen),
                'length_entire_ref': float(tlen),
                'length_query': length_query,
                'length_ref': length_ref,
            })

            qname_secondary = "*" + qname_primary
            if qname_secondary not in high_identity_alignments:
                high_identity_alignments[qname_secondary] = []

            if is_forward:
                # flip signs / coords similarly to your BAM code
                ref_id2 = ref_id if is_reverse else '-' + ref_id
                ref_start_pct2 = 100.0 * ((tlen - tend) / tlen)
                ref_end_pct2 = 100.0 * ((tlen - tstart) / tlen)
                ref_start_cood2 = float(tlen - tend)
                ref_end_cood2 = float(tlen - tstart)

                query_start_pct2 = 100.0 * ((qlen - qend) / qlen)
                query_end_pct2 = 100.0 * ((qlen - qstart) / qlen)
                query_start_cood2 = float(qlen - qend)
                query_end_cood2 = float(qlen - qstart)
            else:
                ref_id2 = ref_id if is_reverse else '-' + ref_id
                ref_start_pct2 = 100.0 * (tstart / tlen)
                ref_end_pct2 = 100.0 * (tend / tlen)
                ref_start_cood2 = float(tstart)
                ref_end_cood2 = float(tend)

                query_start_pct2 = 100.0 * (qstart / qlen)
                query_end_pct2 = 100.0 * (qend / qlen)
                query_start_cood2 = float(qstart)
                query_end_cood2 = float(qend)

            high_identity_alignments[qname_secondary].append({
                'identity': identity,
                'identity_nogap': identity_nogap,
                'ref_id': ref_id2,
                'ref_start': ref_start_pct2,
                'ref_end': ref_end_pct2,
                'ref_start_cood': ref_start_cood2,
                'ref_end_cood': ref_end_cood2,
                'query_start': query_start_pct2,
                'query_end': query_end_pct2,
                'query_start_cood': query_start_cood2,
                'query_end_cood': query_end_cood2,
                'reverse': is_reverse,
                'length_entire_query': float(qlen),
                'length_entire_ref': float(tlen),
                'length_query': length_query,
                'length_ref': length_ref,
            })

    # dedup & merge per query + ref_id (same as your BAM version)
    for query_name, alignments in high_identity_alignments.items():
        new_alignments = []

        # Group by ref_id
        ref_groups = defaultdict(list)
        for aln in alignments:
            ref_groups[aln["ref_id"]].append(aln)

        for ref_id, group in ref_groups.items():
            deduped = dedup_by_query_span(group)
            merged = merge_alignments(deduped)
            final_alignments = dedup_by_query_span(merged)
            new_alignments.extend(final_alignments)

        high_identity_alignments[query_name] = sorted(
            new_alignments,
            key=lambda x: (x["query_start"], x["query_end"],
                           -x["identity"], -x["identity_nogap"])
        )

    return high_identity_alignments


def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract high-identity alignments from a PAF file."
    )
    parser.add_argument("paf_file", help="Path to the input PAF file.")
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.5,
        help="Identity threshold (default: 0)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output file."
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    paf_file_path = args.paf_file
    threshold = args.threshold
    output_path = args.output

    high_identity_alignments = filter_alignments_with_identity(
        paf_file_path,
        threshold=threshold
    )

    # Prepare output (same style as your original output_lines)
    output_lines = []
    for query_name, alignments in high_identity_alignments.items():
        for aln in alignments:
            if ("tig" not in aln['ref_id'] and
                ((aln['query_end'] - aln['query_start'] > 10.0) or
                 (aln['length_query'] > 1_000_000))):
                output_lines.append(
                    f"{query_name}\t{aln['ref_id']}"
                    f"\tQ:{int(aln['length_query']):,}"
                    f"({aln['query_start']:.0f}-{aln['query_end']:.0f})"
                    f"\tR:{int(aln['length_ref']):,}"
                    f"({aln['ref_start']:.2f}-{aln['ref_end']:.2f})"
                    f"\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}"
                )

    try:
        with open(output_path, 'w') as outfile:
            outfile.write("\n".join(output_lines))
    except IOError:
        print(f"Error: Cannot write to file '{output_path}'.", file=sys.stderr)


if __name__ == "__main__":
    main()
