#!/usr/bin/env python
import pysam
import argparse
import sys
import re
from collections import defaultdict
import os
import subprocess
from pathlib import Path

EDGE_RE = re.compile(r'label="[^"]*\s+(\d+)\(')


def read_fasta(fasta_path):
    """Yield (contig_name, sequence) from a FASTA file."""
    name, seq = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            yield name, "".join(seq)


def write_single_fasta(contig_name, seq, out_dir):
    """Write one contig to a FASTA file and return the path."""
    fasta_path = Path(out_dir) / f"{contig_name}.fa"
    with open(fasta_path, "w") as f:
        f.write(f">{contig_name}\n{seq}\n")
    return fasta_path


def run_jumbodbg(contig_name, fasta_path, out_dir, threads=10):
    """Run jumboDBG on a single FASTA file and return path to graph.dot."""
    contig_out = Path(out_dir) / contig_name
    if not os.path.exists(contig_out):
        contig_out.mkdir(parents=True, exist_ok=True)
        cmd = [
            "/Poppy/zmzhang/Consensus_Assembly/cLJA/lib/LJA/bin/jumboDBG",
            "-k", "101",
            "--reads", str(fasta_path),
            "-t", str(threads),
            "--coverage",
            "-o", str(contig_out)
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return contig_out / "graph.dot"


def parse_graph_dot(dot_path):
    """Return the sum of edge lengths from a graph.dot file."""
    total = 0
    with open(dot_path) as f:
        for line in f:
            m = EDGE_RE.search(line)
            if m:
                total += int(m.group(1))
    return total/2


def compute_dbg_ratios(fasta_path, out_dir="jumbodbg_tmp", threads=10):
    """
    Run jumboDBG for each contig in a FASTA file and compute
    (sum of edge lengths / contig length) ratio.

    Returns:
        dict mapping name parts -> ratio, e.g.:
        {
            'contigA': 1.02,
            '001': 1.02,
            ...
        }
    """
    os.makedirs(out_dir, exist_ok=True)
    results = {}

    for name, seq in read_fasta(fasta_path):
        print(f"Processing {name} ...", flush=True)
        contig_len = len(seq)
        contig_fa = write_single_fasta(name, seq, out_dir)
        dot_file = run_jumbodbg(name, contig_fa, out_dir, threads)
        edge_sum = parse_graph_dot(dot_file)
        ratio = edge_sum / contig_len if contig_len > 0 else 0

        parts = name.split("_")
        if len(parts) > 0:
            results[parts[0]] = ratio
        if len(parts) > 1:
            results[parts[1]] = ratio
    
    print(results)

    return results

def classify_dot_edges(dot_file):
    with open(dot_file, 'r') as f:
        dot_text = f.read()
    # Pattern to extract edge info
    edge_pattern = re.compile(r'"?([^"]+)"?\s*->\s*"?(.*?)"?\s*\[label="([^"]+)"')

    edges = []
    node_counts = defaultdict(int)

    # First pass: parse all edges and count node usage
    for line in dot_text.strip().splitlines():
        match = edge_pattern.search(line)
        if match:
            src, tgt, label = match.groups()
            label = label[:label.find(" ")]
            edges.append((src, tgt, label))
            node_counts[src] += 1
            if tgt != src:
                node_counts[tgt] += 1

    # Second pass: classify edges
    results = defaultdict(str)
    for src, tgt, label in edges:
        if node_counts[src] > 1 or node_counts[tgt] > 1:
            classification = "multi-edge"
        else:
            classification = "isolated"
        results[label] = classification

    return results

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
    overlap_len = max(0, overlap_end - overlap_start)

    # Unique lengths
    unique1_len = max(0, ref_end1 - ref_start1 - overlap_len)
    unique2_len = max(0, ref_end2 - ref_start2 - overlap_len)

    # For query lengths (could adjust similarly if you prefer)
    len1 = aln1["length_query"]
    len2 = aln2["length_query"]

    # For weighting, we can proportionally approximate via reference length in percent:
    weight_A = unique1_len
    weight_B = overlap_len
    weight_C = unique2_len

    # For degenerate case (complete overlap)
    if weight_A == 0 and weight_B == 0 and weight_C == 0:
        weight_B = 1

    total_weight = weight_A + weight_B + weight_C

    # Compute weighted identity and identity_nogap
    identity = (
        aln1["identity"] * weight_A +
        ((aln1["identity"] + aln2["identity"]) / 2) * weight_B +
        aln2["identity"] * weight_C
    ) / total_weight

    identity_nogap = (
        aln1["identity_nogap"] * weight_A +
        ((aln1["identity_nogap"] + aln2["identity_nogap"]) / 2) * weight_B +
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

def merge_alignments(alignments):
    alignments = sorted(alignments, key=lambda x: (x["ref_start"], -x["identity_nogap"], -x["identity"]))
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
                    # Keep current, skip nxt
                    j += 1
                    continue
                else:
                    # Replace current with nxt, skip old current
                    current = nxt
                    j += 1
                    continue

            # Else: check if gap is acceptable
            gap = nxt["ref_start"] - current["ref_end"]
            if gap <= 5:  # still allow small gap or overlap
                merge_two_alignments(current, nxt)
                j += 1
                continue
            break  # too far to merge

        merged.append(current)
        i = j
    return merged

def query_overlap(a, b):
    # Compute overlap fraction over the shorter span
    start = max(a["query_start"], b["query_start"])
    end = min(a["query_end"], b["query_end"])
    overlap = max(0, end - start)
    len_a = a["query_end"] - a["query_start"]
    len_b = b["query_end"] - b["query_start"]
    return overlap / min(len_a, len_b) if min(len_a, len_b) > 0 else 0

def ref_overlap(aln1, aln2):
    """
    Return overlap proportion of shorter alignment over the longer one
    based on reference coordinates (percent scale).
    """
    start1, end1 = aln1["ref_start"], aln1["ref_end"]
    start2, end2 = aln2["ref_start"], aln2["ref_end"]

    inter_start = max(start1, start2)
    inter_end = min(end1, end2)
    overlap_len = max(0, inter_end - inter_start)

    len1 = end1 - start1
    len2 = end2 - start2
    shorter = min(len1, len2)

    return overlap_len / shorter if shorter > 0 else 0

def true_query_start_end(cigar):
    # Parse the CIGAR string into list of (length, op)
    tokens = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    tokens = [(int(length), op) for length, op in tokens]

    # Compute the total length of the read (excluding hard clips)
    read_length = sum(length for length, op in tokens if op in ('M', 'I', 'S', '=', 'X', 'H'))

    # Count left and right clip (soft + hard)
    left_clip = 0
    right_clip = 0

    if tokens:
        if tokens[0][1] in ('S', 'H'):
            left_clip = tokens[0][0]
        if tokens[-1][1] in ('S', 'H'):
            right_clip = tokens[-1][0]

    # true query start is where alignment begins in full read
    true_query_start = left_clip
    # true query end is where alignment ends in full read
    true_query_end = read_length - right_clip

    return true_query_start, true_query_end, read_length

def calculate_identities(alignment, gap_threshold=10):
    # Normal identity
    nm = alignment.get_tag('NM')

    # Second identity: matches / query_alignment_length, ignoring gaps >= gap_threshold
    matches = 0
    insertions = 0  # gaps in query (I operations)
    deletions = 0  # gaps in reference (D operations)
    long_gaps = 0

    if hasattr(alignment, 'cigartuples') and alignment.cigartuples is not None:
        for op, length in alignment.cigartuples:
            if op == 0:  # M (match or mismatch)
                matches += length
            elif op == 1:  # I (insertion in query)
                insertions += length
                if length >= gap_threshold:
                    long_gaps += length
            elif op == 7:  # = (match)
                matches += length
            elif op == 2:  # D (deletion in reference)
                deletions += length
                if length >= gap_threshold:
                    long_gaps += length

    mismatches = nm - insertions - deletions
    true_matches = matches - mismatches
   
    denominator = true_matches + mismatches + insertions + deletions
    if denominator <= 0:
        identity = 0
    else:
        identity = true_matches / (true_matches + mismatches + insertions + deletions)

    denominator = true_matches + mismatches + insertions + deletions - long_gaps
    if denominator <= 0:
        identity_nogap = 0.0
    else:
        identity_nogap = true_matches / denominator

    return identity*100,identity_nogap*100

def dedup_by_query_span(alignments):
    alignments = sorted(alignments, key=lambda x: (-(x["query_end"] - x["query_start"]), -x["identity_nogap"], -x["identity"]))
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
    
    deduped = [aln for aln in deduped if aln['query_end']-aln['query_start'] > 10 or aln['length_query']> 1000000]
    return deduped

def filter_alignments_with_identity(bam_file_path, threshold=0):
    high_identity_alignments = {}
    total_alignments = 0
    processed_alignments = 0

    with pysam.AlignmentFile(bam_file_path, "r") as bamfile:
        for alignment in bamfile:
            total_alignments += 1

            if alignment.is_unmapped:
                continue

            identity, identity_nogap = calculate_identities(alignment, gap_threshold=10)
            
            if identity_nogap >= threshold*100:
                processed_alignments += 1
                
                ref_id = alignment.reference_name[alignment.reference_name.find('_') + 1:]
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
                
                query_name = alignment.query_name.split('_')[0]
                # query_name = alignment.reference_name[alignment.query_name.find('_') + 1:]
                if query_name not in high_identity_alignments:
                    high_identity_alignments[query_name] = []
                
                ref_len = bamfile.get_reference_length(alignment.reference_name)

                if alignment.reference_name in id_map:
                    ref_id = id_map[alignment.reference_name]
                
                # if "A" not in ref_id and "B" not in ref_id and "M" not in ref_id and "P" not in ref_id and "Chr" not in ref_id and "X" not in ref_id and "Y" not in ref_id:
                #     continue

                query_alignment_start, query_alignment_end, query_len = true_query_start_end(alignment.cigarstring)

                high_identity_alignments[query_name].append({
                    'identity': identity,
                    'identity_nogap': identity_nogap,
                    'ref_id': ref_id if alignment.is_forward else '-' + ref_id,
                    'ref_start': 100 * (alignment.reference_start / ref_len) if alignment.is_forward else 100 * ((ref_len - alignment.reference_end) / ref_len),
                    'ref_end': 100 * (alignment.reference_end / ref_len) if alignment.is_forward else 100 * ((ref_len - alignment.reference_start) / ref_len),
                    'ref_start_cood': alignment.reference_start if alignment.is_forward else ref_len - alignment.reference_end,
                    'ref_end_cood': alignment.reference_end if alignment.is_forward else ref_len - alignment.reference_start,
                    'query_start': 100 * (query_alignment_start / query_len) if alignment.is_forward else 100 * ((query_len - query_alignment_end) / query_len),
                    'query_end': 100 * (query_alignment_end / query_len) if alignment.is_forward else 100 * ((query_len - query_alignment_start) / query_len),
                    'query_start_cood': query_alignment_start if alignment.is_forward else query_len - query_alignment_end,
                    'query_end_cood': query_alignment_end if alignment.is_forward else query_len - query_alignment_start,
                    'reverse': alignment.is_reverse,
                    'length_entire_query': query_len,
                    'length_entire_ref': ref_len,
                    'length_query': alignment.query_alignment_end - alignment.query_alignment_start,
                    'length_ref': alignment.reference_end - alignment.reference_start,
                    'alignment': alignment
                })

                # query_name = "*" + alignment.query_name.split('_')[0]
                query_name = alignment.query_name.split('_')[1]
                if query_name not in high_identity_alignments:
                    high_identity_alignments[query_name] = []
                
                high_identity_alignments[query_name].append({
                    'identity': identity,
                    'identity_nogap': identity_nogap,
                    'ref_id': ref_id if alignment.is_reverse else '-' + ref_id,
                    'ref_start': 100 * (alignment.reference_start / ref_len) if alignment.is_reverse else 100 * ((ref_len - alignment.reference_end) / ref_len),
                    'ref_end': 100 * (alignment.reference_end / ref_len) if alignment.is_reverse else 100 * ((ref_len - alignment.reference_start) / ref_len),
                    'ref_start_cood': alignment.reference_start if alignment.is_reverse else ref_len - alignment.reference_end,
                    'ref_end_cood': alignment.reference_end if alignment.is_reverse else ref_len - alignment.reference_start,
                    'query_start': 100 * ((query_len - query_alignment_end) / query_len) if alignment.is_forward else 100 * (query_alignment_start / query_len),
                    'query_end': 100 * ((query_len - query_alignment_start) / query_len) if alignment.is_forward else 100 * (query_alignment_end / query_len),
                    'query_start_cood': query_len - query_alignment_end if alignment.is_forward else query_alignment_start,
                    'query_end_cood': query_len - query_alignment_start if alignment.is_forward else query_alignment_end,
                    'reverse': alignment.is_reverse,
                    'length_entire_query': query_len,
                    'length_entire_ref': ref_len,
                    'length_query': alignment.query_alignment_end - alignment.query_alignment_start,
                    'length_ref': alignment.reference_end - alignment.reference_start,
                    'alignment': alignment
                })

        # keep the largest span if overlap
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

            # Sort retained alignments by query coordinates
            high_identity_alignments[query_name] = sorted(new_alignments, key=lambda x: (x["query_start"], x["query_end"], -x["identity"], -x["identity_nogap"]))
            # high_identity_alignments[query_name] = sorted(new_alignments, key=lambda x: (x["query_start"]-x["query_end"], -x["identity"], -x["identity_nogap"]))
    
    
    # print(f"Edges with label: {processed_alignments} of {total_alignments}")
    return high_identity_alignments

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Extract high-identity alignments from a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("fasta_file", default="", help="Path to graph.fasta file.")
    parser.add_argument("-t", "--threshold", type=float, default=0.9,
                        help="Identity threshold (default: 0.9).")
    parser.add_argument("-o", "--output", help="Path to the output file. If not specified, prints to stdout.")
    parser.add_argument("-d", "--dotfile", help="")
    return parser.parse_args()

def main():
    args = parse_arguments()
    bam_file_path = args.bam_file
    threshold = args.threshold
    output_path = args.output
    dotfile = args.dotfile

    high_identity_alignments = filter_alignments_with_identity(bam_file_path, threshold=threshold)
    
    # Prepare output
    output_lines = []
    for query_name, alignments in high_identity_alignments.items():
        for aln in alignments:
            if "tig" not in aln['ref_id'] and aln['query_end']-aln['query_start'] > 10 or aln['length_query']> 1000000:
                # output_lines.append(f"{query_name}\t{aln['ref_id']}\tQ:{aln['length_query']:,}({aln['query_start']:.0f}-{aln['query_end']:.0f})\tQ:{aln['length_query']:,}({aln['query_start_cood']:.0f}-{aln['query_end_cood']:.0f})\tR:{aln['length_ref']:,}({aln['ref_start']:.2f}-{aln['ref_end']:.2f})\tR:{aln['length_ref']:,}({aln['ref_start_cood']:.0f}-{aln['ref_end_cood']:.0f})\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}")
                output_lines.append(f"{query_name}\t{aln['ref_id']}\tQ:{aln['length_query']:,}({aln['query_start']:.0f}-{aln['query_end']:.0f})\tR:{aln['length_ref']:,}({aln['ref_start']:.2f}-{aln['ref_end']:.2f})\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}")
    
    ref_alignments = defaultdict(list)

    # Group ref alignments by chromosome/haplotype structure
    haplo_groups = defaultdict(lambda: defaultdict(list))  # autosomes: {chrom: {'M': [...], 'P': [...]}}
    special_refs = defaultdict(list)  # X, Y, etc.

    # Step 1: Group alignments
    for query_name, alignments in high_identity_alignments.items():
        for aln in alignments:
            ref_id = aln['ref_id']
            ref_span = aln['ref_end'] - aln['ref_start']
            query_span = aln['query_end'] - aln['query_start']

            if "tig" not in ref_id and (query_span > 10 or aln['length_query'] > 1_000_000):
                m = re.match(r"^(\d+)([MPAB])$", ref_id)  # match 1M, 2P, ...
                if m:
                    chrom = m.group(1)
                    hap = m.group(2)
                    haplo_groups[chrom][hap].append((query_name, aln, ref_span, query_span))
                else:
                    # e.g., X, Y, etc.
                    special_refs[ref_id].append((query_name, aln, ref_span, query_span))

    # Step 2: Merge logic
    def merge_spans(spans):
        if not spans:
            return []
        spans.sort()
        merged = []
        current_start, current_end = spans[0]
        for start, end in spans[1:]:
            if start - current_end <= 2:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))
        return merged

    def print_merged_spans(label, spans):
        if not spans:
            return
        merged = merge_spans(spans)
        span_strs = [f"{start:.2f}-{end:.2f}" for start, end in merged]
        span_str = " and ".join(span_strs)
        print(f"--- Merged span positions for {label} (%): {span_str} ---")
    
    edge2comp = {}
    if dotfile:
        edge2comp = classify_dot_edges(dotfile)
    
    uniq_ratios = compute_dbg_ratios(args.fasta_file, output_path + "_jumboDBG", 80) if args.fasta_file != "" else {}

    # Step 3: Print autosomal haplotype groups
    print(f"Edge ID\tLength of edge\tRef ID\tLength of ref\tQuery span\tRef span\tPI\tNon-repetitiveness")
    for chrom in sorted(haplo_groups.keys(), key=lambda x: int(x)):
        all_spans = []
        for hap in ['M', 'P', 'A', 'B']:
            ref_id = f"{chrom}{hap}"
            if hap not in haplo_groups[chrom]:
                continue
            spans = []
            print(f"\n--- Reference: {ref_id} ---")
            for query_name, aln, ref_span, query_span in sorted(haplo_groups[chrom][hap], key=lambda x: x[1]['ref_start']):
                spans.append((aln['ref_start'], aln['ref_end']))
                line = (
                    f"{query_name}\t{aln['length_entire_query']:,}\t{ref_id}\t{aln['length_entire_ref']:,}"
                    f"\tQ:{aln['length_query']:,} ({aln['query_start']:.0f}-{aln['query_end']:.0f}, span={query_span:.0f})"
                    f"\tR:{aln['length_ref']:,} ({aln['ref_start']:.2f}-{aln['ref_end']:.2f}, span={ref_span:.2f})"
                    f"\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}\tNon-repetitiveness={uniq_ratios.get(query_name, 0):.3f}"
                ) if query_name not in edge2comp else (
                    f"{query_name}\t{aln['length_entire_query']:,}\t{ref_id}\t{aln['length_entire_ref']:,}"
                    f"\tQ:{aln['length_query']:,} ({aln['query_start']:.0f}-{aln['query_end']:.0f}, span={query_span:.0f})"
                    f"\tR:{aln['length_ref']:,} ({aln['ref_start']:.2f}-{aln['ref_end']:.2f}, span={ref_span:.2f})"
                    f"\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}\tNon-repetitiveness={uniq_ratios.get(query_name, 0):.3f}"
                    f"\t{edge2comp[query_name]}"
                )
                # output_lines.append(line)
                print(line)
            print_merged_spans(ref_id, spans)
            all_spans.extend(spans)

        if all_spans:
            print_merged_spans(f"{chrom}M+{chrom}P", all_spans)

    # Step 4: Print special chromosomes (X, Y, etc.)
    for ref_id in sorted(special_refs.keys(), key=lambda x: (x != 'X', x != 'Y', 100+int(x[4:]) if "-Chr" in x else (int (x[3:]) if "Chr" in x else -1000000))):
        spans = []
        print(f"\n--- Reference: {ref_id} ---")
        for query_name, aln, ref_span, query_span in sorted(special_refs[ref_id], key=lambda x: x[1]['ref_start']):
            spans.append((aln['ref_start'], aln['ref_end']))
            line = (
                f"{query_name}\t{aln['length_entire_query']:,}\t{ref_id}\t{aln['length_entire_ref']:,}"
                f"\tQ:{aln['length_query']:,} ({aln['query_start']:.0f}-{aln['query_end']:.0f}, span={query_span:.0f})"
                f"\tR:{aln['length_ref']:,} ({aln['ref_start']:.2f}-{aln['ref_end']:.2f}, span={ref_span:.2f})"
                f"\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}\tNon-repetitiveness={uniq_ratios.get(query_name, 0):.3f}"
            )if query_name not in edge2comp else (
                    f"{query_name}\t{aln['length_entire_query']:,}\t{ref_id}\t{aln['length_entire_ref']:,}"
                    f"\tQ:{aln['length_query']:,} ({aln['query_start']:.0f}-{aln['query_end']:.0f}, span={query_span:.0f})"
                    f"\tR:{aln['length_ref']:,} ({aln['ref_start']:.2f}-{aln['ref_end']:.2f}, span={ref_span:.2f})"
                    f"\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}\tNon-repetitiveness={uniq_ratios.get(query_name, 0):.3f}"
                    f"\t{edge2comp[query_name]}"
                )
            # output_lines.append(line)
            print(line)
        print_merged_spans(ref_id, spans)
    # Write to file or stdout
    if output_path:
        try:
            with open(output_path, 'w') as outfile:
                outfile.write("\n".join(output_lines))
            # print(f"Results written to {output_path}")
        except IOError:
            print(f"Error: Cannot write to file '{output_path}'.", file=sys.stderr)
    else:
        print("\n".join(output_lines))

if __name__ == "__main__":
    main()
