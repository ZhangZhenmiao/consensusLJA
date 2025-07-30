#!/usr/bin/env python
import pysam
import argparse
import sys
import re
from collections import defaultdict

def query_overlap(a, b):
    # Compute overlap fraction over the shorter span
    start = max(a["query_start"], b["query_start"])
    end = min(a["query_end"], b["query_end"])
    overlap = max(0, end - start)
    len_a = a["query_end"] - a["query_start"]
    len_b = b["query_end"] - b["query_start"]
    return overlap / min(len_a, len_b) if min(len_a, len_b) > 0 else 0

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
    alignment_length = alignment.query_alignment_length
    if alignment_length == 0 or alignment_length - nm <=0:
        identity = 0
    else:
        identity = (alignment_length - nm) / alignment_length

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
            elif op == 2:  # D (deletion in reference)
                deletions += length

    mismatches = nm - insertions - deletions
    true_matches = matches - mismatches

    denominator = alignment_length - long_gaps
    if denominator <= 0:
        identity_nogap = 0.0
    else:
        identity_nogap = true_matches / denominator

    return identity*100,identity_nogap*100

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
                    "CM075343.1": "Chr14"   
                }
                
                query_name = alignment.query_name.split('_')[0]
                if query_name not in high_identity_alignments:
                    high_identity_alignments[query_name] = []
                
                query_len = contig_lengths[query_name]
                ref_len = bamfile.get_reference_length(alignment.reference_name)

                if alignment.reference_name in id_map:
                    ref_id = id_map[alignment.reference_name]

                query_alignment_start, query_alignment_end, read_length = true_query_start_end(alignment.cigarstring)
                # print(query_name, read_length, query_len)
                if (read_length != query_len):
                    print(alignment.cigarstring)
                assert(read_length == query_len)

                high_identity_alignments[query_name].append({
                    'identity': identity,
                    'identity_nogap': identity_nogap,
                    'ref_id': ref_id if alignment.is_forward else '-' + ref_id,
                    'ref_start': 100 * (alignment.reference_start / ref_len) if alignment.is_forward else 100 * ((ref_len - alignment.reference_end) / ref_len),
                    'ref_end': 100 * (alignment.reference_end / ref_len) if alignment.is_forward else 100 * ((ref_len - alignment.reference_start) / ref_len),
                    'query_start': 100 * (query_alignment_start / query_len) if alignment.is_forward else 100 * ((query_len - query_alignment_end) / query_len),
                    'query_end': 100 * (query_alignment_end / query_len) if alignment.is_forward else 100 * ((query_len - query_alignment_start) / query_len),
                    'reverse': alignment.is_reverse,
                    'length_query': alignment.query_alignment_end - alignment.query_alignment_start,
                    'length_ref': alignment.reference_end - alignment.reference_start,
                    'alignment': alignment
                })

                query_name = alignment.query_name.split('_')[1]
                if query_name not in high_identity_alignments:
                    high_identity_alignments[query_name] = []
                
                high_identity_alignments[query_name].append({
                    'identity': identity,
                    'identity_nogap': identity_nogap,
                    'ref_id': ref_id if alignment.is_reverse else '-' + ref_id,
                    'ref_start': 100 * (alignment.reference_start / ref_len) if alignment.is_reverse else 100 * ((ref_len - alignment.reference_end) / ref_len),
                    'ref_end': 100 * (alignment.reference_end / ref_len) if alignment.is_reverse else 100 * ((ref_len - alignment.reference_start) / ref_len),
                    'query_start': 100 * ((query_len - query_alignment_end) / query_len) if alignment.is_forward else 100 * (query_alignment_start / query_len),
                    'query_end': 100 * ((query_len - query_alignment_start) / query_len) if alignment.is_forward else 100 * (query_alignment_end / query_len),
                    'reverse': alignment.is_reverse,
                    'length_query': alignment.query_alignment_end - alignment.query_alignment_start,
                    'length_ref': alignment.reference_end - alignment.reference_start,
                    'alignment': alignment
                })
        
        ## keep every alignment
        # for query_name, alignments in high_identity_alignments.items():
        #     alignments = sorted(alignments, key=lambda x: (x["ref_id"], x["ref_start"], x["ref_end"], x["query_start"], x["query_end"], -x["identity"]))
        #     # for x in alignments:
        #     #     print(query_name, x['ref_id'], x['ref_start'], x['ref_end'], x['length_query'], x['identity'], sep='\t')
        #     best_map = {}
        #     for x in alignments:
        #         key = (x['ref_id'], x['query_start'], x['query_end'], x['ref_start'], x['ref_end'])
        #         if key not in best_map or x['identity'] > best_map[key]['identity']:
        #             best_map[key] = x
        #     new_list = list(best_map.values())
        #     high_identity_alignments[query_name] =  sorted(new_list, key=lambda x: (x["query_start"], x["query_end"]))

        # keep the largest span if overlap
        for query_name, alignments in high_identity_alignments.items():
            new_alignments = []

            # Group by ref_id
            ref_groups = defaultdict(list)
            for aln in alignments:
                ref_groups[aln["ref_id"]].append(aln)

            for ref_id, group in ref_groups.items():
                group = sorted(group, key=lambda x: -(x["query_end"] - x["query_start"]))  # largest span first
                used = [False] * len(group)

                for i, aln_i in enumerate(group):
                    if used[i]:
                        continue
                    # This alignment will be kept
                    new_alignments.append(aln_i)
                    for j in range(i + 1, len(group)):
                        if used[j]:
                            continue
                        aln_j = group[j]
                        if query_overlap(aln_i, aln_j) >= 0.9:
                            used[j] = True  # Drop overlapping one with smaller span

            # Sort retained alignments by query coordinates
            high_identity_alignments[query_name] = sorted(new_alignments, key=lambda x: (x["query_start"], x["query_end"]))
    
    
    # print(f"Edges with label: {processed_alignments} of {total_alignments}")
    return high_identity_alignments

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Extract high-identity alignments from a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("fasta_file", help="Path to graph.fasta file.")
    parser.add_argument("-t", "--threshold", type=float, default=0.9,
                        help="Identity threshold (default: 0.9).")
    parser.add_argument("-o", "--output", help="Path to the output file. If not specified, prints to stdout.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    bam_file_path = args.bam_file
    threshold = args.threshold
    output_path = args.output
    fasta_file = pysam.FastaFile(args.fasta_file)
    global contig_lengths
    contig_lengths = {}
    # Iterate through each contig in the FASTA file
    for contig in fasta_file.references:
        # Get the length of the contig
        length = fasta_file.get_reference_length(contig)
        contig_lengths[contig.split('_')[0]] = length
        contig_lengths[contig.split('_')[1]] = length
        # print(contig)

    # Close the FASTA file
    fasta_file.close()

    high_identity_alignments = filter_alignments_with_identity(bam_file_path, threshold=threshold)
    
    # Prepare output
    output_lines = []
    for query_name, alignments in high_identity_alignments.items():
        for aln in alignments:
            if "tig" not in aln['ref_id'] and aln['query_end']-aln['query_start'] > 10 or aln['length_query']> 1000000:
                output_lines.append(f"{query_name}\t{aln['ref_id']}\tQ:{aln['length_query']:,}({aln['query_start']:.0f}-{aln['query_end']:.0f})\tR:{aln['length_ref']:,}({aln['ref_start']:.2f}-{aln['ref_end']:.2f})\tPI={aln['identity']:.0f}/{aln['identity_nogap']:.0f}")
    
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
