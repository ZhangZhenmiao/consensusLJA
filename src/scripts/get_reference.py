#!/usr/bin/env python
import pysam
import argparse
import sys

def calculate_identity(alignment):
    try:
        if alignment.has_tag('NM'):
            nm = alignment.get_tag('NM')
        else:
            nm = sum([count for (operation, count) in alignment.cigartuples if operation == 8])
    except:
        nm = 0
    
    alignment_length = alignment.query_alignment_length
    if alignment_length == 0:
        return 0
    
    identity = (alignment_length - nm) / alignment_length
    return identity

def filter_alignments_with_identity(bam_file_path, threshold=0):
    high_identity_alignments = {}
    total_alignments = 0
    processed_alignments = 0

    with pysam.AlignmentFile(bam_file_path, "r") as bamfile:
        for alignment in bamfile:
            total_alignments += 1

            if alignment.is_unmapped:
                continue

            identity = calculate_identity(alignment)
            
            if identity >= threshold:
                processed_alignments += 1
                query_name = alignment.query_name.split('_')[0]
                if query_name not in high_identity_alignments:
                    high_identity_alignments[query_name] = []
                
                ref_id = alignment.reference_name[alignment.reference_name.find('_') + 1:] if alignment.is_forward else '-' + alignment.reference_name[alignment.reference_name.find('_') + 1:]
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
                    "chrY_pat_hsaY": "Y"
                }
                if alignment.reference_name in id_map:
                    ref_id = id_map[alignment.reference_name]
                high_identity_alignments[query_name].append({
                    'identity': identity,
                    'ref_id': ref_id if alignment.is_forward else '-' + ref_id,
                    'ref_start': alignment.reference_start if alignment.is_forward else bamfile.get_reference_length(alignment.reference_name) - alignment.reference_end,
                    'ref_end': alignment.reference_end if alignment.is_forward else bamfile.get_reference_length(alignment.reference_name) - alignment.reference_start,
                    'reverse': alignment.is_reverse,
                    'length': alignment.query_alignment_length,
                    'alignment': alignment
                })

                query_name = alignment.query_name.split('_')[1]
                if query_name not in high_identity_alignments:
                    high_identity_alignments[query_name] = []
                
                high_identity_alignments[query_name].append({
                    'identity': identity,
                    'ref_id': ref_id if alignment.is_reverse else '-' + ref_id,
                    'ref_start': alignment.reference_start if alignment.is_reverse else bamfile.get_reference_length(alignment.reference_name) - alignment.reference_end,
                    'ref_end': alignment.reference_end if alignment.is_reverse else bamfile.get_reference_length(alignment.reference_name) - alignment.reference_start,
                    'reverse': alignment.is_reverse,
                    'length': alignment.query_alignment_length,
                    'alignment': alignment
                })
        
        for query_name, alignments in high_identity_alignments.items():
            alignments = sorted(alignments, key=lambda x: (x["ref_id"], x["ref_start"], -x["ref_end"], -x["identity"]))
            # for x in alignments:
            #     print(query_name, x['ref_id'], x['ref_start'], x['ref_end'], x['length'], x['identity'], sep='\t')
            new_list = []
            for aln in alignments:
                if len(new_list) == 0 or new_list[-1]["ref_id"] != aln["ref_id"]:
                    new_list.append(aln)
                elif aln["ref_end"] > new_list[-1]["ref_end"]:
                    if aln["ref_start"] < new_list[-1]["ref_end"] + 100000:
                        new_list[-1]["ref_end"] = aln["ref_end"]
                        new_list[-1]["identity"] = (new_list[-1]["identity"]*new_list[-1]["length"] + aln["identity"]*aln["length"])/(new_list[-1]["length"] + aln["length"])
                        new_list[-1]["length"] = new_list[-1]["length"] + aln["length"]
                        # new_list[-1]["length"] = new_list[-1]["length"] - (new_list[-1]["ref_end"] - aln["ref_start"]) if new_list[-1]["ref_end"] > aln["ref_start"] else new_list[-1]["length"]
                    else:
                        new_list.append(aln)
                
            high_identity_alignments[query_name] =  sorted(new_list, key=lambda x: x["ref_start"]-x["ref_end"])
    
    
    print(f"Edges with label: {processed_alignments} of {total_alignments}")
    return high_identity_alignments

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Extract high-identity alignments from a BAM file.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("fasta_file", help="Path to graph.fasta file.")
    parser.add_argument("-t", "--threshold", type=float, default=0,
                        help="Identity threshold (default: 0).")
    parser.add_argument("-o", "--output", help="Path to the output file. If not specified, prints to stdout.")
    return parser.parse_args()

def main():
    args = parse_arguments()
    bam_file_path = args.bam_file
    threshold = args.threshold
    output_path = args.output
    fasta_file = pysam.FastaFile(args.fasta_file)
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
            if "tig" not in aln['ref_id'] and aln['ref_end']-aln['ref_start'] > contig_lengths[query_name]/10:
                output_lines.append(f"{query_name}\t{aln['ref_id']}\t{aln['length']}({aln['ref_start']}-{aln['ref_end']})\t{aln['identity']:.2f}")
    
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
