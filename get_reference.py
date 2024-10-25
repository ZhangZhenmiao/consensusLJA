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
                
                high_identity_alignments[query_name].append({
                    'identity': identity,
                    # 'ref_id': alignment.reference_name[:alignment.reference_name.rfind('_')] if alignment.is_forward else '-' + alignment.reference_name[:alignment.reference_name.rfind('_')],
                    'ref_id': alignment.reference_name[alignment.reference_name.find('_') + 1:] if alignment.is_forward else '-' + alignment.reference_name[alignment.reference_name.find('_') + 1:],
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
                    # 'ref_id': alignment.reference_name[:alignment.reference_name.rfind('_')] if alignment.is_reverse else '-' + alignment.reference_name[:alignment.reference_name.rfind('_')],
                    'ref_id': alignment.reference_name[alignment.reference_name.find('_') + 1:] if alignment.is_reverse else '-' + alignment.reference_name[alignment.reference_name.find('_') + 1:],
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
