#!/usr/bin/env python
import argparse
import sys
import os

def classify_relationship(aln):
    """Return containment/overlap classification for query and ref."""
    identity = aln["identity"]

    # query contained in ref
    query_contained = (
        aln["aln_length_query"] / aln["length_query"] >= 0.9 and identity >= 50 and aln["length_query"] <= aln["length_ref"]
    )

    # ref contained in query
    ref_contained = (
        aln["aln_length_ref"] / aln["length_ref"] >= 0.9 and identity >= 50 and aln["length_ref"] <= aln["length_query"]
    )

    # overlap definition
    overlap = (
        aln["aln_length_query"] / aln["length_query"] >= 0.2
        and aln["aln_length_ref"] / aln["length_ref"] >= 0.2
        and identity >= 50
    )

    if query_contained:
        return "Contained"
    elif overlap:
        return "Overlap"
    else:
        return "None"

def generate_paf(query_fasta, ref_fasta, output_paf, threads=8, minimap2_preset="asm20"):
    """
    Generate paf file from query vs reference using minimap2.
    """
    cmd = (
        f"minimap2 -t {threads} -x {minimap2_preset} --eqx -Y -p 0.1 {ref_fasta} {query_fasta} > {output_paf}"
    )
    if os.system(cmd) != 0:
        raise RuntimeError("paf generation failed")

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
        tag, typ, val = f.split(":")
        paf[tag] = val
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


def calculate_identity(paf_entry, gap_threshold=10):
    """
    Identity = matches / aln_len.
    Optionally could parse 'cg:Z:' tag for CIGAR to refine.
    """
    matches = paf_entry["matches"]
    aln_len = paf_entry["aln_len"]
    if aln_len == 0:
        return 0.0, 0.0
    identity = matches / aln_len

    # if we want no-gap identity, need CIGAR string (cg:Z:), otherwise fall back to same value
    identity_nogap = identity
    return identity * 100, identity_nogap * 100

def filter_alignments_with_identity(paf_file_path, threshold=0):
    high_identity_alignments = {}
    total_alignments = 0
    processed_alignments = 0
    best_hits = {}

    with open(paf_file_path) as fh:
        for line in fh:
            paf = parse_paf_line(line)
            if paf is None:
                continue
            total_alignments += 1

            identity, identity_nogap = calculate_identity(paf, gap_threshold=10)
            if identity_nogap < threshold * 100:
                continue
            processed_alignments += 1

            if paf["query_name"] == paf["ref_name"]:
                continue  # skip self-alignments

            query_name1, query_name2 = paf["query_name"].split('_')
            ref_name1, ref_name2 = paf["ref_name"].split('_')

            if query_name1 not in high_identity_alignments:
                high_identity_alignments[query_name1] = []
            
            aln = {
                'identity': identity,
                'identity_nogap': identity_nogap,
                'query_id': query_name1,
                'ref_id': ref_name1,
                'ref_start': paf["ref_start"]/lengths[paf["ref_name"]]*100 if paf["strand"] == "+" else (lengths[paf["ref_name"]] - paf["ref_end"])/lengths[paf["ref_name"]]*100,
                'ref_end': paf["ref_end"]/lengths[paf["ref_name"]]*100 if paf["strand"] == "+" else (lengths[paf["ref_name"]] - paf["ref_start"])/lengths[paf["ref_name"]]*100,
                'query_start': paf["query_start"]/lengths[paf["query_name"]]*100,
                'query_end': paf["query_end"]/lengths[paf["query_name"]]*100,
                'strand': "direct" if paf["strand"] == "+" else "reverse",
                'aln_length_query': paf["query_end"] - paf["query_start"],
                'aln_length_ref': paf["ref_end"] - paf["ref_start"],
                'length_query': lengths[paf["query_name"]],
                'length_ref': lengths[paf["ref_name"]],
                'alignment': paf
            }

            high_identity_alignments[query_name1].append(aln)

            status = classify_relationship(aln)
            key = (aln["query_id"], aln["ref_id"], status)

            # keep only the best identity per (query, ref, status)
            if key not in best_hits or aln["identity"] > best_hits[key]["identity"]:
                best_hits[key] = aln

    # output summary
    print("Query\tQuery Len\tQuery Aln Len\tQuery Start\tQuery End\tRef\tRef strand\tRef Len\tRef Aln Len\tRef Start\tRef End\tIdentity\tStatus")
    for (query, ref, status), aln in best_hits.items():
        print(
            f'{query}\t{aln["length_query"]}\t{aln["aln_length_query"]}\t'
            f'{aln["query_start"]:.0f}\t{aln["query_end"]:.0f}\t'
            f'{ref}\t{aln["strand"]}\t{aln["length_ref"]}\t{aln["aln_length_ref"]}\t'
            f'{aln["ref_start"]:.0f}\t{aln["ref_end"]:.0f}\t'
            f'{aln["identity"]:.2f}\t{status}'
        )

    return high_identity_alignments

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True, help="Query FASTA file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--paf", required=True, help="Output paf file")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()

    if not os.path.exists(args.paf):
        generate_paf(args.query, args.ref, args.paf, threads=args.threads)
    
    global lengths 
    lengths = read_fasta_lengths(args.ref)

    filter_alignments_with_identity(args.paf, threshold=0)
