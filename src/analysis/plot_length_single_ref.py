#!/usr/bin/env python3
import sys
import re
import matplotlib.pyplot as plt
import numpy as np
import argparse

def parse_span(field):
    """Parse fields like Q:3,386,989(0-35) → (span_bp, ref_start_pct, ref_end_pct)"""
    m = re.search(r"Q:([\d,]+)\((\d+)-(\d+)\)", field)
    if not m:
        return 0, None, None
    span = int(m.group(1).replace(",", ""))
    ref_start = int(m.group(2))
    ref_end = int(m.group(3))
    return span, ref_start, ref_end

def parse_ref(raw_ref):
    """Extract (ref_name, is_reverse)"""
    is_rev = raw_ref.startswith("-")
    core = raw_ref[1:] if is_rev else raw_ref
    core = core.rstrip(")")
    m = re.search(r"(\d+)", core)
    ref_name = m.group(1) if m else core
    return ref_name, is_rev

id_map = {
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
    "NC_091271.1": "ChrY",
    "NC_091727.1": "ChrX"
}

def read_fai(fai_file):
    lengths = {}
    with open(fai_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                name, length = parts[0], int(parts[1])
                if name in id_map:
                    name = id_map[name]
                lengths[name] = length
    return lengths

def load_exclude_list(path):
    """Load contigs to exclude (first column only)."""
    exclude = set()
    if path is None:
        return exclude
    with open(path) as f:
        for line in f:
            if line.strip() and "len" in line:
                exclude.add(line.split()[0])
    return exclude

def read_stats(stats_file, contig_lengths, exclude):
    best_refs = {}
    with open(stats_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            contig, raw_ref, span_field = parts[0], parts[1], parts[2]

            if contig in exclude:
                continue  # skip excluded contigs

            span, ref_start, ref_end = parse_span(span_field)
            if ref_start is None:
                continue

            ref_name, is_rev = parse_ref(raw_ref)
            if 'mtDNA' in ref_name or 'chrM' in ref_name:
                continue  # skip mitochondrial mappings
            if 'X' not in ref_name and 'Y' not in ref_name and int(ref_name) > 50:
                continue  # skip non-chromosomal mappings
            if is_rev:
                ref_start, ref_end = (100 - ref_end), (100 - ref_start)
            
            if ref_end - ref_start < 2:
                continue

            if contig not in best_refs or span > best_refs[contig][1]:
                best_refs[contig] = (ref_name, span, ref_start, ref_end)

    print(f"\nBest refs for {stats_file}:")
    for contig, (ref, span, start, end) in best_refs.items():
        print(f"{contig}\t{ref}\t{span}\t({start}-{end})")

    chr_contigs = {}
    for contig, (ref, span, ref_start, ref_end) in best_refs.items():
        contig_len = contig_lengths.get(contig, 0)
        chr_contigs.setdefault(ref, []).append((ref_start, ref_end, contig_len))
    return chr_contigs

def plot_single(ax, stats_file, contigs_fai, stats2_file, contigs2_fai, ref_fai, exclude_file1, exclude_file2, add_dots=False):
    contig_lengths1 = read_fai(contigs_fai)
    contig_lengths2 = read_fai(contigs2_fai)
    ref_lengths = read_fai(ref_fai)
    exclude1 = load_exclude_list(exclude_file1)
    exclude2 = load_exclude_list(exclude_file2)

    chr_contigs1 = read_stats(stats_file, contig_lengths1, exclude1)
    chr_contigs2 = read_stats(stats2_file, contig_lengths2, exclude2)

    def chr_sort_key(c):
        m = re.search(r'\d+', c)
        if c == "X": return 1000
        if c == "Y": return 2000
        return int(m.group()) if m else float('inf')

    chromosomes_sorted = sorted(set(list(chr_contigs1.keys()) + list(chr_contigs2.keys())),
                                key=chr_sort_key)
    x = np.arange(len(chromosomes_sorted))

    hapA_bp = []
    hapB_bp = []
    for c in chromosomes_sorted:
        valA = ref_lengths.get(f"chromosome_{c}A") or next((ref_lengths[k] for k in ref_lengths if f"Chr{c}" in k), 0)
        hapA_bp.append(valA)
    hapA_mb = [v/1e6 for v in hapA_bp]
    hapB_mb = [v/1e6 for v in hapB_bp]

    total_width = 0.8
    hapA_added = False
    hapB_added = False

    bar_height_max = 0
    for i, c in enumerate(chromosomes_sorted):
        slots = []
        slots.append(('hapA', hapA_mb[i], "#4C72B0"))
        slots.append(('cons1', chr_contigs1.get(c, []), "#C44E52"))
        slots.append(('cons2', chr_contigs2.get(c, []), "#8172B3"))

        n_slots = len(slots)
        slot_width = total_width / n_slots
        start_pos = x[i] - total_width / 2

        for j, (name, val, color) in enumerate(slots):
            if 'cons' in name:
                bottom = 0.0
                for rs, re_, clen in sorted(val, key=lambda t: (t[0], t[1])):
                    ax.bar(start_pos + j*slot_width, clen/1e6, slot_width, bottom=bottom,
                           color=color, edgecolor='white', linewidth=0.3,
                           label=name if (i==0 and j==0 and name=='cons1') else None)
                    bottom += clen/1e6

                # --- Add dots for contig counts (only if add_dots=True) ---
                if add_dots and len(val) > 0:
                    # count only contigs longer than 1 Mb
                    long_contigs = [clen for rs, re_, clen in val if clen > 1e6]
                    if long_contigs:
                        # get the consensus bar total height (Mb)
                        bar_height = sum(clen/1e6 for rs, re_, clen in val)
                        if bar_height > bar_height_max:
                            bar_height_max = bar_height

    for i, c in enumerate(chromosomes_sorted):
        slots = []
        slots.append(('hapA', hapA_mb[i], "#4C72B0"))
        slots.append(('cons1', chr_contigs1.get(c, []), "#C44E52"))
        slots.append(('cons2', chr_contigs2.get(c, []), "#8172B3"))

        n_slots = len(slots)
        slot_width = total_width / n_slots
        start_pos = x[i] - total_width / 2

        for j, (name, val, color) in enumerate(slots):
            if 'cons' in name:
                bottom = 0.0
                for rs, re_, clen in sorted(val, key=lambda t: (t[0], t[1])):
                    ax.bar(start_pos + j*slot_width, clen/1e6, slot_width, bottom=bottom,
                           color=color, edgecolor='white', linewidth=0.3,
                           label=name if (i==0 and j==0 and name=='cons1') else None)
                    bottom += clen/1e6

                # --- Add dots for contig counts (only if add_dots=True) ---
                if add_dots and len(val) > 0:
                    # count only contigs longer than 1 Mb
                    long_contigs = [clen for rs, re_, clen in val if clen > 1e6]
                    if long_contigs:
                        # get the consensus bar total height (Mb)
                        # bar_height = sum(clen/1e6 for rs, re_, clen in val)
                        dot_start = bar_height_max + 5.0  # start 5 Mb above bar top
                        dot_spacing = 5.0             # vertical spacing (Mb)

                        if 'cons1' in name:
                            for k in range(len(long_contigs)):
                                ax.scatter(
                                    start_pos + j*slot_width,
                                    dot_start + k*dot_spacing,
                                    s=15, c="#C44E52", marker="o", zorder=5
                                )
                        else:  # 'cons2'
                            for k in range(len(long_contigs)):
                                ax.scatter(
                                    start_pos + j*slot_width,
                                    dot_start + k*dot_spacing,
                                    s=15, c="#8172B3", marker="o", zorder=5
                                )

            else:
                label = None
                if name=='hapA' and not hapA_added:
                    label='Reference Chromosome'; hapA_added=True
                    # label='Haplome A'; hapA_added=True
                if name=='hapB' and not hapB_added:
                    label='Paternal Haplome'; hapB_added=True
                    # label='Haplome B'; hapB_added=True
                ax.bar(start_pos + j*slot_width, val, slot_width, color=color, edgecolor='white', linewidth=0.3, label=label)

    # Dummy bars for consensus legend
    ax.bar(0,0,color="#C44E52", label='MGA', edgecolor='white')
    ax.bar(0,0,color="#8172B3", label='hifiasm', edgecolor='white')
    ax.set_xticks(x)
    ax.set_xticklabels([f"chr{c}" for c in chromosomes_sorted], fontsize=10)
    ax.grid(False)
    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='x', labelsize=10)
    return ax

def main():
    parser = argparse.ArgumentParser(description="Integrate three contig comparison plots into one figure")
    # parser.add_argument('--stats', nargs=3, required=True)
    # parser.add_argument('--fai', nargs=3, required=True)
    # parser.add_argument('--stats2', nargs=3, required=True)
    # parser.add_argument('--fai2', nargs=3, required=True)
    # parser.add_argument('--ref', nargs=3, required=True)
    # parser.add_argument('--exclude1', nargs=3, required=True, help="Files listing contigs to exclude (3)")
    # parser.add_argument('--exclude2', nargs=3, required=True, help="Files listing contigs to exclude (3)")
    # parser.add_argument('--output', required=True)
    parser.add_argument('--stats', required=True)
    parser.add_argument('--fai', required=True)
    parser.add_argument('--stats2', required=True)
    parser.add_argument('--fai2', required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--exclude1', required=True, help="Files listing contigs to exclude (3)")
    parser.add_argument('--exclude2', required=True, help="Files listing contigs to exclude (3)")
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # fig, axes = plt.subplots(3, 1, figsize=(16, 18))

    # for i in range(3):
    #     plot_single(
    #         axes[i],
    #         args.stats[i], args.fai[i],
    #         args.stats2[i], args.fai2[i],
    #         args.ref[i], args.exclude1[i], args.exclude2[i],
    #         add_dots=(i == 2)
    #     )
    #     axes[i].set_ylabel("Length (Mb)", fontsize=14)
    #     axes[i].set_xlabel("Chromosome", fontsize=14)

    # handles, labels = axes[0].get_legend_handles_labels()
    # fig.legend(handles, labels, fontsize=12, frameon=False, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.01))
    # fig.tight_layout(rect=[0, 0.05, 1, 1])
    # plt.savefig(args.output, dpi=300)
    fig, ax = plt.subplots(figsize=(16, 6))   # Only one axis now

    plot_single(
        ax,
        args.stats, args.fai,
        args.stats2, args.fai2,
        args.ref, args.exclude1, args.exclude2,
        add_dots=True        # or False — whatever you want
    )

    ax.set_ylabel("Length (Mb)", fontsize=14)
    ax.set_xlabel("Chromosome", fontsize=14)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=12, frameon=False,
            loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.01))

    fig.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(args.output, dpi=300)

if __name__=="__main__":
    main()
