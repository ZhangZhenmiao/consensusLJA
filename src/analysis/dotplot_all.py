import argparse
import os
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from collections import defaultdict
from itertools import groupby
import numpy as np

from utils.os_utils import expandpath, smart_makedirs
from get_conserved_regions import get_conserved_traces

SCRIPT_FN = os.path.realpath(__file__)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("folders", nargs="+", help="Input folders")
    parser.add_argument("-o", "--output", default="combined_dotplots.html", help="Output HTML file")
    return parser.parse_args()

def parse_cigar(cigar_fn):
    X, Y = [0], [0]
    with open(cigar_fn) as f:
        cigar = f.readline().strip()
    cig_iter = groupby(cigar, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int("".join(length_digits))
        op = next(next(cig_iter)[1])
        if op == "M" or op == "X":
            X.append(X[-1] + length)
            Y.append(Y[-1] + length)
        elif op == "D":
            X.append(X[-1] + length)
            Y.append(Y[-1])
        else:
            assert op == "I"
            X.append(X[-1])
            Y.append(Y[-1] + length)
    return X, Y

def load_unialigner_line(folder):
    cigar_fn = os.path.join(folder, "cigar.txt")
    if not os.path.exists(cigar_fn):
        return None, None
    return parse_cigar(cigar_fn)

def main():
    args = parse_args()

    nplots = len(args.folders)
    ncols = int(nplots ** 0.5)
    nrows = (nplots + ncols - 1) // ncols

    fig = make_subplots(rows=nrows, cols=ncols,
                        shared_xaxes=False, shared_yaxes=False,
                        horizontal_spacing=0.02, vertical_spacing=0.02)

    for i, folder in enumerate(args.folders):
        row = i // ncols + 1
        col = i % ncols + 1

        al_X, al_Y = load_unialigner_line(folder)
        if al_X is None:
            continue

        trace = go.Scattergl(
            x=al_X,
            y=al_Y,
            mode="lines",
            line=dict(color="black", width=2),
            showlegend=False
        )
        fig.add_trace(trace, row=row, col=col)

        # Optional: label subplot
        fig.update_xaxes(title_text=os.path.basename(folder), row=row, col=col)
        fig.update_yaxes(title_text="Y", row=row, col=col)

    fig.update_layout(
        height=250*nrows, width=250*ncols,
        title=dict(text="Combined UniAligner Dotplots", font=dict(size=24))
    )

    fig.write_html(args.output)
    print(f"Saved combined interactive HTML: {args.output}")

if __name__ == "__main__":
    main()
