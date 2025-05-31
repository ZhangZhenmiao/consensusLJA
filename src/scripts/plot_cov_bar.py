import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Parse command line arguments
parser = argparse.ArgumentParser(description="Plot stacked bar chart of edge count by coverage.")
parser.add_argument("input_file", type=str, help="Path to the input file (two columns: length and coverage)")
parser.add_argument("output_file", type=str, help="Path to save the output plot (e.g. output.png)")
args = parser.parse_args()

# Load data
data = np.loadtxt(args.input_file)
df = pd.DataFrame(data, columns=['length', 'coverage'])
df['group'] = np.where(df['length'] <= 10000, 'Length ≤ 10,000', 'Length > 10,000')

# Bin coverage
bins = np.arange(0, 101, 1)
df['coverage_bin'] = pd.cut(df['coverage'], bins=bins, right=False, labels=bins[:-1])

# Count edges in each bin by group
coverage_counts = df.groupby(['coverage_bin', 'group']).size().unstack(fill_value=0)

# Plot settings
plt.figure(figsize=(9, 6))

# Prepare stacked bar chart
x = coverage_counts.index.astype(int)
bottom = np.zeros(len(x))
colors = ["#3A6AB8", "#D56F34"]  # Professional blue and orange

for idx, group in enumerate(coverage_counts.columns):
    plt.bar(
        x,
        coverage_counts[group],
        bottom=bottom,
        label=group,
        color=colors[idx],
        width=0.8,
        edgecolor='black',  # Border color
        linewidth=0.3        # Thin border
    )
    bottom += coverage_counts[group].values

# Axis formatting
plt.xlim(0, 100)
plt.xticks(np.arange(0, 101, 10))
plt.xlabel('Coverage')
plt.ylabel('Number of Edges')
plt.legend(frameon=False, title=None)
plt.xticks(np.arange(0, 101, 10))
plt.yticks(fontsize=10)  # Optional: control font size

# Layout and save
plt.tight_layout()
plt.grid(False)
plt.savefig(args.output_file)
