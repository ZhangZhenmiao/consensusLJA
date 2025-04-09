import matplotlib.pyplot as plt
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description="Plot stacked bar chart of edge count by coverage.")
parser.add_argument("input_file", type=str, help="Path to the input file (two columns: length and coverage)")
parser.add_argument("output_file", type=str, help="Path to save the output plot (e.g. output.png)")
args = parser.parse_args()

# Load data
data = np.loadtxt(args.input_file)
lengths = data[:, 0]
coverages = data[:, 1]

# Define coverage bins (customize if needed)
bins = np.arange(0, max(coverages) + 1, 1)

# Split data into two groups
short_mask = lengths <= 10000
long_mask = lengths > 10000

# Histogram counts for each group
hist_short, _ = np.histogram(coverages[short_mask], bins=bins)
hist_long, _ = np.histogram(coverages[long_mask], bins=bins)

# X positions for each bar group
x = bins[:-1]

# Plotting
plt.figure(figsize=(12, 6))
plt.bar(x, hist_short, width=np.diff(bins), align='edge', color='#4C72B0', label='Length ≤ 10,000')  # Professional blue
plt.bar(x, hist_long, width=np.diff(bins), align='edge', bottom=hist_short, color='#DD8452', label='Length > 10,000')  # Professional orange
plt.xlim(0, 100)
plt.xticks(np.arange(0, 101, 10))

plt.xlabel('Coverage')
plt.ylabel('Number of Edges')
plt.title('Stacked Coverage Distribution of Edges by Length Group')
plt.legend()
plt.tight_layout()
plt.savefig(args.output_file)
