import re
import sys
from collections import defaultdict

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <input_file>")
    sys.exit(1)

filename = sys.argv[1]

# Store results
chrom_data = defaultdict(dict)

with open(filename) as f:
    for line in f:
        if "finish" not in line:
            continue
        parts = line.strip().split()
        if len(parts) != 5:
            continue
        chrom = parts[1]              # e.g., chromosome_15A
        match = re.search(r'_(\d+)', chrom)
        if not match:
            continue
        chrom_num = int(match.group(1))
        group = chrom[-1]             # A or B
        second_last = float(parts[-2])  # the second last column
        last = float(parts[-1])       # last column

        # convert
        second_last_val = round(second_last * 100)
        last_val = round(last * 100)

        chrom_data[chrom_num][group] = (second_last_val, last_val)

print(f"Chr A:")
sum_val1 = 0
sum_val2 = 0
for num in range(1, 19):
    if num not in chrom_data:
        continue
    groups = chrom_data[num]
    val1, val2 = groups['A']
    sum_val1 += val1
    sum_val2 += val2
    print(f"{val1}/{val2}", end='\t')
    
avg1 = round(sum_val1 / 18)
avg2 = round(sum_val2 / 18)
print(f"{avg1}/{avg2}")

print(f"Chr B:")
sum_val1 = 0
sum_val2 = 0
for num in range(1, 19):
    if num not in chrom_data:
        continue
    groups = chrom_data[num]
    val1, val2 = groups['B']
    sum_val1 += val1
    sum_val2 += val2
    print(f"{val1}/{val2}", end='\t')
    
avg1 = round(sum_val1 / 18)
avg2 = round(sum_val2 / 18)
print(f"{avg1}/{avg2}")
