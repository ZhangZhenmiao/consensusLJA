import re
import os
import sys

def process_content(text):
    """Finds (float) and replaces with (rounded_int)."""
    pattern = r"\(([-+]?\d*\.\d+)\)"
    
    def round_match(match):
        num_str = match.group(1)
        # Using int(float(x) + 0.5) if you want traditional rounding, 
        # or just round() for standard Python behavior.
        rounded_val = int(round(float(num_str)))
        return f"({rounded_val})"
    
    return re.sub(pattern, round_match, text)

def update_file(file_path):
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    # 1. Read the data
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # 2. Process the data
    updated_content = process_content(content)

    # 3. Write it back
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(updated_content)
    
    print(f"Success: {file_path} has been updated.")

# Usage
update_file(sys.argv[1])