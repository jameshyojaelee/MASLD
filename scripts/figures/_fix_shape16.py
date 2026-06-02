#!/usr/bin/env python3
"""Add shape = 16 to all geom_point() calls that don't already have a shape parameter."""
import re
import glob
import os

FIG_DIR = os.path.dirname(os.path.abspath(__file__))

# Skip these specific patterns (annotation markers, UpSet dots, etc.)
SKIP_SHAPES = {"shape = 18", "shape = 21", "shape=18", "shape=21"}

def find_matching_paren(text, start):
    """Find the closing ) that matches the ( at position start."""
    depth = 0
    i = start
    while i < len(text):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1

def process_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    original = content
    edits = 0

    # Find all geom_point( occurrences
    pattern = re.compile(r'geom_point\(')
    offset = 0

    while True:
        match = pattern.search(content, offset)
        if not match:
            break

        start_paren = match.end() - 1  # position of (
        end_paren = find_matching_paren(content, start_paren)

        if end_paren == -1:
            offset = match.end()
            continue

        # Extract the full geom_point(...) call
        call_content = content[start_paren + 1:end_paren]

        # Skip if already has shape parameter (fixed or in aes)
        if re.search(r'shape\s*=', call_content):
            offset = end_paren + 1
            continue

        # Skip if it has aes(shape ...) — shape as unmapped aesthetic
        if re.search(r'shape', call_content):
            offset = end_paren + 1
            continue

        # Add shape = 16 before the closing paren
        stripped = call_content.rstrip()
        if stripped == '':
            # geom_point() -> geom_point(shape = 16)
            new_call = 'shape = 16'
        else:
            # Add as last parameter
            new_call = call_content.rstrip() + ', shape = 16'
            # Preserve trailing whitespace/newline
            trailing = call_content[len(call_content.rstrip()):]
            new_call = new_call + trailing

        content = content[:start_paren + 1] + new_call + content[end_paren:]
        edits += 1

        # Adjust offset for the new content length
        offset = start_paren + 1 + len(new_call) + 1

    if edits > 0:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"  {os.path.basename(filepath)}: {edits} geom_point calls updated")

    return edits

# Process all R files in figures directory
r_files = sorted(glob.glob(os.path.join(FIG_DIR, "*.R")))
total = 0
for f in r_files:
    # Skip utility files
    basename = os.path.basename(f)
    if basename in ("publication_theme.R", "load_figure_data.R", "_fix_shape16.py"):
        continue
    n = process_file(f)
    total += n

print(f"\nTotal: {total} geom_point calls updated across {len(r_files)} files")
