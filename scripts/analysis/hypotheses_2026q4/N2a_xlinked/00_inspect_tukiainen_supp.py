"""Print each sheet's first rows of the already-extracted Tukiainen 2017 Supplementary Information file 3.
Usage: python -I 00_inspect_tukiainen_supp.py EXTRACTED_DIR/nature24265-s3"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xlsx_min import read_sheets
base = sys.argv[1]
for f in sorted(os.listdir(base), key=lambda s: (len(s), s)):
    if not f.endswith(".xlsx"):
        continue
    for sh, rows in read_sheets(os.path.join(base, f), max_rows=4).items():
        print("==", f, "|", sh)
        for r in rows[:3]:
            print("   ", [x[:28] for x in r[:9]])
