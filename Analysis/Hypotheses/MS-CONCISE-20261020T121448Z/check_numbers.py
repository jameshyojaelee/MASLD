#!/usr/bin/env python3
"""Check that every number in each baseline file survives in its staged file or Supplementary Note."""

import re
import sys
from collections import Counter
from pathlib import Path

RUN = Path(__file__).resolve().parent
SECTION = {
    "title.md": "front", "abstract.md": "front", "01_intro.md": "front",
    "fig1.md": "res12", "fig2.md": "res12", "fig3.md": "res3", "fig4.md": "res4",
    "fig5.md": "res56disc", "fig5_discussion.md": "res56disc", "fig7.md": "res7",
    "05_figure_legends.md": "legends", "supplementary_notes.md": "notes_methods",
    "METHODS.md": "notes_methods",
}
NUM = re.compile(r"(?<![A-Za-z0-9])\d[\d,]*(?:\.\d+)?(?:e-?\d+)?")
BANNED = re.compile(r"—|\bfrozen\b|\bcascade\b|\bprogression\b|causal gene", re.I)


def numbers(text: str) -> Counter:
    return Counter(t.replace(",", "") for t in NUM.findall(text))


def main() -> None:
    failed = False
    for name, key in SECTION.items():
        base = (RUN / "baseline" / name).read_text()
        staged_path = RUN / "staged" / name
        if not staged_path.exists():
            print(f"{name}\tMISSING_STAGED")
            failed = True
            continue
        staged = staged_path.read_text()
        note = RUN / "staged" / "supplementary_notes.md"
        pool = numbers(staged) + (numbers(note.read_text()) if note.exists() else Counter())
        lost = sorted(set(numbers(base)) - set(pool))
        banned = sorted({m.group(0) for m in BANNED.finditer(staged)})
        print(f"{name}\twords {len(base.split())}->{len(staged.split())}\t"
              f"distinct_numbers_lost {len(lost)}\tbanned {banned}")
        if lost:
            print(f"  lost: {' '.join(lost[:60])}")
            failed = True
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
