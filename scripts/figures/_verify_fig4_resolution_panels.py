#!/usr/bin/env python3
"""Verify the two new Figure 4 panels.

Two checks, using the tools that actually work on this cluster:
  1. Text content via ~/.local/bin/pdftotext (NOT `strings` -- PDF text lives in
     compressed streams, so `strings` cannot read it and proves nothing).
  2. Bold glyph USE via pypdf: resolve every /Font resource, map the Tf operator
     to its BaseFont, and report whether any text-showing operator (Tj/TJ/'/")
     follows a Bold selection. (`strings | grep Helvetica-Bold` gives FALSE
     POSITIVES -- cairo_pdf registers the whole family whether or not bold is
     ever drawn.)
"""
import re
import subprocess
import sys
from pathlib import Path

from pypdf import PdfReader
from pypdf.generic import ContentStream

PDFTOTEXT = Path.home() / ".local/bin/pdftotext"
BANNED = [
    (r"\bTREAT\b", "bare TREAT"),
    (r"tested[ _]negative", "tested negative"),
    (r"\bmega\b", "mega"),
    (r"\bno effect\b", "no effect"),
    (r"\bno difference\b", "no difference"),
    (r"prospectively validated", "prospectively validated"),
]


def extract_text(pdf: Path) -> str:
    out = subprocess.run(
        [str(PDFTOTEXT), "-layout", str(pdf), "-"],
        capture_output=True, text=True, check=True,
    )
    return out.stdout


def bold_glyphs_drawn(pdf: Path):
    """Return list of (page, basefont) where a text-showing op follows a Bold Tf."""
    reader = PdfReader(str(pdf))
    hits = []
    registered = set()
    for pno, page in enumerate(reader.pages, start=1):
        res = page.get("/Resources")
        fonts = {}
        if res is not None:
            fdict = res.get_object().get("/Font")
            if fdict is not None:
                for key, ref in fdict.get_object().items():
                    base = ref.get_object().get("/BaseFont")
                    base = str(base) if base is not None else "?"
                    fonts[str(key)] = base
                    registered.add(base)
        cs = ContentStream(page.get_contents(), reader)
        current = None
        for operands, op in cs.operations:
            o = op.decode("utf-8", "replace") if isinstance(op, bytes) else str(op)
            if o == "Tf" and operands:
                current = fonts.get(str(operands[0]), str(operands[0]))
            elif o in ("Tj", "TJ", "'", '"'):
                if current and "bold" in current.lower():
                    hits.append((pno, current))
    return hits, sorted(registered)


def main(paths):
    rc = 0
    for p in paths:
        pdf = Path(p)
        print(f"\n=== {pdf} ===")
        if not pdf.exists():
            print("  MISSING")
            rc = 1
            continue
        print(f"  size: {pdf.stat().st_size} bytes")

        txt = extract_text(pdf)
        print("  --- pdftotext content ---")
        for line in txt.splitlines():
            if line.strip():
                print("   |", line.rstrip())

        for pat, label in BANNED:
            if re.search(pat, txt, flags=re.IGNORECASE if label == "mega" else 0):
                print(f"  BANNED STRING PRESENT: {label}")
                rc = 1
        print("  banned-string sweep: clean" if rc == 0 else "  banned-string sweep: FAILED")

        hits, registered = bold_glyphs_drawn(pdf)
        print(f"  fonts registered: {registered}")
        if hits:
            print(f"  BOLD GLYPHS DRAWN: {hits}")
            rc = 1
        else:
            print("  bold glyphs drawn: none")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
