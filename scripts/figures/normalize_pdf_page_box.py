#!/usr/bin/env python3
"""Set an exact one-page PDF box without scaling or rewriting its content stream."""

from __future__ import annotations

import argparse
import os
import tempfile
from decimal import Decimal
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.generic import RectangleObject


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("width_in", type=Decimal)
    parser.add_argument("height_in", type=Decimal)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pdf = args.pdf.resolve()
    reader = PdfReader(str(pdf))
    if len(reader.pages) != 1:
        raise SystemExit(f"expected one page, found {len(reader.pages)}: {pdf}")

    width_pt = args.width_in * Decimal(72)
    height_pt = args.height_in * Decimal(72)
    page = reader.pages[0]
    exact = RectangleObject([Decimal(0), Decimal(0), width_pt, height_pt])
    page.mediabox = exact
    page.cropbox = RectangleObject(exact)
    page.trimbox = RectangleObject(exact)
    page.bleedbox = RectangleObject(exact)
    page.artbox = RectangleObject(exact)

    writer = PdfWriter()
    writer.add_page(page)
    if reader.metadata:
        writer.add_metadata({str(k): str(v) for k, v in reader.metadata.items() if v is not None})

    original_mode = pdf.stat().st_mode & 0o7777
    fd, tmp_name = tempfile.mkstemp(prefix=f".{pdf.name}.", suffix=".tmp", dir=pdf.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    # mkstemp creates 0600; os.replace would carry that onto the published panel
    # and make it unreadable to the nslab group.
    os.chmod(tmp, original_mode)
    try:
        with tmp.open("wb") as handle:
            writer.write(handle)
        check = PdfReader(str(tmp))
        box = check.pages[0].mediabox
        if Decimal(str(float(box.width))) != width_pt or Decimal(str(float(box.height))) != height_pt:
            raise RuntimeError(
                f"page-box normalization failed: {float(box.width)} x {float(box.height)} pt"
            )
        os.replace(tmp, pdf)
    finally:
        if tmp.exists():
            tmp.unlink()

    print(f"[page-box] {pdf.name}: {width_pt} x {height_pt} pt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
