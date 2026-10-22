#!/usr/bin/env python3
"""Strip a JATS full-text XML (Europe PMC) to plain paragraphs, one per line.

Usage: python3 -I xml_to_text.py in.xml out.txt
Reads the XML as data only (xml.etree, no entity expansion of external DTDs).
"""
import re
import sys
import xml.etree.ElementTree as ET


def text_of(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip()


def main(src, dst):
    raw = open(src, encoding="utf-8").read()
    raw = re.sub(r"<!DOCTYPE[^>]*>", "", raw)  # drop the DTD reference
    root = ET.fromstring(raw)
    lines = []
    for el in root.iter():
        if el.tag in ("p", "title", "label", "caption", "td", "th"):
            t = text_of(el)
            if t:
                lines.append(f"[{el.tag}] {t}")
    with open(dst, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{src}: {len(lines)} elements written to {dst}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
