#!/usr/bin/env python3
"""Acquire and safely extract GSE238219 replicate methods from the source DOCX."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import urllib.request
import xml.etree.ElementTree as ET
import zipfile


SOURCE_URL = "https://cdn-links.lww.com/permalink/hep/i/hep_2024_08_20_knowles_hep-24-0583_sdc1.docx"
KEYWORDS = (
    "perturb-seq",
    "perturb seq",
    "perturbseq",
    "replicate",
    "crispr",
    "10x",
    "chromium",
)


class ExtendedMethodsError(RuntimeError):
    """Raised when the source document differs from the bounded DOCX contract."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def download(target: Path) -> None:
    request = urllib.request.Request(
        SOURCE_URL,
        headers={"User-Agent": "MASLD-Model-Benchmark/1.0 source-audit"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:  # noqa: S310
        payload = response.read(2 * 1024 * 1024 + 1)
        content_type = response.headers.get_content_type()
    if len(payload) > 2 * 1024 * 1024:
        raise ExtendedMethodsError("supplement exceeds 2 MiB audit ceiling")
    if content_type != "application/vnd.openxmlformats-officedocument.wordprocessingml.document":
        raise ExtendedMethodsError(f"unexpected supplement content type: {content_type}")
    target.write_bytes(payload)


def read_paragraphs(source: Path) -> tuple[list[str], list[str]]:
    with zipfile.ZipFile(source) as archive:
        names = archive.namelist()
        for name in names:
            pure = PurePosixPath(name)
            if pure.is_absolute() or ".." in pure.parts:
                raise ExtendedMethodsError(f"unsafe DOCX member: {name}")
            if name.lower().endswith((".bin", ".exe", ".dll", ".so")):
                raise ExtendedMethodsError(f"executable DOCX member: {name}")
        if "word/document.xml" not in names:
            raise ExtendedMethodsError("DOCX has no word/document.xml")
        info = archive.getinfo("word/document.xml")
        if info.file_size > 10 * 1024 * 1024 or info.compress_size == 0:
            raise ExtendedMethodsError("document XML exceeds bounded loading contract")
        if info.file_size / info.compress_size > 100:
            raise ExtendedMethodsError("document XML compression ratio exceeds ceiling")
        xml_payload = archive.read(info)
    root = ET.fromstring(xml_payload)
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    paragraphs = []
    for paragraph in root.findall(".//w:p", namespace):
        text = "".join(node.text or "" for node in paragraph.findall(".//w:t", namespace)).strip()
        if text:
            paragraphs.append(text)
    relevant = [
        text
        for text in paragraphs
        if any(keyword in text.casefold() for keyword in KEYWORDS)
    ]
    if not relevant:
        raise ExtendedMethodsError("no perturbation or replicate methods found")
    return paragraphs, relevant


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ExtendedMethodsError("refusing to overwrite extended-methods audit")
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.output / "hep-24-0583-sdc1.docx"
    download(source)
    paragraphs, relevant = read_paragraphs(source)
    record = {
        "schema_version": "masld-bench-gse238219-extended-methods-audit-v1",
        "source_url": SOURCE_URL,
        "source_sha256": digest(source),
        "source_size_bytes": source.stat().st_size,
        "document_paragraph_count": len(paragraphs),
        "relevant_paragraph_count": len(relevant),
        "relevant_paragraphs": relevant,
        "purpose": "resolve culture, transduction, capture, sequencing, and biological-replicate ancestry for GSE238219",
        "gse313774_accessed": False,
    }
    (args.output / "extended_methods_audit.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
