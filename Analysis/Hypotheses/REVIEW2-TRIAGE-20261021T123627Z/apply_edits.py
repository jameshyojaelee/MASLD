"""Apply checked old/new edit lists to copies of the manuscript sources.

Usage: python apply_edits.py <RUN> [--drop-decision D3 ...] [--install]

Reads RUN/staged/final_*.json (and RUN/staged/main_*.json for main-loop edits),
checks that each `old` occurs exactly once in the current source and that no two
edits in one file overlap, then writes edited copies to RUN/staged/files/<path>
and a unified diff to RUN/staged/diff.patch. With --install, copies the edited
files over the sources after archiving the originals to RUN/archive_before/.
"""
import argparse
import difflib
import json
import os
import shutil
import sys
from pathlib import Path

REPO = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")


def load_edits(run, drop, globs):
    edits = []
    files = []
    for g in globs:
        files += sorted((run / "staged").glob(g))
    for f in files:
        d = json.loads(f.read_text())
        rows = d.get("final_edits", d.get("edits", []))
        for e in rows:
            if e.get("decision") in drop:
                continue
            e["_src"] = f.name
            edits.append(e)
    return edits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--drop-decision", nargs="*", default=[])
    ap.add_argument("--glob", nargs="*", default=["final_*.json"])
    ap.add_argument("--base", default=None, help="dir holding the files to edit (default: repo)")
    ap.add_argument("--out", default="files")
    ap.add_argument("--install", action="store_true")
    a = ap.parse_args()
    run = REPO / a.run
    edits = load_edits(run, set(a.drop_decision), a.glob)
    base = Path(a.base) if a.base else REPO

    by_file = {}
    for e in edits:
        by_file.setdefault(e["file"], []).append(e)

    errors, patches = [], []
    out_root = run / "staged" / a.out
    for rel, es in sorted(by_file.items()):
        src = base / rel if (base / rel).exists() else REPO / rel
        text = src.read_text()
        spans = []
        for e in es:
            n = text.count(e["old"])
            if n != 1:
                errors.append(f"{e['_src']}:{e['id']}: old occurs {n} times in {rel}")
                continue
            i = text.index(e["old"])
            spans.append((i, i + len(e["old"]), e))
        spans.sort(key=lambda s: s[0])
        for (a0, a1, ea), (b0, b1, eb) in zip(spans, spans[1:]):
            if b0 < a1:
                errors.append(f"overlap in {rel}: {ea['_src']}:{ea['id']} and {eb['_src']}:{eb['id']}")
        new = text
        for s0, s1, e in sorted(spans, key=lambda s: s[0], reverse=True):
            new = new[:s0] + e["new"] + new[s1:]
        dst = out_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(new)
        patches.extend(difflib.unified_diff(text.splitlines(True), new.splitlines(True), f"a/{rel}", f"b/{rel}"))

    (run / "staged" / f"diff_{a.out}.patch").write_text("".join(patches))
    print(f"{len(edits)} edits over {len(by_file)} files; {len(errors)} errors")
    for err in errors:
        print("ERROR", err)
    if errors:
        sys.exit(1)
    if a.install:
        arch = run / "archive_before"
        for rel in by_file:
            (arch / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO / rel, arch / rel)
            target = REPO / rel
            mode = target.stat().st_mode
            os.chmod(target, mode | 0o200)
            shutil.copyfile(out_root / rel, target)
            os.chmod(target, mode)
        print("installed; originals in", arch)


if __name__ == "__main__":
    main()
