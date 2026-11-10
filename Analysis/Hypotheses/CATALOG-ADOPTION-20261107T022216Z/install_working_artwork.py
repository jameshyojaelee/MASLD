#!/usr/bin/env python3
"""Install the explicitly approved Catalog artwork; retain exact prior working files."""
import csv
import hashlib
import json
from pathlib import Path
import os
import shutil

RUN = Path(__file__).resolve().parent
ROOT = RUN.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


rows = list(csv.DictReader((RUN / "figures/INSTALL_PLAN.tsv").open(), delimiter="\t"))
for row in rows:
    src, dst = Path(row["candidate"]), Path(row["canonical"])
    assert src.is_relative_to(RUN / "figures")
    assert dst.is_relative_to(ROOT / "figures/main")
    assert sha(src) == row["sha256"], src
backup = RUN / "working_artwork_before_adoption"
backup.mkdir(exist_ok=False)
records = []
for row in rows:
    src, dst = Path(row["candidate"]), Path(row["canonical"])
    previous = sha(dst) if dst.exists() else None
    if previous:
        saved = backup / dst.relative_to(ROOT)
        saved.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(dst, saved)
        assert sha(saved) == previous
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".catalog-adoption-tmp")
    assert not tmp.exists()
    shutil.copy2(src, tmp)
    os.replace(tmp, dst)
    assert sha(dst) == row["sha256"]
    records.append({**row, "previous_sha256": previous, "installed_sha256": sha(dst)})
(RUN / "working_artwork_install.json").write_text(json.dumps({"status": "PASS", "files": records}, indent=2) + "\n")
print(f"PASS: installed {len(records)} working artwork/source/index files with exact backups")
