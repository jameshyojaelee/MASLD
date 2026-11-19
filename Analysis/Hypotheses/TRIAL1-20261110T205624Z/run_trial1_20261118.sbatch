#!/bin/bash
#SBATCH --job-name=trial1_261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/TRIAL1-20261110T205624Z/logs/trial1_261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/TRIAL1-20261110T205624Z/logs/trial1_261118_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
HERE=$ROOT/Analysis/Hypotheses/TRIAL1-20261110T205624Z
OUT=$HERE/out_20261118
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=20261118

mkdir -p "$HERE/logs" "$OUT"

# Snapshot exact producer scripts into out_20261118 so historical top-level outputs are not overwritten
cp -p "$HERE/build_tables.py" "$OUT/build_tables.py"
cp -p "$HERE/build_tables_v2.py" "$OUT/build_tables_v2.py"

{
  "$PY" --version
  hostname
  date -u
} > "$OUT/versions.txt"

echo "=== Running build_tables.py in $OUT ==="
"$PY" "$OUT/build_tables.py"

echo "=== Running build_tables_v2.py in $OUT ==="
"$PY" "$OUT/build_tables_v2.py"

echo "=== Verifying TRIAL1 N11 repairs ==="
"$PY" - "$OUT" <<'PYEOF'
import csv, os, sys
from collections import Counter

out = sys.argv[1]
for fname in ("trial_outcomes.tsv", "trial_outcomes_v2.tsv"):
    rows = list(csv.DictReader(open(os.path.join(out, fname)), delimiter="\t"))
    assert len(rows) == 72, f"{fname}: expected 72 rows, got {len(rows)}"
    rep = [r for r in rows if r["representative"] == "yes"]
    assert len(rep) == 66, f"{fname}: expected 66 rep rows, got {len(rep)}"
    mechs = sorted({r["mechanism"] for r in rows})
    assert len(mechs) == 28, f"{fname}: expected 28 mechanisms, got {len(mechs)}"

    # Check FGF19 analogue repair and absence of FGF1 in target_genes
    fgf19_rows = [r for r in rows if r["mechanism"] == "FGF19 analogue"]
    assert len(fgf19_rows) == 1, f"{fname}: expected 1 FGF19 analogue row"
    assert fgf19_rows[0]["target_genes"] == "FGF19;FGFR4;KLB", f"{fname}: bad FGF19 targets: {fgf19_rows[0]['target_genes']}"
    all_target_tokens = {tok.strip() for r in rows for tok in r["target_genes"].split(";")}
    assert "FGF1" not in all_target_tokens, f"{fname}: FGF1 still present in target_genes!"
    assert {"FGF19", "FGFR4", "KLB"}.issubset(all_target_tokens), f"{fname}: missing FGF19/FGFR4/KLB"

    # Check MSDC-0602K mechanism label
    msdc_rows = [r for r in rows if r["drug"] == "MSDC-0602K"]
    assert len(msdc_rows) == 3, f"{fname}: expected 3 MSDC-0602K rows"
    assert all(r["mechanism"] == "insulin sensitizer (MPC modulator)" for r in msdc_rows), (
        f"{fname}: bad MSDC-0602K mechanism: {[r['mechanism'] for r in msdc_rows]}"
    )
    assert all("PPARG-sparing" not in r["mechanism"] and "per candidate list" not in r["mechanism"] for r in rows)

    if fname == "trial_outcomes_v2.tsv":
        assert all(r["target_genes"] == "MPC1;MPC2" for r in msdc_rows)
        nv_exact = sum(v == "not verified" for r in rows for v in r.values())
        nv_sub = sum("not verified" in v for r in rows for v in r.values())
        nv_rep_rows = sum(any("not verified" in v for v in r.values()) for r in rep)
        print(f"{fname}: PASS | 72 rows (66 rep), 28 mechs")
        print("  result_vs_control (all):", dict(Counter(r["result_vs_control"] for r in rows)))
        print("  result_vs_control (rep):", dict(Counter(r["result_vs_control"] for r in rep)))
        print(f"  not_verified exact fields (all): {nv_exact}, substring fields (all): {nv_sub}, rep rows flagged: {nv_rep_rows}")
    else:
        print(f"{fname}: PASS | 72 rows (66 rep), 28 mechs")
PYEOF

cd "$OUT"
sha256sum build_tables.py build_tables_v2.py trial_outcomes.tsv excluded.tsv trial_outcomes_v2.tsv excluded_v2.tsv CHANGES_v2.tsv > output_sha256.txt
cat output_sha256.txt
date -u >> "$OUT/versions.txt"
echo "TRIAL1_20261118_DONE"
