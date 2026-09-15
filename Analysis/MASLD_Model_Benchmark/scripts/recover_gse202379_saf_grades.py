#!/usr/bin/env python
"""
recover_gse202379_saf_grades.py

ADDITIVE recovery of the GSE202379 (Gribben/Vallier 2024 Nature, PMID 38778114)
donor-level SAF grades -- Steatosis, Activity, Fibrosis -- from the cached GEO
series matrix.

WHY THIS EXISTS
  343_harvest_documented_fstage.py captures ONLY the F group of the
  "saf score: S<n>A<n>F<n>" characteristic (its regex is
  r"saf\\s*score[:\\s]*S\\d+A\\d+F(\\d)"). The S (steatosis) and A (activity)
  grades are present in the deposited metadata and are discarded.

WHAT THIS DOES NOT DO
  It does NOT modify 343_harvest_documented_fstage.py and does NOT write to
  donor_fstage_documented.tsv. Both remain byte-identical. This is an overlay.

PARSING NOTE (supersedes an earlier mischaracterisation)
  There is no "field misalignment" and nothing to repair. The series matrix is
  rectangular: all 7 !Sample_characteristics_ch1 rows carry exactly 59 fields.
  Two samples (GSM6112243/P30, GSM6112244/P98) simply LACK the 'liver lobe'
  characteristic, so every field after position 3 shifts up one row and the
  last row is blank for them. Parsing each value by its own "key: value" prefix
  PER SAMPLE -- rather than by row position -- recovers all 40 SAF strings
  directly, with no repair step.

343 F-AGREEMENT AUDIT
  extract_fstage() runs SAF_LABEL_MAP over all characteristic columns BEFORE
  the regex pass and returns on first hit. 'disease status: Healthy control'
  therefore short-circuits a present SAF string. One donor (P98, S1A1F1) has a
  real F=1 overwritten to F=0. P30 hits the identical short-circuit but agrees
  by luck (its SAF is F0 anyway), so an audit that only inspects disagreements
  understates the defect. The same mechanism invents 7 stages outright:
  5 'end stage' -> F=4 and 2 'healthy control' -> F=0 donors carry NO SAF
  string at all.

Outputs (into --outdir):
  donor_saf_grades.tsv     one row per donor, 47 rows
  gsm_saf_records.tsv      one row per GSM, 59 rows, with provenance
  marginals.json           realised marginals at n=40 and n=38
  f_agreement_audit.tsv    per-donor F vs 343's F, with override flags
  label_map_defect.json    the mechanism record
  claim_boundary.md        limits that travel with the output file
"""
from __future__ import annotations
import argparse, gzip, hashlib, importlib.util, json, os, re, sys
from pathlib import Path
import pandas as pd

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
STAGE = ROOT/"Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
SERIES = STAGE/"geo_cache/GSE202379_series_matrix.txt.gz"
H343 = ROOT/"Analysis/SingleCell/scripts/343_harvest_documented_fstage.py"
LINEAGE = STAGE/"lineage_cell_counts_donorcollapsed.tsv"
LIANA_DIR = STAGE/"per_donor_lr_donorcollapsed"

SAF_RE = re.compile(r"^S(\d)A(\d)F(\d)$")
CANONICAL_KEYS = ["tissue","disease status","patient id","liver lobe",
                  "gender","saf score","age"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_343():
    """Import 343's REAL functions. The module has a __main__ guard, so
    importing it executes no work and mutates nothing."""
    spec = importlib.util.spec_from_file_location("h343", H343)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_by_key(series_path: Path) -> pd.DataFrame:
    """Key-based parse: each value keyed by its OWN 'key: value' prefix, per sample.
    Row-position parsing is what produced the earlier phantom 'misalignment'."""
    lines = [l.rstrip("\n") for l in gzip.open(series_path, "rt", errors="ignore")]
    acc = next(l for l in lines if l.startswith("!Sample_geo_accession"))
    gsms = [p.strip('"') for p in acc.split("\t")[1:]]
    rows = [[p.strip('"') for p in l.split("\t")[1:]]
            for l in lines if l.startswith("!Sample_characteristics_ch1")]
    # structural invariant: rectangular matrix
    widths = sorted({len(r) for r in rows})
    recs = []
    for j, g in enumerate(gsms):
        d = {"gsm": g, "_keys_present": []}
        for r in rows:
            v = r[j]
            if ": " in v:
                k, val = v.split(": ", 1)
                d[k.strip()] = val.strip()
                d["_keys_present"].append(k.strip())
        recs.append(d)
    df = pd.DataFrame(recs)
    df.attrs["n_gsm"] = len(gsms)
    df.attrs["n_char_rows"] = len(rows)
    df.attrs["char_row_widths"] = widths
    return df


# Source files whose exact bytes this output file depends on. Recorded in
# source.sha256 so any upstream drift is detectable.
SOURCE_FILES = [
    "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/GSE202379_series_matrix.txt.gz",
    "Analysis/SingleCell/scripts/343_harvest_documented_fstage.py",
    "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv",
    "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts_donorcollapsed.tsv",
    "Analysis/MASLD_Model_Benchmark/scripts/recover_gse202379_saf_grades.py",
    "Analysis/MASLD_Model_Benchmark/tests/unit/test_gse202379_saf_grades.py",
    "Analysis/MASLD_Model_Benchmark/slurm/recover_gse202379_saf_grades_cpu.sbatch",
]


def finalize(exec_root: Path) -> int:
    """Write source.sha256, ARTIFACTS.json and COMPLETE for a built execution.
    The metadata block below is hand-written; only the digests are computed."""
    # source.sha256 -- inputs and code, not outputs
    lines = []
    for rel in SOURCE_FILES:
        pth = ROOT/rel
        lines.append(f"{sha256(pth)}  {rel}")
    (exec_root/"source.sha256").write_text("\n".join(lines) + "\n")

    summary = json.loads((exec_root/"validation"/"recovery_stdout.json").read_text())
    arts = []
    for pth in sorted(exec_root.rglob("*")):
        if not pth.is_file() or pth.name in ("ARTIFACTS.json", "COMPLETE"):
            continue
        arts.append({"path": str(pth.relative_to(exec_root)),
                     "sha256": sha256(pth), "size_bytes": pth.stat().st_size})

    manifest = {
        "artifacts": arts,
        "metadata": {
            "artifact_class": "gse202379_saf_grade_recovery",
            "dataset_id": "gse202379",
            "status": "saf_grades_recovered_additively",
            "supersedes_nothing": True,
            "upstream_untouched": [
                "Analysis/SingleCell/scripts/343_harvest_documented_fstage.py",
                "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv",
            ],
            "n_gsm": 59,
            "n_donors": 47,
            "n_saf_graded": 40,
            "n_usable_for_expression_join": 39,
            "n_usable_for_ccc_join": 37,
            "f_agreement_exception_set": ["P98"],
            "f_agreement_gate": ("SAF-derived F equals 343 F on every SAF-graded donor "
                                 "EXCEPT where the label-map pass overrides a present SAF "
                                 "string; exception set enumerated and pinned"),
            "agrees_by_luck": ["P30"],
            "n_stages_invented_by_343_with_no_saf_string": 7,
            "steatosis_levels_at_n40": 4,
            "single_cohort": True,
            "champion_eligible": False,
            "claim_boundary": ("Substrate recovery only. Axes, not power. A SAF-axis x "
                               "per-donor-CCC analysis is n=37 on a single cohort. No "
                               "disease-direction claim may be derived from this artifact."),
        },
        "schema_version": "masld-bench-artifacts-v1",
    }
    mtxt = json.dumps(manifest, separators=(",", ":"), sort_keys=True)
    (exec_root/"ARTIFACTS.json").write_text(mtxt)
    (exec_root/"COMPLETE").write_text(json.dumps({
        "artifact_count": len(arts),
        "manifest_sha256": hashlib.sha256(mtxt.encode()).hexdigest(),
        "schema_version": "masld-bench-complete-v1",
    }, separators=(",", ":"), sort_keys=True))
    print(f"[finalize] {len(arts)} artifacts, manifest written")
    assert summary["n_saf_graded"] == 40
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir")
    ap.add_argument("--finalize", help="execution root to seal")
    args = ap.parse_args()
    if args.finalize:
        return finalize(Path(args.finalize))
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    h343 = load_343()
    gsm = parse_by_key(SERIES)
    # pandas drops .attrs across merge(); capture the structural facts now
    n_gsm_total = int(gsm.attrs["n_gsm"])
    n_char_rows = int(gsm.attrs["n_char_rows"])
    char_row_widths = list(gsm.attrs["char_row_widths"])

    # ---- structural assertions (replace the phantom 38->40 "repair" check)
    assert n_char_rows == 7, n_char_rows
    assert char_row_widths == [n_gsm_total], \
        f"series matrix is ragged: {char_row_widths}"
    missing_lobe = sorted(gsm.loc[
        ~gsm["_keys_present"].apply(lambda ks: "liver lobe" in ks), "gsm"])
    assert len(missing_lobe) == 2, f"expected 2 samples lacking 'liver lobe', got {missing_lobe}"

    # ---- SAF recovery by key
    m = gsm["saf score"].fillna("").str.extract(SAF_RE)
    gsm["steatosis_S"] = pd.to_numeric(m[0], errors="coerce").astype("Int64")
    gsm["activity_A"]  = pd.to_numeric(m[1], errors="coerce").astype("Int64")
    gsm["fibrosis_F"]  = pd.to_numeric(m[2], errors="coerce").astype("Int64")
    gsm["saf_recovered"] = gsm["fibrosis_F"].notna()
    gsm["lacks_liver_lobe_field"] = gsm["gsm"].isin(missing_lobe)
    assert int(gsm["saf_recovered"].sum()) == 40, int(gsm["saf_recovered"].sum())

    # ---- 343's own extractor on 343's own parse of the same file
    h = h343.parse_geo_characteristics(SERIES)
    h["F_343"] = h.apply(h343.extract_fstage, axis=1)
    gsm = gsm.merge(h[["sample", "F_343"]], left_on="gsm", right_on="sample",
                    how="left").drop(columns=["sample"])

    # did 343's label-map pass fire for this sample?
    def label_hit(ds: str):
        s = str(ds).lower()
        for lab, f in h343.SAF_LABEL_MAP.items():
            if lab in s:
                return lab, f
        return None, None
    lab = gsm["disease status"].apply(label_hit)
    gsm["label_map_label"] = [x[0] for x in lab]
    gsm["label_map_F"]     = pd.array([x[1] for x in lab], dtype="Int64")
    gsm["label_map_fired"] = gsm["label_map_label"].notna()

    # ---- donor collapse (patient id is the biological unit)
    gsm = gsm.sort_values("gsm")
    agg = {"gsm": ("gsm", lambda s: ";".join(sorted(s))),
           "n_gsm": ("gsm", "size")}
    donor = gsm.groupby("patient id").agg(**agg).reset_index()
    first = gsm.drop_duplicates("patient id").set_index("patient id")
    for c in ["disease status","gender","age","saf score","steatosis_S",
              "activity_A","fibrosis_F","saf_recovered","lacks_liver_lobe_field",
              "F_343","label_map_label","label_map_F","label_map_fired"]:
        donor[c] = donor["patient id"].map(first[c])
    donor = donor.rename(columns={"patient id":"donor_id","gsm":"gsm_ids",
                                  "saf score":"saf_score_raw",
                                  "disease status":"disease_status"})
    donor["dataset"] = "GSE202379"
    assert len(donor) == 47, len(donor)
    assert int(donor["saf_recovered"].sum()) == 40, int(donor["saf_recovered"].sum())

    # ---- F agreement audit
    g = donor[donor["saf_recovered"]].copy()
    donor["f_agrees_with_343"] = pd.NA
    ok = donor["saf_recovered"]
    donor.loc[ok, "f_agrees_with_343"] = (
        donor.loc[ok, "fibrosis_F"].astype("Int64") ==
        donor.loc[ok, "F_343"].astype("Int64"))
    disagree = sorted(donor.loc[ok & (donor["f_agrees_with_343"] == False), "donor_id"])
    # agrees-by-luck: label map fired AND a real SAF string exists AND F happens to match
    luck = sorted(donor.loc[ok & donor["label_map_fired"] &
                            (donor["f_agrees_with_343"] == True), "donor_id"])
    # invented: label map fired AND NO SAF string
    invented = sorted(donor.loc[~ok & donor["label_map_fired"], "donor_id"])

    # ---- join usability (two separate columns: they differ by two donors)
    lc = pd.read_csv(LINEAGE, sep="\t")
    expr = {s.replace("GSE202379_","") for s in lc["sample"] if s.startswith("GSE202379")}
    pq = {f.replace("GSE202379_","").replace("_lr_scores.parquet","")
          for f in os.listdir(LIANA_DIR) if f.startswith("GSE202379")}
    donor["usable_for_expression_join"] = donor["donor_id"].isin(expr)
    donor["usable_for_ccc_join"] = donor["donor_id"].isin(pq)
    def why(r):
        if r["usable_for_ccc_join"]:
            return ""
        if not r["usable_for_expression_join"]:
            return "absent_from_expression_substrate"
        # in expression but no LIANA parquet -> depth, not phenotype
        return "failed_liana_gate_min2_lineages_at_30_cells__DEPTH_not_missing_phenotype"
    donor["ccc_exclusion_reason"] = donor.apply(why, axis=1)

    # ---- marginals, both sets
    def marg(sub):
        o = {}
        for c, lab_ in [("steatosis_S","steatosis"),("activity_A","activity"),
                        ("fibrosis_F","fibrosis")]:
            vc = sub[c].dropna().astype(int).value_counts().sort_index()
            o[lab_] = {"n": int(vc.sum()), "levels": {str(k): int(v) for k, v in vc.items()}}
        return o
    graded = donor[donor["saf_recovered"]]
    marginals = {
        "n40_all_saf_graded": marg(graded),
        "n38_excluding_P30_P98": marg(graded[~graded["donor_id"].isin(["P30","P98"])]),
        "note": ("The n38 set is the earlier-reported distribution. At n=40 steatosis "
                 "has FOUR levels (an S0 appears); any downstream code assuming a 1-3 "
                 "steatosis scale is wrong."),
    }

    defect = {
        "defect": "SAF_LABEL_MAP pass precedes the regex pass and returns on first hit",
        "location": "Analysis/SingleCell/scripts/343_harvest_documented_fstage.py:169-176",
        "one_mechanism_three_symptoms": {
            "overwritten_real_grade": {"donors": disagree,
                "detail": "present SAF string discarded in favour of a disease-status label"},
            "agrees_by_luck": {"donors": luck,
                "detail": "identical short-circuit; SAF F happens to equal the label-map F, "
                          "so a disagreement-only audit would miss it"},
            "invented_stage_no_saf_string": {"donors": invented,
                "detail": "fibrosis stage imputed from disease-status text with no measured grade"},
        },
        "invented_breakdown": {
            lab_: sorted(donor.loc[~ok & (donor["label_map_label"] == lab_), "donor_id"])
            for lab_ in h343.SAF_LABEL_MAP},
        "gate": ("SAF-derived F must equal 343's F on every SAF-graded donor EXCEPT where "
                 "the label-map pass overrides a present SAF string. Exception set is "
                 "enumerated and pinned."),
        "exception_set": disagree,
    }

    # ---- write
    donor_cols = ["donor_id","dataset","gsm_ids","n_gsm","disease_status","gender","age",
                  "saf_score_raw","steatosis_S","activity_A","fibrosis_F","saf_recovered",
                  "lacks_liver_lobe_field","F_343","f_agrees_with_343","label_map_label",
                  "label_map_F","label_map_fired","usable_for_expression_join",
                  "usable_for_ccc_join","ccc_exclusion_reason"]
    donor[donor_cols].sort_values("donor_id").to_csv(out/"donor_saf_grades.tsv", sep="\t", index=False)
    gsm.drop(columns=["_keys_present"]).to_csv(out/"gsm_saf_records.tsv", sep="\t", index=False)
    (out/"marginals.json").write_text(json.dumps(marginals, indent=2, sort_keys=True) + "\n")
    (out/"label_map_defect.json").write_text(json.dumps(defect, indent=2, sort_keys=True) + "\n")
    donor.loc[ok, ["donor_id","saf_score_raw","fibrosis_F","F_343","f_agrees_with_343",
                   "label_map_label","label_map_fired","disease_status"]] \
         .sort_values("donor_id").to_csv(out/"f_agreement_audit.tsv", sep="\t", index=False)

    (out/"claim_boundary.md").write_text(f"""# Claim boundary -- GSE202379 SAF recovery

## What this is
Donor-level Steatosis / Activity / Fibrosis grades recovered from the deposited
`saf score` characteristic. Substrate recovery only. **No disease-direction
claim is derived here and none may be derived from this artifact alone.**

## Hard limits
- **40 graded donors, ONE cohort.** These are axes, not power. Never write this
  up as a multi-cohort result.
- **A SAF-axis x per-donor-CCC analysis is n=37 on a single cohort.**
  40 graded -> 39 after the expression join -> 37 after the CCC join.
- **At n=40 steatosis has FOUR levels** (an S0 appears). Code assuming a 1-3
  steatosis scale is wrong.
- `usable_for_expression_join` and `usable_for_ccc_join` are **separate
  columns** and differ by two donors. P62 and P67 drop from CCC because they
  fail LIANA's >=2-lineages-at-30-cells gate -- that is a **depth** property of
  the libraries, **not** missing phenotype, and must not be read as such.
- P70 is SAF-graded but absent from the expression substrate entirely.

## Not a valid comparator
`{STAGE.name}/donor_fstage_documented.tsv` is **not** a valid comparator and was
not reconciled against. It is GSM-keyed (67 rows / 58 graded for this cohort,
against 47 donors / 40 graded), and it carries {len(invented)} fibrosis stages
imputed from disease-status text with no measured grade behind them, with false
provenance. It also drops P70. Do not use it to validate this table.

## Upstream defect this artifact records
See `label_map_defect.json`. One mechanism, three symptoms: {len(disagree)}
overwritten real grade(s), {len(luck)} agreeing by luck, {len(invented)}
invented outright.
""")

    print(json.dumps({
        "n_gsm": n_gsm_total, "n_donors": int(len(donor)),
        "n_saf_graded": int(donor["saf_recovered"].sum()),
        "samples_lacking_liver_lobe": missing_lobe,
        "f_disagreements": disagree, "agrees_by_luck": luck, "invented": invented,
        "usable_for_expression_join": int((donor["saf_recovered"] & donor["usable_for_expression_join"]).sum()),
        "usable_for_ccc_join": int((donor["saf_recovered"] & donor["usable_for_ccc_join"]).sum()),
        "marginals": marginals,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
