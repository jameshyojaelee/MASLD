#!/usr/bin/env python3
"""87_run_manifest.py -- Phase-6 WS-6 reproducibility manifest for the seqfunc layer.

Pure-Python stdlib only (no pandas/numpy) so it runs in any environment. It hashes
the key seqfunc INPUT and OUTPUT files and emits two artifacts into results/seqfunc/:

  (A) run_manifest.json          -- per-file path / SHA-256 / byte size / mtime,
                                    env-version file fingerprints, AlphaGenome query
                                    date (if discoverable), and the RNG seeds used by
                                    src/84 and src/85.
  (B) source_verification_table.md -- documents the exact statistic/unit/threshold/
                                    allele-convention/official-vs-reconstruction status
                                    of the two external QTL truth sources.

Determinism note: wall-clock generation timestamps are intentionally omitted. File
mtimes are read from stat() so re-running on unchanged inputs yields identical hashes.
"""
import hashlib
import json
import os
import re
import sys

# --- Paths -------------------------------------------------------------------
# This script lives in GWAS/finemapping/src/ ; the finemapping root is one up.
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
FINEMAP_ROOT = os.path.dirname(SRC_DIR)                       # GWAS/finemapping
PROJECT_ROOT = os.path.abspath(os.path.join(FINEMAP_ROOT, "..", ".."))
SEQFUNC_DIR = os.path.join(FINEMAP_ROOT, "results", "seqfunc")
GWAS_ATAC_DIR = os.path.join(FINEMAP_ROOT, "results", "gwas_atac")
ENVS_DIR = os.path.join(PROJECT_ROOT, "envs")

OUT_MANIFEST = os.path.join(SEQFUNC_DIR, "run_manifest.json")
OUT_TABLE = os.path.join(SEQFUNC_DIR, "source_verification_table.md")

# Key seqfunc INPUT files (feature matrices, truth sets, raw model scores).
INPUT_FILES = [
    os.path.join(SEQFUNC_DIR, "direction_features_caqtl.tsv"),
    os.path.join(SEQFUNC_DIR, "direction_labels_eqtl.tsv"),
    os.path.join(SEQFUNC_DIR, "broadaway_benchmark_truth.tsv"),
    os.path.join(SEQFUNC_DIR, "borzoi_eqtl_logsed_scores.tsv"),
    os.path.join(SEQFUNC_DIR, "alphagenome_eqtl_scores.tsv"),
    os.path.join(SEQFUNC_DIR, "decima_celltype.tsv"),
    os.path.join(SEQFUNC_DIR, "coding_hardening_v2.tsv"),
    os.path.join(GWAS_ATAC_DIR, "caqtl_variant_credset.csv"),
    os.path.join(SEQFUNC_DIR, "mpra_benchmark", "mpra_substrate_truth.tsv"),
]

# Key seqfunc OUTPUT files (verdicts, per-model tables, summaries).
OUTPUT_FILES = [
    os.path.join(SEQFUNC_DIR, "caqtl_zeroshot_verdict.json"),
    os.path.join(SEQFUNC_DIR, "caqtl_zeroshot_permodel.tsv"),
    os.path.join(SEQFUNC_DIR, "eqtl_signal_diagnose_verdict.json"),
    os.path.join(SEQFUNC_DIR, "eqtl_signal_diagnose_permodel.tsv"),
    os.path.join(SEQFUNC_DIR, "decima_celltype_summary.json"),
    os.path.join(SEQFUNC_DIR, "gof_lof_benchmark.json"),
    os.path.join(SEQFUNC_DIR, "mpra_benchmark", "mpra_benchmark_verdict.json"),
    os.path.join(SEQFUNC_DIR, "hardening", "seqfunc_contract.json"),
    os.path.join(SEQFUNC_DIR, "hardening", "decima_celltype_calibrated.tsv"),
    os.path.join(SEQFUNC_DIR, "activity_by_contact", "contract.json"),
    os.path.join(SEQFUNC_DIR, "activity_by_contact", "measured_variant_gene_links.tsv"),
    os.path.join(SEQFUNC_DIR, "coding_protein_bridge", "summary.json"),
    os.path.join(SEQFUNC_DIR, "coding_protein_bridge", "coding_variant_protein_evidence.tsv"),
    os.path.join(SEQFUNC_DIR, "heavy_preflight", "preflight_contract.json"),
    os.path.join(SEQFUNC_DIR, "adult_liver_chrombpnet", "model_manifest.tsv"),
    os.path.join(SEQFUNC_DIR, "disease_splicing", "cohort_feasibility.tsv"),
    os.path.join(SEQFUNC_DIR, "haplotype_saturation", "contract.json"),
    os.path.join(SEQFUNC_DIR, "functional_prior_sensitivity", "contract.json"),
]

# Environment / model-version fingerprint files (recorded, not deeply parsed).
ENV_FILES = [
    os.path.join(SEQFUNC_DIR, "borzoi_env.yml"),
    os.path.join(ENVS_DIR, "alphagenome.yml"),
    os.path.join(ENVS_DIR, "finemapping.yml"),
    os.path.join(ENVS_DIR, "rnaseq.yml"),
    os.path.join(ENVS_DIR, "seqfunc_chrombpnet.yml"),
    os.path.join(ENVS_DIR, "seqfunc_leafcutter.yml"),
]

# Provenance JSONs to scan for an AlphaGenome API query date.
AG_PROVENANCE_FILES = [
    os.path.join(SEQFUNC_DIR, "AG_GATE_VERDICT.json"),
    os.path.join(SEQFUNC_DIR, "alphagenome_eqtl_scores.tsv"),
]

# Scripts to grep for RNG seeds.
SEED_SCRIPTS = [
    os.path.join(SRC_DIR, "84_caqtl_zeroshot.py"),
    os.path.join(SRC_DIR, "85_eqtl_signal_diagnose.py"),
]


def rel(path):
    """Path relative to project root for stable, portable manifest keys."""
    try:
        return os.path.relpath(path, PROJECT_ROOT)
    except ValueError:
        return path


def file_entry(path):
    """Return a manifest record for one file, or a MISSING stub if absent."""
    if not os.path.exists(path):
        return {"path": rel(path), "status": "MISSING"}
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    st = os.stat(path)
    return {
        "path": rel(path),
        "status": "present",
        "sha256": h.hexdigest(),
        "bytes": st.st_size,
        "mtime": int(st.st_mtime),
    }


def scan_ag_query_date():
    """Search AG provenance files for an explicit query/API date. Return record."""
    date_re = re.compile(
        r"(query[_-]?date|api[_-]?date|queried[_-]?(?:at|on)|"
        r"retriev\w*[_-]?date|access[_-]?date)\D{0,3}"
        r"(\d{4}-\d{2}-\d{2})",
        re.IGNORECASE,
    )
    for path in AG_PROVENANCE_FILES:
        if not os.path.exists(path):
            continue
        try:
            with open(path, "r", errors="replace") as fh:
                # Only scan the head; provenance sits in headers/metadata.
                head = fh.read(200000)
        except OSError:
            continue
        m = date_re.search(head)
        if m:
            return {"date": m.group(2), "source": rel(path), "field": m.group(1)}
    return {"date": "unknown", "source": None,
            "note": "no query/api/access date field found in AG provenance files"}


def grep_seeds():
    """Grep src/84 and src/85 for RNG seeds; record file:line:text for each hit."""
    seed_re = re.compile(
        r"^\s*SEED\s*=\s*\d+|(?:^|[^A-Za-z_])seed\s*=\s*\d+|"
        r"set\.seed\(\s*\d+|default_rng\(\s*\d+",
    )
    results = {}
    for path in SEED_SCRIPTS:
        key = rel(path)
        hits = []
        if not os.path.exists(path):
            results[key] = "MISSING"
            continue
        with open(path, "r", errors="replace") as fh:
            for ln, line in enumerate(fh, 1):
                if seed_re.search(line):
                    hits.append({"line": ln, "text": line.strip(),
                                 "ref": "%s:%d" % (key, ln)})
        results[key] = hits if hits else "no explicit numeric seed found"
    return results


def build_manifest():
    ag = scan_ag_query_date()
    manifest = {
        "generated_by": "87_run_manifest.py",
        "note": ("Wall-clock generation timestamps are intentionally omitted for "
                 "determinism; file mtimes are read from stat() instead. Re-running "
                 "on unchanged inputs yields byte-identical hashes."),
        "phase": "Phase-6 WS-6 (seqfunc reproducibility manifest)",
        "project_root": PROJECT_ROOT,
        "hash_algorithm": "sha256",
        "inputs": [file_entry(p) for p in INPUT_FILES],
        "outputs": [file_entry(p) for p in OUTPUT_FILES],
        "env_version_files": [file_entry(p) for p in ENV_FILES],
        "alphagenome_query_date": ag,
        "seeds": grep_seeds(),
        "seed_notes": ("src/84 defines module-level SEED=42 (line 71) and per-routine "
                       "defaults seed=123/321/99; src/85 has no module-level SEED "
                       "constant and relies on per-routine defaults seed=123/99."),
    }
    return manifest


SOURCE_TABLE_MD = """\
# Seqfunc external QTL truth sources -- verification table

**Header note.** AlphaGenome outputs are used **ZERO-SHOT only** (never used to
train, fit, or calibrate any model) and only for non-commercial research, per the
AlphaGenome FAQ / terms of use. Every statistic below is the *exact* quantity the
seqfunc layer consumes, with its unit, threshold, allele convention, and whether it
is an **official source call** or an **internal reconstruction** built on top of the
source's raw summary statistics.

| Source | Assay / cohort | Statistic used | Unit | Threshold | Allele convention | Official call vs internal reconstruction |
|---|---|---|---|---|---|---|
| Broadaway et al. 2024 (AJHG, PMC11393674) | Independent **liver microarray eQTL** meta-analysis, n=1,183 | Marginal **min-p eQTL** per (gene, locus), flagged `is_signal_lead` | Signal-level (169 signals after dedup) | Per-(gene,locus) marginal minimum p-value (lead SNP); NOT a fine-mapped PIP/credible-set cutoff | ref->alt orientation propagated to features; strand-ambiguous A/T and C/G pairs dropped | **INTERNAL reconstruction** -- a marginal lead pick, NOT the source's own fine-mapped causal-variant call |
| Currin et al. 2025 (Genome Research, PMC12212085) | **Bulk liver caQTL**, n=138 donors (bulk liver tissue, NOT primary hepatocyte) | One variant per **peak-lead** (nominal lead); signal-level distance-clump used as a sensitivity arm | Nominal peak-lead variants (11,896); distance-clump signal-level sensitivity (1,470) | Nominal lead per caPeak; **source's official unit is FDR<5% caPeaks (35,361) / 2,126 LD-clumped signals** | ref->alt orientation propagated to accessibility features | **INTERNAL reconstruction** -- our nominal per-peak-lead unit is NOT the source's FDR<5% caPeak (35,361) / 2,126-signal official unit |
| Hu et al. 2025 preprint (PMC12633503; GSE281364) | Episomal MPRA in HepG2 and LX-2, control and PAOA/TGF-beta | Official significant DAV tables S2-S5; barcode-mapped oligo intervals define tested negatives | Variant x cellular-context | DAV FDR < 0.01 in the source; benchmark negatives are barcode-mapped intervals without a DAV call in that context | interval center matched to hg38 substrate; effect direction is not benchmarked until the source logFC allele convention is independently verified | **OFFICIAL DAV calls + source GEO barcode-map experimental universe**; MPRA is cell-line/episomal and not endogenous adult-liver validation |

**Reconstruction caveat.** For both sources the seqfunc layer's benchmark unit
(marginal lead / nominal peak-lead) is an internal convenience reconstruction from
the published summary statistics, not the source's primary fine-mapped or
FDR-thresholded call. Concordance metrics computed against these units are
signal-lead-level, not causal-variant-level.
"""


def main():
    os.makedirs(SEQFUNC_DIR, exist_ok=True)

    manifest = build_manifest()
    with open(OUT_MANIFEST, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
        fh.write("\n")

    with open(OUT_TABLE, "w") as fh:
        fh.write(SOURCE_TABLE_MD)

    n_in = sum(1 for e in manifest["inputs"] if e.get("status") == "present")
    n_out = sum(1 for e in manifest["outputs"] if e.get("status") == "present")
    sys.stderr.write(
        "[87_run_manifest] wrote %s (%d/%d inputs, %d/%d outputs present)\n"
        % (rel(OUT_MANIFEST), n_in, len(INPUT_FILES), n_out, len(OUTPUT_FILES))
    )
    sys.stderr.write("[87_run_manifest] wrote %s\n" % rel(OUT_TABLE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
