#!/usr/bin/env python
"""Freeze the pooled graded-fibrosis substrate for leave-one-cohort-out analysis.

Every column is resolved BY NAME and then VALIDATED BY CONTENT.
Every metadata / proportions file in the Deconvolution tree carries an unnamed
pandas/R rowname at column 0 and is read with index_col=0; the N+1 field count
is asserted before the read.
A zero-row join RAISES.

Unit of inference is the bulk RNA SAMPLE. No donor key exists. Never a donor count.
"""
import hashlib
import json
import os
import re
import sys

import numpy as np
import pandas as pd

SEED = 20260828
N_MC = 5000

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BULK = f"{BASE}/Analysis/Deconvolution/bulk"
DECON = f"{BASE}/Analysis/Deconvolution/results"
KALLISTO = f"{BASE}/RNA-seq/results/kallisto_sensitivity/all_cohorts_gene_counts.tsv.gz"

OUT = sys.argv[1]

# ---------------------------------------------------------------------------
# Cohort declarations. fib_col is resolved by NAME; fib_domain validates CONTENT.
# ---------------------------------------------------------------------------
COHORTS = {
    "GSE135251": dict(
        fib_col="Fibrosis_stage",
        fib_kind="bare_integer",
        control_col="disease",
        control_values={"Control"},
    ),
    "GSE162694": dict(
        fib_col="condition",
        fib_kind="prefixed_string",
        control_col="condition",
        control_values={"Control"},
    ),
    "GSE240729": dict(
        fib_col="fibrosisscore",
        fib_kind="F_string",
        control_col=None,
        control_values=set(),
    ),
    "GSE130970": dict(
        fib_col="fibrosis_stage",
        fib_kind="bare_integer",
        control_col=None,
        control_values=set(),
    ),
    "GSE174478": dict(
        fib_col="condition",
        fib_kind="prefixed_string",
        control_col=None,
        control_values=set(),
    ),
}

# BayesPrism kallisto rerun roster, read off Analysis/Deconvolution/scripts/
# 21_bayesprism_kallisto.sbatch. Only these cohorts have proportions derived
# from the kallisto quantification.
BAYESPRISM_KALLISTO_ROSTER = {"GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"}

LINEAGE_AXIS = [
    "Endothelial cells", "Hepatocytes", "Plasma cells", "T cells", "Cholangiocytes",
    "Fibroblasts", "Macrophages", "Circulating NK/NKT", "Resident NK",
    "Mono+mono derived cells", "Basophils", "B cells", "cDC1s", "cDC2s", "pDCs",
    "Neutrophils",
]
DETECTION_FLOOR = 1e-4
MIN_FRACTION_ABOVE_FLOOR = 0.50


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_ragged(path, expected_header_fields):
    """Assert the N+1 unnamed-rowname offset. Fail closed if it is absent."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
    if len(row) != len(hdr) + 1:
        raise RuntimeError(
            f"{path}: expected the N+1 unnamed-rowname offset "
            f"(header {len(hdr)}, data {len(row)}); got header={len(hdr)} data={len(row)}"
        )
    if expected_header_fields is not None and len(hdr) != expected_header_fields:
        raise RuntimeError(f"{path}: header field count {len(hdr)} != {expected_header_fields}")
    return len(hdr), len(row)


def parse_grade(value, kind):
    v = str(value).strip()
    if kind == "bare_integer":
        m = re.fullmatch(r"([0-4])", v)
        return int(m.group(1)) if m else None
    if kind == "F_string":
        m = re.fullmatch(r"F([0-4])", v)
        return int(m.group(1)) if m else None
    if kind == "prefixed_string":
        m = re.fullmatch(r"(?:NASH|NAFLD|NAFL)_F([0-4])", v)
        return int(m.group(1)) if m else None
    raise ValueError(kind)


def midrank(values):
    s = pd.Series(values)
    return s.rank(method="average").to_numpy(dtype=float)


def pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a - a.mean()
    b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    if d == 0:
        return np.nan
    return float((a * b).sum() / d)


def tie_ceiling(grades):
    """Largest Spearman an UNTIED predictor can reach given this outcome marginal.

    Spearman with ties is Pearson on midranks. The optimum is any predictor whose
    ordering is consistent with the grade ordering; ordering within a tie group is
    irrelevant to the value.
    """
    g = np.asarray(grades, dtype=float)
    order = np.argsort(g, kind="stable")
    pred_rank = np.empty(len(g), dtype=float)
    pred_rank[order] = np.arange(1, len(g) + 1, dtype=float)
    return pearson(pred_rank, midrank(g))


def binary_presence_absence_score(grades, n_draws, rng):
    """Spearman of a PERFECT F0-vs-F1+ call with RANDOM ordering inside each group."""
    g = np.asarray(grades, dtype=float)
    n = len(g)
    yr = midrank(g)
    idx0 = np.flatnonzero(g == 0)
    idx1 = np.flatnonzero(g > 0)
    n0 = len(idx0)
    out = np.empty(n_draws, dtype=float)
    for d in range(n_draws):
        pred = np.empty(n, dtype=float)
        pred[rng.permutation(idx0)] = np.arange(1, n0 + 1, dtype=float)
        pred[rng.permutation(idx1)] = np.arange(n0 + 1, n + 1, dtype=float)
        out[d] = pearson(pred, yr)
    return out


def clr(mat):
    """Centered log-ratio over the columns supplied (the closure is whatever is passed)."""
    a = np.asarray(mat, dtype=float)
    if not np.all(a > 0):
        raise RuntimeError("CLR requires strictly positive parts; exact zeros present")
    lg = np.log(a)
    return lg - lg.mean(axis=1, keepdims=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(SEED)
    report = {
        "seed": SEED,
        "n_monte_carlo_draws": N_MC,
        "unit_of_inference": "bulk_rna_sample",
        "donor_key_exists": False,
        "cohorts": {},
        "excluded_cohorts": {},
        "inputs_sha256": {},
    }

    # -----------------------------------------------------------------
    # 1. Per-cohort labels: verify scale, exclude controls, compute geometry
    # -----------------------------------------------------------------
    frames = []
    for gse, spec in COHORTS.items():
        mpath = f"{BULK}/{gse}/{gse}_metadata.tsv"
        nh, nr = assert_ragged(mpath, None)
        meta = pd.read_csv(mpath, sep="\t", index_col=0, dtype=str)
        if meta.shape[1] != nh:
            raise RuntimeError(f"{gse}: index_col=0 read gave {meta.shape[1]} cols, header had {nh}")
        report["inputs_sha256"][os.path.basename(mpath)] = sha256(mpath)

        col = spec["fib_col"]
        if col not in meta.columns:
            raise RuntimeError(f"{gse}: fibrosis column {col!r} absent from {list(meta.columns)}")
        raw = meta[col].astype(str)

        # CONTENT VALIDATION: the resolved column must actually carry a stage
        # domain, not a run accession or anything else.
        if raw.str.fullmatch(r"[SED]RR\d+").any():
            raise RuntimeError(f"{gse}: column {col!r} carries run accessions, not stages (offset bug)")
        grades = raw.map(lambda v: parse_grade(v, spec["fib_kind"]))
        n_parsed = int(grades.notna().sum())
        if n_parsed == 0:
            raise RuntimeError(f"{gse}: column {col!r} parsed 0 stages under kind {spec['fib_kind']!r}")

        # Control-arm exclusion. A healthy control is a different source
        # population, not the low end of a within-cohort fibrosis gradient.
        if spec["control_col"]:
            cc = meta[spec["control_col"]].astype(str)
            is_control = cc.isin(spec["control_values"])
        else:
            is_control = pd.Series(False, index=meta.index)
            for c in meta.columns:
                if meta[c].astype(str).str.fullmatch("Control", case=False).any():
                    raise RuntimeError(f"{gse}: undeclared 'Control' values in column {c!r}")
        n_controls = int(is_control.sum())

        keep = grades.notna() & (~is_control)
        n_ungraded = int((~grades.notna() & ~is_control).sum())
        g = grades[keep].astype(int)
        if len(g) == 0:
            raise RuntimeError(f"{gse}: zero graded samples after control exclusion")

        marg = {k: int(v) for k, v in g.value_counts().sort_index().items()}
        full_scale = sorted(marg) == [0, 1, 2, 3, 4]
        if not full_scale:
            raise RuntimeError(f"{gse}: scale is NOT full 0-4: {marg}")

        n = len(g)
        ceil = tie_ceiling(g.to_numpy())
        draws = binary_presence_absence_score(g.to_numpy(), N_MC, rng)
        binary_mean = float(draws.mean())
        null_sd = 1.0 / np.sqrt(n - 1)   # permutation Var(Pearson-on-ranks) = 1/(n-1), ties or not
        headroom = ceil - binary_mean

        report["cohorts"][gse] = dict(
            metadata_path=mpath,
            metadata_rows=int(meta.shape[0]),
            fibrosis_column=col,
            fibrosis_column_kind=spec["fib_kind"],
            raw_value_counts={str(k): int(v) for k, v in raw.value_counts().items()},
            controls_excluded=n_controls,
            ungraded_excluded=n_ungraded,
            n_graded_samples=n,
            fibrosis_marginal=marg,
            scale_verified_full_0_4=bool(full_scale),
            tie_ceiling=round(float(ceil), 6),
            binary_only_score_mean=round(binary_mean, 6),
            binary_only_score_sd=round(float(draws.std(ddof=1)), 6),
            binary_only_score_q025=round(float(np.quantile(draws, 0.025)), 6),
            binary_only_score_q975=round(float(np.quantile(draws, 0.975)), 6),
            analytic_null_sd=round(float(null_sd), 6),
            grading_headroom=round(float(headroom), 6),
            headroom_in_null_sd=round(float(headroom / null_sd), 4),
        )
        np.save(os.path.join(OUT, f"binary_mc_draws_{gse}.npy"), draws)

        frames.append(pd.DataFrame({
            "sample_id": g.index,
            "cohort": gse,
            "fibrosis_grade": g.to_numpy(),
        }))

    samples = pd.concat(frames, ignore_index=True)
    if samples["sample_id"].duplicated().any():
        raise RuntimeError("duplicate sample ids across cohorts")

    # -----------------------------------------------------------------
    # 2. Expression on ONE quantification: the kallisto master matrix
    # -----------------------------------------------------------------
    report["inputs_sha256"]["all_cohorts_gene_counts.tsv.gz"] = sha256(KALLISTO)
    kc = pd.read_csv(KALLISTO, sep="\t")
    gene_col = kc.columns[0]
    genes_versioned = kc[gene_col].astype(str).to_numpy()
    kc = kc.drop(columns=[gene_col])
    kc.index = pd.Index([g.split(".")[0] for g in genes_versioned], name="gene_id")

    missing = [s for s in samples["sample_id"] if s not in kc.columns]
    if missing:
        raise RuntimeError(f"{len(missing)} graded samples absent from the kallisto matrix: {missing[:5]}")
    X = kc.loc[:, samples["sample_id"].tolist()]
    if X.shape[1] == 0:
        raise RuntimeError("zero-column join against the kallisto matrix")
    # PAR_Y / duplicate base ids: sum
    n_before = X.shape[0]
    X = X.groupby(level=0).sum()
    counts = X.to_numpy(dtype=np.float64)

    # Shared gene space: EXPRESSION-ONLY, label-free. A gene is retained if in
    # EVERY included cohort at least 20% of that cohort's graded samples carry
    # at least 1 estimated count. This uses no outcome and no held-out label,
    # but it is transductive over expression and is declared as such.
    cohort_vec = samples["cohort"].to_numpy()
    detected = counts >= 1.0
    keep_mask = np.ones(counts.shape[0], dtype=bool)
    per_cohort_frac = {}
    for gse in COHORTS:
        sel = cohort_vec == gse
        frac = detected[:, sel].mean(axis=1)
        per_cohort_frac[gse] = frac
        keep_mask &= frac >= 0.20
    genes_kept = X.index.to_numpy()[keep_mask]
    counts = counts[keep_mask, :]

    lib = counts.sum(axis=0)
    cpm = counts / lib[None, :] * 1e6
    logcpm = np.log2(cpm + 1.0).astype(np.float32).T   # samples x genes

    report["gene_space"] = dict(
        quantification="kallisto tximport gene-level estimated counts, GENCODE v49",
        source=KALLISTO,
        genes_in_source=int(n_before),
        genes_after_version_strip_and_dup_sum=int(X.shape[0]),
        rule="retained if >=20% of the graded samples of EVERY included cohort have count >= 1",
        rule_is_label_free=True,
        rule_is_transductive_over_expression=True,
        genes_retained=int(len(genes_kept)),
        matrix_layout="samples x genes, log2(CPM+1), float32",
    )

    np.save(os.path.join(OUT, "expression_log2cpm.npy"), logcpm)
    with open(os.path.join(OUT, "expression_genes.txt"), "w") as fh:
        fh.write("\n".join(genes_kept.tolist()) + "\n")

    samples = samples.reset_index(drop=True)
    samples["expression_row_index"] = np.arange(len(samples), dtype=int)

    # -----------------------------------------------------------------
    # 3. BayesPrism lineage proportions
    # -----------------------------------------------------------------
    prop_frames = {}
    prop_meta = {}
    for gse in COHORTS:
        ppath = f"{DECON}/{gse}/{gse}_bayesprism_proportions.tsv"
        if not os.path.exists(ppath):
            prop_meta[gse] = dict(available=False, reason="no BayesPrism proportions file")
            continue
        nh, nr = assert_ragged(ppath, len(LINEAGE_AXIS))
        pr = pd.read_csv(ppath, sep="\t", index_col=0)
        if list(pr.columns) != LINEAGE_AXIS:
            raise RuntimeError(f"{gse}: lineage axis mismatch {list(pr.columns)}")
        rs = pr.sum(axis=1)
        if not np.allclose(rs.to_numpy(), 1.0, atol=1e-6):
            raise RuntimeError(f"{gse}: proportion rows do not close to 1 (min {rs.min()}, max {rs.max()})")
        variant = "kallisto" if gse in BAYESPRISM_KALLISTO_ROSTER else "star_era_pre_kallisto_rerun"
        star_backup = os.path.exists(f"{DECON}/{gse}/{gse}_bayesprism_proportions_star_backup.tsv")
        report["inputs_sha256"][os.path.basename(ppath)] = sha256(ppath)
        prop_frames[gse] = pr
        prop_meta[gse] = dict(
            available=True, path=ppath, rows=int(pr.shape[0]),
            deconvolution_input_quantification=variant,
            star_backup_present=bool(star_backup),
        )

    # Eligible lineage family derived on ONE quantification: KALLISTO.
    kall_cohorts = [g for g in COHORTS if prop_meta.get(g, {}).get(
        "deconvolution_input_quantification") == "kallisto"]
    frac_above = {}
    for gse in kall_cohorts:
        pr = prop_frames[gse]
        ids = [s for s in samples.loc[samples.cohort == gse, "sample_id"] if s in pr.index]
        if not ids:
            raise RuntimeError(f"{gse}: zero-row join between graded samples and proportions")
        sub = pr.loc[ids]
        frac_above[gse] = (sub > DETECTION_FLOOR).mean(axis=0)
    ftab = pd.DataFrame(frac_above)
    eligible = sorted(ftab.index[(ftab >= MIN_FRACTION_ABOVE_FLOOR).all(axis=1)].tolist())

    report["lineages"] = dict(
        lineage_axis=LINEAGE_AXIS,
        n_lineages=len(LINEAGE_AXIS),
        detection_floor=DETECTION_FLOOR,
        min_fraction_above_floor=MIN_FRACTION_ABOVE_FLOOR,
        eligibility_derived_on_quantification="kallisto",
        eligibility_derivation_cohorts=kall_cohorts,
        eligible_lineages=eligible,
        n_eligible=len(eligible),
        fraction_above_floor=json.loads(ftab.round(4).to_json()),
        proportions_availability=prop_meta,
    )

    # -----------------------------------------------------------------
    # 4. Per-sample deposit
    # -----------------------------------------------------------------
    for lin in LINEAGE_AXIS:
        samples[f"prop__{lin}"] = np.nan
    for gse, pr in prop_frames.items():
        sel = samples.cohort == gse
        ids = samples.loc[sel, "sample_id"]
        hit = ids.isin(pr.index)
        if hit.sum() == 0:
            raise RuntimeError(f"{gse}: zero-row proportions join")
        vals = pr.reindex(ids[hit]).to_numpy()
        rows = samples.index[sel][hit.to_numpy()]
        for j, lin in enumerate(LINEAGE_AXIS):
            samples.loc[rows, f"prop__{lin}"] = vals[:, j]
        samples.loc[samples.index[sel][~hit.to_numpy()], "proportions_missing"] = True

    samples["proportions_available"] = samples[[f"prop__{l}" for l in LINEAGE_AXIS]].notna().all(axis=1)
    samples["proportions_variant"] = samples["cohort"].map(
        lambda g: prop_meta.get(g, {}).get("deconvolution_input_quantification", "none"))

    # CLR is reported BOTH ways: over the full 16-part closure and over the
    # eligible-lineage subcomposition. Closure can reverse an association.
    have = samples["proportions_available"].to_numpy()
    full16 = samples.loc[have, [f"prop__{l}" for l in LINEAGE_AXIS]].to_numpy(dtype=float)
    clr16 = clr(full16)
    for j, lin in enumerate(LINEAGE_AXIS):
        samples.loc[have, f"clr16__{lin}"] = clr16[:, j]
    subm = samples.loc[have, [f"prop__{l}" for l in eligible]].to_numpy(dtype=float)
    subm = subm / subm.sum(axis=1, keepdims=True)
    clrsub = clr(subm)
    for j, lin in enumerate(eligible):
        samples.loc[have, f"clrsub__{lin}"] = clrsub[:, j]

    # -----------------------------------------------------------------
    # 5. Leave-one-cohort-out folds
    # -----------------------------------------------------------------
    fold_cohorts = sorted(COHORTS)
    samples["loco_fold"] = samples["cohort"].map({c: i for i, c in enumerate(fold_cohorts)})
    folds = []
    for i, c in enumerate(fold_cohorts):
        te = samples.loc[samples.cohort == c, "sample_id"].tolist()
        tr = samples.loc[samples.cohort != c, "sample_id"].tolist()
        folds.append(dict(fold_index=i, held_out_cohort=c,
                          n_test=len(te), n_train=len(tr),
                          test_sample_ids=te, train_sample_ids=tr))
    with open(os.path.join(OUT, "loco_folds.json"), "w") as fh:
        json.dump(dict(seed=SEED, n_folds=len(folds),
                       fold_definition="leave_one_cohort_out",
                       cohort_order=fold_cohorts, folds=folds), fh, indent=2)

    samples.to_csv(os.path.join(OUT, "per_sample_substrate.tsv"), sep="\t", index=False)
    report["total_pooled_samples"] = int(len(samples))
    report["per_cohort_sample_counts"] = {k: int(v) for k, v in samples.cohort.value_counts().items()}

    with open(os.path.join(OUT, "substrate_freeze_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    # digests of every frozen file
    digests = {}
    for fn in sorted(os.listdir(OUT)):
        p = os.path.join(OUT, fn)
        if os.path.isfile(p) and fn != "SHA256SUMS":
            digests[fn] = sha256(p)
    with open(os.path.join(OUT, "SHA256SUMS"), "w") as fh:
        for fn, d in digests.items():
            fh.write(f"{d}  {fn}\n")

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
