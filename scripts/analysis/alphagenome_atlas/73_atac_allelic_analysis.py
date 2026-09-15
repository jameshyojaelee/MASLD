#!/usr/bin/env python3
"""Step 73 (P3): same-variant allelic concordance between GSE281367 snATAC and Atlas liver ATAC predictions.

This is the only same-variant test of Atlas liver predictions in this paper: the variant is inside the peak,
the comparison is within donor (allelic imbalance cancels trans and donor effects), and the assay is adult
human MASLD liver. Genotypes are called from the reads; no donor genotypes exist anywhere in the repository.

Prespecification (filters, seeds, families, written predictions): 73_allelic_prespec.json, copied into the
result root with its sha256 before any number here is read.

Outputs (tables/): allelic_sites.tsv, allelic_concordance.json, allelic_donor_counts.tsv.gz
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import pathlib
import shutil
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
# cohort selected by argv[1]: the counts directory and the output prefix. The target list, the filters and
# the prediction source are identical across cohorts, so GSE244832 is a straight replication of GSE281367.
COHORTS = {"gse281367": ("atac_counts", ""), "gse244832": ("atac_counts_gse244832", "gse244832_")}
TARGETS = TABLES / "atac_targets_in_peaks.tsv"
PRESPEC = la.SCRIPT_DIR / "73_allelic_prespec.json"
# AGA_LIVER_SUMMARIES switches the prediction source. The first pass scored against the preview
# (direct-only, 13,690 variants) because the full run was still building; the full archive is a stated
# re-run, and the path used lands in the summary.
LIVER_SUMMARIES = pathlib.Path(os.environ["AGA_LIVER_SUMMARIES"]).resolve() if os.environ.get("AGA_LIVER_SUMMARIES", "").strip() \
    else la.PROJECT / "GWAS/finemapping/results/alphagenome_atlas/preview-direct-only-20260909T191144Z/tables/liver_summaries.tsv.gz"
STAR = "<*>"
BLOCK_BP = 1_000_000
PERM_SEED = 20260909
BOOT_SEED = 123
SHUFFLE_SEED = 20260910


# ------------------------------------------------------------------ tested helpers

def parse_ad(ref: str, alt_field: str, ad_field: str, target_alt: str) -> tuple[int, int, int]:
    """Split a bcftools AD field into (ref, target-alt, other) reads.

    AD is ordered over [REF] + ALT, and bcftools mpileup always appends the <*> catch-all to ALT. A pileup
    that saw only the reference lists ALT as "<*>" alone, so the target allele is simply absent (0 reads) —
    it must not be read off the catch-all. A third allele at the site belongs in `other`, never in `alt`.
    """
    alts = alt_field.split(",")
    counts = [int(x) for x in ad_field.split(",")]
    alleles = [ref] + alts
    if len(counts) != len(alleles):
        raise la.ContractError(f"AD length {len(counts)} != allele count {len(alleles)} ({ref} {alt_field} {ad_field})")
    n_ref = counts[0]
    n_alt = 0
    n_other = 0
    for allele, c in zip(alleles[1:], counts[1:]):
        if allele == target_alt:
            n_alt += c
        else:
            n_other += c
    return n_ref, n_alt, n_other


def is_indel_record(ref: str, alt_field: str) -> bool:
    """bcftools emits a separate multi-base record wherever it also sees an indel. That record is not our SNV."""
    if len(ref) != 1:
        return True
    return any(len(a) != 1 for a in alt_field.split(",") if a != STAR)


def call_het(n_ref: int, n_alt: int, n_other: int, min_total: int = 10, min_allele: int = 3,
             frac_window: tuple[float, float] = (0.15, 0.85), max_other_frac: float = 0.10) -> bool:
    """Heterozygous call from reads alone (prespec filters). A site with many third-allele reads is refused."""
    total = n_ref + n_alt
    if total < min_total or min(n_ref, n_alt) < min_allele:
        return False
    if n_other > max_other_frac * (total + n_other):
        return False
    frac = n_alt / total
    return frac_window[0] <= frac <= frac_window[1]


def site_statistic(donor_counts: list[tuple[int, int]]) -> dict:
    """Site-level allelic imbalance with the DONOR as the unit: mean over donors of log2((alt+.5)/(ref+.5)).

    Pooling reads instead would let one deep donor outvote the rest, which is pseudoreplication.
    """
    ratios = [math.log2((a + 0.5) / (r + 0.5)) for r, a in donor_counts]
    arr = np.asarray(ratios, dtype=float)
    n = len(arr)
    sd = float(arr.std(ddof=1)) if n > 1 else float("nan")
    return {"mean_log2": float(arr.mean()), "sd_log2": sd, "n_donors": n,
            "se_log2": (sd / math.sqrt(n)) if n > 1 else float("nan"),
            "donor_log2": ratios}


def stratified_permutation_concordance(pred: np.ndarray, meas: np.ndarray, strata: np.ndarray,
                                       draws: int, seed: int) -> tuple[float, np.ndarray]:
    """Null that reassigns WHICH site carries which prediction, within |quantile| strata.

    The statistic is regenerated end to end on every draw; the predictor's magnitude structure is preserved,
    so only the site-to-prediction pairing is destroyed.
    """
    rng = np.random.default_rng(seed)
    observed = float(np.mean(np.sign(pred) == np.sign(meas)))
    out = np.empty(draws, dtype=float)
    idx_by_stratum = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    for d in range(draws):
        shuffled = pred.copy()
        for idx in idx_by_stratum:
            shuffled[idx] = pred[rng.permutation(idx)]
        out[d] = float(np.mean(np.sign(shuffled) == np.sign(meas)))
    p = float((np.sum(out >= observed) + 1) / (draws + 1))
    return p, out


def marginal_expected_concordance(pred: np.ndarray, meas: np.ndarray) -> float:
    """Concordance expected from the two marginal sign skews alone (shared with the benchmark arms)."""
    return la.marginal_expected_concordance(pred, meas)


def block_level_concordance(pred: np.ndarray, meas: np.ndarray, blocks: np.ndarray) -> dict:
    """Concordance with the 1-Mb BLOCK as the independent unit.

    Sites inside one block are in LD: their genotypes and their predictions are correlated, so 400 sites at
    36 loci are not 400 independent observations. Reports the per-block mean, an exact binomial sign test
    over blocks, and the site-level number beside it so the difference is visible.
    """
    uniq = np.unique(blocks)
    per = np.array([float(np.mean(np.sign(pred[blocks == b]) == np.sign(meas[blocks == b]))) for b in uniq])
    above = int(np.sum(per > 0.5))
    decided = int(np.sum(per != 0.5))
    p = float(sum(math.comb(decided, k) for k in range(above, decided + 1)) / (2 ** decided)) if decided else float("nan")
    return {"n_blocks": int(len(uniq)), "mean_per_block": float(per.mean()),
            "n_blocks_above_half": above, "n_blocks_decided": decided, "sign_test_p": p,
            "site_level": float(np.mean(np.sign(pred) == np.sign(meas))),
            "per_block": {str(b): float(v) for b, v in zip(uniq, per)}}


def block_bootstrap_concordance(pred: np.ndarray, meas: np.ndarray, blocks: np.ndarray,
                                draws: int, seed: int) -> dict:
    """1-Mb block bootstrap, one site drawn per sampled block, so nearby sites cannot inflate the interval."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(blocks)
    members = {b: np.flatnonzero(blocks == b) for b in uniq}
    vals = np.empty(draws, dtype=float)
    for d in range(draws):
        picked = rng.choice(uniq, size=len(uniq), replace=True)
        take = np.array([rng.choice(members[b]) for b in picked])
        vals[d] = float(np.mean(np.sign(pred[take]) == np.sign(meas[take])))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"mean": float(vals.mean()), "ci95": [float(lo), float(hi)], "n_blocks": int(len(uniq))}


# ------------------------------------------------------------------ pipeline

def resolve_cohort(argv: list) -> tuple[str, str, str]:
    """Pick the cohort from the command line. Importing the module must never depend on argv (the tests do)."""
    name = (argv[1] if len(argv) > 1 else "gse281367").lower()
    if name not in COHORTS:
        raise la.ContractError(f"unknown cohort {name!r}; expected one of {sorted(COHORTS)}")
    counts_dir, prefix = COHORTS[name]
    return name, counts_dir, prefix


def load_targets() -> dict:
    tgt = {}
    with open(TARGETS) as h:
        for line in h:
            c, p, ref, alt = line.rstrip("\n").split("\t")
            tgt[(c, int(p))] = {"ref": ref, "alt": alt, "uid": f"{c}:{p}:{ref}:{alt}"}
    return tgt


def load_counts(targets: dict, counts: pathlib.Path) -> pd.DataFrame:
    rows = []
    for f in sorted(counts.glob("*.tsv.gz")):
        with gzip.open(f, "rt") as h:
            for r in csv.DictReader(h, delimiter="\t"):
                key = (r["chrom"], int(r["pos"]))
                t = targets.get(key)
                if t is None or is_indel_record(r["ref"], r["alt"]):
                    continue
                n_ref, n_alt, n_other = parse_ad(r["ref"], r["alt"], r["ad"], t["alt"])
                if r["ref"] != t["ref"]:
                    raise la.ContractError(f"pileup REF {r['ref']} != target REF {t['ref']} at {key}")
                rows.append({"donor": r["donor"], "chrom": key[0], "pos": key[1], "uid": t["uid"],
                             "n_ref": n_ref, "n_alt": n_alt, "n_other": n_other,
                             "is_het": call_het(n_ref, n_alt, n_other)})
    return pd.DataFrame(rows)


def load_predictions(uids: set) -> dict:
    """Atlas liver ATAC signed quantile per uid; DNase kept beside it as a second accessibility channel."""
    keep = {"ATAC", "DNASE"}
    out = defaultdict(dict)
    with gzip.open(LIVER_SUMMARIES, "rt") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r["scorer"] not in keep or r["variant_uid"] not in uids:
                continue
            if r["liver_median_quantile"] in ("", "nan"):
                continue
            out[r["variant_uid"]][r["scorer"]] = float(r["liver_median_quantile"])
    return out


def main() -> None:
    COHORT, counts_dir, PREFIX = resolve_cohort(sys.argv)
    COUNTS = ROOT / counts_dir
    TABLES.mkdir(parents=True, exist_ok=True)
    spec = json.loads(PRESPEC.read_text())
    digest = hashlib.sha256(PRESPEC.read_bytes()).hexdigest()
    shutil.copy(PRESPEC, ROOT / "73_allelic_prespec.json")
    (ROOT / "73_allelic_prespec.sha256").write_text(digest + "\n")

    targets = load_targets()
    counts = load_counts(targets, COUNTS)
    if counts.empty:
        raise la.ContractError("no allele counts matched the target list")
    counts.to_csv(TABLES / f"{PREFIX}allelic_donor_counts.tsv.gz", sep="\t", index=False)
    n_donors_total = counts["donor"].nunique()

    het = counts[counts["is_het"]]
    min_het = spec["filters"]["min_het_donors"]
    sites = []
    for uid, g in het.groupby("uid"):
        if len(g) < min_het:
            continue
        st = site_statistic(list(zip(g["n_ref"].tolist(), g["n_alt"].tolist())))
        c, p, ref, alt = uid.split(":")
        sites.append({"uid": uid, "chrom": c, "pos": int(p), "ref": ref, "alt": alt,
                      "n_het_donors": st["n_donors"], "mean_log2_alt_over_ref": st["mean_log2"],
                      "sd_log2": st["sd_log2"], "se_log2": st["se_log2"],
                      "total_reads": int(g["n_ref"].sum() + g["n_alt"].sum())})
    sdf = pd.DataFrame(sites)
    summary = {"cohort": COHORT, "prespec_sha256": digest, "liver_summaries_source": str(LIVER_SUMMARIES), "n_targets": len(targets), "n_donors": int(n_donors_total),
               "n_sites_with_any_counts": int(counts["uid"].nunique()),
               "n_sites_reliable": int(len(sdf)),
               "deviation_from_plan": spec["deviation_from_plan"]}

    if len(sdf) == 0:
        summary["verdict"] = "cannot conclude: no site passed the reliability filter"
        json.dump(summary, (TABLES / f"{PREFIX}allelic_concordance.json").open("w"), indent=1, default=float)
        la.log("P3 allelic: no reliable site")
        return

    preds = load_predictions(set(sdf["uid"]))
    sdf["atac_quantile"] = [preds.get(u, {}).get("ATAC", float("nan")) for u in sdf["uid"]]
    sdf["dnase_quantile"] = [preds.get(u, {}).get("DNASE", float("nan")) for u in sdf["uid"]]
    sdf["block"] = sdf["chrom"] + "~" + (sdf["pos"] // BLOCK_BP).astype(str)

    # reference bias, measured before any concordance number
    bias = float(sdf["mean_log2_alt_over_ref"].mean())
    summary["reference_bias_mean_log2"] = bias
    sdf["mean_log2_centred"] = sdf["mean_log2_alt_over_ref"] - bias
    sdf.to_csv(TABLES / f"{PREFIX}allelic_sites.tsv", sep="\t", index=False)

    for channel in ("atac_quantile", "dnase_quantile"):
        ok = sdf[sdf[channel].notna() & (sdf[channel] != 0)].copy()
        block = f"{channel}_arm"
        if len(ok) < 20:
            summary[block] = {"n": int(len(ok)), "verdict": "cannot conclude: fewer than 20 scored sites"}
            continue
        res = {"n_sites": int(len(ok))}
        for label, mcol in (("raw", "mean_log2_alt_over_ref"), ("bias_centred", "mean_log2_centred")):
            pred = ok[channel].to_numpy(float)
            meas = ok[mcol].to_numpy(float)
            keep = np.sign(meas) != 0
            pred, meas = pred[keep], meas[keep]
            blocks = ok["block"].to_numpy()[keep]
            conc = float(np.mean(np.sign(pred) == np.sign(meas)))
            deciles = pd.qcut(np.abs(pred), q=min(10, max(2, len(pred) // 10)), labels=False, duplicates="drop")
            perm_p, null = stratified_permutation_concordance(pred, meas, np.asarray(deciles), draws=1000, seed=PERM_SEED)
            boot = block_bootstrap_concordance(pred, meas, blocks, draws=2000, seed=BOOT_SEED)
            rng = np.random.default_rng(SHUFFLE_SEED)
            shuffled_sign = np.abs(pred) * rng.choice([-1.0, 1.0], size=len(pred))
            strong = np.abs(meas) > 1.0
            res[label] = {
                "n": int(len(pred)),
                "concordance": conc,
                "block_bootstrap": boot,
                "permutation_p": perm_p,
                "permutation_null_mean": float(null.mean()),
                "sign_shuffle_control": float(np.mean(np.sign(shuffled_sign) == np.sign(meas))),
                "concordance_strong_imbalance": float(np.mean(np.sign(pred[strong]) == np.sign(meas[strong]))) if strong.sum() else None,
                "n_strong": int(strong.sum()),
                "concordance_weak_imbalance": float(np.mean(np.sign(pred[~strong]) == np.sign(meas[~strong]))) if (~strong).sum() else None,
                "n_weak": int((~strong).sum()),
                "spearman": float(pd.Series(pred).corr(pd.Series(meas), method="spearman")),
                "marginal_expected_concordance": marginal_expected_concordance(pred, meas),
                "predicted_positive_fraction": float(np.mean(pred > 0)),
                "measured_positive_fraction": float(np.mean(meas > 0)),
                "block_level": block_level_concordance(pred, meas, blocks),
            }
        summary[block] = res

    diagnostics = []
    if len(sdf) < 50:
        diagnostics.append("fewer than 50 sites pass the reliability filter")
    if abs(bias) > 0.25:
        diagnostics.append(f"reference bias |{bias:.3f}| exceeds 0.25")
    for channel in ("atac_quantile", "dnase_quantile"):
        arm = summary.get(f"{channel}_arm", {})
        raw = arm.get("raw") if isinstance(arm, dict) else None
        if raw and abs(raw["permutation_null_mean"] - 0.5) > 0.02:
            diagnostics.append(f"{channel} permutation null centres at {raw['permutation_null_mean']:.3f}, not 0.5")
        if raw and raw["concordance"] - raw["marginal_expected_concordance"] < 0.02:
            diagnostics.append(f"{channel} concordance {raw['concordance']:.3f} is not above the marginal-skew "
                               f"expectation {raw['marginal_expected_concordance']:.3f}")
    summary["note_on_permutation"] = ("the site-level permutation ignores LD between sites in the same locus and is "
                                     "anti-conservative; the block-level sign test and the 1-Mb block bootstrap are the "
                                     "LD-respecting statements")
    summary["diagnostics_triggered"] = diagnostics
    summary["verdict"] = "cannot conclude" if diagnostics else "concordance reported"
    json.dump(summary, (TABLES / f"{PREFIX}allelic_concordance.json").open("w"), indent=1, default=float)
    la.log(f"P3 allelic [{COHORT}]: {len(sdf)} reliable sites of {len(targets)} targets, {n_donors_total} donors, bias {bias:+.3f}")


if __name__ == "__main__":
    main()
