#!/usr/bin/env python3
"""F5 program profiles and F6 genetic overlap, the two families the original P4 could not run.

F5 was deferred because its size-matched random gene-set null has to draw from the
eligible background gene universe, and 18,570 of that universe's promoter regions
were never in the original region universe, so every draw would have silently
dropped them. Those regions are now retrieved and reduced, so the null can run on
the population the original design named instead of on the scored subset.

F6 is the genetic-overlap family. Step 41 left it empty because Track 0's
variant-peak overlap table did not exist when the universe was defined, eight
hours earlier. The original universe is not reopened; F6 runs over its own
deposit under its own prespecification (data_p4_u2c_prespec.json), whose four
predictions were written before any of its scores was requested.

Feature definitions, track grouping, matching strata, minimum group size and the
matched-draw null are imported from the original producers, unchanged, so an F6
number and an F4 number remain the same quantity. `--self-test` exercises the
statistics on synthetic fixtures with known answers and writes no family result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control

import data_p4 as p4
import data_p4_families as fam

OVERLAPS = p4.ROOT / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/variant_peak_overlaps.tsv.gz"
U2C_PRESPEC = Path(__file__).with_name("data_p4_u2c_prespec.json")
TRACK_GROUPS = ("liver_atac", "liver_atac_adult_verified",
                "liver_dnase_adult_verified", "liver_h3k27ac_adult_verified")
TF_GROUPS = ("chip_tf_liver", "chip_tf_liver_adult_verified")
F5_DRAWS = 1000          # the precedent the F4 TF-composition draws already set
F5_SEED = 20260914       # the original P4 seed; F5 was always part of that prespecification
F6_SEED = 20260917       # the seed registered in data_p4_u2c_prespec.json
FASTA = "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"


def read_features(runs: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Concatenate reduced features from several deposits, refusing a repeated coordinate."""
    feats, tfs = [], []
    for run in runs:
        feats.append(pd.read_csv(run / "saturation_region_features.tsv.gz", sep="\t",
                                 keep_default_na=False, low_memory=False))
        tfs.append(pd.read_csv(run / "saturation_region_tf_top5.tsv.gz", sep="\t",
                               keep_default_na=False))
    feats = pd.concat(feats, ignore_index=True)
    tfs = pd.concat(tfs, ignore_index=True)
    if feats.duplicated(["universe", "region_key"]).any():
        raise ValueError("Repeated (universe, region_key) across deposits")
    if (feats.archive_state != "ok").any():
        raise ValueError("A reduced region did not read cleanly; family tests stay closed")
    return feats, tfs


def promoter_links(da: pd.DataFrame) -> dict:
    links = {}
    for row in da.itertuples():
        for gene in str(row.promoter_genes).replace(";", ",").split(","):
            if gene and gene not in ("NA", "nan", "none", "None"):
                links.setdefault(gene, set()).add(row.peak_coordinate)
    return links


def eligible_da() -> pd.DataFrame:
    da = pd.read_csv(fam.CTX / "da/da_peak_results.tsv.gz", sep="\t", keep_default_na=False)
    manifest = pd.read_csv(fam.CTX / "consensus_peak_manifest.tsv", sep="\t", dtype=str)
    manifest = manifest.loc[~manifest.blacklist_overlap.str.upper().eq("TRUE")]
    keys = set(manifest.lineage + "|" + manifest.chrom + ":" + manifest.start0 + "-" + manifest.end)
    return da.loc[(da.lineage + "|" + da.peak_coordinate).isin(keys)].copy()


def peak_covariates(frame: pd.DataFrame, da: pd.DataFrame, lineage: str, cohort: str) -> pd.DataFrame:
    """GC from the reference, promoter from the source annotation, signal from the cohort counts."""
    import pysam
    genome = pysam.FastaFile(FASTA)
    try:
        gc = [(s.count("G") + s.count("C")) / max(len(s), 1) for s in
              (genome.fetch(r.chrom, int(r.start0), int(r.end)).upper() for r in frame.itertuples())]
    finally:
        genome.close()
    out = frame.copy()
    out["gc"] = gc
    # The reduced feature tables carry chrom/start0/end but not width; the original
    # universe defined it as end - start0 (41_saturation_universe.py:75), and the
    # matching strata need the same quantity for an F6 number to equal an F4 number.
    out["width"] = out["end"].astype(int) - out["start0"].astype(int)
    annotation = da.loc[da.lineage.eq(lineage), ["peak_coordinate", "promoter_genes"]]
    out = out.merge(annotation, left_on="region_key", right_on="peak_coordinate", how="left")
    out["promoter"] = out.promoter_genes.fillna("").apply(
        lambda v: int(v not in ("", "NA", "none", "None", "nan")))
    signals = fam.u2_signal_covariates(cohort, lineage, set(out.region_key))
    out = out.merge(signals, on="region_key", validate="many_to_one")
    return out.reset_index(drop=True)


def f6_family(feats: pd.DataFrame, tfs: pd.DataFrame, out: Path) -> dict:
    """Overlap peaks against matched non-overlap background peaks, per the F6 prespecification."""
    prespec = json.loads(U2C_PRESPEC.read_text())
    if prespec["family"] != "F6_genetic_overlap" or prespec["seed"] != F6_SEED:
        raise ValueError("F6 prespecification identity or seed changed")
    overlaps = pd.read_csv(OVERLAPS, sep="\t", usecols=["lineage", "peak", "queried"], dtype=str)
    overlap_peaks = set(overlaps.peak)
    queried_peaks = set(overlaps.loc[overlaps.queried.str.lower().eq("true"), "peak"])
    da = eligible_da()

    pool = feats.loc[feats.universe.isin(("U2c_genetic_overlap", "U3_f5_background", "U2_snatac"))].copy()
    pool = pool.drop_duplicates("region_key")
    results, compositions, covariates = [], [], []
    for lineage in ("hepatocyte", "stellate"):
        in_lineage = set(da.loc[da.lineage.eq(lineage), "peak_coordinate"])
        frame = pool.loc[pool.region_key.isin(in_lineage)].copy()
        if len(frame) < 2 * fam.MIN_GROUP:
            results.append({"lineage": lineage, "status": "indeterminate_minimum_group",
                            "n_group": int(frame.region_key.isin(overlap_peaks).sum()),
                            "p_nominal": None})
            continue
        for cohort in ("GSE244832", "GSE281367"):
            d = peak_covariates(frame, da, lineage, cohort)
            covariates.append(d[["region_key", "universe", "gc", "width", "promoter",
                                 "signal_mean", "signal_sd", "n_donors"]]
                              .assign(cohort=cohort, lineage=lineage))
            group = d.region_key.isin(overlap_peaks).to_numpy()
            for trackgroup in TRACK_GROUPS:
                for scale in ("rel", "q"):
                    column = f"{trackgroup}_{scale}_mean"
                    if column not in d:
                        continue
                    result = fam.measured_contrast(d, column, group)
                    results.append({"cohort": cohort, "lineage": lineage, "contrast": "overlap_vs_background",
                                    "track_group": trackgroup, "scale": scale,
                                    "prediction": "P6a_primary" if scale == "rel" else "P6a_saturating_scale",
                                    **result})
            # P6c: within the overlap group, does being a queried variant's peak matter?
            sub = d.loc[group].reset_index(drop=True)
            queried = sub.region_key.isin(queried_peaks).to_numpy()
            for trackgroup in TRACK_GROUPS:
                column = f"{trackgroup}_rel_mean"
                if column not in sub:
                    continue
                result = fam.measured_contrast(sub, column, queried)
                results.append({"cohort": cohort, "lineage": lineage, "contrast": "queried_vs_unqueried_within_overlap",
                               "track_group": trackgroup, "scale": "rel", "prediction": "P6c", **result})
            # P6d: TF composition of overlap peaks against matched controls.
            strata = fam.old.stratum_labels(d)
            if min(fam.matching_support(strata, group)) < fam.MIN_GROUP:
                compositions.append({"cohort": cohort, "lineage": lineage, "prediction": "P6d",
                                     "status": "indeterminate_minimum_after_matching", "p_nominal": None})
                continue
            index = {k: i for i, k in enumerate(d.region_key)}
            relevant = tfs.loc[tfs.region_key.isin(index) & tfs.group.isin(TF_GROUPS)]
            for (trackgroup, track), rows in relevant.groupby(["group", "track_name"]):
                present = np.zeros(len(d))
                present[[index[k] for k in rows.region_key]] = 1
                if present.sum() < 50:
                    continue
                cmh = fam.old.cmh_test(present, group, strata)
                draw = fam.old.matched_draw_null(present, strata, group, F5_DRAWS, F6_SEED, statistic="mean")
                compositions.append({"cohort": cohort, "lineage": lineage, "prediction": "P6d",
                                     "group": trackgroup, "track": track,
                                     "transcription_factor": rows.transcription_factor.iloc[0],
                                     "status": "computed_development_comparison",
                                     "share_overlap": draw["observed_group"],
                                     "share_matched_null": draw["null_median"],
                                     "p_nominal": cmh["p"]})
    fam.corrected(results, "F6_all_lineage_cohort_track_contrasts")
    fam.corrected(compositions, "F6_all_lineage_cohort_TF_contrasts")
    pd.DataFrame(results).to_csv(out / "F6_genetic_overlap.tsv", sep="\t", index=False)
    pd.DataFrame(compositions).to_csv(out / "F6_TF_composition.tsv", sep="\t", index=False)
    if covariates:
        pd.concat(covariates, ignore_index=True).to_csv(
            out / "F6_peak_covariates.tsv.gz", sep="\t", index=False)
    computed = [r for r in results if r.get("status") == "computed_development_comparison"]
    return {"contrasts": len(results), "computed": len(computed),
            "indeterminate": sum(str(r.get("status", "")).startswith("indeterminate") for r in results),
            "overlap_coordinates": len(overlap_peaks), "queried_coordinates": len(queried_peaks),
            "seed": F6_SEED, "draws": F5_DRAWS,
            "prespecification_sha256": hashlib.sha256(U2C_PRESPEC.read_bytes()).hexdigest(),
            "control_pool": "non-overlap promoter background and original U2 peaks in the same lineage",
            "admission_used_no_model_score": True}


def gene_set_null(gene_regions: sparse.csr_matrix, region_tf: sparse.csr_matrix,
                  chosen: np.ndarray, n_genes: int, n_draws: int, rng) -> tuple:
    """Share of the gene set's regions carrying each TF, against size-matched random gene sets."""
    def profile(rows):
        mask = np.asarray(gene_regions[rows].max(axis=0).todense()).ravel() > 0
        n = int(mask.sum())
        if not n:
            return None, 0
        return (np.asarray(mask[None, :] @ region_tf).ravel() / n), n
    observed, n_observed = profile(chosen)
    if observed is None:
        return None, 0, None
    total = gene_regions.shape[0]
    draws = np.empty((n_draws, region_tf.shape[1]))
    sizes = np.empty(n_draws, dtype=int)
    for i in range(n_draws):
        pick = rng.choice(total, size=n_genes, replace=False)
        p, n = profile(pick)
        draws[i] = p if p is not None else 0.0
        sizes[i] = n
    return observed, n_observed, (draws, sizes)


def f5_family(feats: pd.DataFrame, tfs: pd.DataFrame, out: Path) -> dict:
    """Per-program TF profile against a size-matched RANDOM GENE-SET null, not a global null."""
    da = eligible_da()
    links = promoter_links(da)
    available = set(feats.region_key)
    missing = {g: sorted(r - available) for g, r in links.items() if r - available}
    background_genes = sorted(links)
    regions = sorted(set().union(*links.values()) & available)
    region_index = {k: i for i, k in enumerate(regions)}
    gene_index = {g: i for i, g in enumerate(background_genes)}

    gr = sparse.lil_matrix((len(background_genes), len(regions)), dtype=bool)
    for gene, rs in links.items():
        for r in rs & available:
            gr[gene_index[gene], region_index[r]] = True
    gene_regions = gr.tocsr()

    relevant = tfs.loc[tfs.region_key.isin(region_index) & tfs.group.isin(TF_GROUPS)].copy()
    relevant["track_key"] = relevant.group + "|" + relevant.track_name
    tracks = sorted(relevant.track_key.unique())
    track_index = {t: i for i, t in enumerate(tracks)}
    factor = dict(zip(relevant.track_key, relevant.transcription_factor))
    rt = sparse.lil_matrix((len(regions), len(tracks)), dtype=np.float64)
    for row in relevant.itertuples():
        rt[region_index[row.region_key], track_index[row.track_key]] = 1.0
    region_tf = rt.tocsr()

    programs = pd.read_csv(fam.PROGRAMS, sep="\t")
    rng = np.random.default_rng(F5_SEED)
    rows, dispositions = [], []
    for uid, part in programs.groupby("program_uid", sort=True):
        genes = sorted({g for g in part.canonical_gene.dropna() if g in gene_index})
        program_regions = set().union(*(links.get(g, set()) for g in genes)) if genes else set()
        scored = len(program_regions & available)
        disposition = {"program_uid": uid, "n_source_genes": int(part.canonical_gene.nunique()),
                       "n_genes_in_background_universe": len(genes),
                       "n_promoter_regions": len(program_regions), "n_scored_regions": scored,
                       "n_missing_regions": len(program_regions - available)}
        if scored < fam.MIN_GROUP or len(genes) < 2:
            disposition.update(status="indeterminate",
                               reason="minimum200_profile_regions_not_met" if scored < fam.MIN_GROUP
                               else "fewer_than_two_genes_in_background_universe")
            dispositions.append(disposition)
            continue
        chosen = np.array([gene_index[g] for g in genes])
        observed, n_observed, null = gene_set_null(gene_regions, region_tf, chosen,
                                                   len(genes), F5_DRAWS, rng)
        draws, sizes = null
        disposition.update(status="computed_size_matched_gene_set_null",
                           n_profile_regions=n_observed,
                           null_region_count_median=float(np.median(sizes)))
        dispositions.append(disposition)
        for track, j in track_index.items():
            column = draws[:, j]
            centre = float(np.median(column))
            extreme = int(np.sum(np.abs(column - centre) >= abs(observed[j] - centre)))
            rows.append({"program_uid": uid, "group": track.split("|", 1)[0],
                         "track": track.split("|", 1)[1],
                         "transcription_factor": factor.get(track, ""),
                         "n_genes": len(genes), "n_profile_regions": n_observed,
                         "share_program": float(observed[j]),
                         "share_random_gene_set_median": centre,
                         "null_lo": float(np.quantile(column, 0.025)),
                         "null_hi": float(np.quantile(column, 0.975)),
                         "share_minus_null": float(observed[j] - centre),
                         "n_draws": F5_DRAWS,
                         "p_nominal": float((extreme + 1) / (F5_DRAWS + 1))})
    frame = pd.DataFrame(rows)
    if len(frame):
        frame["BH_q_within_program"] = np.concatenate([
            false_discovery_control(part.p_nominal.to_numpy(), method="bh")
            for _, part in frame.groupby("program_uid", sort=False)])
    frame.to_csv(out / "F5_program_tf_profiles.tsv.gz", sep="\t", index=False)
    pd.DataFrame(dispositions).to_csv(out / "F5_program_dispositions.tsv", sep="\t", index=False)
    computed = [d for d in dispositions if d["status"] == "computed_size_matched_gene_set_null"]
    return {"programs": int(programs.program_uid.nunique()),
            "programs_computed": len(computed),
            "programs_indeterminate": len(dispositions) - len(computed),
            "eligible_background_genes": len(background_genes),
            "background_regions_available": len(regions),
            "background_genes_still_unscored": len(missing),
            "null": "size_matched_random_gene_set_from_the_eligible_background_gene_universe",
            "draws": F5_DRAWS, "seed": F5_SEED,
            "BH": "within each program's own TF family",
            "minimum_regions_each_side": fam.MIN_GROUP,
            "complete": len(missing) == 0}


def self_test() -> dict:
    """Synthetic fixtures with known answers; no family result is produced."""
    checks = {}
    # A gene set whose regions all carry track 0 must sit above a null drawn from a
    # universe where only its own genes carry it.
    n_genes, n_regions = 40, 40
    gr = sparse.csr_matrix(np.eye(n_genes, n_regions, dtype=bool))
    rt = np.zeros((n_regions, 2))
    rt[:10, 0] = 1.0
    rt[:, 1] = 1.0
    region_tf = sparse.csr_matrix(rt)
    rng = np.random.default_rng(0)
    observed, n_obs, (draws, sizes) = gene_set_null(gr, region_tf, np.arange(10), 10, 200, rng)
    checks["enriched_track_share_is_one"] = bool(abs(observed[0] - 1.0) < 1e-12)
    checks["enriched_track_above_null"] = bool(observed[0] > np.median(draws[:, 0]))
    checks["ubiquitous_track_equals_null"] = bool(abs(observed[1] - np.median(draws[:, 1])) < 1e-12)
    checks["null_sizes_match_gene_set_size"] = bool(set(sizes.tolist()) == {10})
    checks["profile_region_count"] = n_obs == 10
    # A group identical to its controls must not reject.
    frame = pd.DataFrame({"gc": rng.random(1200), "width": np.full(1200, 500),
                          "promoter": rng.integers(0, 2, 1200),
                          "signal_mean": rng.normal(size=1200), "signal_sd": rng.random(1200),
                          "value": rng.normal(size=1200)})
    group = np.zeros(1200, dtype=bool)
    group[:600] = True
    result = fam.measured_contrast(frame, "value", group)
    checks["identical_group_not_rejected"] = bool(result["status"] == "computed_development_comparison"
                                                  and result["p_nominal"] > 0.05)
    checks["small_group_is_indeterminate"] = bool(
        fam.measured_contrast(frame, "value", np.arange(1200) < 50)["status"]
        == "indeterminate_minimum_group")
    if not all(v is True or v == 1 for v in checks.values()):
        raise SystemExit(f"self-test failed: {checks}")
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reduced", type=Path, nargs="*", default=[],
                        help="Reduction output directories to read features from")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--skip-f5", action="store_true")
    parser.add_argument("--skip-f6", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    if args.self_test:
        checks = self_test()
        (args.out / "self_test.json").write_text(json.dumps(
            {"status": "statistics_self_test_passed", "checks": checks,
             "family_results_written": False}, indent=2) + "\n")
        print(json.dumps(checks, indent=2), flush=True)
        return
    if not args.reduced:
        raise SystemExit("--reduced is required unless --self-test")
    feats, tfs = read_features(list(args.reduced))
    summary = {"status": "f5_f6_families_complete",
               "deposits": [str(r) for r in args.reduced],
               "regions_read": int(len(feats)),
               "universes": {k: int(v) for k, v in feats.universe.value_counts().items()},
               "feature_definitions": "imported unchanged from the original P4 producers"}
    if not args.skip_f5:
        summary["F5"] = f5_family(feats, tfs, args.out)
    if not args.skip_f6:
        summary["F6"] = f6_family(feats, tfs, args.out)
    summary["claim_boundary"] = ("in silico saturation in a reference context; every value is a predicted "
                                 "effect for a substitution that mostly exists in no person, and is never "
                                 "presented as measured chromatin or as a variant effect in a donor")
    summary["nominates_nothing_for_followup"] = True
    summary["adopted"] = False
    (args.out / "f5_f6_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
