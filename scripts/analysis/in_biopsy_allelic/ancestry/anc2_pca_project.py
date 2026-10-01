#!/usr/bin/env python3
"""Ancestry step 2: 1kGP reference PCs at the identity SNVs and projection of every MASLD library.

Inputs
  --ref-vcf          1kGP high-coverage phased hardcalls, founders only, at the identity/tag sites
                     (anc1 extraction; GRCh38, same assembly as the BAMs).
  --identity-sites   identity_sites.tsv from anc0 (autosomal A1 targets).
  --exclude-bed      MHC / IG / TR spans (0-based BED); identity SNVs inside are dropped.
  --tag-positions    tags_final.tsv and tags_mappability.tsv (chrom, pos): identity SNVs at any
                     candidate tag position are dropped, so no MASLD genotype at a tag enters the PCs.
  --masld-vcf        A1 identity_panel.merged.vcf.gz (per-library calls, GT missing at DP < 10).
  --crosswalk        A1 sample_to_individual.tsv (run, cohort, individual_id).
  --frozen-crosswalk A2 seal frozen_crosswalk.tsv (run_role, unit_library, exclusion_reason only).
  --library-qc       B-QC library_qc.tsv (e_i = other-allele read fraction at homozygous sites).
  --ped              1kGP ped file (population and superpopulation labels).
Outputs
  --work (restricted): plink2 files; library_ancestry_pcs.tsv (filtered projection, flags, technical
      joins); library_ancestry_pcs_unfiltered.tsv (projection before the RNA-genotype site filter).
  --summary (not per-person): reference_pca_summary.json, reference_centroids.tsv,
      projection_validation_by_sites.tsv, pc_reliability_by_sites.tsv,
      heldout_distance_calibration.tsv, masld_pc_implied_reliability.tsv, rna_genotype_site_qc.tsv, rna_genotype_failing_sites.tsv,
      filter_effect_on_calls.tsv, cohort_superpop_composition.tsv,
      cohort_superpop_composition_sealed_unit.tsv, contamination_vs_ancestry_call.tsv,
      same_individual_consistency.tsv.
Rules
  * Reference = 1kGP founders (no parent listed) minus one of each KING-robust kinship >= 0.0884 pair
    (plink2 --king-cutoff), computed on the identity SNVs.
  * PCA SNVs, filters applied in this order and each counted: outside the exclusion spans; not at a
    candidate tag position; panel HWE_<POP> and ExcHet_<POP> >= 1e-6 in every superpopulation;
    MAF >= 0.05 in the unrelated reference; LD-pruned (plink2 --indep-pairwise 500kb 0.2).
    plink2 --pca 10 allele-wts on the unrelated reference.
  * Projection: plink2 --score on the allele weights with variance-standardize, reference allele
    frequencies (--read-freq) and no-mean-imputation, so a library's score is the average over the
    sites it has called. n_sites_used = ALLELE_CT / 2.
  * RNA-genotype site filter (MASLD side only; the reference PCA is unchanged). Pass 1 projects with
    every PCA SNV and calls the nearest superpopulation (>= MIN_SITES_FLOOR sites). At each identity
    SNV, MASLD units (one library per identity individual; libraries flagged possibly_mixed or from a
    contamination-excluded cohort left out) are compared with unrelated 1kGP founders: EUR-called
    units with EUR, EAS-called units with JPT (they are mostly Japanese-cohort libraries). The site
    fails if, where >= SITEQC_MIN_UNITS units are called, |ALT freq MASLD - ALT freq 1kGP| >
    SITEQC_MAX_FREQ_DIFF, or (where units x expected het >= SITEQC_MIN_EXP_HETS) observed / expected
    heterozygote fraction < SITEQC_MIN_HET_RATIO. AFR-called units are admixed, so the same test
    against 1kGP AFR is reported (diag_ columns) but not used.
    Pass 2 projects with the failing PCA SNVs treated as missing. Pass 2 is the primary result.
  * Superpopulation = nearest reference centroid (Euclidean, PC1-PC4 of the projected scores).
  * Minimum sites: prespecified rule on held-out 1kGP founders (10% per population, seed) masked
    with real MASLD call patterns: the lowest site-count bin from which every higher bin keeps the
    unmasked nearest-centroid call in >= 95% of draws. Libraries below it get no superpopulation.
  * Continuous-PC reliability: held-out 1kGP individuals of superpopulation P are projected with
    masks drawn from MASLD libraries called P (all passing libraries if fewer than
    RELIAB_MIN_POOL are called P). Per PC and site bin: r(masked, unmasked) across individual x draw.
  * Intermediate flag (set after review, 2026-09-29): between_centroids = dist_nearest /
    dist_second > BETWEEN_RATIO; far_from_centroid = dist_nearest above the 99th percentile of
    masked held-out 1kGP distances to their own centroid for the same superpopulation and site bin
    (>= CAL_TARGET draws per bin, masks from passing libraries in that bin). intermediate = either. Intermediate
    libraries keep their PCs but get ancestry_call 'intermediate' and are left out of the
    superpopulation counts.
  * Reads only identity-site genotypes (tag positions removed); no tag allele counts, no stage or
    other outcome. Technical joins: B-QC e_i, frozen-crosswalk run_role / unit_library / exclusion
    reason.
"""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

PLINK2 = "plink"  # the plink/2.0a5.13 module names its plink2 binary "plink"
PLINK_MEM_MB = "16000"  # plink2 otherwise reserves half the node's RAM, above the job's cgroup limit
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
N_PCS, N_ASSIGN_PCS = 10, 4
PCS = [f"PC{k}" for k in range(1, N_PCS + 1)]
APCS = PCS[:N_ASSIGN_PCS]
KING_CUTOFF, MIN_MAF, PRUNE, HWE_MIN_P = 0.0884, 0.05, ["500kb", "0.2"], 1e-6
HOLDOUT_FRAC, MASK_FRACTIONS, MASK_DRAWS = 0.10, [1.0, 0.5, 0.25, 0.1, 0.05], 10
AGREE_MIN, MIN_SITES_FLOOR = 0.95, 100
SITE_BINS = [0, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 10**9]
RELIAB_BINS = [100, 500, 1000, 1500, 2000, 3000, 10**9]
RELIAB_DRAWS, RELIAB_MIN_POOL, CAL_TARGET, BETWEEN_RATIO = 50, 30, 1000, 0.8
SITEQC_MIN_UNITS, SITEQC_MAX_FREQ_DIFF, SITEQC_MIN_HET_RATIO, SITEQC_MIN_EXP_HETS = 40, 0.15, 0.5, 10
SITEQC_REFERENCE = {"EUR": "EUR", "EAS": "JPT"}  # called superpopulation -> 1kGP group the site filter compares with
CONTAM_E = 0.02  # e_i above which a library is reported as possibly carrying a second person's RNA
JAPANESE_COHORTS = ["GSE174478", "GSE193066", "GSE167523"]
EAS_POPS = ["CDX", "CHB", "CHS", "JPT", "KHV"]
MAX_NUMPY_PLINK_DIFF = 1e-4  # the re-derivation must reproduce plink2's scores to this tolerance


def run(cmd):
    subprocess.run([str(c) for c in cmd], check=True)


def plink(*args, threads):
    run([PLINK2, *args, "--threads", threads, "--memory", PLINK_MEM_MB])


def read_iid_column(path):
    """First sample-id column of a plink2 id file ('#IID' or '#FID IID' header) as strings."""
    d = pd.read_csv(path, sep=r"\s+", dtype=str)
    return d["IID" if "IID" in d.columns else "#IID"]


def read_pvar(path):
    """plink2 .pvar with '##' meta lines -> DataFrame of strings (CHROM, POS, ID, REF, ALT, INFO...)."""
    with open(path) as fh:
        n_meta = sum(1 for line in fh if line.startswith("##"))
    d = pd.read_csv(path, sep="\t", skiprows=n_meta, dtype=str, keep_default_na=False)
    return d.rename(columns={"#CHROM": "CHROM"})


def info_float(info, key):
    return pd.to_numeric(info.str.extract(rf"(?:^|;){key}=([^;]+)")[0], errors="coerce")


def read_allele_weights(path):
    """eigenvec.allele -> per-variant ALT weight minus REF weight (m x K), with variant ids."""
    w = pd.read_csv(path, sep="\t")
    ref = w[w["A1"] == w["REF"]].set_index("ID")[PCS]
    alt = w[w["A1"] == w["ALT"]].set_index("ID")[PCS]
    ids = alt.index.intersection(ref.index)
    return ids, alt.loc[ids].to_numpy() - ref.loc[ids].to_numpy()


def read_alt_freq(path, ids):
    f = pd.read_csv(path, sep="\t").set_index("ID")
    return (f.loc[ids, "ALT_CTS"] / f.loc[ids, "OBS_CT"]).to_numpy()


def read_raw(path, ids):
    """plink2 --export A -> (IIDs, n x m ALT-count matrix, NaN = missing) in the order of ids.

    An id absent from the file (site not in that dataset, e.g. ALT '.' in the MASLD merge) is an
    all-missing column, never an error and never a zero.
    """
    raw = pd.read_csv(path, sep=r"\s+")
    geno_cols = [c for c in raw.columns if c.count(":") == 3]  # ids are chr:pos:ref:alt, suffixed _<counted allele>
    col_id = {c.rsplit("_", 1)[0]: c for c in geno_cols}
    g = np.full((len(raw), len(ids)), np.nan)
    for j, i in enumerate(ids):
        if i not in col_id:
            continue
        counted, (ref, alt) = col_id[i].rsplit("_", 1)[1], i.split(":")[2:4]
        assert counted in (ref, alt), f"--export A counted allele {counted} is neither REF nor ALT of {i}"
        x = raw[col_id[i]].to_numpy(dtype=float)
        g[:, j] = x if counted == alt else 2 - x  # always ALT counts
    return raw["IID"].astype(str).to_numpy(), g


def project(g, p_alt, dw):
    """Re-derivation of plink2 --score variance-standardize no-mean-imputation (score averages)."""
    z = (g - 2 * p_alt) / np.sqrt(2 * p_alt * (1 - p_alt))
    called = ~np.isnan(z)
    s = np.where(called, z, 0.0) @ dw
    n = called.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return s / (2 * n)[:, None], n


def read_sscore(path):
    s = pd.read_csv(path, sep="\t")
    s = s.rename(columns={"#IID": "IID"}).drop(columns=["#FID"], errors="ignore")
    s["IID"] = s["IID"].astype(str)
    s["n_sites_used"] = s["ALLELE_CT"] // 2
    s = s.rename(columns={f"PC{k}_AVG": f"PC{k}" for k in range(1, N_PCS + 1)})
    return s[["IID", "n_sites_used"] + PCS]


def nearest(x, centroids):
    """x: n x d; centroids: DataFrame label x d -> nearest, second, distances."""
    c = centroids.to_numpy()
    d = np.sqrt(((x[:, None, :] - c[None, :, :]) ** 2).sum(axis=2))
    order = np.argsort(d, axis=1)
    rows = np.arange(len(x))
    lab = centroids.index.to_numpy()
    return lab[order[:, 0]], lab[order[:, 1]], d[rows, order[:, 0]], d[rows, order[:, 1]]


def axis_coord(x, cent, a, b):
    """Position along the centroid axis a -> b in PC1-4 (0 at centroid a, 1 at centroid b)."""
    v = cent.loc[b].to_numpy() - cent.loc[a].to_numpy()
    return (x - cent.loc[a].to_numpy()) @ v / (v @ v)


def robust_sd(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def main():
    ap = argparse.ArgumentParser()
    for k in ["ref-vcf", "identity-sites", "exclude-bed", "masld-vcf", "crosswalk", "frozen-crosswalk",
              "library-qc", "ped", "work", "summary"]:
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--tag-positions", nargs="+", required=True)
    ap.add_argument("--threads", default="8")
    ap.add_argument("--seed", type=int, default=20260926)
    a = ap.parse_args()
    work, summ = Path(a.work), Path(a.summary)
    work.mkdir(parents=True, exist_ok=True)
    summ.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    t = a.threads
    info = {"seed": a.seed}

    ped = pd.read_csv(a.ped, sep=r"\s+", dtype=str, keep_default_na=False).set_index("SampleID")
    ids = pd.read_csv(a.identity_sites, sep="\t", keep_default_na=False)["id"]
    ids.to_csv(work / "identity_ids.txt", index=False, header=False)
    info["identity_sites_input"] = int(len(ids))

    # 1. reference pgen with allele-aware ids (identity targets fix REF,ALT; ids must match alleles)
    ref = work / "kgp_ref"
    plink("--vcf", a.ref_vcf, "--snps-only", "just-acgt", "--max-alleles", "2", "--chr", "1-22",
          "--set-all-var-ids", "@:#:$r:$a", "--make-pgen", "--out", ref, threads=t)
    kid = work / "kgp_id"
    plink("--pfile", ref, "--extract", work / "identity_ids.txt", "--make-pgen", "--out", kid, threads=t)
    pvar = read_pvar(f"{kid}.pvar")
    info["identity_sites_in_1kgp_allele_matched"] = int(len(pvar))

    # 2. unrelated founders (own --out prefix so each step keeps its plink log)
    king = work / "kgp_king"
    plink("--pfile", kid, "--king-cutoff", KING_CUTOFF, "--out", king, threads=t)
    unrel = read_iid_column(f"{king}.king.cutoff.in.id")
    pd.DataFrame({"#IID": unrel}).to_csv(work / "kgp_unrelated.keep", sep="\t", index=False)
    info["kgp_founders"] = int(sum(1 for l in open(f"{kid}.psam") if not l.startswith("#")))
    info["kgp_unrelated_founders"] = int(len(unrel))
    info["kgp_unrelated_by_superpop"] = ped.loc[unrel, "Superpopulation"].value_counts().sort_index().to_dict()
    # comparison: the panel's own unrelated set (2,504 phase-3 samples), from AN_<POP>_unrel / 2 of the first record
    first = subprocess.run(["bcftools", "query", "-f", "\t".join(f"%INFO/AN_{p}_unrel" for p in SUPERPOPS) + "\n",
                            a.ref_vcf], capture_output=True, text=True, check=True).stdout.split("\n", 1)[0].split("\t")
    info["panel_own_unrelated_by_superpop"] = {p: int(x) // 2 for p, x in zip(SUPERPOPS, first)}

    # 3. PCA SNV set, one filter at a time
    chrom, pos = pvar["CHROM"].str.replace("chr", "", regex=False), pvar["POS"].astype(int)
    bed = pd.read_csv(a.exclude_bed, sep="\t", header=None, dtype={0: str}).iloc[:, :3]
    bed.columns = ["chrom", "start", "end"]
    bed["chrom"] = bed["chrom"].str.replace("chr", "", regex=False)
    in_span = np.zeros(len(pvar), bool)
    for _, b in bed.iterrows():
        in_span |= ((chrom == b["chrom"]) & (pos - 1 >= b["start"]) & (pos - 1 < b["end"])).to_numpy()
    tag_pos = set()
    for f in a.tag_positions:
        tp = pd.read_csv(f, sep="\t", usecols=["chrom", "pos"], dtype=str, keep_default_na=False)
        tag_pos |= set(zip(tp["chrom"].str.replace("chr", "", regex=False), tp["pos"].astype(int)))
    at_tag = np.array([(c, p) in tag_pos for c, p in zip(chrom, pos)])
    hwe_fail = np.zeros(len(pvar), bool)
    for p in SUPERPOPS:
        for key in (f"HWE_{p}", f"ExcHet_{p}"):
            hwe_fail |= (info_float(pvar["INFO"], key) < HWE_MIN_P).to_numpy()
    fq = work / "kgp_freq"
    plink("--pfile", kid, "--keep", work / "kgp_unrelated.keep", "--freq", "--out", fq, threads=t)
    af = pd.read_csv(f"{fq}.afreq", sep="\t").set_index("ID")["ALT_FREQS"].reindex(pvar["ID"]).to_numpy()
    low_maf = np.minimum(af, 1 - af) < MIN_MAF
    keep, steps = np.ones(len(pvar), bool), {}
    for name, bad in [("in_mhc_ig_tr_span", in_span), ("at_candidate_tag_position", at_tag),
                      ("hwe_or_exchet_lt_1e-6_any_superpop", hwe_fail), ("maf_lt_0.05", low_maf)]:
        steps[f"removed_{name}"] = int((keep & bad).sum())
        keep &= ~bad
    steps["candidates_before_prune"] = int(keep.sum())
    cand = work / "pca_candidates.txt"
    pvar.loc[keep, "ID"].to_csv(cand, index=False, header=False)
    prune = work / "kgp_prune"
    plink("--pfile", kid, "--keep", work / "kgp_unrelated.keep", "--extract", cand,
          "--indep-pairwise", *PRUNE, "--out", prune, threads=t)
    prune_in = f"{prune}.prune.in"
    n_prune_in = sum(1 for _ in open(prune_in))
    steps["removed_ld_prune"] = steps["candidates_before_prune"] - n_prune_in
    steps["pca_snvs"] = n_prune_in
    info["pca_snv_filter_steps"] = steps
    info["identity_snvs_at_candidate_tag_positions_total"] = int(at_tag.sum())

    # 4. primary PCA on all unrelated founders
    pca = work / "kgp_pca"
    plink("--pfile", kid, "--keep", work / "kgp_unrelated.keep", "--extract", prune_in,
          "--freq", "counts", "--pca", N_PCS, "allele-wts", "--out", pca, threads=t)
    hdr = open(f"{pca}.eigenvec.allele").readline().rstrip("\n").lstrip("#").split("\t")
    id_col, a1_col = hdr.index("ID") + 1, hdr.index("A1") + 1
    pc_cols = f"{hdr.index('PC1') + 1}-{hdr.index(f'PC{N_PCS}') + 1}"
    info["eigenvalues"] = [float(x) for x in open(f"{pca}.eigenval").read().split()]

    def score(pfile, out, keep=None, weights=pca, exclude=None):
        extra = ["--keep", keep] if keep is not None else []
        extra += ["--exclude", exclude] if exclude is not None and Path(exclude).stat().st_size > 0 else []
        plink("--pfile", pfile, *extra, "--extract", prune_in, "--read-freq", f"{weights}.acount",
              "--score", f"{weights}.eigenvec.allele", id_col, a1_col, "header-read", "no-mean-imputation",
              "variance-standardize", "--score-col-nums", pc_cols, "--out", out, threads=t)
        return read_sscore(f"{out}.sscore")

    ref_proj = score(kid, work / "kgp_proj", keep=work / "kgp_unrelated.keep")
    ref_proj["superpop"] = ped.loc[ref_proj["IID"], "Superpopulation"].to_numpy()
    ref_proj["population"] = ped.loc[ref_proj["IID"], "Population"].to_numpy()
    ref_proj.to_csv(work / "kgp_reference_projection.tsv", sep="\t", index=False)
    centroids = ref_proj.groupby("superpop")[APCS].mean().loc[SUPERPOPS]
    ref_sd = ref_proj.groupby("superpop")[PCS].std()
    ref_proj.groupby("superpop")[PCS].mean().join(ref_sd, rsuffix="_sd").to_csv(summ / "reference_centroids.tsv", sep="\t")
    self_call = nearest(ref_proj[APCS].to_numpy(), centroids)[0]
    info["reference_self_assignment_agreement_by_superpop"] = (
        pd.Series(self_call == ref_proj["superpop"].to_numpy()).groupby(ref_proj["superpop"].to_numpy()).mean().round(4).to_dict())
    eas_ref = ref_proj[ref_proj["superpop"] == "EAS"]
    eas_pop_cent10 = eas_ref.groupby("population")[PCS].mean().loc[EAS_POPS]
    eas_pop_cent4 = eas_ref.groupby("population")[APCS].mean().loc[EAS_POPS]
    info["reference_eas_population_self_assignment"] = {
        "pc1_4": pd.Series(nearest(eas_ref[APCS].to_numpy(), eas_pop_cent4)[0] == eas_ref["population"].to_numpy())
        .groupby(eas_ref["population"].to_numpy()).mean().round(4).to_dict(),
        "pc1_10": pd.Series(nearest(eas_ref[PCS].to_numpy(), eas_pop_cent10)[0] == eas_ref["population"].to_numpy())
        .groupby(eas_ref["population"].to_numpy()).mean().round(4).to_dict()}

    # loading concentration: share of squared allele weight in the heaviest 1-Mb window per PC
    w_ids, dw = read_allele_weights(f"{pca}.eigenvec.allele")
    win = pd.Series(w_ids).str.split(":", expand=True)
    win = (win[0] + ":" + (win[1].astype(int) // 1_000_000).astype(str)).to_numpy()
    sq = pd.DataFrame(dw ** 2, columns=range(1, N_PCS + 1)).groupby(win).sum()
    info["pc_max_1mb_window_share_of_sq_weight"] = (sq.max() / sq.sum()).round(4).tolist()

    # 5. MASLD libraries, pass 1 (every PCA SNV)
    ms = work / "masld_id"
    plink("--vcf", a.masld_vcf, "--vcf-half-call", "m", "--snps-only", "just-acgt", "--max-alleles", "2",
          "--chr", "1-22", "--set-all-var-ids", "@:#:$r:$a", "--make-pgen", "--out", ms, threads=t)
    lib1 = score(ms, work / "masld_proj_unfiltered")
    info["masld_libraries_scored"] = int(len(lib1))
    plink("--pfile", ms, "--extract", prune_in, "--export", "A", "--out", work / "masld_prune", threads=t)
    m_iid, m_g = read_raw(f"{work}/masld_prune.raw", list(w_ids))
    p_alt = read_alt_freq(f"{pca}.acount", w_ids)

    def numpy_check(g, sc, tag):
        x, n = project(g, p_alt, dw)
        chk = sc.set_index("IID").loc[m_iid]
        ok = n > 0
        d = float(np.nanmax(np.abs(x[ok] - chk[PCS].to_numpy()[ok])))
        info[f"numpy_vs_plink2_masld_{tag}_max_abs_diff"] = d
        info[f"numpy_vs_plink2_masld_{tag}_sites_equal"] = bool((n == chk["n_sites_used"].to_numpy()).all())
        info[f"numpy_vs_plink2_masld_{tag}_within_tolerance"] = d <= MAX_NUMPY_PLINK_DIFF

    numpy_check(m_g, lib1, "unfiltered")
    info["pca_snvs_with_no_masld_call"] = int(np.isnan(m_g).all(axis=0).sum())
    lib1 = lib1.set_index("IID")
    pass1_ok = lib1["n_sites_used"] >= MIN_SITES_FLOOR
    lib1["call"] = np.where(pass1_ok, nearest(lib1[APCS].to_numpy(), centroids)[0], "")

    # 6. RNA-genotype site QC at every identity SNV not at a tag position (MASLD side)
    qc_ids = pvar.loc[~at_tag, "ID"]
    qc_ids.to_csv(work / "siteqc_ids.txt", index=False, header=False)
    plink("--pfile", ms, "--extract", work / "siteqc_ids.txt", "--export", "A", "--out", work / "masld_siteqc", threads=t)
    plink("--pfile", kid, "--keep", work / "kgp_unrelated.keep", "--extract", work / "siteqc_ids.txt",
          "--export", "A", "--out", work / "kgp_siteqc", threads=t)
    q_iid, q_g = read_raw(f"{work}/masld_siteqc.raw", list(qc_ids))
    k_iid, k_g = read_raw(f"{work}/kgp_siteqc.raw", list(qc_ids))
    k_pop = ped.loc[k_iid, "Superpopulation"].to_numpy()
    k_popl = ped.loc[k_iid, "Population"].to_numpy()
    xw = pd.read_csv(a.crosswalk, sep="\t", dtype=str, keep_default_na=False)[["run", "cohort", "individual_id"]]
    fz = pd.read_csv(a.frozen_crosswalk, sep="\t", dtype=str, keep_default_na=False,
                     usecols=["run", "run_role", "exclusion_reason", "unit_library"])
    xw = xw.merge(fz, on="run", how="left").fillna("")
    xw["possibly_mixed"] = xw["exclusion_reason"].str.contains("possibly_mixed")
    xw["cohort_contaminated"] = xw["exclusion_reason"].str.contains("cohort_excluded_contamination")
    u = xw.set_index("run").loc[q_iid]
    u["n_sites"] = lib1.loc[q_iid, "n_sites_used"].to_numpy()
    u["call"] = lib1.loc[q_iid, "call"].to_numpy()
    u["row"] = np.arange(len(u))
    u = u[(u["call"] != "") & (u["individual_id"] != "") & ~u["possibly_mixed"] & ~u["cohort_contaminated"]]
    u = u.sort_values("n_sites", ascending=False).drop_duplicates("individual_id")
    site = pd.DataFrame({"id": qc_ids.to_numpy(), "in_pca": qc_ids.isin(set(w_ids)).to_numpy()})
    site["ref"], site["alt"] = site["id"].str.split(":").str[2], site["id"].str.split(":").str[3]
    in_pca = site["in_pca"].to_numpy()

    def compare(g, kg):
        """MASLD units g vs 1kGP group kg at every site -> dict of arrays."""
        n = (~np.isnan(g)).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            f_m = np.nansum(g, axis=0) / (2 * n)
            het_m = np.nansum(g == 1, axis=0) / n
            f_k = np.nanmean(kg, axis=0) / 2
            het_e = 2 * f_k * (1 - f_k)
            ratio = het_m / het_e
        ev = n >= SITEQC_MIN_UNITS
        return {"n": n, "f_m": f_m, "f_k": f_k, "het_m": het_m, "het_e": het_e, "ratio": ratio, "ev": ev,
                "f_bad": ev & (np.abs(f_m - f_k) > SITEQC_MAX_FREQ_DIFF),
                "h_bad": ev & (n * het_e >= SITEQC_MIN_EXP_HETS) & (ratio < SITEQC_MIN_HET_RATIO),
                "flip": ev & (np.abs(f_m - (1 - f_k)) < 0.1) & (np.abs(f_m - f_k) > 0.3)}

    # descriptive comparison of every called superpopulation with the whole 1kGP superpopulation
    cmp, conc = {}, {}
    for p in SUPERPOPS:
        c = cmp[p] = compare(q_g[u.loc[u["call"] == p, "row"].to_numpy()], k_g[k_pop == p])
        site[f"units_{p}"], site[f"alt_freq_masld_{p}"], site[f"alt_freq_1kgp_{p}"] = c["n"], c["f_m"], c["f_k"]
        site[f"het_masld_{p}"], site[f"het_expected_{p}"], site[f"het_ratio_{p}"] = c["het_m"], c["het_e"], c["ratio"]
        site[f"diag_fail_freq_{p}"], site[f"diag_fail_het_{p}"] = c["f_bad"], c["h_bad"]
        if c["ev"].sum() >= 10:
            conc[p] = {"units": int((u["call"] == p).sum()), "sites_evaluated": int(c["ev"].sum()),
                       "alt_freq_correlation": round(float(np.corrcoef(c["f_m"][c["ev"]], c["f_k"][c["ev"]])[0, 1]), 4),
                       "sites_fail_freq": int(c["f_bad"].sum()), "sites_fail_het": int(c["h_bad"].sum()),
                       "flip_like_sites": int(c["flip"].sum())}
    # the filter itself: EUR-called units vs 1kGP EUR; EAS-called units vs 1kGP JPT (the EAS units are
    # mostly Japanese-cohort libraries). AFR-called units are admixed, so AFR is descriptive only.
    fail = np.zeros(len(site), bool)
    for p, ref_grp in SITEQC_REFERENCE.items():
        kg = k_g[k_pop == p] if ref_grp == p else k_g[k_popl == ref_grp]
        c = compare(q_g[u.loc[u["call"] == p, "row"].to_numpy()], kg)
        tag_ = f"{p}_vs_{ref_grp}"
        site[f"alt_freq_1kgp_{ref_grp}"] = c["f_k"] if ref_grp != p else site[f"alt_freq_1kgp_{p}"]
        site[f"het_ratio_{tag_}"], site[f"fail_freq_{tag_}"], site[f"fail_het_{tag_}"] = c["ratio"], c["f_bad"], c["h_bad"]
        fail |= c["f_bad"] | c["h_bad"]
        r_vs = float(np.corrcoef(c["f_m"][c["ev"]], c["f_k"][c["ev"]])[0, 1]) if c["ev"].sum() >= 10 else None
        conc.setdefault(p, {}).update({"reference_for_filter": ref_grp,
                        f"alt_freq_correlation_vs_{ref_grp}": round(r_vs, 4) if r_vs is not None else None,
                        f"sites_fail_freq_vs_{ref_grp}": int(c["f_bad"].sum()), f"sites_fail_het_vs_{ref_grp}": int(c["h_bad"].sum()),
                        f"pca_snvs_fail_vs_{ref_grp}": int(((c["f_bad"] | c["h_bad"]) & in_pca).sum()),
                        f"flip_like_sites_vs_{ref_grp}": int(c["flip"].sum())})
    # why AFR is descriptive only: sites failing in AFR but in neither EUR nor EAS, and whether the MASLD
    # AFR frequency moved toward the EUR frequency (the signature of European admixture)
    a_bad = cmp["AFR"]["f_bad"] | cmp["AFR"]["h_bad"]
    a_only = a_bad & ~(cmp["EUR"]["f_bad"] | cmp["EUR"]["h_bad"]) & ~(cmp["EAS"]["f_bad"] | cmp["EAS"]["h_bad"]) \
        & cmp["AFR"]["f_bad"]
    toward = np.sign(cmp["AFR"]["f_m"] - cmp["AFR"]["f_k"]) == np.sign(cmp["EUR"]["f_k"] - cmp["AFR"]["f_k"])
    eas_units = u[u["call"] == "EAS"]
    site["fail"] = fail
    site["change"] = site["ref"] + ">" + site["alt"]
    site.to_csv(summ / "rna_genotype_site_qc.tsv", sep="\t", index=False, float_format="%.4g")
    fails = site[site["fail"]]
    fails[["id", "change", "in_pca"] + [c for c in site.columns if c.startswith(("fail_", "alt_freq_", "het_ratio_", "units_"))]] \
        .to_csv(summ / "rna_genotype_failing_sites.tsv", sep="\t", index=False, float_format="%.4g")
    pca_fail = fails.loc[fails["in_pca"], "id"]
    pca_fail.to_csv(work / "pca_snvs_failing_rna_qc.txt", index=False, header=False)
    info["rna_genotype_site_qc"] = {
        "rule": (f"fail if |ALT freq MASLD - ALT freq 1kGP| > {SITEQC_MAX_FREQ_DIFF}, or (units x expected het >= "
                 f"{SITEQC_MIN_EXP_HETS}) observed/expected het < {SITEQC_MIN_HET_RATIO}, in EUR-called MASLD units vs 1kGP EUR "
                 f"or in EAS-called MASLD units vs 1kGP JPT, where >= {SITEQC_MIN_UNITS} units are called (units: one library per "
                 f"identity individual, pass-1 call, possibly_mixed and contamination-excluded cohorts left out). The "
                 f"diag_fail_* columns apply the same rule against every whole superpopulation and are not used."),
        "flip_like_rule": "|f_MASLD - (1 - f_1kGP)| < 0.1 and |f_MASLD - f_1kGP| > 0.3",
        "identity_snvs_evaluated": int(len(site)), "identity_snvs_failing": int(fail.sum()),
        "pca_snvs_failing": int(len(pca_fail)),
        "failing_by_change": fails["change"].value_counts().to_dict(),
        "by_superpop": conc,
        "why_afr_is_descriptive_only": {
            "afr_freq_failures_not_failing_in_eur_or_eas": int(a_only.sum()),
            "fraction_of_those_moved_toward_eur_frequency": round(float(toward[a_only].mean()), 3) if a_only.any() else None},
        "why_eas_uses_jpt": {
            "eas_units": int(len(eas_units)),
            "eas_units_from_japanese_cohorts": int(eas_units["cohort"].isin(JAPANESE_COHORTS).sum()),
            "eas_freq_failures_vs_all_eas": int(cmp["EAS"]["f_bad"].sum())}}

    # 7. MASLD libraries, pass 2 (failing PCA SNVs treated as missing) = primary result
    lib = score(ms, work / "masld_proj", exclude=work / "pca_snvs_failing_rna_qc.txt")
    drop = np.isin(np.asarray(w_ids), pca_fail.to_numpy())
    m_g2 = m_g.copy()
    m_g2[:, drop] = np.nan
    numpy_check(m_g2, lib, "filtered")

    # 8. held-out validation with real MASLD call patterns (pass-2 masks)
    keep_ids = ref_proj[["IID", "population"]].sort_values("IID")
    hold_ids = []
    for _, d in keep_ids.groupby("population", sort=True):
        n = max(1, int(round(HOLDOUT_FRAC * len(d))))
        hold_ids += list(rng.choice(d["IID"].to_numpy(), size=n, replace=False))
    hold = keep_ids[keep_ids["IID"].isin(hold_ids)]
    train = keep_ids[~keep_ids["IID"].isin(hold_ids)]
    info["heldout_individuals"] = int(len(hold))
    train[["IID"]].rename(columns={"IID": "#IID"}).to_csv(work / "val_train.keep", sep="\t", index=False)
    hold[["IID"]].rename(columns={"IID": "#IID"}).to_csv(work / "val_hold.keep", sep="\t", index=False)
    vp = work / "val_pca"
    plink("--pfile", kid, "--keep", work / "val_train.keep", "--extract", prune_in,
          "--freq", "counts", "--pca", N_PCS, "allele-wts", "--out", vp, threads=t)
    tr_proj = score(kid, work / "val_train_proj", keep=work / "val_train.keep", weights=vp)
    ho_proj = score(kid, work / "val_hold_proj", keep=work / "val_hold.keep", weights=vp)
    tr_proj["superpop"] = ped.loc[tr_proj["IID"], "Superpopulation"].to_numpy()
    tr_proj["population"] = ped.loc[tr_proj["IID"], "Population"].to_numpy()
    v_cent = tr_proj.groupby("superpop")[APCS].mean().loc[SUPERPOPS]
    v_eas10 = tr_proj[tr_proj["superpop"] == "EAS"].groupby("population")[PCS].mean().loc[EAS_POPS]
    v_eas4 = tr_proj[tr_proj["superpop"] == "EAS"].groupby("population")[APCS].mean().loc[EAS_POPS]
    v_ids, v_dw = read_allele_weights(f"{vp}.eigenvec.allele")
    v_p = read_alt_freq(f"{vp}.acount", v_ids)
    plink("--pfile", kid, "--keep", work / "val_hold.keep", "--extract", prune_in, "--export", "A",
          "--out", work / "val_hold", threads=t)
    h_iid, h_g = read_raw(f"{work}/val_hold.raw", list(v_ids))
    h_full, _ = project(h_g, v_p, v_dw)
    ho = ho_proj.set_index("IID").loc[h_iid]
    info["numpy_vs_plink2_heldout_max_abs_diff"] = float(np.abs(h_full - ho[PCS].to_numpy()).max())
    info["numpy_vs_plink2_heldout_within_tolerance"] = info["numpy_vs_plink2_heldout_max_abs_diff"] <= MAX_NUMPY_PLINK_DIFF
    tr_var = tr_proj[PCS].var().to_numpy()
    info["heldout_over_train_variance_by_pc_unmasked"] = (np.var(h_full, axis=0, ddof=1) / tr_var).round(3).tolist()
    h_label = ped.loc[h_iid, "Superpopulation"].to_numpy()
    h_popl = ped.loc[h_iid, "Population"].to_numpy()
    h_call_full = nearest(h_full[:, :N_ASSIGN_PCS], v_cent)[0]
    info["heldout_full_data_label_agreement"] = float((h_call_full == h_label).mean())

    _, m_gv = read_raw(f"{work}/masld_prune.raw", list(v_ids))
    m_gv[:, np.isin(np.asarray(v_ids), pca_fail.to_numpy())] = np.nan
    masks = ~np.isnan(m_gv)
    rows = []
    for frac in MASK_FRACTIONS:
        for draw in range(MASK_DRAWS):
            pick = masks[rng.integers(len(masks), size=len(h_g))].copy()
            if frac < 1.0:
                pick &= rng.random(pick.shape) < frac
            g = np.where(pick, h_g, np.nan)
            x, n = project(g, v_p, v_dw)
            good = n > 0
            call = np.full(len(g), "", dtype=object)
            call[good] = nearest(x[good, :N_ASSIGN_PCS], v_cent)[0]
            rows.append(pd.DataFrame({"frac": frac, "draw": draw, "iid": h_iid, "superpop": h_label, "n_sites": n,
                                      "agree_full": call == h_call_full, "agree_label": call == h_label}))
    val = pd.concat(rows, ignore_index=True)
    val["site_bin"] = pd.cut(val["n_sites"], SITE_BINS, right=False)
    by_bin = val.groupby("site_bin", observed=True).agg(
        draws=("agree_full", "size"), agree_full=("agree_full", "mean"), agree_label=("agree_label", "mean"))
    by_bin_pop = val.pivot_table(index="site_bin", columns="superpop", values="agree_label", aggfunc="mean", observed=True)
    by_bin = by_bin.join(by_bin_pop.add_prefix("agree_label_")).reset_index()
    by_bin["site_bin_low"] = [iv.left for iv in by_bin["site_bin"]]
    by_bin.to_csv(summ / "projection_validation_by_sites.tsv", sep="\t", index=False, float_format="%.4f")

    # prespecified minimum: lowest bin from which all higher populated bins agree >= AGREE_MIN
    ok_bins = (by_bin["agree_full"] >= AGREE_MIN).to_numpy()
    suffix_ok = np.flip(np.cumprod(np.flip(ok_bins)).astype(bool))
    min_sites = int(by_bin.loc[suffix_ok, "site_bin_low"].min()) if suffix_ok.any() else None
    min_sites = max(min_sites, MIN_SITES_FLOOR) if min_sites is not None else None
    info["min_sites_rule"] = (f"lowest site bin from which every higher bin keeps the unmasked call in >= {AGREE_MIN} "
                              f"of masked held-out draws; floor {MIN_SITES_FLOOR}")
    info["min_sites"] = min_sites

    # 9. continuous-PC reliability and distance calibration: held-out individuals of superpop P with
    #    masks from MASLD libraries called P (pass 2)
    lib = lib.set_index("IID").loc[m_iid].reset_index()
    lib_pass = (lib["n_sites_used"] >= min_sites).to_numpy() if min_sites is not None else np.zeros(len(lib), bool)
    lib_call = np.where(lib_pass, nearest(lib[APCS].to_numpy(), centroids)[0], "")
    h_ax = {"axis_EUR_AFR": axis_coord(h_full[:, :4], v_cent, "EUR", "AFR"),
            "axis_EUR_EAS": axis_coord(h_full[:, :4], v_cent, "EUR", "EAS")}
    rel_rows, mask_source = [], {}
    for p in SUPERPOPS:
        sel = np.where(h_label == p)[0]
        pool_idx = np.where(lib_pass & (lib_call == p))[0]
        mask_source[p] = f"libraries called {p} (n={len(pool_idx)})"
        if len(pool_idx) < RELIAB_MIN_POOL:
            pool_idx = np.where(lib_pass)[0]
            mask_source[p] = f"all passing libraries (only {int((lib_pass & (lib_call == p)).sum())} called {p})"
        if len(pool_idx) == 0:  # smoke runs on three chromosomes: no library reaches the minimum
            pool_idx = np.where(masks.any(axis=1))[0]
            mask_source[p] = "all libraries with any call (none passes the minimum)"
        for draw in range(RELIAB_DRAWS):
            pick = masks[pool_idx[rng.integers(len(pool_idx), size=len(sel))]]
            x, n = project(np.where(pick, h_g[sel], np.nan), v_p, v_dw)
            n1, n2, d1, d2 = nearest(x[:, :N_ASSIGN_PCS], v_cent)
            own = np.sqrt(((x[:, :4] - v_cent.loc[p].to_numpy()) ** 2).sum(axis=1))
            r = pd.DataFrame({"superpop": p, "draw": draw, "iid": h_iid[sel], "population": h_popl[sel], "n_sites": n,
                              "dist_own": own, "ratio": d1 / d2, "call": n1})
            for k, pc in enumerate(PCS):
                r[f"m_{pc}"], r[f"u_{pc}"] = x[:, k], h_full[sel, k]
            for axn, (ca, cb) in {"axis_EUR_AFR": ("EUR", "AFR"), "axis_EUR_EAS": ("EUR", "EAS")}.items():
                r[f"m_{axn}"], r[f"u_{axn}"] = axis_coord(x[:, :4], v_cent, ca, cb), h_ax[axn][sel]
            if p == "EAS":
                r["eas_pop_pc1_4"] = nearest(x[:, :4], v_eas4)[0]
                r["eas_pop_pc1_10"] = nearest(x, v_eas10)[0]
            rel_rows.append(r)
    rel = pd.concat(rel_rows, ignore_index=True)
    rel["site_bin"] = pd.cut(rel["n_sites"], RELIAB_BINS, right=False).astype(str)
    measures = PCS + ["axis_EUR_AFR", "axis_EUR_EAS"]
    out = []
    for (p, b), d in list(rel.groupby(["superpop", "site_bin"])) + [((p, "all"), d) for p, d in rel.groupby("superpop")] \
            + [(("ALL", "all"), rel)]:
        for mname in measures:
            mm, uu = d[f"m_{mname}"].to_numpy(), d[f"u_{mname}"].to_numpy()
            ok = ~np.isnan(mm)
            if ok.sum() < 20:
                continue
            rr = float(np.corrcoef(mm[ok], uu[ok])[0, 1])
            out.append({"superpop": p, "site_bin": b, "measure": mname, "pairs": int(ok.sum()),
                        "individuals": int(d["iid"].nunique()), "r_masked_unmasked": rr, "r2": rr * rr,
                        "sd_unmasked": float(np.std(uu[ok], ddof=1)), "sd_masked": float(np.std(mm[ok], ddof=1)),
                        "sd_masked_minus_unmasked": float(np.std(mm[ok] - uu[ok], ddof=1)),
                        "mean_masked_minus_unmasked": float(np.mean(mm[ok] - uu[ok]))})
    reliab = pd.DataFrame(out)
    reliab.to_csv(summ / "pc_reliability_by_sites.tsv", sep="\t", index=False, float_format="%.4g")
    allbin = reliab[reliab["site_bin"] == "all"].pivot(index="superpop", columns="measure", values="r2")[measures]
    info["pc_reliability"] = {
        "definition": ("r2 of masked vs unmasked projection across held-out individual x draw; masks from MASLD libraries "
                       f"called the same superpopulation; {RELIAB_DRAWS} draws; val PCA (train 90%)"),
        "mask_source": mask_source,
        "r2_all_site_bins": {p: {m: round(float(v), 3) for m, v in allbin.loc[p].items()} for p in allbin.index}}
    eas_rel = rel[rel["superpop"] == "EAS"]
    jpt = {}
    for col in ["eas_pop_pc1_4", "eas_pop_pc1_10"]:
        is_j, call_j = eas_rel["population"] == "JPT", eas_rel[col] == "JPT"
        jpt[col] = {"jpt_recall": round(float((is_j & call_j).sum() / is_j.sum()), 4),
                    "jpt_precision": round(float((is_j & call_j).sum() / max(call_j.sum(), 1)), 4),
                    "all_eas_population_agreement": round(float((eas_rel[col] == eas_rel["population"]).mean()), 4)}
    info["heldout_eas_population_assignment_masked"] = jpt

    # distance calibration for the far_from_centroid flag: masks stratified by site bin, so every bin
    # that holds a passing library gets >= CAL_TARGET draws (same-called pool when it has
    # >= RELIAB_MIN_POOL libraries in that bin, else all passing libraries in that bin)
    lib_bin = pd.cut(lib["n_sites_used"], RELIAB_BINS, right=False).astype(str).to_numpy()
    cal_rows = []
    for p in SUPERPOPS:
        sel = np.where(h_label == p)[0]
        for b in sorted(set(lib_bin[lib_pass])):
            in_b = lib_pass & (lib_bin == b)
            pool_idx, src = np.where(in_b & (lib_call == p))[0], "called"
            if len(pool_idx) < RELIAB_MIN_POOL:
                pool_idx, src = np.where(in_b)[0], "all_in_bin"
            dists, ratios, oks = [], [], []
            for draw in range(int(np.ceil(CAL_TARGET / len(sel)))):
                pick = masks[pool_idx[rng.integers(len(pool_idx), size=len(sel))]]
                x, _ = project(np.where(pick, h_g[sel], np.nan), v_p, v_dw)
                n1, _, d1, d2 = nearest(x[:, :N_ASSIGN_PCS], v_cent)
                dists.append(np.sqrt(((x[:, :4] - v_cent.loc[p].to_numpy()) ** 2).sum(axis=1)))
                ratios.append(d1 / d2)
                oks.append(n1 == p)
            dists, ratios, oks = np.concatenate(dists), np.concatenate(ratios), np.concatenate(oks)
            cal_rows.append({"superpop": p, "site_bin": b, "mask_pool": src, "mask_pool_libraries": int(len(pool_idx)),
                             "draws": int(len(dists)), "q99_dist_own": float(np.quantile(dists, 0.99)),
                             "frac_ratio_gt": float((ratios > BETWEEN_RATIO).mean()), "agree_label": float(oks.mean())})
    cal = pd.DataFrame(cal_rows)
    cal.to_csv(summ / "heldout_distance_calibration.tsv", sep="\t", index=False, float_format="%.4g")
    q99 = cal.set_index(["superpop", "site_bin"])["q99_dist_own"]

    # MASLD spread vs what masking alone predicts, and the implied reliability of each coordinate in
    # MASLD itself: 1 - (masking-noise robust SD / MASLD robust SD)^2, by called superpopulation
    lib_ax = {"axis_EUR_AFR": axis_coord(lib[APCS].to_numpy(), centroids, "EUR", "AFR"),
              "axis_EUR_EAS": axis_coord(lib[APCS].to_numpy(), centroids, "EUR", "EAS")}
    spread, imp_rows = {}, []
    for p in ["EUR", "EAS", "AFR"]:
        sel_l = lib_pass & (lib_call == p)
        if sel_l.sum() < RELIAB_MIN_POOL:
            continue
        d = rel[rel["superpop"] == p]
        spread[p] = {"masld_libraries": int(sel_l.sum()),
                     "masld_robust_sd_pc1_4": [round(robust_sd(lib.loc[sel_l, pc]), 4) for pc in APCS],
                     "masked_heldout_robust_sd_pc1_4": [round(robust_sd(d[f"m_{pc}"]), 4) for pc in APCS],
                     "unmasked_1kgp_robust_sd_pc1_4": [round(robust_sd(ref_proj.loc[ref_proj["superpop"] == p, pc]), 4) for pc in APCS]}
        for mname in measures:
            v = lib.loc[sel_l, mname].to_numpy() if mname in PCS else lib_ax[mname][sel_l]
            sd_m, sd_n = robust_sd(v), robust_sd(d[f"m_{mname}"] - d[f"u_{mname}"])
            imp_rows.append({"superpop": p, "measure": mname, "masld_libraries": int(sel_l.sum()),
                             "masld_robust_sd": sd_m, "masking_noise_robust_sd": sd_n,
                             "implied_reliability_masld": max(0.0, 1 - (sd_n / sd_m) ** 2),
                             "heldout_1kgp_r2": float(reliab.loc[(reliab["superpop"] == p) & (reliab["site_bin"] == "all")
                                                               & (reliab["measure"] == mname), "r2"].iloc[0])})
    info["spread_masld_vs_masked_heldout"] = spread
    imp = pd.DataFrame(imp_rows)
    imp.to_csv(summ / "masld_pc_implied_reliability.tsv", sep="\t", index=False, float_format="%.4g")
    info["implied_reliability_masld"] = {
        "definition": ("max(0, 1 - (robust SD of masked-minus-unmasked held-out 1kGP scores, masks from libraries called P)^2 / "
                       "(robust SD of MASLD libraries called P)^2); an upper bound, because RNA genotype errors add noise "
                       "that 1kGP masking does not"),
        "values": {p: {r["measure"]: round(r["implied_reliability_masld"], 3) for _, r in g.iterrows()}
                   for p, g in imp.groupby("superpop")}}

    # 10. per-library table (restricted)
    lib = xw.merge(lib.rename(columns={"IID": "run"}), on="run", how="left")
    lib["n_sites_used"] = lib["n_sites_used"].fillna(0).astype(int)
    lib["passes_min_sites"] = (lib["n_sites_used"] >= min_sites) if min_sites is not None else False
    x = np.nan_to_num(lib[APCS].to_numpy())
    n1, n2, d1, d2 = nearest(x, centroids)
    lib["axis_EUR_AFR"] = axis_coord(lib[APCS].to_numpy(), centroids, "EUR", "AFR")
    lib["axis_EUR_EAS"] = axis_coord(lib[APCS].to_numpy(), centroids, "EUR", "EAS")
    ps = lib["passes_min_sites"].to_numpy()
    lib["nearest_superpop"] = np.where(ps, n1, "")
    lib["second_superpop"] = np.where(ps, n2, "")
    lib["dist_nearest"] = np.where(ps, d1, np.nan)
    lib["dist_second"] = np.where(ps, d2, np.nan)
    lib["dist_ratio"] = lib["dist_nearest"] / lib["dist_second"]
    lib["between_centroids"] = ps & (lib["dist_ratio"] > BETWEEN_RATIO)
    lib["site_bin"] = pd.cut(lib["n_sites_used"], RELIAB_BINS, right=False).astype(str)
    lib["dist_q99_masked_heldout"] = [q99.get((s, b), np.nan) if p else np.nan
                                      for s, b, p in zip(lib["nearest_superpop"], lib["site_bin"], ps)]
    lib["far_from_centroid"] = ps & (lib["dist_nearest"] > lib["dist_q99_masked_heldout"])
    lib["intermediate"] = lib["between_centroids"] | lib["far_from_centroid"]
    lib["ancestry_call"] = np.where(~ps, "", np.where(lib["intermediate"], "intermediate", lib["nearest_superpop"]))
    eas = ps & (lib["nearest_superpop"] == "EAS").to_numpy()
    lib["nearest_eas_population_pc1_10"] = ""
    lib.loc[eas, "nearest_eas_population_pc1_10"] = nearest(lib.loc[eas, PCS].to_numpy(), eas_pop_cent10)[0]
    u1 = lib1.rename(columns={pc: f"{pc}_unfiltered" for pc in PCS}).rename(
        columns={"n_sites_used": "n_sites_unfiltered", "call": "nearest_superpop_unfiltered"})
    lib = lib.merge(u1[["n_sites_unfiltered", "nearest_superpop_unfiltered"] + [f"{pc}_unfiltered" for pc in APCS]],
                    left_on="run", right_index=True, how="left")
    qc = pd.read_csv(a.library_qc, sep="\t", dtype={"run": str}, usecols=["run", "e_i"])
    lib = lib.merge(qc, on="run", how="left")
    lib.drop(columns=["site_bin"]).to_csv(work / "library_ancestry_pcs.tsv", sep="\t", index=False, float_format="%.6g")
    u1.reset_index().rename(columns={"IID": "run"}).to_csv(work / "library_ancestry_pcs_unfiltered.tsv", sep="\t",
                                                           index=False, float_format="%.6g")
    p = lib[lib["passes_min_sites"]]
    info["libraries_passing_min_sites"] = int(len(p))
    info["intermediate_flags"] = {
        "rule": (f"between_centroids: dist_nearest/dist_second > {BETWEEN_RATIO}; far_from_centroid: dist_nearest > 99th "
                 f"percentile of masked held-out 1kGP distance to own centroid, same superpopulation and site bin "
                 f"(>= {CAL_TARGET} draws per bin, masks from passing libraries in that bin); intermediate = either"),
        "between_centroids": int(p["between_centroids"].sum()), "far_from_centroid": int(p["far_from_centroid"].sum()),
        "intermediate": int(p["intermediate"].sum()),
        "far_not_evaluated_no_calibration": int(p["dist_q99_masked_heldout"].isna().sum())}

    # effect of the RNA-genotype site filter on calls and PCs
    both = p[p["nearest_superpop_unfiltered"] != ""]
    pd.crosstab(both["nearest_superpop_unfiltered"], both["nearest_superpop"]).to_csv(summ / "filter_effect_on_calls.tsv", sep="\t")
    info["site_filter_effect"] = {
        "libraries_compared": int(len(both)),
        "nearest_superpop_changed": int((both["nearest_superpop_unfiltered"] != both["nearest_superpop"]).sum()),
        "median_abs_change_pc1_4": [round(float(np.median(np.abs(both[pc] - both[f"{pc}_unfiltered"]))), 5) for pc in APCS],
        "max_abs_change_pc1_4": [round(float(np.max(np.abs(both[pc] - both[f"{pc}_unfiltered"]))), 5) for pc in APCS],
        "median_sites_lost": float(np.median(both["n_sites_unfiltered"] - both["n_sites_used"]))}

    # 11. cohort composition: libraries; all identity individuals; sealed analysis unit
    def counts(d, by):
        c = d.groupby(by + ["ancestry_call"]).size().unstack(fill_value=0)
        return c.reindex(columns=SUPERPOPS + ["intermediate"], fill_value=0)

    tot = lib.groupby("cohort").agg(libraries=("run", "size"), libraries_passing=("passes_min_sites", "sum"),
                                    median_sites=("n_sites_used", "median"))
    comp_lib = counts(p, ["cohort"])
    ind = p[p["individual_id"] != ""].sort_values("n_sites_used", ascending=False).drop_duplicates("individual_id")
    comp_ind = counts(ind, ["cohort"])
    no_id = p[p["individual_id"] == ""].groupby("cohort").size().rename("lib_passing_no_individual_id")
    clean = p[~p["intermediate"]]
    comp = tot.join(comp_lib.add_prefix("lib_")).join(no_id).join(comp_ind.add_prefix("ind_all_identity_")).fillna(0)
    comp = comp.join(clean.groupby("cohort")[APCS].median().add_prefix("median_"))
    comp = comp.join(clean.groupby("cohort")["n_sites_used"].median().rename("median_sites_passing_not_intermediate"))
    comp.reset_index().to_csv(summ / "cohort_superpop_composition.tsv", sep="\t", index=False, float_format="%.4g")
    info["composition_notes"] = {
        "lib_<POP>": "passing libraries called <POP> and not intermediate; lib_intermediate counted separately",
        "lib_passing_no_individual_id": ("passing libraries with an empty individual_id: A1 identity gave them no individual "
                                         "because they have < 1,000 called autosomal identity sites (A1 passes_min_sites False)"),
        "ind_all_identity_<POP>": ("all A1 identity individuals, each counted once via its library with the most sites; includes "
                                   "GSE193066 and other libraries outside the allelic analysis")}
    su = p[p["unit_library"] == "True"]
    comp_su = counts(su, ["cohort", "run_role"])
    unit_all = lib[lib["unit_library"] == "True"].groupby(["cohort", "run_role"]).agg(
        unit_libraries=("run", "size"), unit_libraries_passing=("passes_min_sites", "sum"))
    unit_all.join(comp_su).fillna(0).reset_index().to_csv(summ / "cohort_superpop_composition_sealed_unit.tsv", sep="\t",
                                                           index=False, float_format="%.4g")
    info["sealed_unit_libraries"] = int((lib["unit_library"] == "True").sum())
    info["sealed_unit_libraries_passing"] = int(len(su))
    info["japanese_cohorts"] = {
        c: {"passing": int(len(d)), "EAS_not_intermediate": int((d["ancestry_call"] == "EAS").sum()),
            "nearest_EAS": int((d["nearest_superpop"] == "EAS").sum()),
            "median_PC1_PC2_not_intermediate": [round(float(d.loc[~d["intermediate"], pc].median()), 4) for pc in ["PC1", "PC2"]],
            "median_sites": float(d["n_sites_used"].median()),
            "nearest_eas_population_pc1_10_JPT": int((d["nearest_eas_population_pc1_10"] == "JPT").sum())}
        for c, d in p.groupby("cohort") if c in JAPANESE_COHORTS}
    info["eas_centroid_pc1_pc2"] = [round(float(centroids.loc["EAS", pc]), 4) for pc in ["PC1", "PC2"]]
    info["eas_reference_sd_pc1_pc2"] = [round(float(ref_sd.loc["EAS", pc]), 4) for pc in ["PC1", "PC2"]]

    # 12. contamination cross-tab (technical quantities only): cohort majority vs minority calls
    rows = []
    for coh, d in p.groupby("cohort"):
        maj = d.loc[~d["intermediate"], "nearest_superpop"].mode()
        maj = maj.iloc[0] if len(maj) else ""
        cat = np.where(d["intermediate"], "intermediate", np.where(d["nearest_superpop"] == maj, "majority", "minority"))
        for c in ["majority", "minority", "intermediate"]:
            dd = d[cat == c]
            rows.append({"cohort": coh, "majority_superpop": maj, "category": c, "libraries": len(dd),
                         f"e_i_gt_{CONTAM_E}": int((dd["e_i"] > CONTAM_E).sum()), "e_i_missing": int(dd["e_i"].isna().sum()),
                         "possibly_mixed": int(dd["possibly_mixed"].sum()),
                         "median_e_i": float(dd["e_i"].median()) if len(dd) else np.nan,
                         "calls": ",".join(f"{k}:{v}" for k, v in dd["nearest_superpop"].value_counts().items())})
    pd.DataFrame(rows).to_csv(summ / "contamination_vs_ancestry_call.tsv", sep="\t", index=False, float_format="%.4g")

    # 13. same-individual libraries should project close together
    rows = []
    for coh, d in p[p["individual_id"] != ""].groupby("cohort"):
        xy = d[APCS].to_numpy()
        iid = d["individual_id"].to_numpy()
        dist = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))
        iu = np.triu_indices(len(d), 1)
        same = iid[iu[0]] == iid[iu[1]]
        if same.any():
            rows.append({"cohort": coh, "same_individual_pairs": int(same.sum()),
                         "median_dist_same": float(np.median(dist[iu][same])),
                         "median_dist_different": float(np.median(dist[iu][~same]))})
    pd.DataFrame(rows).to_csv(summ / "same_individual_consistency.tsv", sep="\t", index=False, float_format="%.4g")

    (summ / "reference_pca_summary.json").write_text(json.dumps(info, indent=2, default=str))
    print(json.dumps(info, indent=2, default=str))


if __name__ == "__main__":
    main()
