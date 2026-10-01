"""Collect the per-gene duplicate counts (07) into library covariates and summaries (B-Dup v2).

The definitions below were fixed in dup/PROGRESS.md before any 07 output existed.

One-end rule (same for all libraries)
  End key (chrom, unclipped 5' position, strand) of read 1 (paired-end) or of the read
  (single-end); reads assigned to one tag gene by aligned blocks on the library's
  transcribed strand (07_gene_dup.py). Duplicates in gene g = reads1_g - distinct1_g.
  RPK_g = reads1_g / (exon_bp_g / 1000), with exon_bp from tag_genes.tsv.
Library covariates
  dup1_sub_rpk2to20  PRIMARY. On the 50,000-read subsample: sum of duplicates / sum of
                     reads over genes with 2 <= RPK < 20 (RPK from the subsample). Same
                     rule, same depth and same expression band for every library, so it
                     can be compared across layouts and cohorts. NA if the library has
                     < 50,000 assigned reads, or the band holds < 5,000 reads or < 50 genes.
  dup1_rpk2to20      the same band at full depth (same NA rule). Rises with depth by
                     construction (more reads per molecule).
  dup1_pooled, dup1_sub_pooled
                     all assigned genes, full depth / subsample. Weighted toward the
                     most-read genes; descriptive only.
  pair_dup_fraction  paired-end only: sum(pairs - distinct_pairs) / sum(pairs) over assigned
                     read-1s whose mate is in the gene-body subset (pair rule, mate-absent
                     reads excluded from numerator and denominator). Within paired-end only.
  genebody_pair_dup_fraction
                     paired-end only: markdup DUPLICATE PAIR / PAIRED over the whole
                     gene-body subset (both in reads; single reads excluded).
  exon_dup_fraction  step 3 (markdup): pair rule for paired-end, one-end rule for
                     single-end, pooled over tag-exon reads. Within one layout only.
  top20_share        share of assigned one-end reads (full depth) in the library's 20
                     most-read tag genes.
Units
  *_reads columns count reads; for paired-end libraries, pairs_* and
  estimated_library_size count pairs, genebody_lw_distinct_reads counts reads (two per
  pair), star_* counts pairs (STAR Log.final.out).
Crosswalk
  Only run, run_role and unit_library are read from frozen_crosswalk.tsv (usecols); its
  stage columns are never loaded.
Outputs (refuses to overwrite)
  <restricted>/dup_per_library_v2.tsv, <restricted>/dup_per_library_gene.tsv.gz
  <exec>/dup_per_cohort_role.tsv   unit libraries, cohort x run_role (+ all_paired, all_single)
  <exec>/dup_spearman_within_cohort.tsv
  <exec>/dup_cohort_prep.tsv        instrument, kit, FFPE, capture per cohort, with sources
  <exec>/dup_columns.tsv            column dictionary
  <exec>/dup_missing_libraries_v2.tsv
  <exec>/sha256sums.txt             every file above plus the step-3 tables
Counts only; no allele, genotype or stage value is read or written.
"""
import argparse
import glob
import gzip
import hashlib
import os
import sys

import numpy as np
import pandas as pd

N_SUB = 50_000
BAND = (2.0, 20.0)
MIN_BAND_READS, MIN_BAND_GENES = 5_000, 50
MEASURES = ["dup1_sub_rpk2to20", "dup1_rpk2to20", "dup1_pooled", "dup1_sub_pooled",
            "pair_dup_fraction", "genebody_pair_dup_fraction", "exon_dup_fraction"]

# Library preparation as stated by the submitters. Each quote is checked against its
# source file at run time; the instrument comes from the ENA run report itself.
PREP = {
    "GSE126848": dict(kit="TruSeq Stranded mRNA (NeoPrep)", selection="poly(A)", capture="no", ffpe="not stated",
                      tissue="needle biopsy", quotes=["TruSeq® Stranded mRNA Library Prep for NeoPrep"]),
    "GSE130970": dict(kit="TruSeq Stranded mRNA", selection="poly(A)", capture="no", ffpe="not stated",
                      tissue="liver biopsy", quotes=["Illumina TruSeq Stranded mRNA Sample Preparation kit"]),
    "GSE135251": dict(kit="not stated ('standard Illumina protocols')", selection="not stated", capture="no",
                      ffpe="no (snap-frozen)", tissue="liver biopsy",
                      quotes=["216 snap-frozen biopsies", "standard Illumina protocols"]),
    "GSE162694": dict(kit="not stated", selection="not stated", capture="not stated", ffpe="not stated",
                      tissue="liver biopsy", quotes=["none provided by the submitter"]),
    "GSE167523": dict(kit="TruSeq Stranded mRNA", selection="poly(A)", capture="no", ffpe="no (frozen)",
                      tissue="liver biopsy", quotes=["mRNA was extracted from frozen tissues",
                                                     "TruSeq stranded mRNA Library Prep Kit"]),
    "GSE174478": dict(kit="TruSeq Stranded Total RNA (rRNA depletion)", selection="rRNA depletion", capture="no",
                      ffpe="not stated", tissue="liver biopsy",
                      quotes=["TruSeq Stranded Total RNA LT Sample Prep Kit (Gold)"]),
    "GSE193066": dict(kit="TruSeq RNA Exome", selection="exome capture", capture="yes", ffpe="yes",
                      tissue="liver biopsy", quotes=["FFPE tissue using High Pure RNA Paraffin kit",
                                                     "Truseq RNA Exome kit"]),
    "GSE213621": dict(kit="NEBNext Ultra RNA Library Prep Kit", selection="not stated", capture="no",
                      ffpe="not stated", tissue="liver tissue",
                      quotes=["NEBNext Ultra RNA Library Prep Kit"]),
    "GSE240729": dict(kit="NEBNext Ultra II Directional (rRNA depletion)", selection="rRNA depletion", capture="no",
                      ffpe="yes", tissue="liver biopsy",
                      quotes=["formalin-fixed paraffin-embedded (FFPE)", "NEBNext Ultra II Directional RNA Library Prep Kit",
                              "rRNA depletion kit"]),
    "PRJNA512027": dict(kit="not stated in SRA/BioProject", selection="not stated", capture="not stated",
                        ffpe="not stated", tissue="liver",
                        quotes=["Transcriptomic Profiling of Obesity-Related Nonalcoholic Steatohepatitis"]),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_markdup(path):
    values = {}
    with open(path) as fh:
        for line in fh:
            key, _, value = line.partition(":")
            if key in ("PAIRED", "SINGLE", "DUPLICATE PAIR", "DUPLICATE SINGLE"):
                values[key] = int(value)
    return values


def band_fraction(reads, distinct, kb):
    rpk = reads / kb
    in_band = (rpk >= BAND[0]) & (rpk < BAND[1])
    n_reads, n_genes = int(reads[in_band].sum()), int(in_band.sum())
    if n_reads < MIN_BAND_READS or n_genes < MIN_BAND_GENES:
        return np.nan, n_reads, n_genes
    return float((reads[in_band] - distinct[in_band]).sum() / n_reads), n_reads, n_genes


def library_measures(gene_table, kb):
    r, d = gene_table["reads1"].to_numpy(float), gene_table["distinct1"].to_numpy(float)
    out = {}
    out["dup1_rpk2to20"], out["band_reads"], out["band_genes"] = band_fraction(r, d, kb)
    out["dup1_pooled"] = (r.sum() - d.sum()) / r.sum() if r.sum() else np.nan
    out["top20_share"] = np.sort(r)[::-1][:20].sum() / r.sum() if r.sum() else np.nan
    if gene_table["reads1_sub"].notna().all():
        rs, ds = gene_table["reads1_sub"].to_numpy(float), gene_table["distinct1_sub"].to_numpy(float)
        out["dup1_sub_rpk2to20"], out["sub_band_reads"], out["sub_band_genes"] = band_fraction(rs, ds, kb)
        out["dup1_sub_pooled"] = (rs.sum() - ds.sum()) / rs.sum()
    else:
        out.update(dup1_sub_rpk2to20=np.nan, sub_band_reads=np.nan, sub_band_genes=np.nan, dup1_sub_pooled=np.nan)
    if gene_table["pairs"].notna().all():
        p, dp = gene_table["pairs"].to_numpy(float), gene_table["distinct_pairs"].to_numpy(float)
        out["pairs_assigned"] = int(p.sum())
        out["pair_dup_fraction"] = (p.sum() - dp.sum()) / p.sum() if p.sum() else np.nan
        out["mate_absent_assigned"] = int(gene_table["mate_absent1"].sum())
    else:
        out.update(pairs_assigned=np.nan, pair_dup_fraction=np.nan, mate_absent_assigned=np.nan)
    return out


def spearman(x, y):
    ok = x.notna() & y.notna()
    if ok.sum() < 5:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(x[ok].rank(), y[ok].rank())[0, 1]), int(ok.sum())


def prep_table(meta_dir, lib):
    rows = []
    for cohort, info in PREP.items():
        ena = glob.glob(os.path.join(meta_dir, f"ena_{cohort}_PRJ*.tsv"))
        if len(ena) != 1:
            sys.exit(f"{cohort}: expected one ENA run report in {meta_dir}")
        runs = pd.read_csv(ena[0], sep="\t", dtype=str, keep_default_na=False)
        ours = set(lib.loc[lib["cohort"] == cohort, "run"])
        runs = runs[runs["run_accession"].isin(ours)]
        instruments = runs["instrument_model"].value_counts()
        sources = [ena[0]]
        matrix = os.path.join(meta_dir, f"{cohort}_series_matrix.txt.gz")
        text = "\n".join(runs["library_construction_protocol"].unique())
        if os.path.exists(matrix):
            sources.append(matrix)
            with gzip.open(matrix, "rt", encoding="utf-8", errors="replace") as fh:
                text += "\n" + fh.read()
        bioproject = os.path.join(meta_dir, f"bioproject_{cohort}.xml")
        if os.path.exists(bioproject):
            sources.append(bioproject)
            with open(bioproject) as fh:
                text += "\n" + fh.read()
        for quote in info["quotes"]:
            if quote not in text:
                sys.exit(f"{cohort}: quote not found in its sources: {quote!r}")
        rows.append({"cohort": cohort, "runs_matched_in_ena": len(runs), "libraries": len(ours),
                     "instrument": "; ".join(f"{k} ({v})" for k, v in instruments.items()),
                     "layout_ena": "; ".join(runs["library_layout"].unique()),
                     "kit": info["kit"], "selection": info["selection"], "capture": info["capture"],
                     "ffpe": info["ffpe"], "tissue": info["tissue"],
                     "evidence_quotes": " | ".join(info["quotes"]), "sources": " ; ".join(sources)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--libraries", required=True)
    ap.add_argument("--step3", required=True, help="restricted step-3 dir (dup_per_library.tsv, markdup_stats/)")
    ap.add_argument("--step3-cohort", required=True, help="step-3 dup_per_cohort.tsv (hashed only)")
    ap.add_argument("--gene-dir", required=True, help="restricted 07 output dir (per_library/, per_gene/)")
    ap.add_argument("--tag-genes", required=True)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--prep-dir", required=True, help="ENA run reports and GEO series matrices")
    ap.add_argument("--restricted-out", required=True)
    ap.add_argument("--exec-out", required=True)
    args = ap.parse_args()

    os.makedirs(args.restricted_out, exist_ok=True)
    os.makedirs(args.exec_out, exist_ok=True)
    out = {name: os.path.join(args.restricted_out, name) for name in ("dup_per_library_v2.tsv", "dup_per_library_gene.tsv.gz")}
    out |= {name: os.path.join(args.exec_out, name) for name in (
        "dup_per_cohort_role.tsv", "dup_spearman_within_cohort.tsv", "dup_cohort_prep.tsv", "dup_columns.tsv",
        "dup_missing_libraries_v2.tsv", "sha256sums.txt")}
    for path in out.values():
        if os.path.exists(path):
            sys.exit(f"refusing to overwrite {path}")

    lib = pd.read_csv(args.libraries, sep="\t", dtype=str, keep_default_na=False)
    tag = pd.read_csv(args.tag_genes, sep="\t", dtype={"gene_key": str})
    kb_by_gene = dict(zip(tag["gene_key"], tag["exon_bp"] / 1000))
    step3 = pd.read_csv(os.path.join(args.step3, "dup_per_library.tsv"), sep="\t", keep_default_na=False,
                        na_values=["NA"], dtype={"individual_id": str, "pg_programs": str})
    xw = pd.read_csv(args.crosswalk, sep="\t", usecols=["run", "run_role", "unit_library"], dtype=str,
                     keep_default_na=False)

    done = {os.path.basename(p)[:-4] for p in glob.glob(os.path.join(args.gene_dir, "per_library", "*.tsv"))}
    missing = lib[~lib["run"].isin(done)]
    rows, gene_frames = [], []
    for run in sorted(done & set(lib["run"])):
        row = pd.read_csv(os.path.join(args.gene_dir, "per_library", f"{run}.tsv"), sep="\t", dtype=str,
                          keep_default_na=False).iloc[0].to_dict()
        genes = pd.read_csv(os.path.join(args.gene_dir, "per_gene", f"{run}.tsv.gz"), sep="\t", dtype={"gene_key": str},
                            na_values=["NA"], keep_default_na=False)
        kb = genes["gene_key"].map(kb_by_gene).to_numpy(float)
        if np.isnan(kb).any():
            sys.exit(f"{run}: gene without exon_bp")
        row.update(library_measures(genes, kb))
        if row["layout"] == "paired":
            md = read_markdup(os.path.join(args.step3, "markdup_stats", f"{run}.txt"))
            row.update(markdup_paired=md["PAIRED"], markdup_dup_pair=md["DUPLICATE PAIR"],
                       markdup_single=md["SINGLE"], markdup_dup_single=md["DUPLICATE SINGLE"],
                       genebody_pair_dup_fraction=md["DUPLICATE PAIR"] / md["PAIRED"] if md["PAIRED"] else np.nan)
        rows.append(row)
        genes.insert(0, "run", run)
        gene_frames.append(genes)
    table = pd.DataFrame(rows)
    for c in ("one_end_assigned", "span_one_end_reads", "span_one_end_distinct", "star_dedup_option_in_header"):
        table[c] = table[c].astype(int)
    table["strand_antisense_fraction"] = table["strand_antisense_fraction"].astype(float)
    table["span_one_end_dup_fraction"] = 1 - table["span_one_end_distinct"] / table["span_one_end_reads"]
    keep3 = ["run", "individual_id", "exon_reads", "exon_dup_reads", "exon_dup_fraction", "exon_mate_absent_reads",
             "estimated_library_size", "genebody_lw_distinct_reads", "star_unique_reads"]
    table = table.merge(step3[keep3], on="run", how="left", validate="1:1")
    table = table.merge(xw, on="run", how="left", validate="1:1")
    table["run_role"] = table["run_role"].replace("", np.nan).fillna("not_in_crosswalk")
    table["unit_library"] = table["unit_library"].eq("True")
    if (table["layout"] != step3.set_index("run").loc[table["run"], "layout"].to_numpy()).any():
        sys.exit("layout disagrees with step 3")
    # single-end check: the span rule over one-end reads must see exactly step 3's exon reads
    single = table["layout"] == "single"
    span_mismatch = int((table.loc[single, "span_one_end_reads"] != table.loc[single, "exon_reads"]).sum())

    columns = ["cohort", "run", "individual_id", "layout", "run_role", "unit_library", "strand_call",
               "strand_antisense_fraction", "pg_programs_noPG", "star_dedup_option_in_header",
               "one_end_assigned", "one_end_multi_gene", "top20_share",
               "dup1_sub_rpk2to20", "sub_band_reads", "sub_band_genes",
               "dup1_rpk2to20", "band_reads", "band_genes", "dup1_pooled", "dup1_sub_pooled",
               "pairs_assigned", "mate_absent_assigned", "pair_dup_fraction",
               "markdup_paired", "markdup_dup_pair", "markdup_single", "markdup_dup_single", "genebody_pair_dup_fraction",
               "exon_reads", "exon_dup_reads", "exon_dup_fraction", "exon_mate_absent_reads",
               "span_one_end_reads", "span_one_end_dup_fraction",
               "estimated_library_size", "genebody_lw_distinct_reads", "star_unique_reads"]
    table = table[columns].sort_values(["cohort", "run"])
    int_cols = ["one_end_multi_gene", "sub_band_reads", "sub_band_genes", "band_reads", "band_genes", "pairs_assigned",
                "mate_absent_assigned", "markdup_paired", "markdup_dup_pair", "markdup_single", "markdup_dup_single",
                "estimated_library_size", "genebody_lw_distinct_reads", "star_unique_reads", "exon_reads",
                "exon_dup_reads", "exon_mate_absent_reads"]
    for c in int_cols:
        table[c] = pd.to_numeric(table[c]).astype("Int64")
    table.to_csv(out["dup_per_library_v2.tsv"], sep="\t", index=False, float_format="%.6g", na_rep="NA")
    pd.concat(gene_frames).to_csv(out["dup_per_library_gene.tsv.gz"], sep="\t", index=False, na_rep="NA",
                                  compression={"method": "gzip", "mtime": 0})

    # per cohort x run_role over unit libraries; pooled rows split by layout
    units = table[table["unit_library"]]

    def summarise(g):
        s = {"n_unit_libraries": len(g), "layouts": ",".join(sorted(g["layout"].unique())),
             "strand_calls": ",".join(f"{k}:{v}" for k, v in g["strand_call"].value_counts().items()),
             "one_end_assigned_median": g["one_end_assigned"].median(), "top20_share_median": g["top20_share"].median()}
        for m in MEASURES:
            v = g[m].astype(float)
            s[f"{m}_n"] = int(v.notna().sum())
            s[f"{m}_median"], s[f"{m}_q25"], s[f"{m}_q75"] = v.median(), v.quantile(0.25), v.quantile(0.75)
        return pd.Series(s)

    per = units.groupby(["cohort", "run_role"]).apply(summarise, include_groups=False).reset_index()
    # group on a copy of layout: include_groups=False hides the grouping columns from summarise
    pooled = (units.assign(layout_group=units["layout"]).groupby(["layout_group", "run_role"])
              .apply(summarise, include_groups=False).reset_index())
    pooled["cohort"] = "all_" + pooled.pop("layout_group")
    per = pd.concat([per, pooled[per.columns]], ignore_index=True)
    per.to_csv(out["dup_per_cohort_role.tsv"], sep="\t", index=False, float_format="%.6g", na_rep="NA")

    # within-cohort Spearman (all finished libraries) against depth and top-20 share
    sp = []
    for cohort, g in table.groupby("cohort"):
        for m in MEASURES + ["top20_share"]:
            for against in ("exon_reads", "one_end_assigned", "top20_share"):
                if m == against:
                    continue
                rho, n = spearman(g[m].astype(float), g[against].astype(float))
                sp.append({"cohort": cohort, "measure": m, "against": against, "spearman": rho, "n": n})
    pd.DataFrame(sp).to_csv(out["dup_spearman_within_cohort.tsv"], sep="\t", index=False, float_format="%.4f",
                            na_rep="NA")

    prep = prep_table(args.prep_dir, lib)
    strand = table.groupby("cohort")["strand_call"].agg(lambda s: ",".join(f"{k}:{v}" for k, v in s.value_counts().items()))
    anti = table.groupby("cohort")["strand_antisense_fraction"].agg(["min", "median", "max"])
    prep["strand_call_from_reads"] = prep["cohort"].map(strand)
    for k in ("min", "median", "max"):
        prep[f"antisense_fraction_{k}"] = prep["cohort"].map(anti[k])
    roles = xw.merge(lib[["run", "cohort"]], on="run").groupby("cohort")["run_role"].agg(
        lambda s: ",".join(f"{k}:{v}" for k, v in s.value_counts().items()))
    prep["run_roles"] = prep["cohort"].map(roles)
    prep["note"] = np.where(prep["cohort"] == "GSE193066", "all runs excluded in frozen_crosswalk; does not enter the model", "")
    prep.to_csv(out["dup_cohort_prep.tsv"], sep="\t", index=False, float_format="%.4f", na_rep="NA")

    dictionary = [
        ("dup1_sub_rpk2to20", "USE: duplication covariate for B-QC/B-MODEL. One-end rule, 50,000-read subsample, genes at 2-20 RPK."),
        ("dup1_rpk2to20", "secondary: same band at full depth; rises with depth."),
        ("dup1_pooled / dup1_sub_pooled", "descriptive: all assigned genes; weighted to the most-read genes."),
        ("pair_dup_fraction", "paired-end only, within paired-end: pair rule over assigned read-1s with the mate present."),
        ("genebody_pair_dup_fraction", "paired-end only: markdup DUPLICATE PAIR / PAIRED, gene-body subset, reads."),
        ("exon_dup_fraction", "step 3; pair rule (paired-end) vs one-end rule (single-end); within one layout only."),
        ("span_one_end_dup_fraction", "check only: one-end rule over read-1s whose span overlaps tag exons, either strand."),
        ("top20_share", "share of assigned one-end reads in the 20 most-read tag genes of the library."),
        ("*_reads, one_end_assigned, band_reads", "reads (read 1 only for paired-end in one-end columns)."),
        ("pairs_assigned, estimated_library_size", "pairs."),
        ("genebody_lw_distinct_reads", "reads (two per pair for paired-end)."),
        ("star_unique_reads", "STAR Log.final.out; pairs for paired-end."),
        ("all duplicate fractions", "include optical/patterned-flowcell duplicates: SRA read names carry no tile or"
                                    " coordinates, so they cannot be separated from PCR duplicates."),
    ]
    pd.DataFrame(dictionary, columns=["column", "meaning"]).to_csv(out["dup_columns.tsv"], sep="\t", index=False)
    missing[["task", "cohort", "run"]].to_csv(out["dup_missing_libraries_v2.tsv"], sep="\t", index=False)

    hashed = [out[k] for k in out if k != "sha256sums.txt"] + [os.path.join(args.step3, "dup_per_library.tsv"),
                                                               args.step3_cohort]
    with open(out["sha256sums.txt"], "w") as fh:
        for path in hashed:
            fh.write(f"{sha256(path)}  {path}\n")

    print(f"finished {len(table)}; missing {len(missing)}; single-end span/exon_reads mismatches {span_mismatch}")
    print(table.groupby(["cohort", "strand_call"]).size().to_string())
    print(per[["cohort", "run_role", "n_unit_libraries", "dup1_sub_rpk2to20_median", "dup1_rpk2to20_median",
               "exon_dup_fraction_median"]].to_string())


if __name__ == "__main__":
    main()
