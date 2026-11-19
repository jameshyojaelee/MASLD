"""Build results/susie_coloc_all_gwas_refit_20261118.csv: the COLOC r2 master with the 90 eQTL
duplicate-ID PSD refit genes swapped in (EQTL-DEDUP-PSD-20261118).

Rules
- Rows of the r2 master are copied byte-for-byte except for (gwas_name, gene) pairs
  where the unversioned Ensembl ID is one of the 89 refit genes AND the study is one of
  the 35 refit studies (= the 35 Tier-1/2 MAIN strata). Those rows are replaced IN PLACE
  (same file position, so "first max in file order" tie rules see the same order) by the
  refit row for the same (study, gene).
- Columns stay identical to r2. Columns the refit output does not carry
  (ancestry, ld_panel, ld_panel_n, stratum_*, locus_lambda_s, lambda_s, ld_reliability)
  are taken from the r2 row of the same (study, gene); eligibility_source is set to
  "eqtl_dedup_psd_20261118_refit".
- MPRIP (ENSG00000133030) is untestable by prespecification (PRESPEC_original.md, "MPRIP"):
  in the 35 studies its SuSiE fields (PP.H3.susie, PP.H4.susie, n_cs_pairs,
  n_cs_pairs_eligible) are blanked and method becomes
  "susie_untestable_unresolved_source_event"; its ABF fields are kept (not refit).
- Rows for refit genes in the 15 non-MAIN strata are not refit and stay as r2.
"""
import csv
import glob
import json
import os

import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RUN = f"{ROOT}/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118"
R2DIR = f"{ROOT}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
R2 = f"{R2DIR}/susie_coloc_all_gwas.csv"
OUT = f"{RUN}/results/susie_coloc_all_gwas_refit_20261118.csv"
OUTJ = f"{RUN}/results/r2_summary_refit_20261118.json"
OUTC = f"{RUN}/results/step1_counts_20261118.json"
MPRIP = "ENSG00000133030"
UNTEST = "susie_untestable_insufficient_shared_posterior"
assert not os.path.exists(OUT), f"{OUT} already exists; refusing to overwrite"

cmp = pd.read_csv(f"{RUN}/results/eqtl_dedup_gene_comparison.tsv", sep="\t")
assert len(cmp) == 90, f"expected 90 genes in eqtl_dedup_gene_comparison.tsv, got {len(cmp)}"
refit_ens = set(
    cmp.loc[cmp.in_refit.astype(str).str.strip().str.lower() == "true", "ensembl"]
    .str.split(".")
    .str[0]
)
assert len(refit_ens) == 89 and MPRIP not in refit_ens, (len(refit_ens), MPRIP in refit_ens)
studies = sorted(os.path.basename(p) for p in glob.glob(f"{RUN}/coloc/*") if os.path.isdir(p))
tier = pd.read_csv(f"{ROOT}/GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
MAIN = set(tier.loc[tier.placement == "main", "study_name"])
assert set(studies) == MAIN and len(MAIN) == 35, (len(studies), set(studies) ^ MAIN)

new = pd.concat(
    [
        pd.read_csv(f, dtype=str, keep_default_na=False)
        for s in studies
        for f in sorted(glob.glob(f"{RUN}/coloc/{s}/susie_coloc_chr*.csv"))
    ],
    ignore_index=True,
)
new["ens0"] = new.ensembl.str.split(".").str[0]
new = new[new.ens0.isin(refit_ens)].copy()
assert not new.duplicated(["gwas_name", "ens0"]).any(), "duplicate (gwas_name, ens0) in refit coloc"
assert len(new) == 89 * 35, f"expected {89 * 35} refit rows, got {len(new)}"
newmap = {(r["gwas_name"], r["ens0"]): r for r in new.to_dict("records")}

with open(R2, newline="") as fh:
    header = next(csv.reader(fh))
refit_cols = [c for c in new.columns if c != "ens0"]
missing_cols = [c for c in header if c not in refit_cols]
print("columns filled from the r2 row:", missing_cols)
assert set(refit_cols) <= set(header), set(refit_cols) - set(header)

n_in = n_out = n_rep = n_mprip = 0
seen = set()
lam_diff = 0
with open(R2, newline="") as fin, open(OUT, "w", newline="") as fout:
    reader = csv.reader(fin)
    w = csv.writer(fout, lineterminator="\n")
    hdr = next(reader)
    w.writerow(hdr)
    ix = {c: i for i, c in enumerate(hdr)}
    for row in reader:
        n_in += 1
        g, e = row[ix["gwas_name"]], row[ix["ensembl"]].split(".")[0]
        if g in MAIN and e in refit_ens:
            key = (g, e)
            assert key not in seen, key
            seen.add(key)
            nr = newmap[key]
            assert nr["gene"] == row[ix["gene"]] or nr["gene"] == "", (key, nr["gene"], row[ix["gene"]])
            out = list(row)
            for c in refit_cols:
                out[ix[c]] = nr[c]
            out[ix["eligibility_source"]] = "eqtl_dedup_psd_20261118_refit"
            if (
                row[ix["locus_lambda_s"]]
                and nr["lambda_s_locus"]
                and abs(float(row[ix["locus_lambda_s"]]) - float(nr["lambda_s_locus"])) > 1e-9
            ):
                lam_diff += 1
            row = out
            n_rep += 1
        elif g in MAIN and e == MPRIP:
            for c in ["PP.H3.susie", "PP.H4.susie", "n_cs_pairs", "n_cs_pairs_eligible"]:
                row[ix[c]] = ""
            row[ix["method"]] = "susie_untestable_unresolved_source_event"
            n_mprip += 1
        w.writerow(row)
        n_out += 1

assert seen == set(newmap), (len(seen), len(newmap))
assert n_rep == 89 * 35, f"expected {89 * 35} replaced rows, got {n_rep}"
assert n_mprip == 35, f"expected 35 MPRIP rows blanked, got {n_mprip}"
print(
    f"rows in {n_in} out {n_out}; replaced in place {n_rep}; MPRIP rows blanked {n_mprip}; "
    f"refit rows whose lambda_s_locus differs from r2 locus_lambda_s (r2 value kept in that column): {lam_diff}"
)


# ---------------- counts under the Panel A producer rule (gene symbol, 35 MAIN strata) ----
def counts(path):
    m = pd.read_csv(
        path,
        usecols=["gwas_name", "gene", "ensembl", "PP.H4.abf", "PP.H4.susie", "method"],
        dtype={"gene": str, "method": str},
        low_memory=False,
    )
    m = m[m.gwas_name.isin(MAIN)]
    m["s4"] = pd.to_numeric(m["PP.H4.susie"], errors="coerce")
    m["a4"] = pd.to_numeric(m["PP.H4.abf"], errors="coerce")
    m = m[m.gene.notna() & m.gene.str.strip().ne("")]
    gg = m.groupby("gene").agg(
        s=("s4", "max"),
        a=("a4", "max"),
        unt=("method", lambda x: (x == UNTEST).any()),
        fin=("s4", lambda x: x.notna().any()),
    )
    gg["s"] = gg.s.fillna(0.0)
    gg["a"] = gg.a.fillna(0.0)
    su = gg[gg.s > 0.5]
    ab = gg[(gg.a > 0.5) & ~(gg.s > 0.5)]

    def bands(v):
        return [
            int((v > 0.9).sum()),
            int(((v >= 0.8) & (v <= 0.9)).sum()),
            int(((v > 0.5) & (v < 0.8)).sum()),
        ]

    res = {
        "multi_signal": len(su),
        "single_signal_only": len(ab),
        "union": len(su) + len(ab),
        "multi_bands_gt0.9_0.8-0.9_0.5-0.8": bands(su.s),
        "single_bands_gt0.9_0.8-0.9_0.5-0.8": bands(ab.a),
        "single_only_untestable": sorted(ab[ab.unt].index.tolist()),
        "abf_gt_0.5": int((gg.a > 0.5).sum()),
        "state": {
            "positive": int((gg.s > 0.5).sum()),
            "untestable": int((~(gg.s > 0.5) & gg.unt).sum()),
            "tested_le_0.5": int((gg.fin & ~(gg.s > 0.5) & ~gg.unt).sum()),
            "no_susie_result": int((~gg.fin & ~gg.unt).sum()),
        },
    }
    return res, set(su.index), set(ab.index)


old, old_su, old_ab = counts(R2)
newc, new_su, new_ab = counts(OUT)
lost = sorted(old_su - new_su)
gained = sorted(new_su - old_su)
union_lost = sorted((old_su | old_ab) - (new_su | new_ab))
union_gained = sorted((new_su | new_ab) - (old_su | old_ab))
moved_multi_to_single = sorted(old_su & new_ab)

print("r2   ", json.dumps(old))
print("refit", json.dumps(newc))
print("multi-signal lost:", lost, "gained:", gained)
print("union lost:", union_lost, "gained:", union_gained)
print("moved multi -> single-only:", moved_multi_to_single)

assert old["multi_signal"] == 436 and old["union"] == 996, old
assert newc["multi_signal"] == 430, f"expected 430 primary multi-signal SuSiE genes, got {newc['multi_signal']}"

summ = json.load(open(f"{R2DIR}/r2_summary.json"))
summ["refit_note"] = (
    "r2 release with the 89 eQTL duplicate-ID PSD refit genes replaced in the 35 MAIN strata "
    "(EQTL-DEDUP-PSD-20261118) and MPRIP SuSiE-untestable; built by "
    f"{RUN}/code/build_refit_master_20261118.py. Only fig2, fig2_susie_untestable, tier12 and "
    "tier12_susie_state are recomputed; other keys are r2 values."
)
summ["fig2"] = {
    "multi_signal": newc["multi_signal"],
    "single_signal_only": newc["single_signal_only"],
    "multi_or_single_signal_coloc_gene_union": newc["union"],
}
summ["fig2_susie_untestable"]["drawn_as_single_signal_only"] = len(newc["single_only_untestable"])
summ["fig2_susie_untestable"]["drawn_as_single_signal_only_genes"] = newc["single_only_untestable"]
summ["tier12"] = {
    "tier12_susie_gt_0.5": newc["multi_signal"],
    "tier12_abf_gt_0.5": newc["abf_gt_0.5"],
    "tier12_union_gt_0.5": newc["union"],
}
summ["tier12_susie_state"] = newc["state"]

with open(OUTJ, "w") as fj:
    json.dump(summ, fj, indent=2)

step1_out = {
    "r2": old,
    "refit": newc,
    "in_refit_genes": len(refit_ens),
    "rows_in": n_in,
    "rows_out": n_out,
    "rows_replaced_in_place": n_rep,
    "mprip_rows_blanked": n_mprip,
    "lambda_s_locus_diffs_kept_r2": lam_diff,
    "multi_signal_lost": lost,
    "multi_signal_gained": gained,
    "union_lost": union_lost,
    "union_gained": union_gained,
    "moved_multi_to_single_only": moved_multi_to_single,
}
with open(OUTC, "w") as fc:
    json.dump(step1_out, fc, indent=2)

print("wrote", OUT, OUTJ, OUTC)
