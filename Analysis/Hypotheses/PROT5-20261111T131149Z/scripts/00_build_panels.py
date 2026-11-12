#!/usr/bin/env python3
"""PROT5 step 0: fix the hepatocyte and negative-control protein panels from non-genetic information only.

Inputs: GWAS Catalog study listing for PMID 39789286 (metadata/gwascatalog_PMID39789286_studies.tsv,
fetched from the GWAS Catalog REST API v2 with pubmed_id=39789286), Human Protein Atlas proteinatlas.tsv
(data/external/hpa), NCBI37.3 gene coordinates (tools/magma/ref). No pQTL statistic is read.

Usage: python -I 00_build_panels.py <project_root> <run_dir>
"""
import sys

import pandas as pd

SEED = ["FBP1", "ARG1", "ADH1B", "ADH4", "GOT1", "KRT18", "FABP1", "ASGR1"]
RETRY = ["ALDOB", "GPT", "CPS1", "BHMT", "HPD"]
CTRL_TISSUES = ["skeletal muscle", "heart muscle", "bone marrow", "lymphoid tissue", "intestine", "kidney", "pancreas"]
MAX_PANEL = 25


def main(root, run):
    s = pd.read_csv(f"{run}/metadata/gwascatalog_PMID39789286_studies.tsv", sep="\t")
    s = s[s["initial_sample_size"].str.startswith("47,745")].copy()
    s["gene"] = s["trait"].str.replace(" protein levels", "", regex=False)
    h = pd.read_csv(f"{root}/data/external/hpa/proteinatlas.tsv", sep="\t", usecols=[
        "Gene", "Ensembl", "Protein class", "RNA tissue specificity", "RNA tissue specificity score",
        "RNA tissue specific nTPM", "RNA single cell type specific nCPM", "Secretome location",
        "Subcellular main location"])
    m = s.merge(h, left_on="gene", right_on="Gene", how="left")
    pc = m["Protein class"].fillna("")
    intra = (m["Secretome location"].isna() & pc.str.contains("Predicted intracellular proteins")
             & ~pc.str.contains("Predicted secreted proteins") & ~pc.str.contains("Predicted membrane proteins"))
    single = m["RNA tissue specificity"] == "Tissue enriched"
    ntpm = m["RNA tissue specific nTPM"].fillna("")
    hep_sc = m["RNA single cell type specific nCPM"].fillna("").str.contains("Hepatocytes")
    strict = single & ntpm.str.match("^liver: ") & intra & hep_sc
    m["hpa_strict_hepatocyte_rule"] = strict
    # the cis rule needs GRCh37 gene coordinates; proteins without them are skipped before ranking
    gl = pd.read_csv(f"{root}/tools/magma/ref/NCBI37.3.gene.loc", sep="\t", header=None,
                     names=["entrez", "chr", "start", "end", "strand", "symbol"], dtype={"chr": str})
    gl = gl.drop_duplicates("symbol")
    has_coord = m["gene"].isin(set(gl["symbol"]))
    rows = []
    for g in SEED:
        r = m[m["gene"] == g]
        rows.append(dict(panel="hepatocyte", gene=g, accession=r["accession"].iloc[0] if len(r) else None,
                         found=len(r) > 0, rule="seed list (task brief)",
                         meets_strict_hpa_rule=bool(r["hpa_strict_hepatocyte_rule"].iloc[0]) if len(r) else None))
    ext = m[strict & has_coord & ~m["gene"].isin(SEED)].sort_values("RNA tissue specificity score", ascending=False)
    ext = ext.head(MAX_PANEL - len(SEED))
    for _, r in ext.iterrows():
        rows.append(dict(panel="hepatocyte", gene=r["gene"], accession=r["accession"], found=True,
                         rule="HPA: liver tissue enriched, not secreted, intracellular, hepatocyte single-cell enriched",
                         meets_strict_hpa_rule=True))
    for g in RETRY:
        r = m[m["gene"] == g]
        rows.append(dict(panel="retry_not_in_panel" if len(r) == 0 else "hepatocyte", gene=g,
                         accession=r["accession"].iloc[0] if len(r) else None, found=len(r) > 0,
                         rule="retry by exact trait name '<GENE> protein levels' in the study listing",
                         meets_strict_hpa_rule=None))
    n_hep = sum(1 for x in rows if x["panel"] == "hepatocyte" and x["found"])
    ctrl = m[single & intra & ~hep_sc & has_coord & ~m["gene"].isin(SEED)].copy()
    ctrl["tissue"] = ntpm[ctrl.index].str.split(":").str[0]
    ctrl = ctrl[ctrl["tissue"].isin(CTRL_TISSUES)].sort_values("RNA tissue specificity score", ascending=False)
    ctrl = ctrl.head(n_hep)
    for _, r in ctrl.iterrows():
        rows.append(dict(panel="negative_control", gene=r["gene"], accession=r["accession"], found=True,
                         rule=f"HPA: {r['tissue']} tissue enriched, not secreted, intracellular, not hepatocyte-listed",
                         meets_strict_hpa_rule=False))
    out = pd.DataFrame(rows)
    hcols = m.set_index("gene")[["RNA tissue specificity", "RNA tissue specificity score", "RNA tissue specific nTPM",
                                 "Secretome location", "Subcellular main location", "RNA single cell type specific nCPM"]]
    out = out.merge(hcols, left_on="gene", right_index=True, how="left")
    out = out.merge(gl[["symbol", "chr", "start", "end"]], left_on="gene", right_on="symbol", how="left").drop(columns="symbol")
    out.to_csv(f"{run}/inputs/panels.tsv", sep="\t", index=False)
    print(out[["panel", "gene", "accession", "found", "meets_strict_hpa_rule", "chr", "start"]].to_string())
    assert out.loc[out["panel"] != "retry_not_in_panel", "chr"].notna().all(), "gene without GRCh37 coordinates"


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
