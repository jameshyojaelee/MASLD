#!/usr/bin/env python
"""
GLP-1RA reversal of MASLD atlas targets (Analysis B) -- META-LEVEL, corrected.

Adversarial review (RevB/RevC) established:
  * The reversal is REAL, not a shared-reference regression-to-mean artifact
    (sample-label permutation floor ~50%/r~-0.1 vs observed ~88%/r~-0.8, p<=0.015).
  * The per-DATASET reversal rule INFLATED "refractory" (26% of them reverse at meta).
    => use a META-LEVEL reversal call as primary.
  * The reversed/refractory split is a significance boundary, not biology
    => add a CONTINUOUS reversal score + an honest 4-way class.
  * Human-serum vs mouse-liver gene-level concordance is CHANCE (kappa~0)
    => report kappa honestly; do NOT headline it; serum is a SEPARATE clinical anchor (218f).

Reversal is defined WITHIN mouse (semaglutide meta-treatment opposes mouse meta-disease),
then overlaid on our human liver targets, with a mouse-vs-human disease-direction
concordance filter (only genes the mouse models correctly are informative).

Outputs (RNA-seq/results/glp1ra/):
  glp1ra_target_reversal_partition.csv
  glp1ra_reversal_summary.txt
"""
import os
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
G = os.path.join(ROOT, "RNA-seq/results/glp1ra")
REF = os.path.join(G, "disease_reference.csv")
META = os.path.join(G, "mouse_reversal/mouse_semaglutide_reversal_meta.csv")
LONG = os.path.join(G, "mouse_reversal/mouse_semaglutide_reversal.csv")
SOMA = os.path.join(G, "human_semaglutide_signature/somascan_full_4979_proteins.csv")

DIS_FDR = 0.05        # disease significance (meta)
TRT_FDR = 0.05        # treatment significance (meta)


def cohen_kappa(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n = len(a)
    po = (a == b).mean()
    cats = set(a) | set(b)
    pe = sum((np.mean(a == c)) * (np.mean(b == c)) for c in cats)
    return (po - pe) / (1 - pe) if pe < 1 else np.nan, po, pe


# --------------------------------------------------------------------------- inputs
ref = pd.read_csv(REF).dropna(subset=["human_symbol"]).drop_duplicates("human_symbol")
ref["human_disease_dir"] = np.sign(ref["logFC"])

meta = pd.read_csv(META)
m = meta.rename(columns={"disease_meta_lfc": "mdis", "treatment_meta_lfc": "mtrt",
                         "disease_meta_padj": "mdis_padj", "treatment_meta_padj": "mtrt_padj"}).copy()
m["mouse_disease_dir"] = np.sign(m["mdis"])
m["mouse_disease_sig"] = m["mdis_padj"] < DIS_FDR
m["treat_opposes"] = (np.sign(m["mtrt"]) != np.sign(m["mdis"])) & (m["mtrt"] != 0) & (m["mdis"] != 0)
m["treat_sig"] = m["mtrt_padj"] < TRT_FDR
# continuous reversal score: >0 means semaglutide pushes the gene back toward normal (LFC units)
m["reversal_index"] = -np.sign(m["mdis"]) * m["mtrt"]
m["fraction_corrected"] = np.clip(m["reversal_index"] / m["mdis"].abs().replace(0, np.nan), -1, 2)


def rev_class(r):
    if not r.mouse_disease_sig:
        return "untested"                       # not robustly disease-dysregulated in mouse
    if r.treat_opposes and r.treat_sig:
        return "reversed"                        # significant meta reversal
    if r.treat_opposes and not r.treat_sig:
        return "reversing_ns"                    # directionally reversing, underpowered
    return "refractory"                          # moving WITH disease = genuine non-reverser


m["mouse_reversal"] = m.apply(rev_class, axis=1)
m["weight_independent"] = m.get("reversal_weight_independent", False) == True  # noqa: E712

# >=2-dataset replication (high-confidence reversed) from the long per-dataset table
lg = pd.read_csv(LONG)
lg["has_trt"] = lg["treatment_lfc"].notna()
lg["ds_dis_sig"] = lg.get("disease_padj", 1) < 0.10
lg["ds_rev"] = (lg["ds_dis_sig"] & lg["has_trt"] & (lg.get("treatment_padj", 1) < 0.10) &
                (np.sign(lg["treatment_lfc"]) != np.sign(lg["disease_lfc"])) &
                (lg["disease_lfc"] != 0) & (lg["treatment_lfc"] != 0))
nrev = lg.groupby("human_symbol")["ds_rev"].sum().rename("mouse_n_datasets_reversed")
m = m.merge(nrev, on="human_symbol", how="left")
m["mouse_n_datasets_reversed"] = m["mouse_n_datasets_reversed"].fillna(0).astype(int)
m["reversed_hi_conf"] = (m["mouse_reversal"] == "reversed") & (m["mouse_n_datasets_reversed"] >= 2)

mkeep = m[["human_symbol", "mdis", "mdis_padj", "mtrt", "mtrt_padj", "mouse_disease_dir",
           "mouse_disease_sig", "reversal_index", "fraction_corrected", "mouse_reversal",
           "weight_independent", "mouse_n_datasets_reversed", "reversed_hi_conf"]]

# --------------------------------------------------------------------------- overlay
tab = ref.merge(mkeep, on="human_symbol", how="left")
tab["mouse_reversal"] = tab["mouse_reversal"].fillna("untested")
# mouse models the human disease direction?
tab["mouse_human_concordant"] = np.where(
    tab["mouse_disease_sig"].fillna(False) & tab["human_disease_dir"].notna(),
    tab["mouse_disease_dir"] == tab["human_disease_dir"], np.nan)
tab.to_csv(os.path.join(G, "glp1ra_target_reversal_partition.csv"), index=False)


# --------------------------------------------------------------------------- summary
def summarize(mask, name):
    s = tab[mask]
    n = len(s)
    vc = s["mouse_reversal"].value_counts()
    # among genes the mouse models correctly (concordant disease direction):
    conc = s[s["mouse_human_concordant"] == True]
    out = [f"\n### {name} (n={n})"]
    for k in ["reversed", "reversing_ns", "refractory", "untested"]:
        out.append(f"    {k:13s}: {int(vc.get(k,0)):5d} ({100*vc.get(k,0)/max(n,1):5.1f}%)")
    # CONTINUOUS headline (effect-size based, threshold-free)
    ds = s[s["mouse_disease_sig"] == True]
    out.append(f"    [HEADLINE, continuous] median fraction-corrected (disease-sig, n={len(ds)}): "
               f"{ds['fraction_corrected'].median():.2f}  "
               f"(>0.5 corrected: {int((ds['fraction_corrected'] > 0.5).sum())}, "
               f"moving-with-disease <0: {int((ds['fraction_corrected'] < 0).sum())})")
    # transparency: strict per-dataset vs meta reversal call
    out.append(f"    [transparency] strict per-dataset reversed (>=1 ds): "
               f"{int((s['mouse_n_datasets_reversed'] >= 1).sum())}  |  meta-reversed: {int(vc.get('reversed',0))}")
    out.append(f"    reversed hi-conf (>=2 datasets): {int(s['reversed_hi_conf'].sum())}")
    out.append(f"    reversed & weight-independent (CDA-HFD, single-model): "
               f"{int((s['mouse_reversal'].eq('reversed') & s['weight_independent']).sum())}")
    out.append(f"    [among {len(conc)} mouse-concordant genes] reversed: "
               f"{int(conc['mouse_reversal'].eq('reversed').sum())} "
               f"({100*conc['mouse_reversal'].eq('reversed').mean() if len(conc) else 0:.1f}%), "
               f"genuine-refractory: {int(conc['mouse_reversal'].eq('refractory').sum())}")
    return "\n".join(out)


L = ["GLP-1RA reversal of atlas targets (META-LEVEL, corrected)", "=" * 56,
     f"disease FDR<{DIS_FDR}, treatment FDR<{TRT_FDR}; reversal defined within mouse meta, "
     "overlaid on human targets."]
L.append(summarize(ref["conv_tier1"] == True, "Convergence Tier-1"))
L.append(summarize(ref["conv_ge3_modalities"] == True, ">=3-modality convergence"))
L.append(summarize(ref["is_deg_tier1"] == True, "Tier-1 DEGs"))

# poolability: mouse-vs-human liver disease-direction concordance
dd = tab.loc[tab["is_deg_tier1"] == True, "mouse_human_concordant"].dropna()
L.append(f"\n### Poolability: mouse-vs-human liver disease-direction concordance "
         f"(Tier-1 DEGs): {dd.mean()*100:.1f}% (n={len(dd)})")

# HONEST serum-vs-mouse concordance (NOT a headline; serum is a separate anchor in 218f)
soma = pd.read_csv(SOMA).rename(columns={"gene": "human_symbol"})
soma["s_dis"] = pd.to_numeric(soma["indep_MASHvsHealthy_logFC"], errors="coerce")
soma["s_trt"] = pd.to_numeric(soma["sema_effect_est1"], errors="coerce")
soma["s_dis_sig"] = soma["indep_signif"].astype(str).str.upper().isin(["Y", "TRUE", "1"])
soma["serum_reversed"] = (soma["s_dis_sig"] & (np.sign(soma["s_trt"]) != np.sign(soma["s_dis"])))
j = tab.merge(soma[["human_symbol", "serum_reversed", "s_dis_sig"]], on="human_symbol", how="inner")
j = j[(j["mouse_disease_sig"] == True) & (j["s_dis_sig"] == True)]
if len(j) > 10:
    a = (j["mouse_reversal"] == "reversed").astype(int).values
    b = j["serum_reversed"].astype(int).values
    kappa, po, pe = cohen_kappa(a, b)
    L.append(f"\n### Human-serum vs mouse-liver reversal concordance (HONEST, not a headline)")
    L.append(f"    n={len(j)} genes disease-sig in both compartments")
    L.append(f"    observed agreement {po*100:.1f}% vs chance {pe*100:.1f}%  =>  Cohen's kappa = {kappa:.3f}")
    L.append("    => gene-level cross-compartment reversal is ~independent (near-zero information);")
    L.append("       serum proteomics is reported as a SEPARATE clinical anchor, not fused.")

summary = "\n".join(L)
with open(os.path.join(G, "glp1ra_reversal_summary.txt"), "w") as fh:
    fh.write(summary + "\n")
print(summary)
print("\n[done] wrote partition + summary")
