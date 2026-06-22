"""
MASLD Multi-Evidence Convergence Dashboard
===========================================
Standalone Streamlit app (does NOT touch the archived DEG explorer or the
Next.js portal). Surfaces, per gene and globally:
  - the convergence ranking + per-source evidence (S1-S8 log Bayes factors)
  - druggability (genome-wide DGIdb, fixed 2026-06-10) + DepMap essentiality
  - the Tier-4A evidence-ORTHOGONALITY reframe (modality correlations + complementarity)
  - the Tier-4B/3D tractability summary

Run:  micromamba run -n rnaseq streamlit run streamlit_convergence/app.py
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- paths
APP_DIR = Path(__file__).resolve().parent
ROOT = APP_DIR
for p in [APP_DIR, *APP_DIR.parents]:
    if (p / "RNA-seq").exists() and (p / "GWAS").exists():
        ROOT = p; break
ME    = ROOT / "RNA-seq/results/multi_evidence"
ATLAS = ME / "multi_evidence_atlas.csv"
CONV  = ME / "convergence_evidence.csv"
DGIDB = ROOT / "data/dgidb/dgidb_druggable_genes.csv"

CTRL = "#9E9E9E"
SOURCE_LABELS = {
    "log_BF_S1": "S1 Human bulk RNA-seq",
    "log_BF_S2_coloc": "S2 Genetic causal (COLOC)",
    "log_BF_S3": "S3 Genetic causal (GWAS+eQTL)",
    "log_BF_S4": "S4 Essentiality (DepMap)",
    "log_BF_S5": "S5 Epigenomic (ATAC/SCENIC+)",
    "log_BF_S6": "S6 Spatial",
    "log_BF_S7": "S7 Single-cell",
    "log_BF_S8": "S8 Cross-species",
}

st.set_page_config(page_title="MASLD Convergence Atlas", layout="wide", page_icon="🧬")

# ---------------------------------------------------------------- loaders
@st.cache_data(show_spinner=False)
def load_conv() -> pd.DataFrame:
    df = pd.read_csv(CONV)
    df["human_symbol"] = df["human_symbol"].astype(str)
    return df

@st.cache_data(show_spinner=False)
def load_atlas() -> pd.DataFrame:
    wish = ["human_symbol", "gene_biotype", "bulk_logFC", "bulk_padj", "is_conserved",
            "essentiality_chronos", "is_essential", "dgidb_druggable",
            "coloc_best_susie_pp4_polyfun", "coloc_best_pp4_polyfun", "sceqtl_coloc_best_pp4"]
    have = pd.read_csv(ATLAS, nrows=0).columns.tolist()
    use = [c for c in wish if c in have]
    df = pd.read_csv(ATLAS, usecols=use, low_memory=False)
    df["human_symbol"] = df["human_symbol"].astype(str)
    return df

@st.cache_data(show_spinner=False)
def load_dgidb() -> pd.DataFrame:
    if DGIDB.exists():
        d = pd.read_csv(DGIDB)
        d["gene_name"] = d["gene_name"].astype(str).str.upper()
        return d
    return pd.DataFrame(columns=["gene_name", "n_drugs", "n_approved_drugs"])

@st.cache_data(show_spinner=False)
def load_modality_corr():
    f = ME / "convergence_evidence_modality_correlations.csv"
    return pd.read_csv(f) if f.exists() else None

@st.cache_data(show_spinner=False)
def load_complementarity():
    f = ME / "complementarity_pvalues.csv"
    return pd.read_csv(f) if f.exists() else None

conv = load_conv()
atlas = load_atlas()
dgidb = load_dgidb()
dgidb_map = dict(zip(dgidb["gene_name"], zip(dgidb.get("n_drugs", 0), dgidb.get("n_approved_drugs", 0))))
atlas_map = atlas.set_index(atlas["human_symbol"].str.upper())

# ---------------------------------------------------------------- header
st.title("🧬 MASLD Multi-Evidence Convergence Atlas")
st.caption(
    "Per-gene evidence cards across 8 orthogonal sources, druggability (genome-wide DGIdb), "
    "and the evidence-orthogonality view. Convergence is rare *because* sources are near-independent — "
    "see the Orthogonality tab.")

tab_gene, tab_top, tab_orth, tab_tract = st.tabs(
    ["🔎 Gene evidence card", "🏆 Top convergent targets", "🔗 Evidence orthogonality", "💊 Tractability"])

# ================================================================ TAB 1: gene card
with tab_gene:
    genes = sorted(conv.loc[conv.get("excluded_from_ranking", False) != True, "human_symbol"].unique()) \
        if "excluded_from_ranking" in conv.columns else sorted(conv["human_symbol"].unique())
    default = "THRB" if "THRB" in genes else (genes[0] if genes else "")
    g = st.selectbox("Select a gene", genes, index=genes.index(default) if default in genes else 0)
    row = conv[conv["human_symbol"] == g]
    if row.empty:
        st.warning("Gene not found in convergence table.")
    else:
        r = row.iloc[0]
        a = atlas_map.loc[g.upper()] if g.upper() in atlas_map.index else None
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Convergence rank", int(r["convergence_rank"]) if pd.notna(r.get("convergence_rank")) else "—")
        c2.metric("Convergence score", f"{r.get('convergence_score', float('nan')):.2f}")
        c3.metric("Modalities active", f"{int(r.get('n_modalities_active', 0))} / 8")
        c4.metric("Tier", str(r.get("tier", "—")))
        c5.metric("Concordance", str(r.get("concordance_state", "—")))

        # per-source evidence bars
        st.subheader("Per-source evidence (log Bayes factor)")
        vals = {SOURCE_LABELS[k]: r.get(k, np.nan) for k in SOURCE_LABELS if k in conv.columns}
        ev = pd.Series(vals).fillna(0.0).sort_values()
        fig, ax = plt.subplots(figsize=(6, 3))
        colors = ["#C2185B" if v > 0 else CTRL for v in ev.values]
        ax.barh(ev.index, ev.values, color=colors)
        ax.axvline(0, color=CTRL, lw=0.5); ax.set_xlabel("log BF (>0 = active)")
        ax.spines[["top", "right"]].set_visible(False)
        st.pyplot(fig, use_container_width=True)

        # druggability + essentiality + effect
        st.subheader("Target annotation")
        d1, d2, d3, d4 = st.columns(4)
        nd, na = dgidb_map.get(g.upper(), (0, 0))
        druggable = g.upper() in dgidb_map
        d1.metric("DGIdb druggable", "Yes" if druggable else "No",
                  f"{int(nd)} drugs / {int(na)} approved" if druggable else "RNA-targeting candidate")
        if a is not None:
            ess = str(a.get("is_essential", "")).upper() in ("TRUE", "1")
            d2.metric("DepMap essential", "Yes" if ess else "No",
                      f"chronos {a.get('essentiality_chronos', float('nan')):.2f}" if pd.notna(a.get("essentiality_chronos")) else "")
            d3.metric("Biotype", str(a.get("gene_biotype", "—")))
            pp4 = a.get("coloc_best_susie_pp4_polyfun", a.get("coloc_best_pp4_polyfun", np.nan))
            d4.metric("COLOC PP.H4 (best)", f"{pp4:.2f}" if pd.notna(pp4) else "—")
            lfc, padj = a.get("bulk_logFC", np.nan), a.get("bulk_padj", np.nan)
            if pd.notna(lfc):
                st.caption(f"Bulk disease-vs-control: logFC = {lfc:+.2f}, padj = {padj:.1e} | "
                           f"conserved (cross-species): {a.get('is_conserved','—')}")
        if not druggable:
            st.info(f"**{g}** has no small-molecule drug in DGIdb → a candidate for RNA-targeting (Cas13) "
                    "knockdown, the rationale for the perturbation library.")

# ================================================================ TAB 2: top targets
with tab_top:
    st.subheader("Top convergent targets (filterable, downloadable)")
    cols = ["human_symbol", "convergence_rank", "convergence_score", "n_modalities_active",
            "tier", "concordance_state"]
    cols = [c for c in cols if c in conv.columns]
    t = conv[cols].copy()
    # join druggability + biotype
    t["dgidb_druggable"] = t["human_symbol"].str.upper().isin(dgidb_map)
    bt = atlas.set_index(atlas["human_symbol"].str.upper())["gene_biotype"] if "gene_biotype" in atlas.columns else None
    if bt is not None:
        t["gene_biotype"] = t["human_symbol"].str.upper().map(bt)
    f1, f2, f3 = st.columns(3)
    tiers = ["(all)"] + sorted(t["tier"].dropna().unique().tolist()) if "tier" in t.columns else ["(all)"]
    sel_tier = f1.selectbox("Tier", tiers)
    min_mod = f2.slider("Min modalities active", 0, 8, 3)
    only_drug = f3.selectbox("Druggability", ["(all)", "DGIdb-druggable only", "Undruggable only (RNA-target candidates)"])
    q = t[t["n_modalities_active"] >= min_mod] if "n_modalities_active" in t.columns else t
    if sel_tier != "(all)" and "tier" in q.columns:
        q = q[q["tier"] == sel_tier]
    if only_drug == "DGIdb-druggable only":
        q = q[q["dgidb_druggable"]]
    elif only_drug.startswith("Undruggable"):
        q = q[~q["dgidb_druggable"]]
    q = q.sort_values("convergence_rank") if "convergence_rank" in q.columns else q
    st.caption(f"{len(q):,} genes match")
    st.dataframe(q.head(500), use_container_width=True, height=520)
    st.download_button("Download filtered table (CSV)", q.to_csv(index=False).encode(),
                       "convergent_targets_filtered.csv", "text/csv")

# ================================================================ TAB 3: orthogonality
with tab_orth:
    st.subheader("Evidence sources are near-orthogonal")
    st.caption("All pairwise Spearman |ρ| < 0.21 — sources capture non-redundant biology, so high "
               "multi-modal convergence is *expected* to be rare. That is the finding, not a weakness.")
    mc = load_modality_corr()
    if mc is not None:
        piv = mc.pivot(index="Mod_A", columns="Mod_B", values="Spearman_rho")
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-0.25, vmax=0.25)
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index, fontsize=7)
        for i in range(len(piv.index)):
            for j in range(len(piv.columns)):
                v = piv.values[i, j]
                if pd.notna(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
        fig.colorbar(im, ax=ax, shrink=0.7, label="Spearman ρ")
        st.pyplot(fig, use_container_width=False)
    cp = load_complementarity()
    if cp is not None:
        st.subheader("Each source adds independent validation signal (complementarity)")
        sub = cp[cp["validation"].isin(["V1_DGIdb", "V3_Conserved"])].copy()
        sub["sig"] = np.where(sub["perm_pvalue"] < 0.05, "perm p<0.05", "n.s.")
        st.dataframe(sub[["source", "validation", "improvement", "perm_pvalue", "sig"]]
                     .sort_values(["validation", "improvement"], ascending=[True, False]),
                     use_container_width=True, height=320)

# ================================================================ TAB 4: tractability
with tab_tract:
    st.subheader("Druggability across convergence tiers (genome-wide DGIdb)")
    st.caption("DGIdb druggability fixed 2026-06-10 (was a broken 300-gene column → now 5,012 genes). "
               "Convergent targets are enriched for tractable biology, but a majority still need "
               "non-small-molecule approaches — the rationale for an RNA-targeting (Cas13) library.")
    rows = []
    nmod = conv["n_modalities_active"] if "n_modalities_active" in conv.columns else pd.Series(0, index=conv.index)
    for label, mask in [("All atlas genes", pd.Series(True, index=conv.index)),
                        ("Convergent ≥3/8", nmod >= 3),
                        ("Convergent ≥4/8", nmod >= 4)]:
        sub = conv[mask]
        n = len(sub)
        dr = sub["human_symbol"].str.upper().isin(dgidb_map).sum()
        rows.append({"Set": label, "n": n, "DGIdb-druggable": int(dr),
                     "% druggable": f"{100*dr/max(1,n):.0f}%", "Undruggable": int(n-dr)})
    st.table(pd.DataFrame(rows))
    tf = ME / "tractability_4B.csv"
    if tf.exists():
        st.caption("Full tractability table (incl. Cas13 library):")
        st.dataframe(pd.read_csv(tf), use_container_width=True)

st.caption(f"Data: {CONV.name}, {ATLAS.name}, {DGIDB.name} · ROOT={ROOT}")
