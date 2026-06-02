#!/usr/bin/env python3
"""
296_community_subtype_labels.py
Generate descriptive subtype labels for F0-F1 and F3-F4 Leiden communities
based on MSigDB Hallmark enrichment patterns and transition dynamics.

This version computes its OWN hypergeometric enrichment for every community
against the full 50-pathway MSigDB Hallmark panel. The upstream enrichment
CSVs (292_*) only test 4 hallmarks, so most communities are unlabeled there.

Inputs:
  RNA-seq/results/network/communities_f2/
    communities_F01.csv   (per-gene community assignments)
    communities_F34.csv
    community_transition_table.csv
  Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt

Outputs:
  community_labels.json   (JSON keyed by stage -> community_id)
  community_labels.tsv    (human-readable table)

Deploy:
  Copies JSON + updates portal_export_v2/communities_F01.json, communities_F34.json
  adding `label` field per community.
"""
from __future__ import annotations
import json
import shutil
import subprocess
from pathlib import Path
import pandas as pd
from scipy.stats import hypergeom

# ------------------------------------------------------------------ paths
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
NET_DIR = ROOT / "RNA-seq/results/network/communities_f2"
PORTAL_DIR = ROOT / "masld-atlas-v2/public/data/network/portal_export_v2"
GMT = ROOT / "Analysis/downstream_analysis/pathway_analysis/data/msigdb.v2025.1.Hs.symbols.gmt"

OUT_JSON = NET_DIR / "community_labels.json"
OUT_TSV = NET_DIR / "community_labels.tsv"

# ------------------------------------------------------------------ rules
# Ordered rules: first match wins among SIGNIFICANT hallmarks.
LABEL_RULES = [
    ("inflammatory", {"HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"}, "Inflammatory"),
    ("oxphos", {"HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM"}, "OxPhos"),
    ("fibrogenic", {"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_APICAL_JUNCTION",
                    "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_MYOGENESIS"}, "Fibrogenic"),
    ("proliferation", {"HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
                       "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                       "HALLMARK_MITOTIC_SPINDLE"}, "Proliferation"),
    ("interferon", {"HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"}, "Interferon"),
    ("hypoxic", {"HALLMARK_HYPOXIA"}, "Hypoxic"),
    ("immune_signaling", {"HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                          "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION"},
     "Immune-Signaling"),
    ("hepatic_detox", {"HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM"}, "Hepatic-Detox"),
    ("stress_signaling", {"HALLMARK_COAGULATION", "HALLMARK_KRAS_SIGNALING_UP",
                          "HALLMARK_KRAS_SIGNALING_DN",
                          "HALLMARK_UV_RESPONSE_UP", "HALLMARK_UV_RESPONSE_DN",
                          "HALLMARK_P53_PATHWAY", "HALLMARK_APOPTOSIS"}, "Stress-Signaling"),
    ("adipogenic", {"HALLMARK_ADIPOGENESIS"}, "Adipogenic"),
    ("cholesterol", {"HALLMARK_CHOLESTEROL_HOMEOSTASIS"}, "Cholesterol"),
    ("glycolysis", {"HALLMARK_GLYCOLYSIS"}, "Glycolysis"),
    ("mtor", {"HALLMARK_MTORC1_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING"}, "MTOR"),
]

# Bonferroni-adjusted threshold across 50 hallmarks (strict)
P_BONF_THRESH = 0.01 / 50  # = 2e-4
# Nominal fallback threshold (for communities with no Bonferroni-significant hit)
P_NOMINAL_THRESH = 0.05


def clean_hallmark(name: str) -> str:
    return name.replace("HALLMARK_", "").replace("_", " ").title()


def load_hallmarks(gmt_path: Path) -> dict[str, set[str]]:
    sets = {}
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0].startswith("HALLMARK_"):
                continue
            name = parts[0]
            genes = set(parts[2:])
            sets[name] = genes
    return sets


def compute_enrichment(comm_genes: set[str], universe: set[str],
                       hallmarks: dict[str, set[str]]) -> pd.DataFrame:
    """Hypergeometric enrichment for one community vs all hallmarks.

    K = |universe ∩ hallmark|  (white balls in urn)
    N = |universe|             (total balls)
    n = |community|            (draws)
    k = |community ∩ hallmark| (white drawn)
    P(X >= k) = hypergeom.sf(k-1, N, K, n)
    """
    N = len(universe)
    n = len(comm_genes)
    rows = []
    for hname, hgenes in hallmarks.items():
        hg_in_u = hgenes & universe
        K = len(hg_in_u)
        k = len(comm_genes & hg_in_u)
        if K == 0 or k == 0:
            pval = 1.0
        else:
            pval = float(hypergeom.sf(k - 1, N, K, n))
        rows.append({"gene_set": hname, "overlap": k, "set_size": K,
                     "comm_size_in_universe": n, "universe": N, "pvalue": pval})
    return pd.DataFrame(rows).sort_values("pvalue").reset_index(drop=True)


def assign_label(bonf_sig: list[str], nominal_sig: list[str],
                 top_set: str, top_p: float) -> tuple[str, str]:
    """Two-tier rule matching.

    1. Walk Bonferroni-significant hallmarks (p < 2e-4) in ascending-p order;
       return first curated rule match.
    2. If no curated rule matches among Bonferroni-significant, try again
       among nominally-significant hallmarks (p < 0.05).
    3. If a Bonferroni-significant top hit exists but no curated rule matches,
       use the cleaned hallmark name as label (top_fallback_bonf).
    4. If nothing is even nominally significant, label as Mixed.
    """
    for hname in bonf_sig:
        for rule_name, hset, base in LABEL_RULES:
            if hname in hset:
                return base, rule_name
    for hname in nominal_sig:
        for rule_name, hset, base in LABEL_RULES:
            if hname in hset:
                return base, rule_name
    if bonf_sig:
        return clean_hallmark(bonf_sig[0]), "top_fallback_bonf"
    if nominal_sig:
        return clean_hallmark(nominal_sig[0]), "top_fallback_nominal"
    return "Mixed", "no_enrichment"


# ------------------------------------------------------------------ load data
print(f"[load] hallmarks: {GMT}")
HALLMARKS = load_hallmarks(GMT)
print(f"[load] {len(HALLMARKS)} hallmark sets")

trans = pd.read_csv(NET_DIR / "community_transition_table.csv")
com_f01 = pd.read_csv(NET_DIR / "communities_F01.csv")
com_f34 = pd.read_csv(NET_DIR / "communities_F34.csv")

# Universe: all genes with a community assignment (same across stages)
universe_f01 = set(com_f01["gene"].unique())
universe_f34 = set(com_f34["gene"].unique())
print(f"[univ] F01={len(universe_f01)}  F34={len(universe_f34)}")

# community -> gene set
genes_by_comm_f01 = com_f01.groupby("community_id")["gene"].apply(set).to_dict()
genes_by_comm_f34 = com_f34.groupby("community_id")["gene"].apply(set).to_dict()

# ------------------------------------------------------------------ per-community enrichment
def enrich_all(genes_by_comm: dict, universe: set[str]) -> dict[int, pd.DataFrame]:
    out = {}
    for cid, gs in genes_by_comm.items():
        out[cid] = compute_enrichment(gs & universe, universe, HALLMARKS)
    return out


print("[enrich] F01 ...")
enr_f01 = enrich_all(genes_by_comm_f01, universe_f01)
print("[enrich] F34 ...")
enr_f34 = enrich_all(genes_by_comm_f34, universe_f34)

# ------------------------------------------------------------------ transition matching
size_f01 = {c: len(g) for c, g in genes_by_comm_f01.items()}
size_f34 = {c: len(g) for c, g in genes_by_comm_f34.items()}

best_f01_to_f34 = (trans.sort_values("jaccard", ascending=False)
                        .drop_duplicates("F01_id", keep="first")
                        .set_index("F01_id"))
best_f34_to_f01 = (trans.sort_values("jaccard", ascending=False)
                        .drop_duplicates("F34_id", keep="first")
                        .set_index("F34_id"))

f01_class = (trans.drop_duplicates("F01_id", keep="first")
                  .set_index("F01_id")["F01_class"].to_dict())
f34_class = (trans.drop_duplicates("F34_id", keep="first")
                  .set_index("F34_id")["F34_class"].to_dict())


def build_stage(stage: str, enr_by_comm: dict, size_self: dict,
                size_other: dict, best_map: pd.DataFrame,
                trans_class: dict, other_key: str):
    out = {}
    for cid, df in enr_by_comm.items():
        bonf_sig = df[df["pvalue"] < P_BONF_THRESH].sort_values("pvalue")["gene_set"].tolist()
        nominal_sig = df[df["pvalue"] < P_NOMINAL_THRESH].sort_values("pvalue")["gene_set"].tolist()
        tops = df.head(3)["gene_set"].tolist()
        p_top = float(df["pvalue"].iloc[0]) if len(df) else 1.0
        base, rule = assign_label(bonf_sig, nominal_sig, tops[0] if tops else "", p_top)

        n_self = int(size_self.get(cid, 0))
        matched_other = None
        size_delta = None
        if cid in best_map.index:
            matched_other = int(best_map.loc[cid, other_key])
            n_other = int(size_other.get(matched_other, 0))
            # size_delta = (matched community size on OTHER side) - (this side)
            # For F01 stage: how much does the matched F34 community expand/contract?
            # For F34 stage: how much did this F34 grow vs the matched F01 source?
            size_delta = (n_other - n_self) if stage == "F01" else (n_self - n_other)

        # Expansion/Contraction qualifier based on 10% threshold
        ratio = None
        if size_delta is not None and n_self > 0:
            ratio = size_delta / n_self
        suffix = ""
        curated = rule not in ("top_fallback_bonf", "top_fallback_nominal", "no_enrichment")
        if curated:
            if ratio is not None and ratio > 0.10:
                suffix = "-Expanding"
            elif ratio is not None and ratio < -0.10:
                suffix = "-Contracting"
            else:
                suffix = "-Stable"

        label = f"{base}{suffix}"

        if size_delta is None:
            dshift = "unmatched"
        elif abs(size_delta) < 50:
            dshift = "stable"
        elif size_delta > 0:
            dshift = "expanding"
        else:
            dshift = "contracting"

        tclass = trans_class.get(cid, "unknown")
        out[str(int(cid))] = {
            "label": label,
            "top_hallmarks": tops,
            "p_top": p_top,
            "n_sig_hallmarks_bonf": int(len(bonf_sig)),
            "n_sig_hallmarks_nominal": int(len(nominal_sig)),
            "n_genes": n_self,
            "transition_class": tclass,
            "matched_community_other_side": matched_other,
            "size_delta_to_matched": size_delta,
            "dominant_shift": dshift,
            "rule_matched": rule,
        }
    return out


labels = {
    "F01": build_stage("F01", enr_f01, size_f01, size_f34,
                       best_f01_to_f34, f01_class, "F34_id"),
    "F34": build_stage("F34", enr_f34, size_f34, size_f01,
                       best_f34_to_f01, f34_class, "F01_id"),
}

# ------------------------------------------------------------------ write outputs
OUT_JSON.write_text(json.dumps(labels, indent=2))
print(f"[write] {OUT_JSON}")

tsv_rows = []
for stage, comms in labels.items():
    for cid, rec in comms.items():
        tsv_rows.append({
            "stage": stage,
            "community_id": cid,
            "label": rec["label"],
            "rule_matched": rec["rule_matched"],
            "n_genes": rec["n_genes"],
            "n_sig_hallmarks_bonf": rec["n_sig_hallmarks_bonf"],
            "n_sig_hallmarks_nominal": rec["n_sig_hallmarks_nominal"],
            "transition_class": rec["transition_class"],
            "matched_community_other_side": rec["matched_community_other_side"],
            "size_delta_to_matched": rec["size_delta_to_matched"],
            "dominant_shift": rec["dominant_shift"],
            "top_hallmark_1": rec["top_hallmarks"][0] if len(rec["top_hallmarks"]) > 0 else "",
            "top_hallmark_2": rec["top_hallmarks"][1] if len(rec["top_hallmarks"]) > 1 else "",
            "top_hallmark_3": rec["top_hallmarks"][2] if len(rec["top_hallmarks"]) > 2 else "",
            "p_top": rec["p_top"],
        })
pd.DataFrame(tsv_rows).to_csv(OUT_TSV, sep="\t", index=False)
print(f"[write] {OUT_TSV}")

# ------------------------------------------------------------------ deploy
PORTAL_DIR.mkdir(parents=True, exist_ok=True)
shutil.copy(OUT_JSON, PORTAL_DIR / "community_labels.json")
print(f"[copy]  {PORTAL_DIR / 'community_labels.json'}")

for stage, fname in [("F01", "communities_F01.json"), ("F34", "communities_F34.json")]:
    p = PORTAL_DIR / fname
    if not p.exists():
        print(f"[warn] {p} missing; skipping label injection")
        continue
    data = json.loads(p.read_text())
    stage_labels = labels[stage]
    n_patched = 0
    for c in data.get("communities", []):
        cid = str(int(c["id"]))
        if cid in stage_labels:
            c["label"] = stage_labels[cid]["label"]
            c["label_rule"] = stage_labels[cid]["rule_matched"]
            c["top_hallmarks"] = stage_labels[cid]["top_hallmarks"]
            c["transition_class"] = stage_labels[cid]["transition_class"]
            c["dominant_shift"] = stage_labels[cid]["dominant_shift"]
            n_patched += 1
    p.write_text(json.dumps(data, indent=2))
    print(f"[patch] {p}  ({n_patched}/{len(data.get('communities', []))} communities labeled)")

# ------------------------------------------------------------------ validate JSONs
for f in [OUT_JSON, PORTAL_DIR / "community_labels.json",
          PORTAL_DIR / "communities_F01.json", PORTAL_DIR / "communities_F34.json"]:
    if f.exists():
        subprocess.run(["python", "-m", "json.tool", str(f)], stdout=subprocess.DEVNULL, check=True)
        print(f"[valid] {f}")

# ------------------------------------------------------------------ summary
print("\n=== Label distribution ===")
for stage in ["F01", "F34"]:
    recs = labels[stage]
    total = len(recs)
    fallback_rules = ("top_fallback_bonf", "top_fallback_nominal", "no_enrichment")
    curated_recs = [r for r in recs.values() if r["rule_matched"] not in fallback_rules]
    n_curated = len(curated_recs)
    n_fallback = total - n_curated
    pct = 100.0 * n_curated / total if total else 0.0
    genes_curated = sum(r["n_genes"] for r in curated_recs)
    genes_total = sum(r["n_genes"] for r in recs.values())
    pct_genes = 100.0 * genes_curated / genes_total if genes_total else 0.0
    # Also report any-labeled (curated + fallback with some enrichment)
    n_any = sum(1 for r in recs.values() if r["rule_matched"] != "no_enrichment")
    genes_any = sum(r["n_genes"] for r in recs.values() if r["rule_matched"] != "no_enrichment")
    pct_any_genes = 100.0 * genes_any / genes_total if genes_total else 0.0
    print(f"\n[{stage}] n_communities={total}")
    print(f"  curated:       {n_curated}/{total} ({pct:.1f}%)  |  genes: {genes_curated}/{genes_total} ({pct_genes:.1f}%)")
    print(f"  any enriched:  {n_any}/{total} ({100*n_any/total:.1f}%)  |  genes: {genes_any}/{genes_total} ({pct_any_genes:.1f}%)")
    print(f"  no_enrichment (Mixed): {total - n_any} (all tiny communities, median size ~2 genes)")
    counts = pd.Series([r["label"] for r in recs.values()]).value_counts()
    for lab, n in counts.items():
        print(f"  {n:3d}  {lab}")
