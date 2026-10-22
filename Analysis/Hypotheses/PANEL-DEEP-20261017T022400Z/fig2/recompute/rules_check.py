"""Style-rule checks on every Fig 2 panel dump: font sizes, bold faces, non-black text, text beyond the page,
axis-title / n tokens. usage: rules_check.py <dump dir> <tag>"""
import sys, json, glob, os, pandas as pd
D, tag = sys.argv[1], sys.argv[2]
W = {"fig2A_gwas_cascade": 254.16, "Fig2B_PIP_vs_SuSiE-coloc": 129, "Fig2C_coding_noncoding_split": 93, "Fig2D_pip_architecture_by_trait_directness": 133.2,
     "Fig2D_coding_by_ancestry": 403, "Fig2D_direct_subclass_split": 302, "Fig2E_ancestry_unique_coloc_GWS": 144, "Fig2F_crossancestry_coloc": 185.76,
     "Fig2G_RORA_NAFLD": 234.72, "Fig2H_FABP1_ALT": 234.72, "Fig2I_locus_legend": 475.92}
TOK = {"fig2A_gwas_cascade": ["5 ancestries", "35 GWAS", "Fine-mapped", "COLOC", "PP.H4"],
       "Fig2B_PIP_vs_SuSiE-coloc": ["SuSiE PIP of COLOC lead", "PP.H4", "327 genes"],
       "Fig2C_coding_noncoding_split": ["996 colocalizing genes"],
       "Fig2D_pip_architecture_by_trait_directness": ["Direct", "Liver enzymes"],
       "Fig2D_coding_by_ancestry": ["Direct minus enzyme", "percentage points"],
       "Fig2D_direct_subclass_split": ["Direct subset minus all liver enzymes", "percentage points"],
       "Fig2E_ancestry_unique_coloc_GWS": ["EUR-only", "shared", "GWS"],
       "Fig2F_crossancestry_coloc": ["EUR", "AFR", "AMR", "EAS", "SAS", "PP.H4", "multi-signal", "single-signal only"],
       "Fig2G_RORA_NAFLD": ["chr15", "Mb"], "Fig2H_FABP1_ALT": ["chr2", "Mb"], "Fig2I_locus_legend": ["n.s."]}
out = {}
for p, w in W.items():
    ch = pd.read_csv(f"{D}/{p}.chars.tsv", sep="\t"); wd = pd.read_csv(f"{D}/{p}.words.tsv", sep="\t")
    txt = " ".join(wd.text.astype(str))
    rot = "".join(ch[ch.upright == False].sort_values("top", ascending=False).text) if "upright" in ch else ""
    out[p] = {"font_sizes": sorted(ch["size"].round(2).unique().tolist()), "bold": sorted(ch[ch.font.str.contains("Bold", case=False)].text.unique().tolist())[:10],
              "fonts": sorted(ch.font.unique().tolist()), "text_fills": ch.fill.value_counts().to_dict(),
              "beyond_page": ch[(ch.x1 > w + 0.05) | (ch.x0 < -0.05)].text.tolist(),
              "missing_tokens": [t for t in TOK[p] if t not in txt and t.replace(" ", "") not in rot.replace(" ", "")],
              "rotated_text": rot[:80]}
json.dump(out, open(f"out_geom/rules_{tag}.json", "w"), indent=1)
print(json.dumps(out, indent=1))
