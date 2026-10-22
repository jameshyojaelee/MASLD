#!/bin/bash
# Job 1161 (review-R16, K05): how many plotted Figure 2B points change x when the
# fine-mapping PIP is SuSiE-only instead of max(recommended_pip, susie_pip)?
# Read-only. Streams the two fine-mapping CSVs and keeps only rows whose
# (study, chromosome, position) match a Figure 2B source row.
set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$PWD
SRC=$BASE/figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv
FM1=$BASE/GWAS/finemapping/results/combined_finemapping.csv
FM2=$BASE/GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv

# keys: study \t chr \t pos from source columns pp4_gwas (3) and x_top_snp (4)
awk -F, 'NR>1{split($4,a,":"); print $3"\t"a[1]"\t"a[2]}' "$SRC" | sort -u > "$OUT/keys.tsv"
wc -l "$OUT/keys.tsv"

filter() {
  local f=$1 tag=$2
  awk -F, -v tag="$tag" 'BEGIN{OFS="\t"; while((getline l < "'"$OUT"'/keys.tsv")>0){split(l,k,"\t"); want[k[1]"|"k[2]"|"k[3]]=1}}
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "source","study","chromosome","position","susie_pip","carma_pip","susiex_pip","mesusie_pip","recommended_pip"; next}
    { key=$h["study"]"|"$h["chromosome"]"|"$h["position"]; if(key in want) print tag,$h["study"],$h["chromosome"],$h["position"],$h["susie_pip"],(("carma_pip" in h)?$h["carma_pip"]:""),(("susiex_pip" in h)?$h["susiex_pip"]:""),(("mesusie_pip" in h)?$h["mesusie_pip"]:""),$h["recommended_pip"] }' "$f"
}
head -1 "$FM2" | tr ',' '\n' | grep -n "study\|chromosome\|position\|pip" > "$OUT/fm2_header_cols.txt"
{ filter "$FM1" canonical_finemapping; filter "$FM2" coloc_targeted_finemapping | tail -n +2; } > "$OUT/fm_rows_matching_fig2b.tsv"
wc -l "$OUT/fm_rows_matching_fig2b.tsv"

PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python - "$SRC" "$OUT" <<'PY'
import sys, pandas as pd, numpy as np
src, out = sys.argv[1], sys.argv[2]
s = pd.read_csv(src)
fm = pd.read_csv(f"{out}/fm_rows_matching_fig2b.tsv", sep="\t")
for c in ["susie_pip","carma_pip","susiex_pip","mesusie_pip","recommended_pip"]:
    fm[c] = pd.to_numeric(fm[c], errors="coerce")
# Producer definition: causal_pip = pmax(recommended_pip, susie_pip, na.rm=TRUE), max over rows per key
fm["causal_pip"] = fm[["recommended_pip","susie_pip"]].max(axis=1, skipna=True)
g = fm.groupby(["study","chromosome","position"])
per = g.agg(causal_pip_max=("causal_pip","max"), susie_only_max=("susie_pip","max"),
            carma_max=("carma_pip","max"), susiex_max=("susiex_pip","max"), mesusie_max=("mesusie_pip","max"),
            n_rows=("susie_pip","size")).reset_index()
s[["chromosome","position"]] = s["x_top_snp"].str.split(":", expand=True).astype(int)
m = s.merge(per, left_on=["pp4_gwas","chromosome","position"], right_on=["study","chromosome","position"], how="left")
m["repro_diff"] = (m["causal_pip"] - m["causal_pip_max"]).abs()
p = m[m["plotted"] == True].copy()
lines = []
lines.append(f"source_rows\t{len(s)}")
lines.append(f"plotted_rows\t{len(p)}")
lines.append(f"plotted_rows_matched_to_finemapping\t{p['causal_pip_max'].notna().sum()}")
lines.append(f"plotted_rows_reproducing_causal_pip_within_1e-6\t{(p['repro_diff']<1e-6).sum()}")
d = (p["causal_pip"] - p["susie_only_max"])
lines.append(f"plotted_rows_with_finite_susie_only\t{p['susie_only_max'].notna().sum()}")
lines.append(f"plotted_rows_without_any_susie_pip\t{p['susie_only_max'].isna().sum()}")
for thr in [0.01, 0.05, 0.1, 0.25]:
    lines.append(f"plotted_rows_x_drop_gt_{thr}\t{(d>thr).sum()}")
lines.append(f"plotted_rows_x_ge_0.5_as_plotted\t{(p['causal_pip']>=0.5).sum()}")
lines.append(f"plotted_rows_x_ge_0.5_susie_only\t{(p['susie_only_max']>=0.5).sum()}")
lines.append(f"plotted_rows_cross_0.5_downward\t{((p['causal_pip']>=0.5)&~(p['susie_only_max']>=0.5)).sum()}")
lines.append(f"median_x_drop_among_changed_gt_0.01\t{d[d>0.01].median() if (d>0.01).any() else 'NA'}")
# which method supplied the max when it exceeds SuSiE
def src_of(r):
    vals = {"carma": r.carma_max, "susiex": r.susiex_max, "mesusie": r.mesusie_max}
    best = max(vals, key=lambda k: -1 if pd.isna(vals[k]) else vals[k])
    return best
ch = p[d>0.01]
if len(ch):
    lines.append("changed_rows_max_method_counts\t" + ch.apply(src_of, axis=1).value_counts().to_json())
lab = p[p["is_label"]==True][["gene","pp4_gwas","x_top_snp","causal_pip","susie_only_max","carma_max","susiex_max","mesusie_max","pp4"]]
open(f"{out}/k05_summary.tsv","w").write("\n".join(lines)+"\n")
lab.to_csv(f"{out}/k05_labelled_points.tsv", sep="\t", index=False)
p[["gene","pp4_gwas","x_top_snp","class","causal_pip","susie_only_max","carma_max","susiex_max","mesusie_max","pp4"]].to_csv(f"{out}/k05_plotted_points.tsv", sep="\t", index=False)
print("\n".join(lines)); print(lab.to_string())
PY
PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
echo DONE
