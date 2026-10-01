#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
# T1 v2 founders: EUR and EAS samples with no recorded parents in the 1kGP pedigree
# (trio parents kept, children dropped); KING kinship on common biallelic SNVs from
# chr15-22 (MAF >= 0.05 within the kept samples, thinned to 20%, LD-pruned 200 kb r2 0.2);
# for each pair with kinship > 0.0884 the sample with the lexicographically larger ID is
# written to exclude_related.txt. Public data only.
# Usage: sbatch --job-name=<n> --output=<log> t1_founders_king.sh <out_dir>
set -euo pipefail
out=$1; mkdir -p "$out"
kg=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/allelic_refs/kgp_highcov
module load plink/2.0a5.13 >/dev/null 2>&1
py="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
$py - "$kg/20130606_g1k_3202_samples_ped_population.txt" "$out" <<'PY'
import sys, pandas as pd
p = pd.read_csv(sys.argv[1], sep=r"\s+")
f = p[p.Superpopulation.isin(["EUR", "EAS"]) & (p.FatherID.astype(str) == "0") & (p.MotherID.astype(str) == "0")]
f[["SampleID"]].to_csv(f"{sys.argv[2]}/founders_keep.txt", index=False, header=False)
f[["SampleID", "Population", "Superpopulation", "FamilyID"]].to_csv(f"{sys.argv[2]}/founders.tsv", sep="\t", index=False)
print(f.groupby(["Superpopulation", "Population"]).size().to_string())
PY
ls $kg/1kGP_high_coverage_Illumina.chr{15,16,17,18,19,20,21,22}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz > /dev/null
: > "$out/merge_list.txt"
for c in 15 16 17 18 19 20 21 22; do
  plink --vcf "$kg/1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" --keep "$out/founders_keep.txt" \
    --snps-only just-acgt --max-alleles 2 --maf 0.05 --thin 0.2 --seed 20260929 --set-all-var-ids @:#:\$r:\$a \
    --rm-dup exclude-all --make-pgen --threads 8 --out "$out/chr$c" > /dev/null
  echo "$out/chr$c" >> "$out/merge_list.txt"
done
plink --pmerge-list "$out/merge_list.txt" pfile --make-pgen --threads 8 --out "$out/merged" > /dev/null
plink --pfile "$out/merged" --indep-pairwise 200kb 0.2 --threads 8 --out "$out/prune" > /dev/null
plink --pfile "$out/merged" --extract "$out/prune.prune.in" --make-king-table --king-table-filter 0.0884 --threads 8 --out "$out/king" > /dev/null
echo "pruned SNVs used: $(wc -l < "$out/prune.prune.in")"
$py - "$out" <<'PY'
import sys, pandas as pd
o = sys.argv[1]
k = pd.read_csv(f"{o}/king.kin0", sep=r"\s+")
k = k.rename(columns={"#IID1": "IID1"})
fo = pd.read_csv(f"{o}/founders.tsv", sep="\t").set_index("SampleID")
k["pop1"], k["pop2"] = k.IID1.map(fo.Population), k.IID2.map(fo.Population)
k["dropped"] = [max(a, b) for a, b in zip(k.IID1, k.IID2)]
k.to_csv(f"{o}/related_pairs_kinship_gt_0.0884.tsv", sep="\t", index=False)
drop = sorted(set(k.dropped))
open(f"{o}/exclude_related.txt", "w").write("\n".join(drop) + ("\n" if drop else ""))
kept = fo.drop(index=[d for d in drop if d in fo.index])
print("related pairs", len(k), "| samples dropped", len(drop))
print(kept.groupby(["Superpopulation", "Population"]).size().to_string())
PY
