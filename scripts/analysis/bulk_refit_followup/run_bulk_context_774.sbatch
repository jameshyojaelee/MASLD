#!/bin/bash
#SBATCH --job-name=bulk_context774
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/bulk_context774_%j.out
#SBATCH --error=logs/bulk_context774_%j.err
set -euo pipefail
export MASLD_PROJECT_ROOT="$PWD"
export BULK_FOLLOWUP_ROOT="$PWD/figures/candidates/bulk-context-774-$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}"
export FIGURE_CANDIDATE_ROOT="$BULK_FOLLOWUP_ROOT"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONNOUSERSITE=1
r_env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export BULK_DEG_FILE="$PWD/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
export FIG5B_BULK_DEG="$BULK_DEG_FILE"
mkdir -p "$BULK_FOLLOWUP_ROOT/source" "$BULK_FOLLOWUP_ROOT/class_difference/results" "$BULK_FOLLOWUP_ROOT/class_difference/prespec"
printf '%s\n' "$BULK_FOLLOWUP_ROOT" > "logs/bulk_context774_${SLURM_JOB_ID}.root"
cp scripts/analysis/bulk_refit_followup/restore_constant_sex_cohort.R scripts/analysis/bulk_refit_followup/run_bulk_context_774.sbatch scripts/figures/fig5b_protein_triage.R scripts/figures/composite_mrna_protein.R "$BULK_FOLLOWUP_ROOT/source/"
cp Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R scripts/analysis/bulk_refit_followup/class_pair_label_negative_control.py "$BULK_FOLLOWUP_ROOT/source/"
export FIG5C_BULK_ID_FILE="$BULK_FOLLOWUP_ROOT/panel4c_bulk_gene_ids.tsv"
"$r_env/python" - <<'PY2'
import os,pandas as pd
from pathlib import Path
base=Path(os.environ['MASLD_PROJECT_ROOT']); out=Path(os.environ['BULK_FOLLOWUP_ROOT'])
m=pd.read_csv(base/'Analysis/Multimodal_Program_Projection/config/panel4c_bulk_gene_ids.tsv',sep='\t')
b=pd.read_csv(os.environ['BULK_DEG_FILE']); b['base']=b.gene.str.split('.').str[0]
for i,row in m.iterrows():
    if row.bulk_gene_id.split('.')[0] not in set(b['base']):
        hit=b[b.symbol.eq(row.gene)]
        assert len(hit)==1,(row.gene,len(hit))
        m.loc[i,'bulk_gene_id']=hit.iloc[0].gene
assert m[m.gene.eq('HLA-DQA1')].bulk_gene_id.iloc[0].split('.')[0]=='ENSG00000196735'
m.to_csv(out/'panel4c_bulk_gene_ids.tsv',sep='\t',index=False)
PY2
sha256sum "$BULK_DEG_FILE" "$FIG5C_BULK_ID_FILE" > "$BULK_FOLLOWUP_ROOT/inputs_before.sha256"
"$r_env/Rscript" --vanilla "$BULK_FOLLOWUP_ROOT/source/restore_constant_sex_cohort.R"
PROTEOMICS_OUT_DIR="$BULK_FOLLOWUP_ROOT/proteomics" "$r_env/Rscript" --vanilla "$BULK_FOLLOWUP_ROOT/source/02_proteomics_projection.R"
export FIG5C_PROTEOMICS_DIR="$BULK_FOLLOWUP_ROOT/proteomics"
"$r_env/Rscript" --vanilla "$BULK_FOLLOWUP_ROOT/source/fig5b_protein_triage.R"
"$r_env/Rscript" --vanilla "$BULK_FOLLOWUP_ROOT/source/composite_mrna_protein.R"
old_class="$PWD/RNA-seq/results/class_difference_test/cdt-corrected-20260925T203546Z"
cp "$old_class/prespec/PRESPEC.md" "$BULK_FOLLOWUP_ROOT/class_difference/prespec/"
cp "$old_class/"0[0-4]_*.py "$BULK_FOLLOWUP_ROOT/class_difference/"
export CDT_BULK="$BULK_DEG_FILE"
export CDT_COLOC="$PWD/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
"$r_env/python" - <<'PY3'
import os,pandas as pd
from pathlib import Path
out=Path(os.environ['BULK_FOLLOWUP_ROOT'])/'class_difference'
for p in out.glob('*.py'):
    s=p.read_text().replace('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design',os.environ['MASLD_PROJECT_ROOT']).replace('bulk_background_14920','bulk_background')
    if p.name=='03_run_class_test.py':
        s=s.replace('R.to_csv(OUT/"class_difference_results.tsv", sep="\\t", index=False)', 'R["relabel_rematch_calibrated"] = False\nR["relabel_rematch_interpretation"] = "invalid_anti_conservative_diagnostic"\nR.to_csv(OUT/"class_difference_results.tsv", sep="\\t", index=False)')
    if p.name=='04_signflip_diagnostic.py':
        s=s.replace('result = pd.DataFrame(rows)', 'result = pd.DataFrame(rows)\nresult["signflip_interpretation"] = "post_hoc_nominal_pair_exchangeability"\nresult["gene_dependence_adjusted"] = False')
    p.write_text(s)
b=pd.read_csv(os.environ['CDT_BULK']); c=pd.read_csv(os.environ['CDT_COLOC']);b['g']=b.gene.str.split('.').str[0];c['g']=c.ensembl.str.split('.').str[0]
c=c.sort_values(['coloc_best_susie_pp4','coloc_best_abf_pp4'],ascending=False).drop_duplicates('g')
u=b.merge(c,on='g');co=u.coloc_best_susie_pp4.ge(.5);de=u.padj.lt(.05)&u.logFC.abs().gt(.5)
(out/'expected_counts.txt').write_text(','.join(map(str,[len(u),int(co.sum()),int(de.sum()),int((co&de).sum())])))
PY3
export CDT_EXPECT="$(cat "$BULK_FOLLOWUP_ROOT/class_difference/expected_counts.txt")"
for script in 00_assemble_substrate 01_matching_power 02_caliper_scan; do
  "$r_env/python" "$BULK_FOLLOWUP_ROOT/class_difference/${script}.py" "$BULK_FOLLOWUP_ROOT/class_difference/results"
done
for script in 03_run_class_test 04_signflip_diagnostic; do
  "$r_env/python" "$BULK_FOLLOWUP_ROOT/class_difference/${script}.py" "$BULK_FOLLOWUP_ROOT/class_difference"
done
CLASS_COMPARISON_ROOT="$BULK_FOLLOWUP_ROOT/class_difference" \
CLASS_NEGATIVE_CONTROL_ROOT="$BULK_FOLLOWUP_ROOT/class_difference/negative_control" \
  "$r_env/python" "$BULK_FOLLOWUP_ROOT/source/class_pair_label_negative_control.py"
sha256sum -c "$BULK_FOLLOWUP_ROOT/inputs_before.sha256"
touch "$BULK_FOLLOWUP_ROOT/COMPLETE"
