"""Light n_actual counts for lens-L03 candidates. Reads saved tables only."""
import glob
import pandas as pd

R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
P = R + 'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
g = pd.read_parquet(P + 'passport_gene_index.parquet')
e = pd.read_parquet(P + 'passport_evidence_long.parquet')
e['ens'] = e['ensembl_id'].astype(str).str.split('.').str[0]
g['ens'] = g['ensembl_id'].astype(str).str.split('.').str[0]

testable = e[e['testability_state'] != 'untestable'] if 'testability_state' in e else None
tt = e[(e.evidence_domain == 'transcriptomics') & (e.call_state != 'untestable')].ens
gt = e[(e.evidence_domain == 'genetics') & (e.call_state != 'untestable')].ens
joint = set(tt) & set(gt)
print('C1 joint testable universe:', len(joint))
print('C1 classes:', g.primary_evidence_class.value_counts().to_dict())
state = set(g[g.primary_evidence_class.isin(['established_state_associated', 'concordant'])].ens)
print('C1 established_state+concordant:', len(state), 'in joint:', len(state & joint))
print('C1 next_experiment:', g.next_experiment_rule_id.value_counts().to_dict())

eq = pd.concat([pd.read_csv(f, sep='\t') for f in
                glob.glob(R + 'GWAS/finemapping/results/eqtl_susie/chr*_summary.tsv')])
eq['ens'] = eq['ensg'].astype(str).str.split('.').str[0]
eq = eq.groupby('ens')['n_cs'].max()
eg = set(eq[eq > 0].index)
print('C1 eQTL tested genes:', eq.shape[0], 'eGenes:', len(eg),
      'joint tested:', len(joint & set(eq.index)), 'joint eGenes:', len(joint & eg))
print('C1 state genes that are eGenes:', len(state & joint & eg))

cs = pd.read_csv(R + 'figures/candidates/bulk-context-774-20261014T140904Z_22236311/cohort_summary.tsv', sep='\t')
print('C2 cohort_summary:\n', cs.to_string())
co = pd.read_csv(R + 'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/'
                 'primary_source_controls/coefficient_table.tsv.gz', sep='\t')
print('C2 pooled genes:', co.shape[0], 'calls padj<0.05:', int((co.padj < 0.05).sum()))

comp = pd.read_csv(R + 'Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz', sep='\t')
print('C4 composition shape:', comp.shape)

stress = set(g[g.next_experiment_rule_id == 'EXP_STRESS_PERTURB_V1'].ens)
gn = pd.read_csv(R + 'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv', sep='\t',
                 usecols=['gene', 'gene_id', 'mane_select', 'lof.oe_ci.upper'])
gn = gn[gn.mane_select.astype(str).str.lower().isin(['true','1'])]
gn_ens = set(gn.gene_id.astype(str).str.split('.').str[0])
print('C6 MANE rows:', gn.shape[0], 'stress genes:', len(stress), 'with LOEUF:', len(stress & gn_ens),
      'joint with LOEUF:', len(joint & gn_ens))
hpa = pd.read_csv(R + 'data/external/hpa/proteinatlas.tsv', sep='\t', usecols=['Ensembl'])
print('C6 HPA genes:', hpa.shape[0], 'stress in HPA:', len(stress & set(hpa.Ensembl)))

dt = pd.read_csv(R + 'data/external/drug_targets/drug_target_classification.tsv', sep='\t')
print('C7 drug cols:', list(dt.columns)[:12])
print('C2 coef cols:', list(co.columns))
for thr in [0.0, 0.5, 0.585, 1.0]:
    print('C2 padj<0.05 & |logFC|>', thr, int(((co.padj < 0.05) & (co.logFC.abs() > thr)).sum()))
print('C6 gnomAD unique gene_id:', gn.gene_id.nunique())
dt['ens'] = dt.ensembl_id.astype(str).str.split('.').str[0]
clin = dt[dt.drug_dev_status.astype(str).str.contains('masld', case=False)]
print('C7 drug_dev_status:', dt.drug_dev_status.value_counts().to_dict())
print('C7 masld-status genes:', clin.shape[0], 'in joint:', len(set(clin.ens) & joint),
      'max_phase_masld>=2:', int((dt.max_phase_masld >= 2).sum()))
