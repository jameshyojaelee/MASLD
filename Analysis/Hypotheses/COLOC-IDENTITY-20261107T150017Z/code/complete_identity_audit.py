#!/usr/bin/env python3
"""Count confirmed fit defects separately from reference-position exposure."""
from pathlib import Path
import json
import pandas as pd

RUN = Path(__file__).resolve().parents[1]
ROOT = RUN.parents[2]
OUT = RUN / 'audit'
fits = pd.read_csv(OUT / 'all_saved_eqtl_fit_identity_checks.tsv', sep='\t')
positive = pd.read_csv(OUT / 'positive_gene_study_identity_risk.tsv', sep='\t')
# Panel assignment in the first screen was ancestry-derived, not recovered from
# each historical GWAS fit. Rename it explicitly; do not claim fitted exposure.
positive = positive.rename(columns={'actual_panel': 'reference_panel_for_screen'})
exports = pd.read_csv(OUT / 'posterior_export_inventory.tsv', sep='\t')
rows = []
for r in exports.itertuples(index=False):
    d = pd.read_csv(r.path, sep='\t')
    assert {'signal_pair_index', 'snp'}.issubset(d.columns)
    duplicated = d.duplicated(['signal_pair_index', 'snp'], keep=False)
    rows.append(dict(gwas_name=r.study, ensembl=r.ensembl,
        exported_rows=len(d), duplicate_pair_variant_rows=int(duplicated.sum()),
        duplicate_positions=int(d.loc[duplicated, 'snp'].nunique()), path=r.path))
export_checks = pd.DataFrame(rows)
export_checks.to_csv(OUT / 'posterior_export_identity_checks.tsv', sep='\t', index=False)
positive = positive.merge(export_checks, on=['gwas_name', 'ensembl'], how='left', validate='one_to_one')
positive['confirmed_duplicate_eqtl_fit'] = positive.n_duplicate_ids.gt(0)
positive['confirmed_duplicate_export'] = positive.duplicate_pair_variant_rows.gt(0)
positive.to_csv(OUT / 'confirmed_positive_identity_checks.tsv', sep='\t', index=False)
fm = ROOT / 'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z'
tier = pd.read_csv(fm / 'gene_level_coloc_tier12.csv')
assert {'ensembl', 'coloc_best_susie_pp4'}.issubset(tier.columns), list(tier.columns)
selected = tier[tier.coloc_best_susie_pp4.gt(.5)].copy()
assert len(selected) == 436
selected = selected.merge(fits, on='ensembl', how='left', validate='one_to_one', suffixes=('', '_fit'))
selected['confirmed_duplicate_eqtl_fit'] = selected.n_duplicate_ids.gt(0)
selected.to_csv(OUT / 'primary_436_identity_checks.tsv', sep='\t', index=False)
summary = dict(saved_eqtl_fits=len(fits), duplicate_eqtl_fits=int(fits.n_duplicate_ids.gt(0).sum()),
    positive_gene_studies=len(positive), positive_genes=int(positive.ensembl.nunique()),
    positive_gene_studies_with_duplicate_eqtl_fit=int(positive.confirmed_duplicate_eqtl_fit.sum()),
    positive_genes_with_duplicate_eqtl_fit=int(positive.loc[positive.confirmed_duplicate_eqtl_fit, 'ensembl'].nunique()),
    primary_genes=len(selected), primary_genes_with_duplicate_eqtl_fit=int(selected.confirmed_duplicate_eqtl_fit.sum()),
    primary_genes_missing_eqtl_fit=int(selected.n_variants.isna().sum()),
    posterior_exports_checked=len(export_checks), duplicate_posterior_exports=int(export_checks.duplicate_pair_variant_rows.gt(0).sum()),
    positive_records_with_export=int(positive.exported_rows.notna().sum()),
    positive_records_with_duplicate_export=int(positive.confirmed_duplicate_export.sum()),
    interpretation='Duplicate fitted eQTL columns invalidate the existing signal interpretation; they do not prove a gene is biologically false. Reference-only collisions are not confirmed GWAS-fit defects.')
(OUT / 'confirmed_identity_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
