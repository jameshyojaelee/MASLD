"""Recheck source-native r2 examples and eQTL target gene biotypes."""
import collections
import csv
import hashlib
import json
import pathlib
import sys

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

OUT = pathlib.Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=False)
R2 = pathlib.Path('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z')
LEADS = pathlib.Path('data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv')
GENES = {'GNMT': 'ENSG00000124713', 'MAT1A': 'ENSG00000151224',
         'CYP2C19': 'ENSG00000165841', 'IGFBP7': 'ENSG00000163453',
         'BICC1': 'ENSG00000122870'}
input_hashes = {}


def sha(path):
    value = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    input_hashes[str(path)] = value.hexdigest()


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float('nan')


def write(name, rows):
    with (OUT / name).open('w') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


primary = {r['ensembl']: r for r in csv.DictReader((R2 / 'gene_level_coloc_tier12.csv').open())
           if number(r['coloc_best_susie_pp4']) > .5}
assert len(primary) == 436
types = collections.defaultdict(set)
type_sources = {}
for row in csv.DictReader(LEADS.open(), delimiter='\t'):
    if row['Ensembl'] in primary:
        types[row['Ensembl']].add(row['Gene_BioType'])
        type_sources[row['Ensembl']] = str(LEADS)
missing = set(primary) - types.keys()
tiers = {r['study_name'] for r in csv.DictReader(open('GWAS/finemapping/config/gwas_trait_tier.tsv'), delimiter='\t')
         if r['placement'] == 'main' and r['tier'] in ('1', '2')}
chr_for_missing = {}
examples = []
with (R2 / 'susie_coloc_all_gwas.csv').open() as handle:
    for row in csv.DictReader(handle):
        if row['ensembl'] in missing:
            chr_for_missing[row['ensembl']] = row['chr']
        if row['ensembl'] in GENES.values() and row['gwas_name'] in tiers:
            examples.append(row)
named = []
for symbol, gene in GENES.items():
    candidates = [r for r in examples if r['ensembl'] == gene]
    for model, field in [('ABF', 'PP.H4.abf'), ('SuSiE', 'PP.H4.susie')]:
        finite = [r for r in candidates if number(r[field]) == number(r[field])]
        if finite:
            best = max(finite, key=lambda r: number(r[field]))
            named.append(dict(gene=symbol, ensembl=gene, model=model,
                              study=best['gwas_name'], pp_h4=best[field],
                              native_method=best['method'], shared_variant=best['top_snp'],
                              shared_variant_pp=best['top_snp_PP'],
                              status='positive' if number(best[field]) > .5 else 'tested_below_threshold',
                              source=str(R2 / 'susie_coloc_all_gwas.csv')))
        else:
            named.append(dict(gene=symbol, ensembl=gene, model=model, study='', pp_h4='',
                              native_method='', shared_variant='', shared_variant_pp='',
                              status='no_result', source=str(R2 / 'susie_coloc_all_gwas.csv')))
write('recurring_example_exact_pairs.tsv', named)
gnmt = next(r for r in named if r['gene'] == 'GNMT' and r['model'] == 'ABF')
assert gnmt['study'] == '2021_34841290_NAFLD_EUR'
assert gnmt['shared_variant'] == '6:42931261'
assert abs(float(gnmt['pp_h4']) - .705831502782787) < 1e-12

# Gene-level eQTL biotypes are source-native. They do not identify an isoform.
positions = {42928758: 'rs2296805', 42931261: 'rs2296804',
             42920549: 'rs9471976', 42928017: 'rs11752813'}
eqtl = {}
for chromosome in sorted(set(chr_for_missing.values()) | {'6'}, key=int):
    path = pathlib.Path(f'data/broadaway_eqtl/chr{chromosome}_marginal_summary_results.tsv')
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['ENSG'] in missing:
                types[row['ENSG']].add(row['Gene_Biotype'])
                type_sources[row['ENSG']] = str(path)
            if chromosome == '6' and row['ENSG'] == GENES['GNMT'] and int(row['POS']) in positions:
                position = int(row['POS'])
                assert position not in eqtl, 'Ambiguous native eQTL rows'
                eqtl[position] = row
    sha(path)
assert set(types) == set(primary), 'A primary gene has no source-native biotype'
assert all(len(value) == 1 for value in types.values()), 'Conflicting source biotypes'
biotypes = [dict(ensembl=gene, biotype=next(iter(types[gene])),
                 source=type_sources[gene], eligibility='primary_SuSiE_PP_H4_gt_0.5',
                 isoform_identified='false', selected_coloc_release=R2.name)
            for gene in sorted(primary)]
write('primary_eqtl_target_gene_biotypes.tsv', biotypes)
counts = collections.Counter(r['biotype'] for r in biotypes)
write('target_gene_biotype_counts.tsv', [dict(biotype=key, genes=value, denominator=436)
                                       for key, value in sorted(counts.items())])
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 6,
                     'axes.labelsize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
labels = sorted(counts, key=lambda key: counts[key])
fig, ax = plt.subplots(figsize=(4, 2))
ax.barh(labels, [counts[x] for x in labels], height=.55,
        color=['#1565C0' if x == 'protein_coding' else '#7B1FA2' if x in ('lncRNA', 'lincRNA') else '#9E9E9E' for x in labels])
for index, key in enumerate(labels):
    ax.text(counts[key] + 4, index, str(counts[key]), va='center')
ax.set_xlim(0, max(counts.values()) * 1.13)
ax.set_xlabel('Primary SuSiE genes (n = 436)')
fig.tight_layout()
fig.savefig(OUT / 'Fig2_companion_eqtl_target_gene_biotype.pdf')
plt.close(fig)

# Re-read exact native GWAS effects; nonpalindromic proxies support signed alignment.
# for the signed GNMT bridge; the C/G secondary marker retains a strand caveat.
gwas_path = pathlib.Path('GWAS/finemapping/data/sumstats_af/2021_34841290_NAFLD_EUR_af.tsv')
gwas = {}
with gwas_path.open() as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        if row['chromosome'] == '6' and int(row['position']) in positions:
            position = int(row['position'])
            assert position not in gwas
            gwas[position] = row
assert set(eqtl) == set(positions) == set(gwas)
block = pathlib.Path('GWAS/finemapping/data/ld_ref/polyfun_eur/chr6/42038721.43756169/42038721.43756169')
bim = list(csv.reader(pathlib.Path(str(block) + '.bim').open(), delimiter='\t'))
indexes = {int(row[3]): index for index, row in enumerate(bim) if int(row[3]) in positions}
assert set(indexes) == set(positions)
with pathlib.Path(str(block) + '.ld').open() as handle:
    for index, line in enumerate(handle):
        if index == indexes[42928758]:
            correlation = list(map(float, line.split()))
            break
assert len(correlation) == len(bim)
chain = []
complement = str.maketrans('ACGT', 'TGCA')
for position, rsid in positions.items():
    first, second = gwas[position], eqtl[position]
    if second['EA'] == first['allele1']:
        risk_beta = float(first['beta'])
    elif second['EA'] == first['allele2']:
        risk_beta = -float(first['beta'])
    else:
        raise ValueError('Unaligned native alleles')
    palindromic = first['allele1'] == first['allele2'].translate(complement)
    chain.append(dict(rsid=rsid, chromosome=6, position_grch37=position,
                      is_current_coloc_anchor=(position == 42931261),
                      source_role='current_ABF_shared_variant' if position == 42931261 else 'nonpalindromic_LD_proxy' if position == 42928758 else 'source_reported_functional_variant',
                      gwas_study='2021_34841290_NAFLD_EUR',
                      gwas_effect_allele=first['allele1'], gwas_other_allele=first['allele2'],
                      gwas_beta=first['beta'], gwas_se=first['se'], gwas_nominal_p=first['pval'],
                      gwas_effect_af=first['af'], eqtl_effect_allele=second['EA'],
                      eqtl_other_allele=second['NEA'], eqtl_beta=second['Beta'],
                      eqtl_se=second['SE'], eqtl_nominal_p=second['PVAL'], eqtl_effect_af=second['EAF'],
                      eqtl_n=second['N'], gwas_beta_aligned_to_eqtl_effect_allele=risk_beta,
                      lower_expression_higher_risk='' if palindromic else (risk_beta * float(second['Beta']) < 0),
                      reported_label_risk_expression_product=(risk_beta * float(second['Beta'])),
                      palindromic=palindromic,
                      sign_qualification='reported_alleles_only_palindromic' if palindromic else 'nonpalindromic_native_alignment',
                      ukbb_eur_r2_to_rs2296805=correlation[indexes[position]] ** 2,
                      gnmt_direct_abf_pp_h4=gnmt['pp_h4'], gnmt_direct_abf_anchor=gnmt['shared_variant'],
                      published_gnmt_protein_association=rsid in ('rs9471976', 'rs11752813'),
                      published_liver_biopsy_n=268 if rsid in ('rs9471976', 'rs11752813') else '',
                      external_source='https://pmc.ncbi.nlm.nih.gov/articles/PMC3463826/' if rsid in ('rs9471976', 'rs11752813') else '',
                      gwas_source=str(gwas_path), eqtl_source='data/broadaway_eqtl/chr6_marginal_summary_results.tsv',
                      ld_source=str(pathlib.Path(str(block) + '.ld')), coloc_source=gnmt['source']))
assert chain[0]['lower_expression_higher_risk']
write('gnmt_risk_eqtl_ld_source_chain.tsv', chain)
for path in [R2 / 'susie_coloc_all_gwas.csv', R2 / 'gene_level_coloc_tier12.csv',
             LEADS, gwas_path, pathlib.Path(str(block) + '.bim'), pathlib.Path(str(block) + '.ld')]:
    sha(path)
(OUT / 'input_sha256.json').write_text(json.dumps(input_hashes, indent=2))
(OUT / 'COMPLETE.txt').write_text('Five example model states and exact pairs checked;436 source-native target gene biotypes; native GNMT risk/eQTL/UKBB-EUR LD chain extracted. Palindromic secondary-marker signs remain qualified; no causal gene or isoform claim.\n')
print(json.dumps(dict(target_gene_biotypes=dict(counts), recurring_examples=named, gnmt_chain=chain), indent=2))
