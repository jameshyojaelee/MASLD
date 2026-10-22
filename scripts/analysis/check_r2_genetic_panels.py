"""Check candidate genetic panels against r2 and compare actual rendered pages."""
import collections
import csv
import hashlib
import json
import pathlib
import subprocess
import sys
from PIL import Image, ImageChops

out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
c2 = pathlib.Path('figures/candidates/fig2-coloc-r2-20260924T053244Z')
m2 = pathlib.Path('figures/main/fig2_genetics/panels')
files = ('fig2A_gwas_cascade.pdf', 'Fig2B_PIP_vs_SuSiE-coloc.pdf',
         'Fig2C_coding_noncoding_split.pdf', 'Fig2E_ancestry_unique_coloc_GWS.pdf',
         'Fig2F_crossancestry_coloc.pdf', 'Fig2G_RORA_NAFLD.pdf',
         'Fig2H_FABP1_ALT.pdf', 'Fig2I_locus_legend.pdf')
pairs = [(c2 / name, m2 / name) for name in files]
pairs += [(pathlib.Path('figures/candidates/resource-corrections-20261005-v1/figure1/fig1c_genetic_trait_scope.pdf'),
           pathlib.Path('figures/main/fig1_atlas_overview/fig1c_genetic_trait_scope.pdf')),
          (pathlib.Path('figures/candidates/fig5d-coloc-r2-20260924T053244Z/figure5/panels/fig5d_snatac_accessibility.pdf'),
           pathlib.Path('figures/main/fig5_molecular_context/panels/fig5d_snatac_accessibility.pdf'))]
results = []
for candidate, main in pairs:
    images = []
    for label, path in [('candidate', candidate), ('main', main)]:
        prefix = out / (candidate.stem + '_' + label)
        subprocess.run(['pdftoppm', '-r', '150', '-singlefile', '-png', str(path), str(prefix)], check=True)
        text = subprocess.check_output(['pdftotext', '-layout', str(path), '-']).decode()
        prefix.with_suffix('.txt').write_text(text)
        images.append(Image.open(str(prefix) + '.png').convert('RGB'))
    different = images[0].size != images[1].size or ImageChops.difference(*images).getbbox() is not None
    results.append(dict(panel=candidate.name, candidate=str(candidate), main=str(main), rendered_difference=different))

tiers = {row['study_name']: row['tier_label'] for row in csv.DictReader(open('GWAS/finemapping/config/gwas_trait_tier.tsv'), delimiter='\t') if row['placement'] == 'main' and row['tier'] in ('1', '2')}
supports = collections.defaultdict(set)
with open('GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv') as handle:
    for row in csv.DictReader(handle):
        try:
            positive = float(row['PP.H4.susie']) > .5
        except ValueError:
            positive = False
        if positive and row['gwas_name'] in tiers:
            supports[row['ensembl']].add(tiers[row['gwas_name']])
counts = collections.Counter('both' if len(scope) == 2 else 'direct_disease' if 'direct_MASLD' in scope else 'enzyme' for scope in supports.values())
assert counts == {'enzyme': 414, 'both': 10, 'direct_disease': 12}, counts

hash_checks = []
for file in [c2 / 'input_sha256.tsv', pathlib.Path('figures/candidates/fig5d-coloc-r2-20260924T053244Z/input_output_sha256.tsv')]:
    for line in file.read_text().splitlines():
        expected, path = line.split(maxsplit=1)
        if len(expected) != 64:
            continue
        digest = hashlib.sha256()
        with open(path, 'rb') as handle:
            for block in iter(lambda: handle.read(1024*1024), b''):
                digest.update(block)
        assert digest.hexdigest() == expected, path
        hash_checks.append(path)

(out / 'panel_comparison.json').write_text(json.dumps(dict(panels=results, primary_trait_scope=dict(counts), checked_hashes=hash_checks), indent=2))
(out / 'COMPLETE.txt').write_text('Ten candidate/main page comparisons; Figure1C native SuSiE partition414/10/12; recorded candidate input/output hashes match.\n')
print(json.dumps(results, indent=2))
