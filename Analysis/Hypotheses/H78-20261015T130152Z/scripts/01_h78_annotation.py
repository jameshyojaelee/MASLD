"""H78: per-gene constraint, essentiality and secretion fields for the Gene Catalog.

Implements PRESPEC_H78.md (sha256 faaf7f41...c42855) exactly:
  S1 input hashes -> G1 counts (S2) -> G2 keys (S3) -> G3 release -> join yield ->
  PC-1..PC-3 / NC-1..NC-4 (S4-S6) -> D2, D3 -> spot check (S3) -> table -> D1, D4, D5 ->
  RESULT_H78.md with the verdict computed by the prespec's Reading rule.
No hypothesis test and no p value are computed anywhere (family size 0).

--debug-null replaces the Hart lists and the PC-3 genes with random draws and permutes the route
labels used by D4/D5, so the control and description code can be exercised without computing any
real control value. It is never used for the final run.
"""
import argparse
import bisect
import gzip
import hashlib
import json
import os
import re
import statistics
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--out', required=True)
ap.add_argument('--external', required=True)
ap.add_argument('--seed', type=int, default=20260923)
ap.add_argument('--debug-null', action='store_true')
ap.add_argument('--debug-force-ship', action='store_true', help='debug-null only: ignore S4-S6 so D4/D5 code runs')
args = ap.parse_args()
np.random.seed(args.seed)
RNG = np.random.default_rng(args.seed)
DEBUG = args.debug_null
OUT = args.out
os.makedirs(OUT, exist_ok=True)
pd.set_option('display.width', 250)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 60)

R = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
H78 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RR = R + 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/'
R4 = RR + 'round4-20261015T112148Z/'
W = R4 + 'write/H78/'
D24 = R + 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/myojin_hlf/source/depmap_24q4/'
EXT = args.external.rstrip('/') + '/'

PRESPEC_SHA = 'faaf7f419e797d28aa6f939ba31fa23407902d13c9a9b4ce9fa1c38c1cf42855'
PRESPEC_PATHS = [RR + 'plans/PRESPEC_H78.md', H78 + '/PRESPEC_H78.md']

CAT = R + 'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet'
GN = R + 'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv'
HPA = R + 'data/external/hpa/proteinatlas.tsv'
GE = D24 + 'CRISPRGeneEffect.csv'
MODEL = D24 + 'Model.csv'
CE_SRC = W + 'depmap24q4/CRISPRInferredCommonEssentials.csv'
CTL_SRC = W + 'depmap24q4/AchillesCommonEssentialControls.csv'
CE = EXT + 'CRISPRInferredCommonEssentials.csv'
CTL = EXT + 'AchillesCommonEssentialControls.csv'
CEG = R + 'data/public_r4/Hart2017_CEG2/CEGv2.txt'
NEG = R + 'data/public_r4/Hart2017_CEG2/NEGv1.txt'
GTF = '/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz'

# Inputs table of the prespec (role, path, sha256). The two DepMap lists are checked at the source
# path and at the run's external/ copy, which is the copy that is read.
INPUTS = [
    ('catalog_r8_gene_index', CAT, 'ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b'),
    ('gnomad_v41_constraint', GN, '68d8abdb7fc48f570869b02dfaa74b9fecaece7fcc5f301ddca40ec1ce12da00'),
    ('hpa_proteinatlas', HPA, '5c546331ec5dda4bcec615c1c16eadc5f440e7e7a1fbb3a178f8ecaf983b0c17'),
    ('depmap24q4_gene_effect', GE, '3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b'),
    ('depmap24q4_model', MODEL, 'b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa'),
    ('depmap24q4_inferred_common_essentials_source', CE_SRC, '0352de5335ae3bf72250dd136d1b0df87a6b72f653df34c2560cce17c0e7d928'),
    ('depmap24q4_inferred_common_essentials_run_copy', CE, '0352de5335ae3bf72250dd136d1b0df87a6b72f653df34c2560cce17c0e7d928'),
    ('depmap24q4_achilles_ce_controls_source', CTL_SRC, 'a1a48ad516011c521b1981769e2f78de84d1410371f8517a7428349be6e6b98a'),
    ('depmap24q4_achilles_ce_controls_run_copy', CTL, 'a1a48ad516011c521b1981769e2f78de84d1410371f8517a7428349be6e6b98a'),
    ('hart_cegv2', CEG, '1b1df46f829e2c47dc8a4c164f60f5f1dc67513bafba7d8b9609421e0c75d27e'),
    ('hart_negv1', NEG, '18a6b37a8f182a1bf501548e6b28fe4fcef02fcd296f873ef500a21eeb433739'),
    ('gencode_v49_gtf', GTF, '73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9'),
    ('gate_count_script', R4 + 'gate/NG17/count_ng17.py', 'd07895ff7a91652bb1a5d2a111d8e1cc3ae09c646a48099728245360a983ee3e'),
    ('gate_log', R4 + 'gate/NG17/count_22251965.log', '9df93570eed766a32d940243d54843187a96538f4df259faa1aa0d29bd75b220'),
    ('h78_count_script', W + 'count_h78.py', 'da33c1f7b69d33270ad39d13664be253647c06a3c054d622979062a688c02821'),
    ('h78_count_log', W + 'count_22252308.log', '690b964672f0ddd0e20ad6deebd58ac624ff8ba374f4b853601407cf3df8ef28'),
    ('h78_depmap24q4_count_script', W + 'count_depmap24q4.py', 'c3d94eede113611b4874a0dc56bf94b2f40790b16433f14e3364c95a003171cd'),
    ('h78_depmap24q4_count_log', W + 'count24q4_22252330.log', '32046ba7ed307b9bb572161258f586014c9d20be72f7962802774d6018d048ab'),
]

ROUTES = ['EXP_ALLELE_AWARE_V1', 'EXP_DIRECT_ASSAY_V1', 'EXP_MULTI_LAYER_V1', 'EXP_STRESS_PERTURB_V1', 'EXP_UNRESOLVED_V1']
STRESS = 'EXP_STRESS_PERTURB_V1'
PRIMARY = set([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])
PC3_GENES = ['ALB', 'APOA1', 'APOB', 'FGA', 'FGB', 'FGG', 'SERPINA1', 'TF', 'HP', 'AHSG']
SPOT_GENES = ['PNPLA3', 'TM6SF2', 'GCKR', 'HSD17B13', 'MBOAT7', 'APOB', 'ALB', 'RPL3', 'POLR2A', 'MYC']
HPA_MS = 'Blood concentration - Conc. blood MS [pg/L]'
HPA_IM = 'Blood concentration - Conc. blood IM [pg/L]'
HPA_SEC = 'Secretome location'

DIAG = []  # rows: order, section, item, value, expected, rule, status


def log(*a):
    print(*a, flush=True)


def diag(section, item, value, expected='', rule='', status=''):
    DIAG.append({'order': len(DIAG) + 1, 'section': section, 'item': item, 'value': value,
                 'expected': expected, 'rule': rule, 'status': status})
    log(f'[{section}] {item}: {value}' + (f' | expected {expected}' if expected != '' else '') +
        (f' | rule {rule}' if rule else '') + (f' | {status}' if status else ''))


def write_diag():
    pd.DataFrame(DIAG).to_csv(os.path.join(OUT, 'diagnostics.tsv'), sep='\t', index=False)


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def write_result(verdict, rule_line, body):
    name = 'RESULT_H78_DEBUGNULL.md' if DEBUG else 'RESULT_H78.md'
    with open(os.path.join(OUT, name), 'w') as f:
        f.write(f'# H78 result: Gene Catalog follow-up annotation (gnomAD v4.1, DepMap 24Q4, HPA)\n\n')
        f.write(f'Written by `scripts/01_h78_annotation.py` at {datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}.'
                f' Prespec sha256 `{PRESPEC_SHA}`. Seed {args.seed}.' + (' DEBUG-NULL RUN: control values are on random sets and mean nothing.' if DEBUG else '') + '\n\n')
        f.write(f'**Verdict (computed by the prespec Reading rule): {verdict}**\n\n{rule_line}\n\n')
        f.write(body)
    log(f'VERDICT: {verdict}')


def stop(rule, why):
    diag('stop', rule, why, status='STOP')
    write_diag()
    write_result(f'stopped ({rule})', f'Stop rule {rule} fired: {why}. Nothing ships; the Catalog is unchanged; the paper says nothing.',
                 'See `diagnostics.tsv` for every check run before the stop.\n')
    sys.exit(0)


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')
    if r.returncode not in (0, 1):  # grep returns 1 on no match
        raise RuntimeError(f'command failed ({r.returncode}): {cmd}\n{r.stderr}')
    return r.stdout


log(f'H78 run start {datetime.now(timezone.utc).isoformat()} out={OUT} debug_null={DEBUG} seed={args.seed}')

# ---------------------------------------------------------------- 1. hashes (S1)
for p in PRESPEC_PATHS:
    s = sha256(p)
    diag('1_sha256', f'prespec {p.replace(R, "")}', s, PRESPEC_SHA, 'equal', 'pass' if s == PRESPEC_SHA else 'FAIL')
    if s != PRESPEC_SHA:
        stop('S1', f'prespec sha256 differs at {p}')
rows = []
bad = []
for role, p, exp in INPUTS:
    s = sha256(p)
    ok = s == exp
    rows.append({'role': role, 'path': p, 'bytes': os.path.getsize(p), 'sha256': s, 'expected': exp, 'match': ok})
    diag('1_sha256', role, s, exp, 'equal', 'pass' if ok else 'FAIL')
    if not ok:
        bad.append(role)
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'input_sha256.tsv'), sep='\t', index=False)
if bad:
    stop('S1', 'sha256 differs for ' + ', '.join(bad))
diag('1_sha256', 'S1', 'all inputs match', status='pass')
diag('1_sha256', 'hpa_release', 'proteinatlas.tsv states no release number; identifier = file date 2026-05-25 (zip member timestamp) + sha256 above [release number Unverified]')

# ---------------------------------------------------------------- load
g = pd.read_parquet(CAT, columns=['ensembl_id', 'symbol', 'primary_evidence_class', 'next_experiment_rule_id'])
g['ens'] = g.ensembl_id.astype(str).str.split('.').str[0]

gn = pd.read_csv(GN, sep='\t', usecols=['gene', 'gene_id', 'transcript', 'mane_select', 'lof.oe_ci.upper', 'lof.pLI',
                                        'lof.oe_ci.upper_bin_decile', 'constraint_flags', 'syn.oe'], low_memory=False)
mane = gn[gn.mane_select.astype(str).str.lower().eq('true')]
ensm = mane[mane.gene_id.astype(str).str.startswith('ENSG')].copy()
ezm = mane[mane.gene_id.astype(str).str.fullmatch(r'\d+')].copy()
L = ensm.dropna(subset=['lof.oe_ci.upper']).sort_values('gene_id').copy()
L['dec'] = pd.qcut(L['lof.oe_ci.upper'].rank(method='first'), 10, labels=False)
pair = ensm[['gene_id', 'gene']].merge(ezm[['gene', 'gene_id']].rename(columns={'gene_id': 'entrez'}), on='gene')
pair['entrez'] = pair.entrez.astype(int)

gt_rows = []
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.split('\t', 9)
        if p[2] != 'gene':
            continue
        a = p[8]
        gid = re.search(r'gene_id "([^"]+)"', a).group(1)
        nm = re.search(r'gene_name "([^"]+)"', a)
        hg = re.search(r'hgnc_id "([^"]+)"', a)
        gt = re.search(r'gene_type "([^"]+)"', a)
        gt_rows.append((gid.split('.')[0], gid, nm.group(1) if nm else None, hg.group(1) if hg else None, p[0],
                        gt.group(1) if gt else None))
gc = pd.DataFrame(gt_rows, columns=['ens', 'gid_full', 'sym', 'hgnc', 'chr', 'gtype'])
gc['primary'] = gc.chr.isin(PRIMARY)
gc_any = gc[~gc.gid_full.str.endswith('_PAR_Y')].drop_duplicates('ens')
gc_prim = gc_any[gc_any.primary]
hgnc2ens = gc_prim.dropna(subset=['hgnc']).drop_duplicates('hgnc').set_index('hgnc').ens
name2ens = gc_prim.dropna(subset=['sym']).drop_duplicates('sym').set_index('sym').ens
gtype = gc_any.set_index('ens').gtype

ge_cols = pd.read_csv(GE, nrows=0).columns
dmc = pd.DataFrame({'col': ge_cols[1:]})
dmc['entrez'] = dmc.col.str.extract(r'\((\d+)\)$')[0].astype(int)
ge = pd.read_csv(GE, index_col=0)
ge_med = ge.median(axis=0, skipna=True)
ge_n = ge.notna().sum(axis=0)
dmc['median'] = dmc.col.map(ge_med)
dmc['n_lines'] = dmc.col.map(ge_n).astype(int)
ce = pd.read_csv(CE)
ce['entrez'] = ce.iloc[:, 0].str.extract(r'\((\d+)\)$')[0].astype(int)
ctl = pd.read_csv(CTL)
ctl['entrez'] = ctl.iloc[:, 0].str.extract(r'\((\d+)\)$')[0].astype(int)
ce_set = set(ce.entrez)
dmc['ce'] = dmc.entrez.isin(ce_set)

h = pd.read_csv(HPA, sep='\t', usecols=['Gene', 'Ensembl', HPA_SEC, HPA_IM, HPA_MS])


def secreted_level(txt):
    if pd.isna(txt) or str(txt).strip() == '':
        return 'not_annotated'
    t = str(txt)
    if t.startswith('Secreted'):
        return 'secreted'
    if t == 'Intracellular and membrane':
        return 'intracellular_and_membrane'
    if t == 'Immunoglobulin genes':
        return 'immunoglobulin'
    return 'UNMAPPED'


h['level'] = h[HPA_SEC].map(secreted_level)

# ---------------------------------------------------------------- build table in memory
t = g[['ens', 'symbol', 'next_experiment_rule_id']].copy()
ensm_i = ensm.set_index('gene_id')
t['in_mane'] = t.ens.isin(ensm_i.index)
t['gnomad_v41_mane_transcript'] = t.ens.map(ensm_i.transcript)
t['gnomad_v41_loeuf'] = t.ens.map(ensm_i['lof.oe_ci.upper'])
t['gnomad_v41_pli'] = t.ens.map(ensm_i['lof.pLI'])
t['gnomad_v41_loeuf_decile_mane_ens'] = t.ens.map(L.set_index('gene_id').dec).astype('Int64')
t['gnomad_v41_loeuf_decile_official'] = t.ens.map(ensm_i['lof.oe_ci.upper_bin_decile']).astype('Int64')
flags = ensm_i.constraint_flags.astype(str).str.strip()
flags = flags.where(~flags.isin(['[]', '', 'nan', 'None']), '')
t['gnomad_v41_constraint_flags'] = t.ens.map(flags).fillna('')
t['entrez_id'] = t.ens.map(pair.set_index('gene_id').entrez).astype('Int64')
dmi = dmc.set_index('entrez')
t['in_dm'] = t.entrez_id.isin(dmi.index).astype(bool)
t['depmap24q4_common_essential'] = np.where(t.in_dm, np.where(t.entrez_id.isin(ce_set), 'TRUE', 'FALSE'), '')
t['depmap24q4_median_gene_effect'] = t.entrez_id.map(dmi['median']).astype(float)
t['depmap24q4_n_lines'] = t.entrez_id.map(dmi['n_lines']).astype('Int64')
hi = h.set_index('Ensembl')
t['in_hpa'] = t.ens.isin(hi.index)
t['hpa_secretome_location'] = t.ens.map(hi[HPA_SEC]).fillna('')
t['hpa_secreted'] = np.where(t.in_hpa, t.ens.map(hi.level), '')
t['hpa_blood_ms_pg_per_l'] = t.ens.map(hi[HPA_MS]).astype(float)
t['hpa_blood_im_pg_per_l'] = t.ens.map(hi[HPA_IM]).astype(float)
t['gtype'] = t.ens.map(gtype)
t['has_loeuf'] = t.gnomad_v41_loeuf.notna()
r_gn = np.where(t.has_loeuf, '', np.where(t.in_mane, 'gnomad_mane_no_loeuf',
                np.where(t.gtype.eq('protein_coding'), 'protein_coding_no_gnomad_mane_row', 'noncoding_gene')))
r_dm = np.where(t.entrez_id.isna(), 'no_entrez_partner', np.where(~t.in_dm, 'not_in_depmap24q4', ''))
r_hp = np.where(~t.in_hpa, 'not_in_hpa', '')
t['r_gn'], t['r_dm'], t['r_hp'] = r_gn, r_dm, r_hp

# ---------------------------------------------------------------- 2. G1 counts (S2), G2 keys (S3), G3 release
lev_atlas = h.loc[h[HPA_SEC].notna() & (h[HPA_SEC].astype(str).str.strip() != ''), 'level'].value_counts()
G1 = [
    ('catalog_genes', len(t), 28233),
    ('route EXP_ALLELE_AWARE_V1', int((t.next_experiment_rule_id == 'EXP_ALLELE_AWARE_V1').sum()), 400),
    ('route EXP_DIRECT_ASSAY_V1', int((t.next_experiment_rule_id == 'EXP_DIRECT_ASSAY_V1').sum()), 13551),
    ('route EXP_MULTI_LAYER_V1', int((t.next_experiment_rule_id == 'EXP_MULTI_LAYER_V1').sum()), 26),
    ('route EXP_STRESS_PERTURB_V1', int((t.next_experiment_rule_id == STRESS).sum()), 1129),
    ('route EXP_UNRESOLVED_V1', int((t.next_experiment_rule_id == 'EXP_UNRESOLVED_V1').sum()), 13127),
    ('ensembl_mane_genes_with_loeuf', int(L.gene_id.nunique()), 17325),
    ('depmap24q4_genes', len(dmc), 17916),
    ('depmap24q4_common_essentials', len(ce), 1523),
    ('hpa_genes', int(h.Ensembl.nunique()), 20162),
    ('catalog_with_loeuf', int(t.has_loeuf.sum()), 17152),
    ('catalog_in_depmap24q4', int(t.in_dm.sum()), 16761),
    ('catalog_in_hpa', int(t.in_hpa.sum()), 18675),
    ('ensembl_mane_rows_with_entrez_partner', len(pair), 17468),
    ('catalog_with_entrez_key', int(t.entrez_id.notna().sum()), 17290),
    ('blank noncoding_gene', int((t.r_gn == 'noncoding_gene').sum()), 9579),
    ('blank protein_coding_no_gnomad_mane_row', int((t.r_gn == 'protein_coding_no_gnomad_mane_row').sum()), 1354),
    ('blank gnomad_mane_no_loeuf', int((t.r_gn == 'gnomad_mane_no_loeuf').sum()), 148),
    ('blank no_entrez_partner', int((t.r_dm == 'no_entrez_partner').sum()), 10943),
    ('blank not_in_depmap24q4', int((t.r_dm == 'not_in_depmap24q4').sum()), 529),
    ('blank not_in_hpa', int((t.r_hp == 'not_in_hpa').sum()), 9558),
    ('hpa_level secreted (atlas-wide)', int(lev_atlas.get('secreted', 0)), 1754),
    ('hpa_level intracellular_and_membrane (atlas-wide)', int(lev_atlas.get('intracellular_and_membrane', 0)), 619),
    ('hpa_level immunoglobulin (atlas-wide)', int(lev_atlas.get('immunoglobulin', 0)), 147),
]
g1_bad = []
for item, v, e in G1:
    diag('2_G1', item, v, e, 'equal', 'pass' if v == e else 'FAIL')
    if v != e:
        g1_bad.append(f'{item} {v} vs {e}')
if g1_bad:
    stop('S2', '; '.join(g1_bad))

unmapped = h.loc[h.level == 'UNMAPPED', HPA_SEC].value_counts()
diag('2_G2', 'hpa secretome texts mapping to no level (atlas-wide)', int(unmapped.sum()), 0, 'equal', 'pass' if unmapped.sum() == 0 else 'FAIL')
if unmapped.sum() > 0:
    stop('S3', 'unmapped HPA secretome text: ' + '; '.join(unmapped.index))
cat_pair = pair[pair.gene_id.isin(t.ens)]
G2 = [
    ('output rows', len(t), 28233),
    ('output unique ensembl_id', int(t.ens.nunique()), 28233),
    ('gnomAD Ensembl MANE gene_id with >1 row', int(ensm.gene_id.duplicated().sum()), 0),
    ('Catalog Ensembl ID with >1 Entrez', int(cat_pair.gene_id.duplicated().sum()), 0),
    ('Catalog Entrez ID with >1 Ensembl', int(cat_pair.entrez.duplicated().sum()), 0),
    ('Catalog Ensembl ID with >1 HPA row', int(h[h.Ensembl.isin(t.ens)].Ensembl.duplicated().sum()), 0),
    ('symbols with >1 Entrez MANE row', int(ezm.gene.duplicated().sum()), 0),
    ('DepMap Entrez with >1 column', int(dmc.entrez.duplicated().sum()), 0),
]
g2_bad = []
for item, v, e in G2:
    diag('2_G2', item, v, e, 'equal', 'pass' if v == e else 'FAIL')
    if v != e:
        g2_bad.append(f'{item} {v}')
if g2_bad:
    stop('S3', '; '.join(g2_bad))
n_ce_in = int(ce.entrez.isin(dmc.entrez).sum())
diag('2_G3', 'common-essential Entrez IDs that are 24Q4 matrix columns', n_ce_in, 1523, 'equal', 'pass' if n_ce_in == 1523 else 'FAIL')
if n_ce_in != 1523:
    stop('S1', f'G3 release check: {n_ce_in} of 1523 common essentials in matrix')
diag('2_G3', 'depmap release', 'DepMap Public 24Q4 gene effect and both list sha256 equal the prespec table', status='pass')
diag('2_G3', 'depmap24q4 lines x genes', f'{ge.shape[0]} x {ge.shape[1]}', '1178 x 17916')
diag('2_info', 'Catalog genes absent from GENCODE v49', int(t.gtype.isna().sum()), 0)
diag('2_info', 'Catalog genes on a non-primary scaffold', int((~t.ens.isin(gc_prim.ens)).sum()), 51)
diag('2_info', 'gnomAD Ensembl MANE rows with LOEUF and an Entrez partner', int(pair.gene_id.isin(L.gene_id).sum()), 17314)
pp = L[['gene_id', 'gene', 'lof.oe_ci.upper']].merge(pair[['gene_id', 'entrez']], on='gene_id').merge(
    ezm.assign(entrez=ezm.gene_id.astype(int))[['entrez', 'lof.oe_ci.upper']].rename(columns={'lof.oe_ci.upper': 'loeuf_ez'}), on='entrez')
diag('2_info', 'paired rows with identical LOEUF on Ensembl and Entrez rows', f'{int((pp["lof.oe_ci.upper"] == pp.loeuf_ez).sum())} of {len(pp)}', '17314 of 17314')
diag('2_info', 'gnomAD constraint_flags non-empty (Ensembl MANE rows)', int((flags != '').sum()), 595)

# G4: stress-route coverage (decides D4 fields only)
st = t[t.next_experiment_rule_id == STRESS]
G4 = {'gnomad': st.has_loeuf.mean(), 'depmap': st.in_dm.mean(), 'hpa': st.in_hpa.mean()}
for k, v in G4.items():
    diag('2_G4', f'stress-route coverage {k}', round(float(v), 4), '>=0.90', 'D4 field printed if >=0.90', 'pass' if v >= 0.90 else 'below 0.90: left out of D4')

# ---------------------------------------------------------------- 3. join yield
for k, v in [('catalog genes with Ensembl MANE row', int(t.in_mane.sum())), ('catalog genes with LOEUF', int(t.has_loeuf.sum())),
             ('catalog genes with official decile', int(t.gnomad_v41_loeuf_decile_official.notna().sum())),
             ('catalog genes with constraint flag', int((t.gnomad_v41_constraint_flags != '').sum())),
             ('catalog genes with Entrez key', int(t.entrez_id.notna().sum())), ('catalog genes in 24Q4 matrix', int(t.in_dm.sum())),
             ('catalog genes with an HPA row', int(t.in_hpa.sum())), ('catalog genes with secretome text', int((t.hpa_secretome_location != '').sum())),
             ('catalog genes with blood MS', int(t.hpa_blood_ms_pg_per_l.notna().sum())), ('catalog genes with blood IM', int(t.hpa_blood_im_pg_per_l.notna().sum())),
             ('inferred common essentials reachable from Catalog by Entrez', int(ce.entrez.isin(t.entrez_id.dropna().astype(int)).sum()))]:
    diag('3_join_yield', k, v)
for lev, n in t.loc[t.in_hpa, 'hpa_secreted'].value_counts().items():
    diag('3_join_yield', f'catalog hpa_secreted={lev}', int(n))
mult_all = int((gc_any.dropna(subset=['hgnc']).groupby('hgnc').ens.nunique() > 1).sum())
mult_prim_s = gc_prim.dropna(subset=['hgnc']).groupby('hgnc').ens.nunique()
mult_prim = set(mult_prim_s[mult_prim_s > 1].index)
diag('3_join_yield', 'GENCODE HGNC IDs with >1 Ensembl ID (all scaffolds)', mult_all, 3135)
diag('3_join_yield', 'GENCODE HGNC IDs with >1 Ensembl ID (primary assembly)', len(mult_prim), 242)


def load_hart(path):
    x = pd.read_csv(path, sep='\t')
    x['ENTREZ_ID'] = x.ENTREZ_ID.astype('Int64')
    return x


ceg = load_hart(CEG)
neg = load_hart(NEG)
if DEBUG:
    # random sets of the same sizes: HGNC IDs from GENCODE, Entrez IDs from the 24Q4 matrix
    hg_pool = np.array(sorted(hgnc2ens.index))
    ez_pool = np.array(sorted(dmc.entrez))
    for nm, x in [('CEGv2', ceg), ('NEGv1', neg)]:
        x['HGNC_ID'] = RNG.choice(hg_pool, len(x), replace=False)
        x['ENTREZ_ID'] = pd.array(RNG.choice(ez_pool, len(x), replace=False), dtype='Int64')
    log('DEBUG-NULL: Hart lists replaced by random HGNC/Entrez draws')

Lset = set(L.gene_id)
hart = {}
for nm, x in [('CEGv2', ceg), ('NEGv1', neg)]:
    x = x.copy()
    x['ens'] = x.HGNC_ID.map(hgnc2ens)
    n0 = len(x)
    n_ens = int(x.ens.notna().sum())
    n_L = int(x.ens.isin(Lset).sum())
    n_ez = int(x.ENTREZ_ID.notna().sum())
    n_dm = int(x.ENTREZ_ID.isin(dmc.entrez).sum())
    n_cat = int(x.ens.isin(t.ens).sum())
    n_multi = int(x.HGNC_ID.isin(mult_prim).sum())
    exp = {'CEGv2': (684, 681, 653, 674, 0), 'NEGv1': (927, 919, 723, 803, 3)}[nm]
    diag('3_join_yield', f'{nm} genes', n0, exp[0])
    diag('3_join_yield', f'{nm} HGNC->Ensembl (primary assembly, first line) [lost]', f'{n_ens} [{n0 - n_ens}]', exp[1])
    diag('3_join_yield', f'{nm} Ensembl->Ensembl MANE with LOEUF [lost]', f'{n_L} [{n_ens - n_L}]', exp[2])
    diag('3_join_yield', f'{nm} Entrez non-blank [blank]', f'{n_ez} [{n0 - n_ez}]', 'blank 3')
    diag('3_join_yield', f'{nm} Entrez->24Q4 matrix [lost]', f'{n_dm} [{n0 - n_dm}]', exp[3])
    diag('3_join_yield', f'{nm} Ensembl in Catalog', n_cat)
    diag('3_join_yield', f'{nm} HGNC IDs with two primary-assembly Ensembl IDs', n_multi, exp[4])
    hart[nm] = x

# ---------------------------------------------------------------- 4. controls
pc3_genes = list(PC3_GENES)
if DEBUG:
    pc3_genes = list(RNG.choice(np.array(sorted(t.symbol[t.in_hpa & t.symbol.isin(name2ens.index)])), 10, replace=False))
    log('DEBUG-NULL: PC-3 genes replaced by', pc3_genes)
Li = L.set_index('gene_id')
all_med = float(L['lof.oe_ci.upper'].median())
all_syn = float(L['syn.oe'].median())
CTRL = []


def ctrl(cid, source, value, n, predicted, rule, passed, note=''):
    CTRL.append({'control': cid, 'source': source, 'value': value, 'n': n, 'predicted': predicted, 'pass_rule': rule,
                 'result': 'pass' if passed else 'FAIL', 'note': note})
    diag('4_controls', cid, value, predicted, rule, ('pass' if passed else 'FAIL') + (f' ({note})' if note else ''))


ceg_L = Li.loc[sorted(set(hart['CEGv2'].ens.dropna()) & Lset)]
neg_L = Li.loc[sorted(set(hart['NEGv1'].ens.dropna()) & Lset)]
pc1_frac = float((ceg_L.dec == 0).mean())
pc1_med = float(ceg_L['lof.oe_ci.upper'].median())
ctrl('PC-1a', 'gnomad', round(pc1_frac, 4), len(ceg_L), '0.25 to 0.33', '>= 0.20', pc1_frac >= 0.20, 'CEGv2 bottom-decile fraction (recomputed decile)')
ctrl('PC-1b', 'gnomad', round(pc1_med, 4), len(ceg_L), '0.65 to 0.75', f'< all-gene median {all_med:.4f}', pc1_med < all_med, 'CEGv2 median LOEUF')
cegsyn = float(ceg_L['syn.oe'].median())
nc1 = abs(cegsyn - all_syn)
ctrl('NC-1', 'gnomad', round(nc1, 4), int(ceg_L['syn.oe'].notna().sum()), 'within 0.03', '<= 0.05', nc1 <= 0.05,
     f'CEGv2 syn.oe median {cegsyn:.4f} vs all {all_syn:.4f} (n all {int(L["syn.oe"].notna().sum())})')
nc2 = float((neg_L.dec == 0).mean())
ctrl('NC-2', 'gnomad', round(nc2, 4), len(neg_L), '0.01 to 0.04', '<= 0.06', nc2 <= 0.06, 'NEGv1 bottom-decile fraction')
ceg_dm = hart['CEGv2'][hart['CEGv2'].ENTREZ_ID.isin(dmc.entrez)]
neg_dm = hart['NEGv1'][hart['NEGv1'].ENTREZ_ID.isin(dmc.entrez)]
pc2 = float(ceg_dm.ENTREZ_ID.isin(ce_set).mean())
ctrl('PC-2', 'depmap', round(pc2, 4), len(ceg_dm), '0.80 to 0.95', '>= 0.70', pc2 >= 0.70,
     'CEGv2 flagged common-essential; checks join and release, not DepMap accuracy (489 CEGv2 genes are DepMap controls)')
nc4 = float(neg_dm.ENTREZ_ID.isin(ce_set).mean())
ctrl('NC-4', 'depmap', round(nc4, 4), len(neg_dm), '0 to 0.005', '<= 0.02', nc4 <= 0.02, 'NEGv1 flagged common-essential')
ctl_in = int(ctl.entrez.isin(ce_set).sum())
diag('4_controls', 'Achilles common-essential controls in inferred list (printed, not a rule)', f'{ctl_in} / {len(ctl)} = {ctl_in / len(ctl):.4f}', '973 / 1247 = 0.78')
ti = t.set_index('ens')
pc3_ens = [name2ens.get(s) for s in pc3_genes]
pc3_rows = ti.reindex([e for e in pc3_ens if e is not None])
n_blood = int((pc3_rows.hpa_secretome_location == 'Secreted to blood').sum())
n_ms = int(pc3_rows.hpa_blood_ms_pg_per_l.notna().sum())
diag('4_controls', 'PC-3 genes mapped to a Catalog row', f'{len(pc3_rows)} of 10 ({",".join(f"{s}={e}" for s, e in zip(pc3_genes, pc3_ens))})', '10 of 10')
ctrl('PC-3a', 'hpa', f'{n_blood} of 10', 10, '>= 9 of 10', '>= 8 of 10', n_blood >= 8 and len(pc3_rows) == 10, 'Secreted to blood')
ctrl('PC-3b', 'hpa', f'{n_ms} of 10', 10, '>= 9 of 10', '>= 8 of 10', n_ms >= 8 and len(pc3_rows) == 10, 'non-blank blood MS concentration')
hlev = h.set_index('Ensembl').level
ceg_h = hart['CEGv2'].ens.dropna()
ceg_h = ceg_h[ceg_h.isin(hlev.index)]
nc3 = float((hlev.reindex(ceg_h) == 'secreted').mean())
ctrl('NC-3', 'hpa', round(nc3, 4), len(ceg_h), '<= 0.03', '<= 0.08', nc3 <= 0.08, 'CEGv2 genes with an HPA row whose hpa_secreted is secreted')
cdf = pd.DataFrame(CTRL)
cdf.to_csv(os.path.join(OUT, 'controls.tsv'), sep='\t', index=False)
fails = set(cdf.loc[cdf.result == 'FAIL', 'control'])
withheld = []
if fails & {'PC-1a', 'PC-1b', 'NC-1', 'NC-2'}:
    withheld.append('gnomad')
    diag('4_controls', 'S4', 'gnomAD columns withheld: ' + ','.join(sorted(fails & {'PC-1a', 'PC-1b', 'NC-1', 'NC-2'})), status='WITHHELD')
if fails & {'PC-2', 'NC-4'}:
    withheld.append('depmap')
    diag('4_controls', 'S5', 'DepMap columns withheld: ' + ','.join(sorted(fails & {'PC-2', 'NC-4'})), status='WITHHELD')
if fails & {'PC-3a', 'PC-3b', 'NC-3'}:
    withheld.append('hpa')
    diag('4_controls', 'S6', 'HPA columns withheld: ' + ','.join(sorted(fails & {'PC-3a', 'PC-3b', 'NC-3'})), status='WITHHELD')
if args.debug_force_ship:
    assert DEBUG, '--debug-force-ship is only allowed with --debug-null'
    withheld = []
    log('DEBUG-NULL: S4-S6 ignored so that D4/D5 code paths run')
ship = [s for s in ['gnomad', 'depmap', 'hpa'] if s not in withheld]
diag('4_controls', 'sources shipping', ','.join(ship) if ship else 'none')

# ---------------------------------------------------------------- 5. D2 and D3
d2 = L.dropna(subset=['lof.oe_ci.upper_bin_decile'])
ct = pd.crosstab(d2.dec, d2['lof.oe_ci.upper_bin_decile'].astype(int))
ct.index.name = 'recomputed_decile'
ct.columns.name = 'official_decile'
ct.to_csv(os.path.join(OUT, 'd2_decile_crosstab.tsv'), sep='\t')
diag_frac = float((d2.dec == d2['lof.oe_ci.upper_bin_decile']).mean())
within1 = float(((d2.dec - d2['lof.oe_ci.upper_bin_decile']).abs() <= 1).mean())
diag('5_D2', 'genes with both deciles', len(d2), 17319)
diag('5_D2', 'fraction on diagonal', round(diag_frac, 4))
diag('5_D2', 'fraction within one decile', round(within1, 4))
diag('5_D2', 'recomputed decile bin sizes', L.dec.value_counts().sort_index().tolist())
diag('5_D2', 'recomputed decile-0 max LOEUF', float(L.loc[L.dec == 0, 'lof.oe_ci.upper'].max()))
log(ct.to_string())
frac_lt = (ge < -0.5).mean(axis=0)  # gate rule: NaN counts as not < -0.5, denominator all lines
approx = set(frac_lt[frac_lt >= 0.9].index)
approx_ez = set(dmc.loc[dmc.col.isin(approx), 'entrez'])
both = len(approx_ez & ce_set)
d3 = pd.DataFrame([{'rule': 'Chronos < -0.5 in >= 90% of all 1,178 lines (gate rule, NaN = not below)', 'approx_n': len(approx_ez),
                    'official_n': len(ce_set), 'both': both, 'approx_only': len(approx_ez - ce_set), 'official_only': len(ce_set - approx_ez),
                    'jaccard': round(both / len(approx_ez | ce_set), 4), 'gate_count_on_unrecorded_matrix': 981}])
d3.to_csv(os.path.join(OUT, 'd3_depmap_approx_vs_official.tsv'), sep='\t', index=False)
for c in d3.columns:
    diag('5_D3', c, d3.iloc[0][c])

# ---------------------------------------------------------------- 6. spot check by a second code path (grep/awk on raw files)
gn_hdr = open(GN).readline().rstrip('\n').split('\t')
gi = {c: i for i, c in enumerate(gn_hdr)}
hp_hdr = [c.strip('"') for c in open(HPA).readline().rstrip('\n').split('\t')]
hpi = {c: i for i, c in enumerate(hp_hdr)}
ge_hdr = open(GE).readline().rstrip('\n').split(',')
names_re = '|'.join(SPOT_GENES)
gtf_lines = sh(f"zcat {GTF} | awk -F'\\t' '$3==\"gene\"' | grep -E 'gene_name \"({names_re})\";'").splitlines()


def gtf_attr(line, key):
    for part in line.split('\t')[8].split(';'):
        part = part.strip()
        if part.startswith(key + ' '):
            return part[len(key) + 1:].strip('"')
    return None


spot_ens = {}
spot_gtype = {}
for line in gtf_lines:
    nm = gtf_attr(line, 'gene_name')
    if nm in SPOT_GENES and line.split('\t')[0] in PRIMARY and nm not in spot_ens:
        spot_ens[nm] = gtf_attr(line, 'gene_id').split('.')[0]
        spot_gtype[nm] = gtf_attr(line, 'gene_type')
# second-path decile reference: all Ensembl MANE LOEUF values from awk, ranked by (LOEUF, gene_id)
li, gii, mi = gi['lof.oe_ci.upper'] + 1, gi['gene_id'] + 1, gi['mane_select'] + 1
aw = sh(f"awk -F'\\t' 'NR>1 && ${mi}==\"true\" && ${gii} ~ /^ENSG/ && ${li}!=\"NA\" {{print ${gii}\"\\t\"${li}}}' {GN}").splitlines()
ref = sorted(((float(v), e) for e, v in (x.split('\t') for x in aw)), key=lambda z: (z[0], z[1]))
nref = len(ref)
rank_of = {e: i + 1 for i, (v, e) in enumerate(ref)}


def decile2(r, n):
    # qcut on ranks 1..n into 10 equal-frequency bins: bin k holds e_k < r <= e_{k+1}, e_k = 1 + (n-1)k/10 (exact integer form)
    k = 0
    for j in range(1, 10):
        if (n - 1) * j < 10 * (r - 1):
            k = j
    return k


ce_lines = open(CE).read().splitlines()[1:]
spot = []
cat_pa = pd.read_parquet(CAT, engine='pyarrow', columns=['ensembl_id', 'symbol']).set_index('ensembl_id').symbol


def fnum(s):
    return None if s in ('', 'NA', None) else float(s)


def same(a, b):
    a_blank = a is None or (isinstance(a, float) and np.isnan(a)) or a is pd.NA or (isinstance(a, str) and a == '')
    b_blank = b is None or (isinstance(b, float) and np.isnan(b)) or b is pd.NA or (isinstance(b, str) and b == '')
    if a_blank or b_blank:
        return a_blank and b_blank
    if isinstance(b, float) or isinstance(a, float):
        return bool(np.isclose(float(a), float(b), rtol=1e-9, atol=0))
    return str(a) == str(b)


for s in SPOT_GENES:
    e = spot_ens.get(s)
    exp = {'ensembl_id': e, 'symbol': cat_pa.get(e)}
    raw = [x.split('\t') for x in sh(f"grep -P '\\t{e}\\t' {GN}").splitlines()]
    raw = [x for x in raw if x[gi['gene_id']] == e and x[gi['mane_select']] == 'true']
    if len(raw) > 1:
        stop('S3', f'spot check {s}: {len(raw)} Ensembl MANE rows')
    r0 = raw[0] if raw else None
    exp['gnomad_v41_mane_transcript'] = r0[gi['transcript']] if r0 else None
    exp['gnomad_v41_loeuf'] = fnum(r0[gi['lof.oe_ci.upper']]) if r0 else None
    exp['gnomad_v41_pli'] = fnum(r0[gi['lof.pLI']]) if r0 else None
    exp['gnomad_v41_loeuf_decile_mane_ens'] = decile2(rank_of[e], nref) if e in rank_of else None
    od = fnum(r0[gi['lof.oe_ci.upper_bin_decile']]) if r0 else None
    exp['gnomad_v41_loeuf_decile_official'] = int(od) if od is not None else None
    fl = r0[gi['constraint_flags']].strip() if r0 else ''
    exp['gnomad_v41_constraint_flags'] = '' if fl == '[]' else fl
    ez = None
    if r0:
        sym_gn = r0[gi['gene']]
        zl = sh(f"awk -F'\\t' 'NR>1 && $1==\"{sym_gn}\" && ${mi}==\"true\" && ${gii} ~ /^[0-9]+$/ {{print ${gii}}}' {GN}").split()
        if len(zl) > 1:
            stop('S3', f'spot check {s}: symbol {sym_gn} has {len(zl)} Entrez MANE rows')
        ez = int(zl[0]) if zl else None
    exp['entrez_id'] = ez
    colk = [i for i, c in enumerate(ge_hdr) if ez is not None and c.endswith(f' ({ez})')]
    if colk:
        vals = [v for v in sh(f"cut -d, -f{colk[0] + 1} {GE} | tail -n +2").split('\n')]
        vals = [float(v) for v in vals if v.strip() != '']
        exp['depmap24q4_median_gene_effect'] = statistics.median(vals)
        exp['depmap24q4_n_lines'] = len(vals)
        exp['depmap24q4_common_essential'] = 'TRUE' if any(x.endswith(f'({ez})') and x.split(' (')[-1] == f'{ez})' for x in ce_lines) else 'FALSE'
    else:
        exp['depmap24q4_median_gene_effect'] = None
        exp['depmap24q4_n_lines'] = None
        exp['depmap24q4_common_essential'] = ''
    hl = [x.split('\t') for x in sh(f"grep -P '\\t{e}\\t' {HPA}").splitlines()]
    hl = [x for x in hl if x[hpi['Ensembl']].strip('"') == e]
    if len(hl) > 1:
        stop('S3', f'spot check {s}: {len(hl)} HPA rows')
    if hl:
        sec = hl[0][hpi[HPA_SEC]].strip('"')
        exp['hpa_secretome_location'] = sec
        exp['hpa_secreted'] = ('not_annotated' if sec == '' else 'secreted' if sec.startswith('Secreted') else
                               'intracellular_and_membrane' if sec == 'Intracellular and membrane' else
                               'immunoglobulin' if sec == 'Immunoglobulin genes' else 'UNMAPPED')
        exp['hpa_blood_ms_pg_per_l'] = fnum(hl[0][hpi[HPA_MS]].strip('"'))
        exp['hpa_blood_im_pg_per_l'] = fnum(hl[0][hpi[HPA_IM]].strip('"'))
    else:
        exp['hpa_secretome_location'] = ''
        exp['hpa_secreted'] = ''
        exp['hpa_blood_ms_pg_per_l'] = None
        exp['hpa_blood_im_pg_per_l'] = None
    reasons = []
    if exp['gnomad_v41_loeuf'] is None:
        reasons.append('gnomad_mane_no_loeuf' if r0 else ('protein_coding_no_gnomad_mane_row' if spot_gtype.get(s) == 'protein_coding' else 'noncoding_gene'))
    if ez is None:
        reasons.append('no_entrez_partner')
    elif not colk:
        reasons.append('not_in_depmap24q4')
    if not hl:
        reasons.append('not_in_hpa')
    exp['annotation_blank_reason'] = ';'.join(reasons)
    spot.append((s, e, exp))

OUTCOLS = ['ensembl_id', 'symbol', 'gnomad_v41_mane_transcript', 'gnomad_v41_loeuf', 'gnomad_v41_pli', 'gnomad_v41_loeuf_decile_mane_ens',
           'gnomad_v41_loeuf_decile_official', 'gnomad_v41_constraint_flags', 'entrez_id', 'depmap24q4_common_essential',
           'depmap24q4_median_gene_effect', 'depmap24q4_n_lines', 'hpa_secretome_location', 'hpa_secreted', 'hpa_blood_ms_pg_per_l',
           'hpa_blood_im_pg_per_l', 'annotation_blank_reason']
t['ensembl_id'] = t.ens
t['annotation_blank_reason'] = [';'.join(x for x in (a, b, c) if x) for a, b, c in zip(t.r_gn, t.r_dm, t.r_hp)]
ti = t.set_index('ens')
sp_rows = []
mism = []
for s, e, exp in spot:
    row = ti.loc[e] if e in ti.index else None
    for c in OUTCOLS:
        o = None if row is None else row[c]
        if isinstance(o, (np.integer,)):
            o = int(o)
        if isinstance(o, (np.floating,)):
            o = float(o)
        ok = row is not None and same(o, exp[c])
        sp_rows.append({'gene': s, 'ensembl_id': e, 'column': c, 'output': '' if DEBUG else o, 'second_path': '' if DEBUG else exp[c],
                        'match': ok})
        if not ok:
            mism.append(f'{s}.{c}')
spdf = pd.DataFrame(sp_rows)
spdf.to_csv(os.path.join(OUT, 'spot_check.tsv'), sep='\t', index=False)
diag('6_spot', 'GENCODE primary-assembly Ensembl IDs', json.dumps(spot_ens))
diag('6_spot', 'second-path decile reference size', nref, 17325)
diag('6_spot', 'cells compared', len(spdf), len(SPOT_GENES) * len(OUTCOLS))
diag('6_spot', 'mismatches', len(mism), 0, 'equal', 'pass' if not mism else 'FAIL')
pol = ti.loc[spot_ens['POLR2A']] if 'POLR2A' in spot_ens else None
pol_ok = pol is not None and pd.isna(pol.gnomad_v41_loeuf) and pol.depmap24q4_common_essential == '' and \
    pol.annotation_blank_reason == 'protein_coding_no_gnomad_mane_row;no_entrez_partner'
diag('6_spot', 'POLR2A prediction (blank gnomAD and DepMap; reason protein_coding_no_gnomad_mane_row;no_entrez_partner)',
     pol.annotation_blank_reason if pol is not None else 'missing', status='as predicted' if pol_ok else 'NOT as predicted')
if not DEBUG:
    log(spdf.pivot(index='column', columns='gene', values='output').reindex(OUTCOLS).to_string())
if mism:
    stop('S3', 'spot-check mismatch: ' + ', '.join(mism))

# ---------------------------------------------------------------- withholding and the primary table
S7 = len(ship) == 0
cols = ['ensembl_id', 'symbol']
gcols = ['gnomad_v41_mane_transcript', 'gnomad_v41_loeuf', 'gnomad_v41_pli', 'gnomad_v41_loeuf_decile_mane_ens',
         'gnomad_v41_loeuf_decile_official', 'gnomad_v41_constraint_flags']
dcols = ['depmap24q4_common_essential', 'depmap24q4_median_gene_effect', 'depmap24q4_n_lines']
hcols = ['hpa_secretome_location', 'hpa_secreted', 'hpa_blood_ms_pg_per_l', 'hpa_blood_im_pg_per_l']
if 'gnomad' in ship:
    cols += gcols
cols += ['entrez_id']
if 'depmap' in ship:
    cols += dcols
if 'hpa' in ship:
    cols += hcols
parts = [t.r_gn if 'gnomad' in ship else pd.Series('', index=t.index), t.r_dm if 'depmap' in ship else pd.Series('', index=t.index),
         t.r_hp if 'hpa' in ship else pd.Series('', index=t.index)]
t['annotation_blank_reason'] = [';'.join(x for x in z if x) for z in zip(*parts)]
cols += ['annotation_blank_reason']
write_diag()
tab_name = 'catalog_followup_annotation.tsv' if not S7 else 'catalog_followup_annotation.NOT_SHIPPED_S7.tsv'
outt = t[cols].sort_values('ensembl_id')
assert len(outt) == 28233 and outt.ensembl_id.is_unique
outt.to_csv(os.path.join(OUT, tab_name), sep='\t', index=False, float_format='%.10g')
log(f'wrote {tab_name}: {outt.shape}')

# ---------------------------------------------------------------- D1, D4, D5
route = t.next_experiment_rule_id.copy()
if DEBUG:
    route = pd.Series(RNG.permutation(route.values), index=t.index)
    log('DEBUG-NULL: route labels permuted for D1/D4/D5')
t['route'] = route


def d1row(x):
    return pd.Series({'genes': len(x), 'loeuf': int(x.has_loeuf.sum()), 'depmap24q4': int(x.in_dm.sum()), 'hpa_row': int(x.in_hpa.sum()),
                      'hpa_secretome_text': int((x.hpa_secretome_location != '').sum()), 'blood_ms': int(x.hpa_blood_ms_pg_per_l.notna().sum()),
                      'blood_im': int(x.hpa_blood_im_pg_per_l.notna().sum()),
                      'noncoding_gene': int((x.r_gn == 'noncoding_gene').sum()),
                      'protein_coding_no_gnomad_mane_row': int((x.r_gn == 'protein_coding_no_gnomad_mane_row').sum()),
                      'gnomad_mane_no_loeuf': int((x.r_gn == 'gnomad_mane_no_loeuf').sum()),
                      'no_entrez_partner': int((x.r_dm == 'no_entrez_partner').sum()), 'not_in_depmap24q4': int((x.r_dm == 'not_in_depmap24q4').sum()),
                      'not_in_hpa': int((x.r_hp == 'not_in_hpa').sum())})


d1 = pd.DataFrame({r: d1row(t[t.route == r]) for r in ROUTES}).T
d1.loc['ALL'] = d1row(t)
d1.index.name = 'route'
d1.to_csv(os.path.join(OUT, 'd1_coverage.tsv'), sep='\t')
log('D1 coverage\n' + d1.to_string())

d4_fields = [k for k in ['gnomad', 'depmap', 'hpa'] if k in ship and G4[k] >= 0.90]
d4 = []
for r in ROUTES + ['ALL']:
    x = t if r == 'ALL' else t[t.route == r]
    row = {'route': r, 'genes': len(x)}
    if 'gnomad' in d4_fields:
        nl = int(x.has_loeuf.sum())
        n0 = int((x.gnomad_v41_loeuf_decile_mane_ens == 0).sum())
        row.update({'n_loeuf': nl, 'n_loeuf_decile0': n0, 'frac_loeuf_decile0': round(n0 / nl, 4) if nl else np.nan})
    if 'depmap' in d4_fields:
        nd = int(x.in_dm.sum())
        nc = int((x.depmap24q4_common_essential == 'TRUE').sum())
        row.update({'n_depmap24q4': nd, 'n_common_essential': nc, 'frac_common_essential': round(nc / nd, 4) if nd else np.nan})
    if 'hpa' in d4_fields:
        nh = int(x.in_hpa.sum())
        ns = int((x.hpa_secreted == 'secreted').sum())
        row.update({'n_hpa_row': nh, 'n_secreted': ns, 'frac_secreted': round(ns / nh, 4) if nh else np.nan,
                    'n_blood_ms': int(x.hpa_blood_ms_pg_per_l.notna().sum())})
    d4.append(row)
d4 = pd.DataFrame(d4)
d4.to_csv(os.path.join(OUT, 'd4_route_description.tsv'), sep='\t', index=False)
log('D4 route description (unmatched fractions; no test; LOEUF direction known from the gate)\n' + d4.to_string())

d5 = t[(t.route == STRESS)].copy()
c0 = d5.gnomad_v41_loeuf_decile_mane_ens.eq(0).fillna(False) if 'gnomad' in ship else pd.Series(False, index=d5.index)
c1 = d5.depmap24q4_common_essential.eq('TRUE') if 'depmap' in ship else pd.Series(False, index=d5.index)
d5 = d5[c0 | c1].copy()
d5['criterion'] = np.where(c0[d5.index] & c1[d5.index], 'loeuf_decile0_and_common_essential',
                           np.where(c0[d5.index], 'loeuf_decile0', 'common_essential'))
d5cols = ['ensembl_id', 'symbol', 'criterion'] + (['gnomad_v41_loeuf', 'gnomad_v41_loeuf_decile_mane_ens'] if 'gnomad' in ship else []) + \
         (['depmap24q4_common_essential', 'depmap24q4_median_gene_effect'] if 'depmap' in ship else []) + \
         (['hpa_secreted', 'hpa_blood_ms_pg_per_l'] if 'hpa' in ship else [])
d5 = d5[d5cols].sort_values(['criterion', 'symbol'])
d5.to_csv(os.path.join(OUT, 'd5_stress_route_followup_list.tsv'), sep='\t', index=False, float_format='%.10g')
log(f'D5: {len(d5)} stress-route genes; by criterion {d5.criterion.value_counts().to_dict()}')

# ---------------------------------------------------------------- verdict and RESULT
if S7:
    verdict = 'untestable (S7)'
    rule_line = 'S7: all three sources withheld by S4 to S6. Nothing ships; the Catalog is unchanged; the paper says nothing.'
elif len(ship) == 3:
    verdict = 'pass (all three sources ship)'
    rule_line = 'All seven controls passed their pass rules, so gnomAD, DepMap and HPA columns all ship in the supplementary table.'
else:
    verdict = 'partial (' + ', '.join(ship) + ' ship; ' + ', '.join(withheld) + ' withheld)'
    rule_line = 'One or two sources withheld by S4 to S6. Ship the passing sources; the Methods sentence names only those.'

ctrl_md = cdf[['control', 'source', 'value', 'n', 'predicted', 'pass_rule', 'result', 'note']].to_markdown(index=False)
g1_md = '\n'.join(f'- {i}: {v} (expected {e})' for i, v, e in G1)
src_names = {'gnomad': 'gnomAD v4.1 LOEUF and pLI for the Ensembl MANE Select transcript (deciles over 17,325 genes)',
             'depmap': 'the DepMap Public 24Q4 common-essential call and median Chronos gene effect',
             'hpa': 'the Human Protein Atlas secretome location and blood concentrations'}
named = [src_names[s] for s in ship]
methods = ('For every Catalog gene we added ' + (', '.join(named[:-1]) + ', and ' + named[-1] if len(named) > 2 else ' and '.join(named)) +
           '. These fields describe how a gene might behave in a follow-up experiment. They do not enter any evidence class or route.') if ship else ''
d1a = d1.loc['ALL']
body = f"""## Controls (pass/fail rules on fractions; not tests; family size 0, no correction, no p value)

{ctrl_md}

Achilles common-essential controls in the inferred list: {ctl_in} / {len(ctl)} = {ctl_in / len(ctl):.4f} (printed, not a rule).
Withheld sources: {', '.join(withheld) if withheld else 'none'}. Shipping: {', '.join(ship) if ship else 'none'}.

## Gate G1 (all equal to the prespec)

{g1_md}

G2 key checks: all zero duplicates. G3: {n_ce_in} of 1,523 common essentials are 24Q4 matrix columns.
G4 stress-route coverage: LOEUF {G4['gnomad']:.4f}, DepMap {G4['depmap']:.4f}, HPA {G4['hpa']:.4f} (D4 fields printed: {', '.join(d4_fields)}).

## Spot check (second code path: grep/awk on raw files)

{len(spdf)} cells for 10 genes, {len(mism)} mismatches. POLR2A prediction: {'as predicted' if pol_ok else 'NOT as predicted'} (reason `{pol.annotation_blank_reason if pol is not None else 'missing'}`).

## Primary output

`{tab_name}`: {len(outt)} rows, one per Catalog r8 gene, columns {', '.join(cols)}.
Atlas-wide coverage: LOEUF {int(d1a.loeuf)} ({d1a.loeuf / 28233:.1%}), DepMap 24Q4 {int(d1a.depmap24q4)} ({d1a.depmap24q4 / 28233:.1%}), HPA row {int(d1a.hpa_row)} ({d1a.hpa_row / 28233:.1%}).
Blank reasons: noncoding_gene {int(d1a.noncoding_gene)} (non-applicability); protein_coding_no_gnomad_mane_row {int(d1a.protein_coding_no_gnomad_mane_row)} and gnomad_mane_no_loeuf {int(d1a.gnomad_mane_no_loeuf)} (missing data); no_entrez_partner {int(d1a.no_entrez_partner)}; not_in_depmap24q4 {int(d1a.not_in_depmap24q4)}; not_in_hpa {int(d1a.not_in_hpa)}.
`gnomad_v41_constraint_flags` is the gnomAD text with the empty list `[]` written as blank.

## D2: recomputed vs official LOEUF decile

{len(d2)} genes with both: {diag_frac:.4f} on the diagonal, {within1:.4f} within one decile (`d2_decile_crosstab.tsv`).

## D3: gate approximate DepMap rule on the 24Q4 matrix

Chronos < -0.5 in >= 90% of lines: {len(approx_ez)} genes (gate: 981 on the unrecorded 1,186-line matrix). Overlap with the official 1,523: {both}; approximate only {len(approx_ez - ce_set)}; official only {len(ce_set - approx_ez)}; Jaccard {both / len(approx_ez | ce_set):.4f}.

## D4: route description (unmatched fractions; no test, no interval; not blind for LOEUF)

{d4.to_markdown(index=False)}

The in-data stress-route LOEUF difference was below the minimum detectable difference, so no test was run and these fractions support no statement about direction.

## D5: stress-route genes with LOEUF decile 0 or a common-essential call

{len(d5)} genes (`d5_stress_route_followup_list.tsv`); by criterion {d5.criterion.value_counts().to_dict()}. Tied to Catalog r8 route membership.

## Proposed text (returned as data; not installed)

- Methods (`docs/manuscript/METHODS.md:761-788`): "{methods}"
- Figure 6 legend (`docs/manuscript/05_figure_legends.md:427`): "Constraint, essentiality and secretion fields for each gene are in Supplementary Table [N]; they are annotations and do not change classes."
- `docs/RESULTS.md:2707-2735`: "Follow-up annotation fields: {', '.join({'gnomad': 'gnomAD v4.1', 'depmap': 'DepMap 24Q4', 'hpa': 'HPA'}[s] for s in ship)}; 28,233 genes; coverage {' / '.join(f'{v / 28233:.1%}' for s, v in [('gnomad', d1a.loeuf), ('depmap', d1a.depmap24q4), ('hpa', d1a.hpa_row)] if s in ship)}"

In no reading does the paper state that stress-route genes are more or less constrained, essential or secreted than other genes.
"""
write_result(verdict, rule_line, body)
write_diag()
log('done')
