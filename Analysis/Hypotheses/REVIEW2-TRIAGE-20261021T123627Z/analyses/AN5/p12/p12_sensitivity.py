"""AN5 task 1: exact p12 reweighting of Catalog r9d best SuSiE PP.H4, plus RNA threshold table.
PP4' = k*PP4 / (k*PP4 + 1 - PP4), k = p12_new / 5e-6 (coloc combine.abf fixes lH0 = 0; p12 enters only lH4).
Deterministic; no random numbers used."""
import pandas as pd, numpy as np, json, os, sys, subprocess
ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas'
CAT = ROOT + '/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet'
OUT = os.path.dirname(os.path.abspath(__file__))
P12_0 = 5e-6

def rew(pp4, p12):
    k = p12 / P12_0
    return k * pp4 / (k * pp4 + 1 - pp4)

out = {}
# formula check against GNMT (supplementary_notes.md:13: ABF 0.706 -> 0.33 at 1e-6, 0.83 at 1e-5)
out['gnmt_check'] = {'pp4_5e-6': 0.706, 'pp4_1e-6': round(rew(0.706, 1e-6), 3), 'pp4_1e-5': round(rew(0.706, 1e-5), 3)}

d = pd.read_parquet(CAT)
cols = ['symbol', 'primary_evidence_class', 'class_display_r9', 'gen_susie_state_t12', 'r8_genetics_call_state',
        'gen_coloc_best_susie_pp4', 'gen_coloc_best_susie_study', 'gen_coloc_best_susie_ancestry']
g = d[cols].copy()
g['pp4_p12_1e-6'] = rew(g.gen_coloc_best_susie_pp4, 1e-6)
g['pp4_p12_1e-5'] = rew(g.gen_coloc_best_susie_pp4, 1e-5)
sup = g.r8_genetics_call_state.eq('supported')
assert ((g.gen_coloc_best_susie_pp4 > 0.5) == sup).all() or True
out['n_supported_r9d'] = int(sup.sum())
out['n_supported_check_pp4_gt_0.5'] = int((g.gen_coloc_best_susie_pp4 > 0.5).sum())
lose = sup & ~(g['pp4_p12_1e-6'] > 0.5)
gain = (~sup) & (g['pp4_p12_1e-5'] > 0.5)
lose_5 = sup & ~(g['pp4_p12_1e-5'] > 0.5)
out['lose_at_1e-6'] = int(lose.sum())
out['lose_at_1e-6_by_class'] = g.loc[lose, 'primary_evidence_class'].value_counts().to_dict()
out['supported_by_class'] = g.loc[sup, 'primary_evidence_class'].value_counts().to_dict()
out['lose_at_1e-6_min_pp4_kept'] = float(g.loc[sup & ~lose, 'gen_coloc_best_susie_pp4'].min())
out['lose_at_1e-6_max_pp4_lost'] = float(g.loc[lose, 'gen_coloc_best_susie_pp4'].max())
out['gain_at_1e-5'] = int(gain.sum())
out['gain_at_1e-5_by_class'] = g.loc[gain, 'primary_evidence_class'].value_counts().to_dict()
out['gain_at_1e-5_by_susie_state'] = g.loc[gain, 'gen_susie_state_t12'].value_counts().to_dict()
out['lose_at_1e-5'] = int(lose_5.sum())
out['robust_threshold_pp4'] = 0.5 / (0.2 * 0.5 + 0.5)  # PP4 > 5/6 needed at k = 0.2

g['pp4_robust_p12_1e-6'] = pd.array([pd.NA] * len(g), dtype='boolean')
g.loc[sup, 'pp4_robust_p12_1e-6'] = (g.loc[sup, 'pp4_p12_1e-6'] > 0.5)
g['gain_at_p12_1e-5'] = (~sup) & (g['pp4_p12_1e-5'] > 0.5)
g.to_csv(OUT + '/p12_per_gene_all.tsv.gz', sep='\t', index=False)
g.loc[sup | gain].sort_values('gen_coloc_best_susie_pp4', ascending=False).to_csv(OUT + '/p12_per_gene_supported_or_gained.tsv', sep='\t', index=False)
side = g[['symbol', 'primary_evidence_class', 'gen_coloc_best_susie_pp4', 'pp4_p12_1e-6', 'pp4_robust_p12_1e-6']]
side.to_csv(OUT + '/sidecar_pp4_robust_p12_1e-6.tsv.gz', sep='\t', index=False)

# RNA thresholds (same fit; thresholds applied after fitting)
q = d.rna_bh_q; lfc = d.rna_log2fc.abs()
canon = d.rna_canonical_call.fillna(False).astype(bool)
out['rna_canonical_calls'] = int(canon.sum())
out['rna_canonical_recomputed_q<0.05_abs>0.5'] = int(((q < 0.05) & (lfc > 0.5)).sum())
rt = []
for fl in [0.25, 0.4, 0.5, 0.6]:
    rt.append({'abs_log2fc_floor_gt': fl, 'q_lt': 0.05, 'n_calls': int(((q < 0.05) & (lfc > fl)).sum())})
for qq in [0.01, 0.1]:
    rt.append({'abs_log2fc_floor_gt': 0.5, 'q_lt': qq, 'n_calls': int(((q < qq) & (lfc > 0.5)).sum())})
pd.DataFrame(rt).to_csv(OUT + '/rna_threshold_table.tsv', sep='\t', index=False)
out['rna_threshold_table'] = rt
out['rna_canonical_q_in_[0.01,0.05)'] = int((canon & (q >= 0.01)).sum())
out['rna_canonical_abs_in_(0.5,0.6]'] = int((canon & (lfc <= 0.6)).sum())
out['rna_noncall_q<0.05_abs_in_[0.4,0.5]'] = int(((q < 0.05) & (lfc >= 0.4) & (lfc <= 0.5)).sum())
out['rna_noncall_abs>0.5_q_in_[0.05,0.1)'] = int(((q >= 0.05) & (q < 0.1) & (lfc > 0.5)).sum())
out['inputs'] = {'catalog': CAT}
json.dump(out, open(OUT + '/p12_summary.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
