"""Re-read the producing files behind the no_source rows and stage RESULTS additions.

Prints TSV rows (target_id, source_path, column_or_row, value_reread) to stdout and
writes nosource_edits.json next to this script. Reads only; writes only into work/.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
WORK = Path(__file__).resolve().parent
RESULTS = ROOT / "docs/RESULTS.md"

T4 = ("RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/"
      "workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction")
FM = "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321"
F2RUN = "figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
H58 = "Analysis/Hypotheses/H58-20261015T112423Z/jobs/5800/H58-20261015T113519Z_22252009"
R4 = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
R14 = R4 + "/check/Claims/R14-10"
MPRA = R4 + "/repair/F7a/exec/mpra/mpra_bonferroni.tsv"
H79 = "Analysis/Hypotheses/H79-20261015T125125Z/results/run-20261015T132031Z/s3_ancestry_by_cohort_group.tsv"
GSE83452 = "Analysis/MASLD_Model_Benchmark/executions/gse83452-intervention-response-20260831/results_20260831/results.json"
STUDY = "Analysis/MASLD_Model_Benchmark/executions/encoder-field-score-21314196/evaluation/study_metrics.tsv"
F6 = "figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair/panels_20261015T173321Z"

rows = []


def emit(tid, path, where, value):
    rows.append((tid, path, where, str(value)))


def tsv(rel, **kw):
    return pd.read_csv(ROOT / rel, sep="\t", **kw)


# ---- T1 title.md: modularity nulls and the 15-program community
summ = tsv(f"{T4}/interaction_summary.tsv").set_index("metric")["value"].astype(float)
mod_obs = summ["modularity_observed"]
mod_rank = summ["modularity_rank_matched_null_mean"]
mod_rew = summ["modularity_rewired_mean"]
emit("T1a", f"{T4}/interaction_summary.tsv", "metric=modularity_observed", f"{mod_obs:.6f}")
emit("T1a", f"{T4}/interaction_summary.tsv", "metric=modularity_rank_matched_null_mean", f"{mod_rank:.6f}")
emit("T1a", f"{T4}/interaction_summary.tsv", "metric=modularity_rewired_mean", f"{mod_rew:.6f}")
mn = tsv(f"{T4}/modularity_null.tsv")
n_null = mn.groupby("null").size().to_dict()
emit("T1a", f"{T4}/modularity_null.tsv", "rows per null", json.dumps(n_null))
assert abs(mn[mn.null == "rank_matched_permutation"].modularity.mean() - mod_rank) < 1e-9
assert abs(mn[mn.null == "degree_preserving_rewiring"].modularity.mean() - mod_rew) < 1e-9
com = tsv(f"{T4}/program_communities.tsv")
com["lineage"] = com.feature_id.str.split("_").str[1]
tab = pd.crosstab(com.community, com.lineage)
c15 = [c for c in tab.index if tab.loc[c].sum() == 15]
assert len(c15) == 1
c15 = c15[0]
n_fib, n_chol = int(tab.loc[c15, "fibroblasts"]), int(tab.loc[c15, "cholangiocytes"])
assert tab.loc[c15].sum() == n_fib + n_chol
emit("T1b", f"{T4}/program_communities.tsv", f"community={c15}",
     f"15 programs; fibroblasts {n_fib}; cholangiocytes {n_chol}; other 0")
n_comm = int(summ["n_communities"])

# ---- T2 fig2.md ~124: VEP annotation covers all normalized PIP of gated sets
arch = tsv(f"{FM}/credible_set_pip_architecture.tsv")
gate = arch[arch.pip_sum_gate_passed.astype(str).str.lower() == "true"]
emit("T2", f"{FM}/credible_set_pip_architecture.tsv", "gated rows; max unresolved_pip_mass; min annotated_pip_mass",
     f"{len(gate)} of {len(arch)}; {gate.unresolved_pip_mass.max():.6g}; {gate.annotated_pip_mass.min():.6g}")

# ---- T3 fig2.md ~126: PNPLA3 cluster, 124 direct sets, mean coding+splice 0.117
sys.path.insert(0, str(ROOT / F2RUN / "snapshot/scripts/figures"))
import fig2_multiplicity_correction as m  # noqa: E402

cc = pd.read_csv(ROOT / F2RUN / "multiplicity/corrected_locus_coordinates.tsv", sep="\t", dtype={"locus": str})
cc = cc[cc.trait_scope.isin(m.SCOPES)].copy()
assert cc.pip_sum_gate_passed.astype(str).str.lower().eq("true").all()
cc["coding"] = cc[m.CODING_COLS].sum(axis=1)
cc = m.assign_clusters(cc)
direct = cc[cc.trait_scope == m.SCOPES[0]]
cl = direct.groupby("cluster_id").agg(n=("coding", "size"), mean_coding=("coding", "mean"),
                                      chrom=("chrom", "first"), lo=("pos", "min"), hi=("pos", "max"))
cl["share"] = cl.mean_coding / cl.mean_coding.sum()
pn = cl[(cl.chrom.astype(str) == "22") & (cl.lo <= 44324727) & (cl.hi >= 44324727)]
assert len(pn) == 1
pn_n, pn_mean = int(pn.n.iloc[0]), float(pn.mean_coding.iloc[0])
top5 = cl.share.sort_values(ascending=False).head(5).sum()
emit("T3", f"{F2RUN}/multiplicity/corrected_locus_coordinates.tsv",
     "direct sets in the chr22 cluster holding 44,324,727 (producer assign_clusters, 1 Mb gap); mean protein-altering+canonical-splice mass",
     f"cluster {pn.index[0]}; n {pn_n}; mean {pn_mean:.6f}; span {int(pn.lo.iloc[0])}-{int(pn.hi.iloc[0])}")
emit("T3", f"{F2RUN}/multiplicity/corrected_locus_coordinates.tsv", "direct clusters; top-5 share; largest share",
     f"{len(cl)}; {top5:.4f}; {cl.share.max():.4f}")

# ---- T4 fig5_discussion.md ~207: locus-clustered sign-flip
log = (ROOT / R14 / "r14_10_recompute.log").read_text()
mm = re.search(r"seed 20260907 pair signflip p (\S+) locus signflip p (\S+)", log)
p_locus = float(mm.group(2))
p_locus_sign = float(re.search(r"locus sign test p (\S+)", log).group(1))
pairs = tsv(f"{R14}/pairs_with_locus_recomputed.tsv")
lm = pairs.groupby("locus").d.mean()
n_loci, n_pos_loci, n_pairs = int(lm.size), int((lm > 0).sum()), len(pairs)
emit("T4", f"{R14}/r14_10_recompute.log", "seed 20260907 locus signflip p (100,000 draws)", f"{p_locus:.6f}")
emit("T4", f"{R14}/pairs_with_locus_recomputed.tsv", "n loci; locus means > 0; pairs",
     f"{n_loci}; {n_pos_loci}; {n_pairs}")
emit("T4", f"{R14}/r14_10_recompute.log", "locus sign test p", f"{p_locus_sign:.4f}")

# ---- T5 fig7.md ~27: Bonferroni lower bound for HyenaDNA
mp = tsv(MPRA)
f7 = mp[mp.family == "F7_displayed_rows"]
hy = f7[f7.model == "hyenadna"].iloc[0]
n_pos_bound = int((f7.bonf_pct_lo > 0).sum())
emit("T5", MPRA, "family=F7_displayed_rows, model=hyenadna: m; gain; lo95; hi95; bonf_pct_lo; n_valid",
     f"{int(hy.m)}; {hy.gain:.6f}; {hy.lo95:.6f}; {hy.hi95:.6f}; {hy.bonf_pct_lo:.6f}; {int(hy.n_valid)}")
emit("T5", MPRA, "F7_displayed_rows with bonf_pct_lo > 0; distinct models", f"{n_pos_bound}; {f7.model.nunique()}")

# ---- O1 Fig 2A: 23 of 767 EUR enzyme lead windows within 500 kb
g7 = tsv(f"{H58}/G7_coverage.tsv")
eur = g7[(g7.ancestry == "EUR") & (g7.trait_scope == "tier2_liver_enzyme")].iloc[0]
emit("O1", f"{H58}/G7_coverage.tsv", "EUR tier2_liver_enzyme: lead_windows; windows_with_cs; windows_with_cs_pos500kb",
     f"{int(eur.lead_windows)}; {int(eur.windows_with_cs)}; {int(eur.windows_with_cs_pos500kb)}")

# ---- O2 Fig 2D: all non-European studies pooled
dr = tsv(f"{H58}/descriptive_refits.tsv").set_index("label").loc["nonEUR_pooled"]
emit("O2", f"{H58}/descriptive_refits.tsv", "label=nonEUR_pooled: fraction_direct; fraction_enzyme; difference [lo, hi]; clusters d/e",
     f"{dr.fraction_direct:.4f}; {dr.fraction_enzyme:.4f}; {dr.difference:.4f} [{dr.diff_boot_lo:.4f}, {dr.diff_boot_hi:.4f}]; "
     f"{int(dr.n_clusters_direct)}/{int(dr.n_clusters_enzyme)}")

# ---- O3 Fig 3C: GSE135251 disease units by genetic similarity
anc = tsv(H79).set_index(["dataset", "group_binary"])
dis = anc.loc[("GSE135251", "Disease")]
ctl = anc.loc[("GSE135251", "Control")]
emit("O3", H79, "GSE135251 Disease EUR of total; Control AFR/EUR of total",
     f"{int(dis.EUR)} of {int(dis.sum())}; {int(ctl.AFR)}/{int(ctl.EUR)} of {int(ctl.sum())}")

# ---- O4 Fig 7E: GSE83452 resolved vs persisted by arm
res = json.loads((ROOT / GSE83452).read_text())["histology"]["P2_by_arm"]
emit("O4", GSE83452, "histology.P2_by_arm", json.dumps(res))

# ---- O5 Fig 7F: smallest held-out study
sm = tsv(STUDY).drop_duplicates("held_out_study")[["held_out_study", "rows", "donors"]]
g189 = sm[sm.held_out_study == "GSE189600"].iloc[0]
tot_cells = int(sm.rows.sum())
emit("O5", STUDY, "GSE189600 rows; donors; all-study rows", f"{int(g189.rows)}; {int(g189.donors)}; {tot_cells}")

# ---- O6 Fig 6B: row and column totals
f6s = json.loads((ROOT / F6 / "f6_repair_summary.json").read_text())["A_gen_untestable_by_reason"]
mx = tsv(f"{F6}/fig6b_call_matrix_3x3_source.tsv")
rna_unt = int(mx[mx.rna_call == "untestable"].n_genes.sum())
gen_unt = int(mx[mx.genetics_call == "untestable"].n_genes.sum())
emit("O6", f"{F6}/f6_repair_summary.json", "A_gen_untestable_by_reason", json.dumps(f6s))
emit("O6", f"{F6}/fig6b_call_matrix_3x3_source.tsv", "sum n_genes rna_call=untestable; genetics_call=untestable",
     f"{rna_unt}; {gen_unt}")
assert gen_unt == f6s["underpowered_source"] + f6s["insufficient_shared_posterior"]

for r in rows:
    print("\t".join(r))

# ---------------------------------------------------------------- staged RESULTS edits
text = RESULTS.read_text()
edits = []


def stage(eid, anchor, addition, reason, source, numbers, reread, position="after"):
    count = text.count(anchor)
    new = anchor + addition if position == "after" else addition + anchor
    edits.append({"id": eid, "file": "docs/RESULTS.md", "old": anchor, "new": new, "reason": reason,
                  "source": source, "numbers_in_new": numbers, "value_reread": reread, "old_count": count})


def pct(x):
    return f"{100 * x:.2f}%"


stage("N-01",
      "testing family, 2,426 coordinated beyond the leading residual component (650\nsame-lineage, 1,776 cross-lineage).",
      f" Modularity of the {n_comm} fixed communities is {mod_obs:.3f} against null means of"
      f" {mod_rank:.3f} ({n_null['rank_matched_permutation']} rank-matched label permutations) and"
      f" {mod_rew:.3f} ({n_null['degree_preserving_rewiring']} degree-preserving rewirings); one community"
      f" holds 15 programs, {n_fib} fibroblast and {n_chol} cholangiocyte (`interaction_summary.tsv`,"
      f" `modularity_null.tsv`, `program_communities.tsv` in the same directory).",
      "title.md ~56-57 quotes 0.138 / 0.089 nulls and a 15-program community; RESULTS has only the 0.358 discovery value",
      f"{T4}/interaction_summary.tsv (modularity_observed, modularity_rank_matched_null_mean, modularity_rewired_mean); "
      f"{T4}/modularity_null.tsv (rows per null); {T4}/program_communities.tsv (community {c15})",
      [f"{n_comm}", f"{mod_obs:.3f}", f"{mod_rank:.3f}", f"{mod_rew:.3f}", str(n_null['rank_matched_permutation']),
       str(n_null['degree_preserving_rewiring']), "15", str(n_fib), str(n_chol)],
      f"{mod_obs:.6f}; {mod_rank:.6f}; {mod_rew:.6f}; nulls {n_null}; community {c15}: fib {n_fib}, chol {n_chol}")

stage("N-02",
      "cluster fractions and one PDFF credible set at chr11:823,586 carries 22.5%;",
      f" the PNPLA3 cluster holds {pn_n} direct credible sets with a mean protein-altering plus"
      f" canonical-splice fraction of {pn_mean:.3f} (recomputed from"
      f" `{F2RUN}/multiplicity/corrected_locus_coordinates.tsv` with the producer's 1-Mb clustering);",
      "fig2.md ~126 quotes the PNPLA3 cluster (124 direct sets, mean 0.117); no table on disk holds it",
      f"{F2RUN}/multiplicity/corrected_locus_coordinates.tsv, clustered with snapshot fig2_multiplicity_correction.assign_clusters, "
      f"cluster {pn.index[0]} (chr22 {int(pn.lo.iloc[0])}-{int(pn.hi.iloc[0])})",
      [str(pn_n), f"{pn_mean:.3f}"], f"n {pn_n}; mean {pn_mean:.6f}; top-5 share {top5:.4f}")

stage("N-03",
      "post hoc and does not correct dependence among genes, shared loci or gene\nsets; it is not confirmatory evidence of biological separation or opposition.",
      f" Flipping signs jointly within the {n_loci} H10 loci that hold the {n_pairs} matched genetic genes gives"
      f" nominal p={p_locus:.5f} (100,000 draws, seed 20260907), and {n_pos_loci} of {n_loci} locus means"
      f" are positive (sign test p={p_locus_sign:.3f}); both are post hoc"
      f" (`{R14}/r14_10_recompute.log`).",
      "fig5_discussion.md ~207 quotes the locus-clustered sign-flip (216 loci, P=0.01606; 99 of 216, P=0.247); RESULTS lacks it",
      f"{R14}/r14_10_recompute.log (seed 20260907 locus signflip p; locus sign test p); {R14}/pairs_with_locus_recomputed.tsv (locus, d)",
      [str(n_loci), str(n_pairs), f"{p_locus:.5f}", "100,000", "20260907", str(n_pos_loci), f"{p_locus_sign:.3f}"],
      f"p {p_locus:.6f}; loci {n_loci}; positive {n_pos_loci}; pairs {n_pairs}; sign test {p_locus_sign:.4f}")

g189_pct = 100 * int(g189.rows) / tot_cells
stage("N-04",
      "[Source tables and checks](../Analysis/MASLD_Model_Benchmark/executions/fig7-source-exposure-correction-20261014/README.md).",
      f"\n\n**Figure 7B, 7E and 7F legend values.** 7B: across the {int(hy.m)} displayed reporter-assay"
      f" rows ({f7.model.nunique()} models), HyenaDNA gains {hy.gain:.6f} Spearman over the allele-identity"
      f" ridge (95% interval {hy.lo95:.6f} to {hy.hi95:.6f}; Bonferroni-corrected percentile lower bound"
      f" {hy.bonf_pct_lo:.6f} from the same {int(hy.n_valid):,} draws); it is the only row whose corrected"
      f" lower bound is above zero (`{MPRA}`). 7E: in GSE83452 the diet arm has {res['Diet_resolved']}"
      f" resolved and {res['Diet_persisted']} persisted participants and the surgery arm"
      f" {res['BS_resolved']} and {res['BS_persisted']} (`{GSE83452}`, `histology.P2_by_arm`). 7F: the"
      f" smallest held-out study, GSE189600, contributes {int(g189.donors)} donors and {int(g189.rows)} of"
      f" {tot_cells:,} cells ({g189_pct:.1f}%) (`{STUDY}`).",
      "fig7.md ~27 and legend 7B quote the corrected lower bound 0.009; legend 7E/7F quote 6 vs 16 and 149 cells; RESULTS lacks all three",
      f"{MPRA} (family F7_displayed_rows, model hyenadna: m, gain, lo95, hi95, bonf_pct_lo, n_valid); "
      f"{GSE83452} (histology.P2_by_arm); {STUDY} (held_out_study GSE189600 rows, donors; sum of rows)",
      [str(int(hy.m)), str(f7.model.nunique()), f"{hy.gain:.6f}", f"{hy.lo95:.6f}", f"{hy.hi95:.6f}",
       f"{hy.bonf_pct_lo:.6f}", f"{int(hy.n_valid):,}", str(res['Diet_resolved']), str(res['Diet_persisted']),
       str(res['BS_resolved']), str(res['BS_persisted']), str(int(g189.donors)), str(int(g189.rows)),
       f"{tot_cells:,}", f"{g189_pct:.1f}%"],
      f"bonf_pct_lo {hy.bonf_pct_lo:.6f}; rows above zero {n_pos_bound}; diet {res['Diet_resolved']}/{res['Diet_persisted']}; "
      f"surgery {res['BS_resolved']}/{res['BS_persisted']}; GSE189600 {int(g189.rows)} cells / {int(g189.donors)} donors of {tot_cells}")

stage("N-05",
      "coverage (enzyme 8 of 767 lead windows, direct 47 of 47)",
      "",
      "legend Fig 2A ~39 quotes 23 of 767 within 500 kb; RESULTS gives only 8 of 767",
      f"{H58}/G7_coverage.tsv (EUR tier2_liver_enzyme: lead_windows, windows_with_cs, windows_with_cs_pos500kb)",
      [str(int(eur.windows_with_cs)), str(int(eur.lead_windows)), str(int(eur.windows_with_cs_pos500kb))],
      f"{int(eur.windows_with_cs)} / {int(eur.windows_with_cs_pos500kb)} of {int(eur.lead_windows)}")
edits[-1]["new"] = (f"coverage (enzyme {int(eur.windows_with_cs)} of {int(eur.lead_windows)} lead windows by locus identifier,"
                    f" {int(eur.windows_with_cs_pos500kb)} of {int(eur.lead_windows)} within 500 kb by position"
                    f" (`{H58}/G7_coverage.tsv`); direct 47 of 47)")

stage("N-06",
      "[-0.0489, 0.1065] (95% intervals shown); no call for any member.",
      f" All non-European studies pooled ({int(dr.n_clusters_direct)} direct and {int(dr.n_clusters_enzyme)}"
      f" enzyme clusters) give {pct(dr.fraction_direct)} against {pct(dr.fraction_enzyme)}, difference"
      f" {dr.difference:.4f} [{dr.diff_boot_lo:.4f}, {dr.diff_boot_hi:.4f}], a descriptive value with no"
      f" calibration control (`{H58}/descriptive_refits.tsv`, row nonEUR_pooled).",
      "legend Fig 2D ~92-93 and fig2.md quote non-European pooled 7.0% / 3.7%, 0.0331 [-0.0383, 0.1274]; RESULTS lacks it",
      f"{H58}/descriptive_refits.tsv (label nonEUR_pooled: n_clusters_direct, n_clusters_enzyme, fraction_direct, fraction_enzyme, difference, diff_boot_lo, diff_boot_hi)",
      [str(int(dr.n_clusters_direct)), str(int(dr.n_clusters_enzyme)), pct(dr.fraction_direct), pct(dr.fraction_enzyme),
       f"{dr.difference:.4f}", f"{dr.diff_boot_lo:.4f}", f"{dr.diff_boot_hi:.4f}"],
      f"{dr.fraction_direct:.6f}; {dr.fraction_enzyme:.6f}; {dr.difference:.6f} [{dr.diff_boot_lo:.6f}, {dr.diff_boot_hi:.6f}]")

stage("N-07",
      "EUR controls; 6 of its 10 controls are of African genetic similarity).",
      "",
      "legend Fig 3C ~190-191 quotes 171 of 197 disease units of European similarity; RESULTS gives only the controls",
      f"{H79} (dataset GSE135251, group Disease: EUR and row total)",
      [str(int(ctl.EUR)), str(int(ctl.AFR)), str(int(ctl.sum())), str(int(dis.EUR)), str(int(dis.sum()))],
      f"control AFR {int(ctl.AFR)} EUR {int(ctl.EUR)} of {int(ctl.sum())}; disease EUR {int(dis.EUR)} of {int(dis.sum())}")
edits[-1]["new"] = (f"EUR controls; {int(ctl.AFR)} of its {int(ctl.sum())} controls are of African genetic similarity,"
                    f" against {int(dis.EUR)} of its {int(dis.sum())} disease units of European similarity;"
                    f" `{H79}`).")
assert int(ctl.EUR) == 2

stage("N-08",
      "indeterminate and genetics untestable 8,621. Figure 6D bars:",
      "",
      "legend Fig 6B ~590-592 and 6D ~610 quote 9,312 and 4,232; RESULTS gives the cells but not these totals",
      f"{F6}/f6_repair_summary.json (A_gen_untestable_by_reason); {F6}/fig6b_call_matrix_3x3_source.tsv (n_genes by rna_call / genetics_call)",
      [f"{rna_unt:,}", f"{gen_unt:,}", f"{f6s['underpowered_source']:,}", str(f6s["insufficient_shared_posterior"])],
      f"RNA untestable {rna_unt}; genetics untestable {gen_unt} = {f6s['underpowered_source']} + {f6s['insufficient_shared_posterior']}")
edits[-1]["new"] = (f"indeterminate and genetics untestable 8,621. Totals: RNA untestable {rna_unt:,}; genetics"
                    f" untestable {gen_unt:,}, of which {f6s['underpowered_source']:,} have no Tier-1/2 COLOC record and"
                    f" {f6s['insufficient_shared_posterior']} are SuSiE-untestable for insufficient shared posterior"
                    f" (`{F6}/f6_repair_summary.json`, `fig6b_call_matrix_3x3_source.tsv`). Figure 6D bars:")

# anchors must be unique and must not overlap one another
for e in edits:
    assert e["old_count"] == 1, (e["id"], e["old_count"])
spans = sorted((text.index(e["old"]), text.index(e["old"]) + len(e["old"]), e["id"]) for e in edits)
for a, b in zip(spans, spans[1:]):
    assert a[1] <= b[0], (a, b)
(WORK / "nosource_edits.json").write_text(json.dumps(edits, indent=1, ensure_ascii=False) + "\n")
