"""Fork C: needs_user items 9, 28, 29, 30, 31. Builds part_C.json against the current live files."""
import json, os, re
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/open_items-20261017T105826Z")
LEG = "docs/manuscript/05_figure_legends.md"; RES = "docs/RESULTS.md"
TAB = "Analysis/MASLD_Model_Benchmark/executions/tabicl-nonlinear-comparator-20260905T000000Z/results"
F6 = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F6/repair_20261015T172841Z"
E = []
def live(p): return open(os.path.join(ROOT, p), encoding="utf-8").read()
def add(id_, f, old, new, reason, decision, src, prod):
    pat = r"\s+".join(re.escape(w) for w in old.split())
    hits = [m.group(0) for m in re.finditer(pat, live(f))]
    assert len(hits) == 1, (id_, len(hits))
    E.append(dict(id=id_, file=f, old=hits[0], new=new, reason=reason, decision=decision, source_ids=src, producing_file=prod))

add("C-09-legend", LEG,
    "differs from ridge by −0.014 to +0.001 Spearman across four comparisons in 109 GSE268273 participants, against a minimum detectable effect of 0.067; none exceeds it,",
    "differs from ridge by −0.014 to +0.001 Spearman across the four prespecified comparisons in 109 GSE268273 participants, against a minimum detectable effect of 0.067 (not drawn in this panel); none exceeds it,",
    "Values match the producer (4 TabICLv2 vs ridge deltas on the same 30 or 100 PCs, GSE268273: -0.0085, -0.0039, -0.0135, +0.0014). The comparison is not drawn in 7D; the legend now says so.",
    "item 9", ["needs_user:9"], f"{TAB}/results.json (deltas: tabicl30/100 vs ridge30/100 rawy and em)")
add("C-09-results", RES,
    "| Generic nonlinearity on the severity task (TabICLv2 vs ridge, 8 comparisons) | mean **+0.0031**, all below MDE 0.0667 | 653 participants, three cohorts; `MODELS.md` |",
    "| Generic nonlinearity on the severity task (TabICLv2 vs ridge on the same 30 or 100 training-fold PCs) | Prespecified GSE268273 test (109 participants): four comparisons, **−0.014 to +0.001**. Eight comparisons against ridge on raw stage (two in GSE268273, six post hoc in GSE213621 under 3-level and 4-level coding and in GSE276114): mean **+0.0031**, range −0.014 to +0.027. All below MDE 0.0667 | GSE268273 109, GSE213621 299 and 367 records (two codings of one cohort), GSE276114 177 participants; `executions/tabicl-nonlinear-comparator-20260905T000000Z/results/results.json` (deltas) and `secondary_deltas.json` (post hoc) |",
    "The 8-comparison mean +0.0031 reproduces from results.json + secondary_deltas.json; it mixes 2 prespecified and 6 post hoc deltas, and '653 participants, three cohorts' counts GSE213621 only at 367 records. The legend's four comparisons are the prespecified test. Both are now stated with their producer.",
    "item 9", ["needs_user:9"], f"{TAB}/results.json; {TAB}/secondary_deltas.json")

add("C-28-results", RES,
    "two programs (Stromal ECM and Ductular injury) in seven datasets from the 2026-08-07 candidate.",
    "two programs (Stromal ECM and Ductular injury) in seven datasets from the 2026-08-07 candidate. Liver protein (PXD051911) is untestable for both programs under the Figure 5E scoring rule (at least eight measured members and at least 20% of L1 weight retained): Stromal ECM (IGFBP7) has 8 measured members carrying 16.6% of its weight, and Ductular injury (BICC1) has 2 carrying 9.0% (`figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv`).",
    "Adds the RESULTS line behind the 6A legend clause 'liver protein (PXD051911) was untestable for both'.",
    "item 28", ["needs_user:28"], "figures/main/fig5_molecular_context/panels/data/panel_deep_20261017/5E.tsv (n_measured, retained_l1_weight, testable)")

add("C-29-results", RES,
    "(two POLR1HASP, two ABCC6P2 and one CES1P1), with no quarantine record.",
    "(two POLR1HASP, two ABCC6P2 and one CES1P1), with no quarantine record. Counted by symbol instead, 103 of the 114 share their symbol with another Ensembl ID tested in the bulk fit (the count in the Figure 6B legend).",
    "Both splits come from one producer: 109/5 by quarantine reason and 103 of 114 by shared symbol (C_bulk_absent_symbol_shared_in_bulk).",
    "item 29", ["needs_user:29"], f"{F6}/f6_repair_summary.json (C_bulk_absent_by_reason, C_bulk_absent_symbol_shared_in_bulk); {F6}/C_R07_05_bulk_ids_absent.tsv")
add("C-29-legend", LEG,
    "103 of which share a symbol with another ID.",
    "103 of which share their symbol with another Ensembl ID tested in the bulk fit.",
    "Names the comparison set for the 103.", "item 29", ["needs_user:29"], f"{F6}/f6_repair_summary.json")

add("C-30-results", RES,
    "membership or weight changes through these joins.",
    "membership or weight changes through these joins.\n\nFigure 6C worked entries (r8 classes, Catalog r9d values; `figures/main/fig6_gene_catalog/panels/data/panel_deep_20261017/6C.tsv`). GLP1R is untested: no Tier-1/2 SuSiE result (best ABF PP.H4 0.022, MVP NAFLD EAS) and not detected in the bulk fit. HKDC1 is supported in both domains: SuSiE PP.H4 0.988 (UK Biobank AST; > 0.5 in 3 of 5 Tier-1/2 studies with a SuSiE result) and log2FC +1.07 (95% CI 0.68 to 1.47; BH q 8.0e-7). MTARC1 is unresolved: no Tier-1/2 SuSiE result (best ABF PP.H4 0.046, FinnGen NASH) and log2FC −0.15 (−0.22 to −0.08; q 1.2e-4); its fine-mapped missense variant (PIP 0.76 in a PDFF GWAS; PIP ≥ 0.5 in 3 studies) does not enter the genetic call. THRB is genetically anchored: SuSiE PP.H4 > 0.999 (UK Biobank GGT; > 0.5 in 1 of 4 studies) and log2FC −0.20 (−0.26 to −0.14; q 4.2e-10).",
    "Adds RESULTS lines for the 6C worked entries. The installed 6C.tsv text field prints '> 0.999' for THRB (SuSiE 0.999967), so the old sidecar '1.000 SuSiE' no longer reaches the panel.",
    "item 30", ["needs_user:30"], "figures/main/fig6_gene_catalog/panels/data/panel_deep_20261017/6C.tsv; catalog_r9d_gene_table.parquet (rna_bh_q)")

add("C-31-legend", LEG,
    "In GSE276114, 9 or 10 of the 81 MASLD and 24 of the 82 chronic viral hepatitis participants are patients with HCC sampled from peritumoral liver, with no per-sample label.",
    "In GSE276114, 9 or 10 of the 81 MASLD and 24 of the 82 chronic viral hepatitis participants are patients with HCC sampled from peritumoral liver at least 1 cm from the tumor, with no per-sample label. These counts come from the source study (34 participants with HCC: 24 with chronic viral hepatitis and 10 with MASLD); one MASLD sample is not deposited.",
    "Counts checked: GEO deposit has 81 MASLD and 82 CVH participants; the source study reports 34 HCC split 70.6%/29.4% = 24 CVH / 10 MASLD; one of 82 MASLD samples is withheld, so 9 or 10. No per-sample label exists.",
    "item 31", ["needs_user:31"], "Analysis/Proteomics/results/gse276114_disease_metadata.csv (disease counts); Analysis/MASLD_Model_Benchmark/MODEL_RESULTS_SUMMARY.md:3135-3176 (source-study HCC counts)")
add("C-31-results", RES,
    "No demonstrated advantage over the recorded transductive in-cohort ridge; intervals including zero do not establish equivalence |",
    "No demonstrated advantage over the recorded transductive in-cohort ridge; intervals including zero do not establish equivalence. In GSE276114, 9 or 10 of the 81 MASLD and 24 of the 82 chronic viral hepatitis participants have HCC (peritumoral tissue at least 1 cm from the tumor; source-study counts, no per-sample label; `MODEL_RESULTS_SUMMARY.md:3135-3176`) |",
    "Adds the RESULTS line behind the legend HCC counts.", "item 31", ["needs_user:31"],
    "Analysis/MASLD_Model_Benchmark/MODEL_RESULTS_SUMMARY.md:3135-3176; Analysis/Proteomics/results/gse276114_disease_metadata.csv")

# order-independent checks: unique, non-overlapping, no banned terms added
BAN = re.compile(r"—|\bfrozen\b|\bprogression\b|\bcascade\b", re.I)
for f in {e["file"] for e in E}:
    t = live(f); spans = sorted((t.find(e["old"]), t.find(e["old"]) + len(e["old"]), e["id"]) for e in E if e["file"] == f)
    for a, b in zip(spans, spans[1:]): assert b[0] >= a[1], (a, b)
for e in E:
    assert len(BAN.findall(e["new"])) <= len(BAN.findall(e["old"])), e["id"]
json.dump(E, open(os.path.join(OUT, "part_C.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print(len(E), "edits")
