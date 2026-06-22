// H10 — Deep Literature & Competitive Benchmarking workflow.
// Per method: multi-modal source sweep (PubMed + bioRxiv + WebSearch) -> defensibility memo with citations
// -> citation-accuracy check -> competitor head-to-heads -> synthesis (dossier + actionable findings).
export const meta = {
  name: 'megareview-H10-literature',
  description: 'H10: deep external-literature benchmarking — per-method defensibility + better-method scan + paper-ready citations + competitor head-to-heads',
  phases: [
    { title: 'Research', detail: 'per-method multi-modal literature sweep' },
    { title: 'Competitors', detail: 'competitor paper head-to-heads' },
    { title: 'CiteCheck', detail: 'verify each cited PMID/DOI supports its claim' },
    { title: 'Synthesize', detail: 'literature dossier + actionable findings' },
  ],
}

const ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'
const CTX = ROOT + '/docs/audit/2026-06-13_megareview/REVIEW_CONTEXT.md'
const PRE = `You are a methods-literature analyst on the MASLD mega-review, TEAM H10. FIRST Read ${CTX} for project context (methods, the v3 thesis, the known vulnerabilities). You do DEEP EXTERNAL LITERATURE research. Use ToolSearch to load: the PubMed MCP tools (mcp__claude_ai_PubMed__search_articles, mcp__claude_ai_PubMed__get_article_metadata, mcp__claude_ai_PubMed__find_related_articles), the bioRxiv MCP tools (mcp__claude_ai_bioRxiv__search_preprints, mcp__claude_ai_bioRxiv__search_published_preprints), and WebSearch/WebFetch. Prefer primary methods papers (2018-2026) and recent (2024-2026) alternatives. Every citation MUST be a real PMID or DOI you actually retrieved — NEVER fabricate. READ-ONLY; no repo edits.`

const A = typeof args === 'string' ? JSON.parse(args) : (args || {})
const _agent = agent
const tryAgent = async (p, o) => { let r = null; for (let i = 0; i < 4 && r == null; i++) { r = await _agent(p, o) } return r }

const STR = { type:'string' }; const NUM = { type:'number' }
const ENUM = (e) => ({ type:'string', enum:e }); const ARR = (items) => ({ type:'array', items })
const S = (props, req) => ({ type:'object', additionalProperties:false, required:req, properties:props })

const CITE = S({ ref:STR, type:ENUM(['PMID','DOI','URL']), supports:STR }, ['ref','supports'])
const MEMO_SCHEMA = S({
  method:STR,
  current_choice:STR,
  verdict:ENUM(['defensible','defensible_with_caveats','switch_recommended','uncertain']),
  summary:STR,
  better_or_newer_methods:ARR(STR),
  switch_cost:STR,
  reviewer_risk:STR,
  positioning:STR,
  citations:ARR(CITE)
}, ['method','verdict','summary','citations'])
const COMP_SCHEMA = S({
  paper:STR,
  our_advantages:ARR(STR),
  their_advantages:ARR(STR),
  parity_or_overlap:ARR(STR),
  novelty_verdict:STR,
  citations:ARR(CITE)
}, ['paper','our_advantages','their_advantages','novelty_verdict'])
const CITECHECK_SCHEMA = S({ checks: ARR(S({ ref:STR, exists:ENUM(['verified','not_found','uncertain']), supports_claim:ENUM(['yes','no','partial','uncertain']), note:STR }, ['ref','exists','supports_claim'])) }, ['checks'])
const FINDING_ITEM = S({ title:STR, severity:ENUM(['Critical','High','Medium','Low']), type:ENUM(['bug','rigor','repro','literature']), file:STR, line:STR, claim:STR, evidence:STR, recommendation:STR, lit_question:STR }, ['title','severity','type','file','claim','evidence','recommendation'])
const REPORT_SCHEMA = S({ executive_summary:STR, tally:S({critical:NUM,high:NUM,medium:NUM,low:NUM},['critical','high','medium','low']), markdown_report:STR, findings:ARR(FINDING_ITEM) }, ['executive_summary','tally','markdown_report','findings'])

// ---- Method clusters to benchmark ----
const METHODS = A.methods || [
  {k:'bulk-de', label:'Multi-cohort bulk DE: limma-voom-qw fixed-effect vs dream/DESeq2/metafor', q:'Is pooled limma-voom quality-weighted with a FIXED dataset effect (vs dream random-effects or metafor random-effects meta-analysis) the right canonical for a 9-cohort MASLD mega-analysis? Is ashr lfsr<0.05 & |shrunk_logFC|>0.5 standard? Is a "rank-of-ranks over LOCO + external metrics" pre-registered method-selection a recognized paradigm? What do recent (2024-26) RNA-seq mega-analysis best-practice papers recommend?'},
  {k:'coloc', label:'Colocalization: coloc.abf vs coloc.susie + cross-ancestry LD', q:'Best practice for multi-ancestry colocalization: coloc.susie vs coloc.abf, SuSiEx/MESuSiE, and the effect of MISMATCHED cross-ancestry LD (EUR eQTL LD vs EAS/AFR GWAS LD) on PP.H4 calibration (Wallace 2020/2021). Is colocalizing an EUR-LD eQTL SuSiE fit against an EAS-LD GWAS SuSiE fit defensible? PolyFun LD reference appropriateness.'},
  {k:'twas', label:'TWAS: OTTERS vs S-PrediXcan/FUSION + INTACT/cTWAS/HyPrColoc', q:'Is OTTERS (multi-method ACAT) appropriate vs S-PrediXcan/FUSION for liver TWAS? Is single-method P+T (lassosum empty) defensible as "multi-method"? INTACT and cTWAS calibration. Cross-ancestry TWAS validity when eQTL is EUR and GWAS is non-EUR.'},
  {k:'sc-integration', label:'scRNA integration: scVI vs Harmony/ScaleSC + batch granularity', q:'scVI vs Harmony vs ScaleSC for a 7-dataset 1.2M-cell liver atlas; sample-level vs dataset-level batch key; reproducibility/seed requirements. Recent (2024-26) atlas-integration benchmarks (e.g. scIB).'},
  {k:'hotspot', label:'Hotspot autocorrelation modules (DANB Z sample-size inflation)', q:'DeTomaso & Yosef Hotspot DANB autocorrelation Z scales with cell count; at 650K cells the Z>=7 FDR threshold is saturated (100% pass). Is a top-N-by-Z cap defensible? Should an autocorrelation effect-size (C) floor be used instead? What do module-discovery best-practice papers say?'},
  {k:'cnmf-dialogue', label:'cNMF k-selection + DIALOGUE multicellular programs', q:'cNMF (Kotliar 2019) k-selection: cophenetic vs biological-resolution; is k=6 over k=4 (cophenetic-optimal) defensible? DIALOGUE (Jerby-Arnon 2022): is computing multicellular programs on duplicated pseudobulk (not single cells) valid, or does it pseudoreplicate?'},
  {k:'ccc', label:'Differential cell-cell communication: LIANA vs CellChat/CellPhoneDB/NicheNet', q:'Best practice for DIFFERENTIAL (disease-vs-control) cell-cell communication with donor replication. Is a single pooled-cell point estimate with no donor replicates / permutation / p-value valid? How do LIANA/CellChat/CrossTalkeR handle donor-level inference and the pseudoreplication problem (Squair 2021 pseudobulk)?'},
  {k:'deconvolution', label:'Bulk + spatial deconvolution: BayesPrism/MuSiC/cell2location', q:'BayesPrism vs MuSiC for liver bulk deconvolution; cell2location for Visium. Resolution limits (only 2/17,512 spots >=80% hepatocyte) — is Visium deconvolution adequate to claim within-cell-type sub-states, or is Xenium/MERFISH needed?'},
  {k:'pseudotime', label:'Pseudotime: Palantir/DPT/CytoTRACE/PAGA consensus', q:'Validity of 4-method consensus pseudotime; root-cell selection; donor-level vs cell-level trajectory metrics; sensitivity of pseudotime-vs-stage correlations to atlas contamination.'},
  {k:'nmf-ordinal', label:'NMF k-program subtyping + ordinal staging (CORN/CORAL)', q:'NMF reproducibility standards (seed stability, ARI/cophenetic); is a k=6 decomposition with ARI=0.295 across seeds publishable as canonical subtypes? Ordinal CORN/CORAL for fibrosis staging; concept-bottleneck models.'},
  {k:'spatial-svg', label:'Spatial SVG: Moran I vs SPARK/SpatialDE + spot pseudoreplication', q:'Moran-I SVG detection vs SPARK-X/SpatialDE; spot-level vs sample/donor-level inference (6,546 spots from 2-3 donors = pseudoreplication); MT/RP confounder removal before SVG. Best-practice 2024-26.'},
  {k:'regulon-heritability', label:'chromVAR/SCENIC+ regulons + LDSC cell-type heritability', q:'chromVAR per-cell motif-deviation testing (donor aggregation vs pseudoreplication); SCENIC+ regulon-activity DE power at n=18; LDSC-SEG cell-type partitioned heritability interpretation (NK vs hepatocyte enrichment).'},
  {k:'drug-reversal', label:'Drug repurposing: LINCS/CGP connectivity reversal + network proximity', q:'Connectivity-map (LINCS L1000) signature-reversal validity for drug repurposing; network-proximity (Guney 2016); is a stage-specific reversal claim robust? Best-practice for transcriptomic drug repurposing 2024-26.'},
  {k:'fstage-transfer', label:'Ordinal F-stage label transfer from scRNA + leakage-free eval', q:'Transferring ordinal Kleiner fibrosis stage to unlabeled donors via scVI-latent kNN/anchors; leakage-free evaluation (anchors must not leak into held-out); is QWK 0.76 from endpoint-anchored training honest if jackknife gives 0.286 and held-out gives 0.0? Best-practice for ordinal label transfer.'},
  {k:'convergence', label:'Multi-evidence convergence scoring (heuristic vs Open Targets/Bayesian)', q:'Is an evidence-weighted HEURISTIC convergence score (explicitly not a posterior) a recognized target-prioritization paradigm vs Open Targets locus-to-gene, PoPS, or Bayesian integration? Per-gene permutation FDR validity; the Mountjoy 2021 / Schipper 2025 held-out paradigm.'},
]

// ---- Competitor papers ----
const COMPETITORS = A.competitors || [
  {k:'feng', label:'Feng et al. 2026 (bioRxiv) — MASLD meta-analysis', q:'Compare our pooled cohort-adjusted mega-analysis + COLOC + deconvolution + multi-ancestry to Feng 2026 (per-study meta-analysis, no COLOC/deconvolution/multi-ancestry). Our advantages, their advantages, novelty.'},
  {k:'tzouanas', label:'Tzouanas et al. 2026 (Cell) — mouse-longitudinal + Xenium multiome', q:'Compare to Tzouanas 2026 (scRNA+snATAC multiome + Xenium, mouse-first, MATCHA, 4 hep stress programs). What do they have that we lack (Xenium spatial resolution, multiome) and vice versa (mega-analysis, COLOC, drug repurposing)?'},
  {k:'li', label:'Li et al. 2025 (Nat Genet) — 61-liver scRNA+Visium+MALDI', q:'Compare to Li 2025 (61 livers, MITF LAM master TF, RSPO3-LGR6 axis, no GWAS/drug/cross-species). EPHA2/TM6SF2 overlap; our scale + causal + drug advantages.'},
  {k:'govaere', label:'Govaere et al. 2020 STM + 2023 Nat Metab — signature + proteo-transcriptomic', q:'Compare to Govaere 2020 (25-gene signature, single-cohort) and 2023 (SomaScan proteo-transcriptomic, NMF 3 subtypes). We integrate their GeoMx/CosMx — positioning.'},
  {k:'broadaway', label:'Broadaway et al. 2024 (AJHG) — source liver eQTL panel', q:'Our COLOC uses their N=1,183 eQTL panel. Their colocalization = 747 eGenes (lipid-dominated, EUR-only, no PIP); ours 751 share 36% (269). Is our 482 novel-to-us defensible (68% need data they lack: EUR AST / non-EUR / MASLD GWAS)? Apples-to-apples ABF overlap 51%.'},
  {k:'seagle-kuchenhoff', label:'Seagle 2025 (AJHG) multi-ancestry GWAS + Kuchenhoff 2026 cross-organ fibrosis', q:'Seagle 2025 (multi-ancestry GWAS + MR drug repurposing, FADS1/S1PR2); Kuchenhoff 2026 (953 samples 4-organ fibrosis scRNA; liver has FEWEST consensus genes, AUROC 0.60). Positioning vs our multi-ancestry COLOC + cross-species.'},
  {k:'others', label:'Piras 2024, Suppli 2019, Wang 2021, Kozumi 2021, Hoang 2019, LITMUS 2024', q:'Brief head-to-heads: Piras (GeneMeta 10-dataset meta + MEGENA), Suppli (SHH+ hep, male-only controls), Wang (NPC scRNA, RUNX1/CREB3L1, CCl4 not metabolic), Kozumi (THBS2 Japanese), Hoang (ordinal gNAS/gFib), Vacca/LITMUS (mouse-model ranking). Our differentiators vs each.'},
]

log(`H10: ${METHODS.length} method clusters + ${COMPETITORS.length} competitor head-to-heads`)

// ---- Phase: Research (per-method defensibility memos) ----
phase('Research')
const memos = (await parallel(METHODS.map(m => () =>
  tryAgent(`${PRE}\nMETHOD CLUSTER: "${m.label}".\nResearch question: ${m.q}\nProduce a defensibility memo: (1) the project's current choice; (2) verdict (defensible / defensible_with_caveats / switch_recommended / uncertain); (3) a concise evidence-based summary; (4) better/newer methods if any, with the switch cost; (5) the reviewer risk if we keep the current choice; (6) competitive positioning; (7) >=5 REAL citations (PMID/DOI you actually retrieved via PubMed/bioRxiv/WebSearch) each with the specific claim it supports. Do genuine multi-source searches before writing.`,
    {schema:MEMO_SCHEMA, label:`research:${m.k}`, phase:'Research'})
))).filter(Boolean)

// ---- Phase: Competitors ----
phase('Competitors')
const comps = (await parallel(COMPETITORS.map(c => () =>
  tryAgent(`${PRE}\nCOMPETITOR HEAD-TO-HEAD: "${c.label}".\n${c.q}\nProduce: our_advantages, their_advantages, parity/overlap, a novelty verdict (what we genuinely add to the field given this competitor), and REAL citations (PMID/DOI). Search to confirm the competitor's actual scope/claims; do not rely on memory alone.`,
    {schema:COMP_SCHEMA, label:`competitor:${c.k}`, phase:'Competitors'})
))).filter(Boolean)

// ---- Phase: Citation accuracy check ----
phase('CiteCheck')
const allCitedRefs = []
for (const m of memos) for (const c of (m.citations||[])) allCitedRefs.push({src:`method:${m.method}`, ...c})
for (const c of comps) for (const ci of (c.citations||[])) allCitedRefs.push({src:`competitor:${c.paper}`, ...ci})
// chunk refs into ~12 groups for parallel verification
const CHUNKS = 12
const groups = Array.from({length:CHUNKS}, () => [])
allCitedRefs.forEach((r,i) => groups[i % CHUNKS].push(r))
const citeChecks = (await parallel(groups.filter(g=>g.length).map((g,i) => () =>
  tryAgent(`${PRE}\nCITATION-ACCURACY CHECK. For EACH reference below, use the PubMed/bioRxiv MCP tools (or WebSearch for DOIs/URLs) to verify (a) the reference EXISTS (verified/not_found/uncertain) and (b) whether it actually SUPPORTS the stated claim (yes/no/partial/uncertain). Flag any fabricated/hallucinated citation as not_found. References:\n${JSON.stringify(g)}`,
    {schema:CITECHECK_SCHEMA, label:`citecheck:${i}`, phase:'CiteCheck'})
))).filter(Boolean).flatMap(r=>r.checks)

// ---- Phase: Synthesize ----
phase('Synthesize')
const fabricated = citeChecks.filter(c=>c.exists==='not_found')
const report = await tryAgent(`${PRE}\nYou are the SYNTHESIZER for team H10 (Deep Literature & Competitive Benchmarking).\nMETHOD MEMOS:\n${JSON.stringify(memos).slice(0,70000)}\nCOMPETITOR BRIEFS:\n${JSON.stringify(comps).slice(0,30000)}\nCITATION-CHECK (fabricated/not_found = ${fabricated.length}):\n${JSON.stringify(citeChecks).slice(0,20000)}\nWrite a markdown_report: (1) Executive summary — which methods are DEFENSIBLE, which need a SWITCH or caveat, and the top reviewer-risks; (2) Per-method defensibility table (method | verdict | better-method | switch cost | key citation); (3) Competitive-positioning section (per competitor: what we add); (4) Citation-integrity note (any fabricated citations from our own memos — flag for correction); (5) Paper-ready citation list. ALSO emit a findings array (type='literature'): one finding per method whose verdict is switch_recommended or defensible_with_caveats (severity High if switch_recommended, else Medium; file = the relevant pipeline script/doc; recommendation = the methodological action), plus one finding per fabricated citation (severity High, type='rigor'). Fill the tally from these findings.`,
  {schema:REPORT_SCHEMA, label:'synthesize', phase:'Synthesize'})

return { team:'H10', findings:(report&&report.findings)||[], memos, competitors:comps, citeChecks, report }
