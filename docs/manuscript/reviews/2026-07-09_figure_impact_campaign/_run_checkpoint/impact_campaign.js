export const meta = {
  name: 'masld-figure-impact-campaign',
  description: 'Adversarial impact+argument review of the 5 MASLD main figures + connective writing for a NatGenet/NatMetab target; delivers a strategy + concrete rewrites',
  phases: [
    { title: 'Spine', detail: 'forge the single defensible headline + arc + venue bar + competitor lines + buried gems' },
    { title: 'Attack', detail: 'per-unit adversarial lenses (impact/reviewer/logic/design/rigor)' },
    { title: 'Defend', detail: 'blue-team rebuts each attack; filter stale/false-positive' },
    { title: 'Adjudicate', detail: 'reconcile -> ranked surviving findings + drafted concrete fixes' },
    { title: 'Iterate', detail: 'completeness critic -> augment flagged units' },
    { title: 'Synthesize', detail: 'capstone strategy: headline verdict, top-10 moves, rebuttal preview' },
  ],
}

// ------------------------------------------------------------------ FACT PACK
const FACTPACK = [
'TARGET VENUE: Nature Genetics / Nature Metabolism (top specialist). Bar = ONE surprising, defensible, mechanistic claim + population scale. A "comprehensive breadth-is-the-value resource" is desk-reject material at these venues. The paper OWN stated defensible value is "multi-modal breadth itself, not any single headline claim" (paper_outline.md:25) -- treat that as an editorial liability to fix, not an asset.',
'',
'CENTRAL THESIS (current, 2026-05-09): MASLD progression is a MULTI-STEP CELLULAR CASCADE (not a single F2 switch); every CRN fibrosis transition F0->F1..F3->F4 shows a substantial transcriptomic shift + a stage-specific cell-type program; inherited causal signal (colocalization) and the disease-state transcriptome are LARGELY ORTHOGONAL (DEG-COLOC Jaccard ~0.01), so multi-evidence convergence -- not any single assay -- is what nominates credible targets. Load-bearing convergent trio recurring in every figure: THRB, RORA, HKDC1.',
'',
'RETRACTION HISTORY (each is a "why trust the current framing" attack line): single "F2 switch" -> "two-transition" -> "multi-step cascade"; "female-enriched progressor" RETIRED (sex-dimorphic genes 290->8); "convergence beats best single channel" RETIRED (now: matches/ties, only Moylan is a significant win); cross-species/LITMUS pillar CUT 2026-07-05 (now human-only; "7 modalities" must be recounted -- convergence SCORE uses only 6 channels); categorical DEG-COLOC Fisher OR=1.67 RETRACTED; "123 genes 7/7 OR=4.85" RETRACTED; PNPLA3 SuSiEx PIP 0.987 was a MYTH (real 0.003-0.007); all plasma Olink AUROCs WITHDRAWN; "3.76x proteo-cap-spatial above-chance" does NOT survive on the 9,882 universe.',
'',
'CANONICAL CURRENT NUMBERS BY FIGURE:',
'Fig1 (atlas): 1,259 QC-passing / 9 cohorts; 846 samples in the 5-control-bearing pooled DE; 1,918 canonical Tier-1 DEGs (1,419 up/499 down; limma-voom-qw C2; effect-size-aware interval-null FDR<0.05 at lfc=0.25 -- DO NOT use the word "TREAT" in manuscript); atlas 27,187 x 444; LOO-CV mean DEG recovery 81.6% / rho 0.952 / direction 99.9%; 7,842 robust genes.',
'Fig2 (genetics): MAIN SuSiE-COLOC 473 effector genes (PP.H4>0.5) across 35 Tier-1/2 liver GWAS; 369>0.8; 267>0.9; ABF 863; union 1,031. Full-50 supp: 736 SuSiE/1,527 union (496 rest on Tier-3/4). Coding-led only 34/1,031 (3.3%). Ancestry: EUR-only 356/shared 81/non-EUR-only 36; gated non-EUR = 14 GWS + 20 suggestive; 9 genes in >=3 ancestries, none in all 5. Marquee PP.H4: THRB 0.9999(GGT), RORA 0.998, HKDC1 0.992(AST), MTTP 0.835; PNPLA3 0.012, GCKR 0.545/0.426, TM6SF2 0.67(abf only). Broadaway overlap 269/751=36% (their eQTL panel is REUSED; novelty = trait-anchoring + fine-mapping + expanded GWAS, not the eQTLs).',
'Fig3 (progression+scRNA): fibrosis DEG burden 14/305/2,117/3,621 at F1/F2/F3/F4; up:down 8.8->2.1; per-study DE rho=0.25 (fragmented); NMF k=6 (hep fraction 92->78%); single-cell 269 QC donors/211 stage-resolved/~796k cells; Hotspot 54/117 stage-changing (36/54 non-hep); CCC 13 LR pairs; F-stage analyses restricted to 58 Andrews donors. Sex-dimorphic = 8 (1F/7M).',
'Fig4 (physical corroboration): prioritized universe 9,882 (transcriptomic 8,088 UNION genetic 3,038; convergent INTERSECTION 1,244); validated proteomics 714/spatial 447/snATAC 560; >=2 assays 141; all-3 = 3. HONEST above-chance stat (seed-pinned): ONLY proteo-cap-spatial beats chance (obs 70 vs exp 51.4, 1.36x, p=2.6e-3); >=2 overall at chance (141 vs 137, p=0.35); genetic-only subset (1,794) enriched in NOTHING. mRNA-protein rho=0.44/66% concordant (indep DIA-MS PXD051911 n=58); triage 145 confirmed/17 discordant/2,882 not-corroborated. Visium GSE192741: 370/415 SVGs; 125 gain/290 keep/80 lose. gsMap OR 1.55-1.78 (ALT/AST/GGT) on convergent set only.',
'Fig5 (convergence+drugs): 6-channel; max 4/6 modalities (16 genes); 371 at exactly 3; 387 at >=3; NO gene at 5-6, NO Tier-2. Tiers: 677 Tier-1 genetic-validated (579 EUR+98 cross-anc)/204 Tier-3/25,668 Tier-4. Drug-dev recovery: 2 approved (THRB/resmetirom, GLP1R/semaglutide)/208 clinical/5 discontinued/1,504 preclinical; RORA->TB-840 Ph1 + HKDC1->KO = post-hoc validations. Benchmark ROC: Govaere 0.93 vs 0.95 (p=0.62 tie), Moylan fibrosis 0.94 vs 0.77 (p=7e-8 WIN -- the ONLY significant win), biomarkers 0.87 vs 0.85 (p=0.60 tie), genetic 0.68 vs 0.85 (LOSE p=0.004), NIDDK 0.69 vs 0.78 (p=0.07). Per-gene FDR = 0 genes at BH<0.05 (permutation-floor-limited); only RANK-level null holds.',
'',
'LOAD-BEARING VULNERABILITIES (soft spots -- but note the writing ALREADY hedges most; the risk is UNDER-claiming):',
'A. F-stage: within-Andrews donor-jackknife QWK 0.638 only; cross-cohort transfer UNTESTABLE (no external Kleiner labels); "healthy control" label heterogeneous across cohorts. Drives cascade thesis/stage-CCC/HKDC1 panel -> restricted to 58 Andrews donors.',
'B. NMF k=6 unsupported by internal metrics (cophenetic/silhouette max at k=3); factor subspace stable (cosine 0.976) but hard partition fuzzy (ARI ~0.465). Present programs as CONTINUOUS axes, not discrete subtypes.',
'C. Spatial batch confound: Vu 2025 + GSE192741 fully study-confounded on disease axis; Visium supports LOCALIZATION/ZONATION only, never a pooled disease fold-change. CosMx cell-level Wilcoxon inflated a count 639->111 at slide resolution.',
'D. Benchmark framing: integration matches/ties single channels; ONLY Moylan is a significant win; LOSES on genetic + clinical panels. Method-selection: frozen pre-registered winner was limma_trend__C2; adopted limma_voom_qw__C2 is rank-8 via a post-hoc size-match -> CANNOT claim "a pre-registered benchmark selected our method." Random-gene ceiling: 500 random genes give F>=3 AUROC 0.860 vs curated 0.853 (n.s.) -- fibrosis prediction is transcriptome-pervasive.',
'E. Convergence per-gene FDR unresolvable (0 survivors; min p~1e-4); integer ranks are a saturated-plateau artifact (THRB moved #90->#3,831 with identical PP.H4); cite marquee genes by tier+PP.H4, not rank. Mountjoy/Schipper SUPERVISED held-out paradigm invoked for an UNSUPERVISED score.',
'F. Orthogonality thesis hinges on a VERSION-UNSTABLE hepatocyte COLOC fold (~1.0; 1.16x old <-> 0.93x current) the authors will NOT stake the thesis on; tensions with WS3 "0/7 transitions hepatocyte-dominant". This is simultaneously the most-novel conceptual claim AND the softest foundation.',
'G. Pseudoreplication residue (mega-review): GSE244832 (dominant disease cohort) = 18 real donors as 117 SRR runs; donor-level re-derivations pending in places.',
'',
'COMPETITORS: Feng 2026 (bioRxiv; 29 datasets/2,640 samples; per-study meta + 6-layer score -> 39 genes; MLIP in-vitro; portal) = CLOSEST SCOOP RISK; differentiator = pooled cohort-adjusted mega-analysis + COLOC + deconvolution + multi-ancestry + spatial/ATAC + drug repurposing; never name directly. Tzouanas 2026 (Cell; mouse Multiome + human Xenium; MATCHA TF-enhancer triads; 4 hepatocyte stress programs; HMGCS2 LiverKO->HCC IN VIVO -- this paper CANNOT match in-vivo); differentiator = human population scale + causal + convergence + drug repurposing. Li 2025 (Nat Genet; 540K cells but only 61 livers, no bulk/GWAS). Vacca/LITMUS (dropped with cross-species). Broadaway 2024 = SOURCE eQTL panel reused. Govaere 2020/2023 = external benchmark (GSE135251 overlap = circularity risk).',
'',
'CURRENT WRITING (source of truth = docs/manuscript/working/*.md; the working files are the CURRENT claims):',
'Fig2 = topic-sentence style, heavy explicit hedging (ancestry porting, PNPLA3 coding-blindness, PDFF no clean coloc, novelty rests on portfolio not expression). Fig3 = rigor-firewall structure (LEADS with per-study non-reproducibility rho=0.25, then 4 trust axes, THEN biology). Fig4 = bold-subheading-as-claim style (physical corroboration); each subheading names its own limit. Fig5+Discussion = 3 lessons + explicit 9-item Limitations; benchmark honesty foregrounded. Intro (01_intro.md) drafted, cross-species cut. NO abstract.md exists. NO fig1.md exists.',
'',
'GROUNDING RULES (mandatory for every agent): (1) READ the actual current file(s) named for your unit before critiquing -- attack the CURRENT text/panels, not a summary or a retired draft. (2) CITE the real file / sidecar CSV a claim lives in. (3) FLAG any critique that may target a stale/superseded number (numbers churned across >=4 DEG defs and >=5 convergence builds). (4) NUMBERS.md "Key Numbers by Figure" uses STALE legacy 6-figure numbering -- do NOT trust it for figure mapping; the working/*.md files + paper_outline.md sections 1/3 are canonical. (5) Method<->modality rule: never propose a fix that transplants a tool across assays (e.g. a methylation diff-variability method onto RNA-seq counts). (6) The prose is deliberately hedged: distinguish JUSTIFIED honesty from IMPACT-KILLING over-hedging; flag places the paper BURIES a strong true result.',
].join('\n')

const ROOT = 'docs/manuscript/'

const UNITS = [
  { key:'fig1', kind:'figure', title:'Figure 1 - Multi-cohort atlas overview',
    files:'NO fig1.md prose exists. Read the Fig1 spec in docs/paper_outline.md (search "Fig 1"/"Figure 1"), the roster in scripts/figures/fig1_atlas_overview_v2.R header, and note figures/main/fig1_atlas_overview/ survives only as raster figure1_JL.png with most panels migrated to the Fig3 supplement.',
    note:'LEAST settled. Core question: WHAT should Fig 1 be at the NG/NatMetab bar -- a claim-carrying opener or a resource-credibility panel? It currently barely exists.' },
  { key:'fig2', kind:'figure', title:'Figure 2 - Fine-mapping & multi-ancestry colocalization',
    files:'Read docs/manuscript/working/fig2.md (full) and the panel roster in figures/main/fig2_genetics/panels/ (A-H: cascade, PIP x PP.H4, coding split, ancestry counts, cross-ancestry, RORA/FABP1 locus zooms, legend).',
    note:'Final. Very heavily hedged. Novelty admitted to rest on portfolio+fine-mapping, not the reused Broadaway eQTLs.' },
  { key:'fig3', kind:'figure', title:'Figure 3 - Integrated bulk + single-cell progression cascade',
    files:'Read docs/manuscript/working/fig3.md (full) and figures/main/fig3_RNAseq/README.md (panel map A-I).',
    note:'Converged; rigor-firewall structure LEADS with non-reproducibility before biology. WS3 multicellular reframe tensions with hepatocyte-centric leg. 4->10->16 CCC triple is UNREPRODUCIBLE (softened to qualitative).' },
  { key:'fig4', kind:'figure', title:'Figure 4 - Physical corroboration in tissue/protein/chromatin',
    files:'Read docs/manuscript/working/fig4.md (full, includes legends) and figures/main/fig4_validation/README.md (panel map 4a-4f).',
    note:'MOST in-flux. Fig4a has 5 unfinalized design candidates. Honest stat: ONLY proteo-cap-spatial beats chance (1.36x, p=2.6e-3); genetic-only enriches in nothing. Spatial = localization only (batch confound).' },
  { key:'fig5', kind:'figure', title:'Figure 5 - Convergence & drug-target prioritization',
    files:'Read the Results portion of docs/manuscript/working/fig5_discussion.md and the panel roster figures/main/fig5_convergence/panels/ (a=6-channel matrix, b=channel histogram, c=drug calibration plane, d=benchmark ROC).',
    note:'Final. Benchmark: only Moylan is a significant win; per-gene FDR=0; 2 approved drugs recovered. Is the payoff strong enough for NG?' },
  { key:'thesis_arc', kind:'connective', title:'Overall thesis & Fig1->Fig5 narrative arc',
    files:'Read docs/paper_outline.md section 1 (thesis, story arc, retired framings).',
    note:'The stated defensible value ("breadth, not a single claim") is editorially fatal at NG/NatMetab. Core question: what single headline should the arc build to, and does the current arc build to it?' },
  { key:'abstract', kind:'connective', title:'Abstract (does not yet exist)',
    files:'NO abstract.md exists. Read docs/paper_outline.md section 1 + the topic sentences of working/fig2.md..fig5_discussion.md to infer what the abstract MUST assert.',
    note:'Task = DRAFT the highest-impact honest abstract AND critique the framing it implies. Must survive a hostile 30-second editor skim.' },
  { key:'intro', kind:'connective', title:'Introduction',
    files:'Read docs/manuscript/working/01_intro.md (full).',
    note:'Cross-species pillar cut. Check the why-now/why-us GAP statement, the positioning vs recent meta-analyses (Feng, un-named), and whether it sets up the chosen headline.' },
  { key:'discussion', kind:'connective', title:'Discussion + Limitations',
    files:'Read the Discussion + Limitations portion of docs/manuscript/working/fig5_discussion.md.',
    note:'3 lessons + explicit 9-item Limitations. Check over- vs under-claiming and whether the closing payoff is strong enough for NG/NatMetab.' },
]

// figures get all 5 lenses; connective units skip panel-design
const LENSES = [
  { id:'impact_editor', effort:'high', desc:'IMPACT / EDITOR lens. You are a hard-nosed Nature Genetics / Nature Metabolism editor deciding send-out vs desk-reject. Judge: what is the single most exciting DEFENSIBLE point in this unit, and is it foregrounded or buried under hedging? Is the claim novel vs Feng/Tzouanas/Li/Govaere/Broadaway or incremental? What is the "so what" for a hepatology/genetics reader? Where does the unit read as a comprehensive-but-unsurprising resource (the death frame)? Name the buried gems that should be promoted.' },
  { id:'reviewer2', effort:'high', desc:'REVIEWER #2 lens. You are a hostile expert MASLD hepatologist/geneticist trying to trigger major-revision or reject. Find the fatal flaw; name the exact analysis, control, or experiment you would demand; identify every place the unit OVER-claims beyond its data AND every place it UNDER-claims a real result. Be specific and merciless but fair.' },
  { id:'argument_logic', effort:'high', desc:'ARGUMENT / LOGIC lens. For each topic sentence / bold subheading, decompose claim -> evidence -> warrant. Does the cited panel/number actually support the sentence? Any inferential leaps, non-sequiturs, or internal contradictions? Check CROSS-FIGURE number consistency (e.g. coloc counts 1,031 vs 1,162 vs 862; DEG-count churn; modality count 6 vs 7). Flag hedges that are so heavy they negate the claim.' },
  { id:'panel_design', effort:'medium', figuresOnly:true, desc:'PANEL-LOGIC / DESIGN lens. Does each panel earn its place, or is it redundant/decorative? Is there a MISSING panel the argument needs? Is the final/capstone panel an actual payoff? Does the figure tell ONE story at a glance? Deliver an explicit verdict where relevant: for Fig4a pick among the 5 candidates (forest/corefield/physlayers/ghostflow/triangle) tied to the honest above-chance stat; for Fig1 state concretely what the figure should contain.' },
  { id:'rigor_spot', effort:'medium', desc:'RIGOR SPOT-CHECK lens (light -- two big audits already ran; find only NEW issues on the CURRENT lineup). Look for: circularity/leakage, pseudoreplication, wrong baseline/null, a number computed on a DIFFERENT universe than the one now shown (e.g. spatial 447 vs the audited 111), batch confound leaking into a disease claim, or a reproducibility bug. Do NOT re-report already-hedged known limitations unless the hedge is now insufficient.' },
]

// -------------------------------------------------------------------- SCHEMAS
const SPINE_SCHEMA = {
  type:'object', additionalProperties:false,
  properties:{
    chosen_headline:{ type:'string', description:'the single sharpest defensible+novel claim this paper can make at NG/NatMetab' },
    headline_rationale:{ type:'string' },
    alternatives:{ type:'array', items:{ type:'string' } },
    arc:{ type:'string', description:'the Fig1->Fig2->Fig3->Fig4->Fig5 arc that builds to the headline' },
    venue_bar:{ type:'string', description:'what NG vs NatMetab each require here + a realistic-target verdict if the headline cannot carry NG' },
    competitor_lines:{ type:'array', items:{ type:'string' }, description:'the exact differentiation line vs Feng, Tzouanas, Li, Broadaway, Govaere' },
    buried_gems:{ type:'array', items:{ type:'string' }, description:'strong true results currently hidden under hedging that should be promoted' },
    rubric:{ type:'string', description:'the yardstick downstream reviewers must apply to every unit' },
  },
  required:['chosen_headline','arc','venue_bar','competitor_lines','buried_gems','rubric'],
}

const ATTACK_SCHEMA = {
  type:'object', additionalProperties:false,
  properties:{
    lens:{ type:'string' }, unit:{ type:'string' },
    findings:{ type:'array', items:{
      type:'object', additionalProperties:false,
      properties:{
        target:{ type:'string', description:'the exact current claim/panel/sentence attacked (quote or precise ref)' },
        axis:{ type:'string', enum:['impact','novelty','logic','rigor','design','positioning','underclaim'] },
        severity:{ type:'string', enum:['Critical','High','Medium','Low'] },
        attack:{ type:'string' },
        why_it_matters:{ type:'string', description:'consequence at the NG/NatMetab bar' },
        cited_source:{ type:'string', description:'real file / sidecar CSV the claim lives in' },
        stale_risk:{ type:'boolean' },
        proposed_fix:{ type:'string', description:'concrete direction: rewrite gist or panel change' },
      },
      required:['target','axis','severity','attack','why_it_matters','proposed_fix'],
    } },
  },
  required:['lens','unit','findings'],
}

const DEFENSE_SCHEMA = {
  type:'object', additionalProperties:false,
  properties:{
    unit:{ type:'string' },
    verdicts:{ type:'array', items:{
      type:'object', additionalProperties:false,
      properties:{
        target:{ type:'string' },
        survives:{ type:'boolean' },
        rationale:{ type:'string', description:'already-addressed / factually-wrong / stale-number / genuine' },
        corrected_severity:{ type:'string', enum:['Critical','High','Medium','Low','Dismissed'] },
      },
      required:['target','survives','rationale','corrected_severity'],
    } },
  },
  required:['unit','verdicts'],
}

const ADJ_SCHEMA = {
  type:'object', additionalProperties:false,
  properties:{
    unit:{ type:'string' },
    headline_verdict:{ type:'string', description:'2-4 sentences: does this unit advance the chosen headline AND clear the bar' },
    surviving_findings:{ type:'array', items:{
      type:'object', additionalProperties:false,
      properties:{
        rank:{ type:'integer' },
        target:{ type:'string' },
        axis:{ type:'string' },
        severity:{ type:'string', enum:['Critical','High','Medium','Low'] },
        leverage:{ type:'string', enum:['high','medium','low'] },
        problem:{ type:'string' },
        concrete_fix:{ type:'string', description:'FOR WRITING: the DRAFTED ready-to-paste replacement sentence/subheading/legend. FOR PANELS: exact panel id + precise change. NO placeholders.' },
        effort:{ type:'string', enum:['low','medium','high'] },
      },
      required:['rank','target','axis','severity','leverage','problem','concrete_fix','effort'],
    } },
    verdict_md:{ type:'string', description:'full markdown for this unit verdict file: headline verdict, then ranked findings each with its drafted concrete fix' },
  },
  required:['unit','headline_verdict','surviving_findings','verdict_md'],
}

const CRITIC_SCHEMA = {
  type:'object', additionalProperties:false,
  properties:{
    gaps:{ type:'array', items:{
      type:'object', additionalProperties:false,
      properties:{
        unit:{ type:'string' },
        gap:{ type:'string', description:'a claim not attacked / competitor not compared / cross-figure inconsistency / place the paper UNDER-claims' },
        priority:{ type:'string', enum:['high','medium','low'] },
      },
      required:['unit','gap','priority'],
    } },
  },
  required:['gaps'],
}

const AUG_SCHEMA = {
  type:'object', additionalProperties:false,
  properties:{
    unit:{ type:'string' },
    extra_findings:{ type:'array', items:ADJ_SCHEMA.properties.surviving_findings.items },
    addendum_md:{ type:'string' },
  },
  required:['unit','extra_findings','addendum_md'],
}

// ------------------------------------------------------------------- HELPERS
async function tryAgent(prompt, opts, attempts) {
  const n = attempts || 2
  for (let i = 0; i < n; i++) {
    const r = await agent(prompt, opts)
    if (r !== null && r !== undefined) return r
    log('retry ' + (opts.label || '') + ' attempt ' + (i + 1))
  }
  return null
}

function compactFindings(atkResults) {
  const flat = []
  for (const a of atkResults) {
    if (!a || !a.findings) continue
    for (const f of a.findings) flat.push('[' + a.lens + '][' + f.axis + '/' + f.severity + '] TARGET: ' + f.target + ' | ATTACK: ' + f.attack + ' | WHY: ' + f.why_it_matters + ' | FIX: ' + f.proposed_fix + (f.stale_risk ? ' | STALE_RISK' : '') + (f.cited_source ? ' | SRC: ' + f.cited_source : ''))
  }
  return flat.join('\n')
}

// STAGE gate: 'attacks' = run spine + ALL attacks then STOP (no defend); 'full' = whole campaign.
const STAGE = (args && args.stage) || 'full'

// ------------------------------------------------------------------- PHASE A
phase('Spine')
const spineHint = 'You are forging the impact spine for a MASLD human multi-omic atlas paper aimed at Nature Genetics / Nature Metabolism.\n\nFACT PACK:\n' + FACTPACK + '\n\n'

const spineAgents = await parallel([
  () => tryAgent(spineHint + 'TASK: Propose exactly 3 candidate HEADLINE CLAIMS the paper can honestly make that SURVIVE every item in the retraction history. For each: state it in one sentence, and score novelty / defensibility / evidence-strength / NG-fit (1-5 each) with a one-line justification. Then pick the single strongest and describe the Fig1->Fig5 arc that would build to it. Be concrete and MASLD-specific; do not default to "breadth". Return a short structured brief as text.', { label:'spine:headlines', phase:'Spine', effort:'high' }),
  () => tryAgent(spineHint + 'TASK: Calibrate the exact bar. What specifically does Nature Genetics require of this paper vs Nature Metabolism? Then write the precise one-line DIFFERENTIATION line the paper must land vs each of: Feng 2026, Tzouanas 2026 (Cell, has in-vivo this paper lacks), Li 2025, Broadaway 2024 (reused eQTL panel), Govaere. End with a blunt verdict: can the strongest honest headline carry NG, or is NatMetab/Nat Commun the realistic ceiling? Return as text.', { label:'spine:bar+competitors', phase:'Spine', effort:'high' }),
  () => tryAgent(spineHint + 'TASK: UNDER-CLAIM SCOUT. This paper has hedged itself into potential mediocrity. Identify the strongest TRUE results currently buried under caveats or demoted to supp that a top editor would find exciting -- e.g. the genetics(non-parenchymal)-vs-expression(hepatocyte) orthogonality, the drug-development-gradient recovery WITHOUT labels (2 approved + RORA/HKDC1 post-hoc), the regulatory-not-coding architecture (96% non-coding credible-set leads), the effector-gene map extending Broadaway. For each buried gem: where it lives now, why it is under-sold, and how to promote it honestly. Return as text.', { label:'spine:buried-gems', phase:'Spine', effort:'high' }),
])

const spine = await tryAgent(
  spineHint +
  'Three specialists produced these briefs.\n\n--- HEADLINES ---\n' + (spineAgents[0] || 'n/a') +
  '\n\n--- BAR + COMPETITORS ---\n' + (spineAgents[1] || 'n/a') +
  '\n\n--- BURIED GEMS ---\n' + (spineAgents[2] || 'n/a') +
  '\n\nTASK: Synthesize into the definitive IMPACT SPINE. Choose the single headline, the arc, the venue verdict, the competitor lines, the buried gems to promote, and a crisp RUBRIC that every downstream figure/writing reviewer must apply (the yardstick: does this unit advance the chosen headline AND clear the bar, without over- or under-claiming).',
  { label:'spine:synthesize', phase:'Spine', schema:SPINE_SCHEMA, effort:'high' })

log('SPINE headline: ' + (spine ? spine.chosen_headline : 'FAILED'))

const spineStr = spine ? [
  'CHOSEN HEADLINE: ' + spine.chosen_headline,
  'RATIONALE: ' + (spine.headline_rationale || ''),
  'ARC: ' + spine.arc,
  'VENUE BAR: ' + spine.venue_bar,
  'COMPETITOR LINES:\n- ' + (spine.competitor_lines || []).join('\n- '),
  'BURIED GEMS TO PROMOTE:\n- ' + (spine.buried_gems || []).join('\n- '),
  'RUBRIC: ' + spine.rubric,
].join('\n') : 'SPINE UNAVAILABLE -- default to: headline = genetics/expression orthogonality makes multi-evidence convergence a necessity, validated by label-free drug-gradient recovery; do not lead with breadth.'

// ------------------------------------------------------------- PHASES B + iterate
const baseCtx = (u) =>
  'You are adversarially reviewing ONE unit of a MASLD atlas paper targeting Nature Genetics / Nature Metabolism.\n\n' +
  'FACT PACK:\n' + FACTPACK + '\n\nIMPACT SPINE (the yardstick):\n' + spineStr + '\n\n' +
  'UNIT: ' + u.title + '\nUNIT NOTE: ' + u.note + '\nFILES TO READ FIRST (mandatory -- attack the CURRENT text/panels): ' + u.files + '\n\n'

async function attackUnit(u) {
  const lenses = LENSES.filter(l => !(l.figuresOnly && u.kind !== 'figure'))
  const results = await parallel(lenses.map(l => () =>
    tryAgent(
      baseCtx(u) +
      'YOUR LENS: ' + l.desc + '\n\nInstructions: Read the named file(s) now. Then produce your findings. Every finding must name the exact current target, cite the real source file, set stale_risk=true if it might target a superseded number, and give a concrete proposed fix. Prioritize IMPACT and ARGUMENT over minor rigor (rigor was twice-audited). Return 3-8 of your sharpest findings; quality over quantity.',
      { label:'atk:' + u.key + ':' + l.id, phase:'Attack', schema:ATTACK_SCHEMA, effort:l.effort }
    ).then(r => r || { lens:l.id, unit:u.key, findings:[] })
  ))
  return { unit:u, attacks:results }
}

async function defendUnit(prev, u) {
  const atkStr = compactFindings(prev.attacks)
  const def = await tryAgent(
    baseCtx(u) +
    'You are the BLUE TEAM defender. Below are attacks from up to five adversarial lenses. For EACH distinct attack, judge whether it SURVIVES: dismiss it if (a) the current text already addresses/hedges it adequately, (b) it is factually wrong (prior audits produced wrong findings -- verify against the fact pack and the actual file), or (c) it targets a stale/superseded/retired number rather than the current claim. Keep it if it is a genuine, current, impact- or argument-relevant weakness. Read the unit file(s) as needed to verify.\n\nATTACKS:\n' + atkStr,
    { label:'def:' + u.key, phase:'Defend', schema:DEFENSE_SCHEMA, effort:'medium' }
  )
  return { unit:u, attacks:prev.attacks, defense:def || { unit:u.key, verdicts:[] } }
}

async function adjudicateUnit(prev, u) {
  const atkStr = compactFindings(prev.attacks)
  const defStr = (prev.defense.verdicts || []).map(v => (v.survives ? 'KEEP' : 'DROP') + ' [' + v.corrected_severity + '] ' + v.target + ' :: ' + v.rationale).join('\n')
  const adj = await tryAgent(
    baseCtx(u) +
    'You are the ADJUDICATOR. Reconcile the attacks with the blue-team defense; keep only genuine, CURRENT, high-leverage weaknesses. Rank them by impact-leverage x severity x how-silently-the-error-travels (impact/argument weighted above minor rigor). For EACH surviving finding, write a CONCRETE FIX: for writing, the DRAFTED, ready-to-paste replacement sentence / bold subheading / legend (match the paper voice: active, finding-first, no colored fonts, no "TREAT" -- say "adjusted P<0.05"); for panels, the exact panel id + precise change spec (e.g. which of the 5 Fig4a candidates + why; what Fig1 should contain). NO placeholders -- every fix must be usable as-is. Also render verdict_md: a clean markdown section for this unit (## verdict, then a ranked list, each finding with problem + drafted fix + effort).\n\nATTACKS:\n' + atkStr + '\n\nDEFENSE:\n' + defStr,
    { label:'adj:' + u.key, phase:'Adjudicate', schema:ADJ_SCHEMA, effort:'high' }
  )
  if (adj) { adj.unit = u.key }  // force canonical key for downstream matching
  return adj || { unit:u.key, headline_verdict:'(adjudication failed)', surviving_findings:[], verdict_md:'## ' + u.key + '\n(adjudication failed)' }
}

phase('Attack')
const attacked = (await parallel(UNITS.map(u => () => attackUnit(u)))).filter(Boolean)
log('all units attacked: ' + attacked.length + '/' + UNITS.length)

if (STAGE === 'attacks') {
  // HARD STOP after attacks -- no defend/adjudicate/iterate/synthesize agent() calls fire.
  return {
    stage: 'attacks_only',
    spine, spineStr,
    attacks: attacked.map(a => ({ unit: a.unit.key, unit_title: a.unit.title, results: a.attacks })),
  }
}

let unitResults = (await pipeline(
  attacked,
  (a) => defendUnit(a, a.unit),
  (d, a0) => adjudicateUnit(d, a0.unit)
)).filter(Boolean)

// -------------------------------------------------------------- PHASE C iterate
phase('Iterate')
const adjDigest = unitResults.map(a => '### ' + a.unit + '\nVERDICT: ' + a.headline_verdict + '\nFINDINGS: ' + (a.surviving_findings || []).map(f => f.target).join(' ; ')).join('\n\n')

const critic = await tryAgent(
  'You are the COMPLETENESS CRITIC for an adversarial review of a MASLD atlas paper (NG/NatMetab target).\n\nIMPACT SPINE:\n' + spineStr + '\n\nHere is what the per-unit adjudications currently cover:\n' + adjDigest + '\n\nTASK: Find what EVERY panel missed. Specifically hunt for: (a) a load-bearing claim that was never attacked; (b) a competitor comparison never made; (c) a CROSS-FIGURE inconsistency (numbers, thesis legs that contradict, e.g. hepatocyte-centric vs WS3 multicellular); (d) a place the paper UNDER-claims a strong true result; (e) an arc-level gap (does Fig1->Fig5 actually build to the headline?). Return only genuinely NEW high/medium gaps.',
  { label:'critic:round1', phase:'Iterate', schema:CRITIC_SCHEMA, effort:'high' }
)

const gaps = (critic && critic.gaps) ? critic.gaps.filter(g => g.priority !== 'low') : []
const byUnit = {}
for (const g of gaps) { (byUnit[g.unit] = byUnit[g.unit] || []).push(g.gap) }
const flaggedUnits = Object.keys(byUnit)
log('critic flagged ' + flaggedUnits.length + ' units for augmentation')

const augmentations = await parallel(flaggedUnits.map(uk => () => {
  const u = UNITS.find(x => x.key === uk) || { key:uk, title:uk, note:'', files:'(cross-cutting)' }
  return tryAgent(
    baseCtx(u) +
    'A completeness critic flagged these GAPS for this unit that the first review round missed:\n- ' + byUnit[uk].join('\n- ') +
    '\n\nTASK: Address each gap. Read the file(s) as needed. Return extra surviving findings (same structure: rank continues from where the unit left off; each with a DRAFTED concrete_fix, no placeholders) and an addendum_md markdown block to append to this unit verdict.',
    { label:'aug:' + uk, phase:'Iterate', schema:AUG_SCHEMA, effort:'high' }
  )
}))

// merge augmentations (match by index -> flaggedUnits, not the agent-filled unit field)
augmentations.forEach((aug, i) => {
  if (!aug) return
  const uk = flaggedUnits[i]
  const tgt = unitResults.find(a => a.unit === uk)
  if (tgt) {
    tgt.surviving_findings = (tgt.surviving_findings || []).concat(aug.extra_findings || [])
    tgt.verdict_md = (tgt.verdict_md || '') + '\n\n### Round-2 additions\n' + (aug.addendum_md || '')
  }
})

// ------------------------------------------------------------- PHASE D capstone
phase('Synthesize')
const fullDigest = unitResults.map(a =>
  '### ' + a.unit + ' -- ' + a.headline_verdict + '\n' +
  (a.surviving_findings || []).slice(0, 12).map(f => '- [' + f.leverage + '/' + f.severity + '/' + f.axis + '] ' + f.target + ' -> ' + f.concrete_fix).join('\n')
).join('\n\n')

const strategy = await tryAgent(
  'You are the CAPSTONE STRATEGIST. Synthesize an adversarial impact review of a MASLD atlas paper (target: Nature Genetics / Nature Metabolism) into the FINAL STRATEGY that will strengthen the main figures and writing.\n\n' +
  'FACT PACK:\n' + FACTPACK + '\n\nIMPACT SPINE:\n' + spineStr + '\n\nCROSS-FIGURE COMPLETENESS GAPS:\n' + (gaps.map(g => '- (' + g.unit + '/' + g.priority + ') ' + g.gap).join('\n') || 'none') + '\n\nPER-UNIT ADJUDICATIONS:\n' + fullDigest + '\n\n' +
  'Produce a single polished markdown strategy document with these sections:\n' +
  '1. HEADLINE VERDICT -- the single sharpest defensible claim, the venue reality (NG vs NatMetab), and the explicit Fig1->Fig5 arc that builds to it; call out where the paper currently under-claims or mis-emphasizes.\n' +
  '2. TOP 10 HIGHEST-LEVERAGE MOVES -- ranked, each = the move + which figure/section + expected impact + effort. These must be drawn from the per-unit findings.\n' +
  '3. CROSS-CUTTING FIXES -- number-consistency reconciliation (coloc 1,031/1,162/862; modality 6 vs 7; DEG churn), the hedging-vs-impact rebalance (which hedges to keep vs which over-hedges to cut), and the competitor-positioning lines (Feng scoop-defense; standing next to Tzouanas in-vivo).\n' +
  '4. REVIEWER REBUTTAL PREVIEW -- the 5 hardest questions an NG/NatMetab reviewer will ask, and for each: does the paper currently answer it, and the one-line strengthening.\n' +
  'Be decisive and concrete. This is the deliverable the author will act on.',
  { label:'capstone:strategy', phase:'Synthesize', effort:'high' }
)

return {
  spine,
  spineStr,
  units: unitResults.map(a => ({ key:a.unit, headline_verdict:a.headline_verdict, surviving_findings:a.surviving_findings || [], verdict_md:a.verdict_md || '' })),
  gaps,
  strategy_md: strategy || '(capstone failed -- reconstruct from per-unit verdicts)',
}
