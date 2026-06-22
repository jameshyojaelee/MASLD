// Capstone — Synthesis & Meta-Review.
// args: { themes:[{key,label,keywords,fix_hint}], reverify:[{id,title,file,severity,claim,evidence}] }
// Phase 1 Themes: one agent per systemic theme -> consolidated brief.
// Phase 2 Reverify: 3-lens cross-team re-verification of the highest-impact still-uncertain Criticals.
// Phase 3 Deliverables: exec summary + remediation roadmap + reviewer-defense dossier + coverage matrix.
export const meta = {
  name: 'megareview-capstone',
  description: 'Capstone: theme consolidation -> cross-team reverification -> executive summary + remediation roadmap + reviewer-defense dossier + coverage matrix',
  phases: [
    { title: 'Themes', detail: 'consolidate systemic themes' },
    { title: 'Reverify', detail: 'cross-team re-verification of top Criticals' },
    { title: 'Deliverables', detail: 'exec summary, roadmap, dossier, coverage' },
  ],
}

const A = typeof args === 'string' ? JSON.parse(args) : (args || {})
const ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'
const DIR = ROOT + '/docs/audit/2026-06-13_megareview'
const CSV = DIR + '/00_master_findings.csv'
const PRE = `You are a synthesis lead on the MASLD mega-review CAPSTONE. The full findings database is ${CSV} (1,900+ findings across teams V1-V7, H8-H10, columns: id,team,severity,type,adversarial_verdict,title,file,line,claim,evidence,recommendation). Per-team reports are ${DIR}/team_*.md. READ-ONLY. Be precise, cite file:line and finding ids.`

const _agent = agent
const tryAgent = async (p, o) => { let r = null; for (let i = 0; i < 4 && r == null; i++) { r = await _agent(p, o) } return r }
const STR = { type:'string' }; const NUM = { type:'number' }
const ENUM = (e) => ({ type:'string', enum:e }); const ARR = (items) => ({ type:'array', items })
const S = (props, req) => ({ type:'object', additionalProperties:false, required:req, properties:props })

const THEME_SCHEMA = S({ theme:STR, canonical_issue:STR, severity:ENUM(['Critical','High','Medium','Low']), n_findings:NUM, affected_files:ARR(STR), root_cause:STR, blast_radius:STR, fix:STR, manuscript_impact:STR }, ['theme','canonical_issue','severity','root_cause','fix','manuscript_impact'])
const VERDICT_SCHEMA = S({ verdict:ENUM(['confirmed','refuted','uncertain']), corrected_severity:ENUM(['Critical','High','Medium','Low','drop']), reasoning:STR }, ['verdict','reasoning'])
const MD_SCHEMA = S({ markdown:STR }, ['markdown'])
const COVERAGE_SCHEMA = S({ markdown:STR, gaps:ARR(STR) }, ['markdown','gaps'])

const THEMES = A.themes || []
const REVERIFY = A.reverify || []
log(`Capstone: ${THEMES.length} themes, ${REVERIFY.length} reverify targets`)

// ---- Phase 1: Theme consolidation ----
phase('Themes')
const themeBriefs = (await parallel(THEMES.map(t => () =>
  tryAgent(`${PRE}\nCONSOLIDATE the systemic theme: "${t.label}".\nKeywords/scope: ${t.keywords}\nFix direction: ${t.fix_hint||''}\nGrep ${CSV} and read the relevant team reports for findings on this theme. Produce ONE canonical-issue brief: the single root cause, the full list of affected files/scripts/columns (deduplicated), the blast radius (which manuscript claims/figures/numbers depend on it), the concrete fix, and the manuscript impact if unfixed. Count how many of the 1,900 findings collapse into this one theme.`,
    {schema:THEME_SCHEMA, label:`theme:${t.key}`, phase:'Themes'})
))).filter(Boolean)

// ---- Phase 2: Cross-team re-verification ----
phase('Reverify')
const LENSES = ['correctness','statistics','reproduce-it']
const reverified = (await parallel(REVERIFY.map(f => () =>
  parallel(LENSES.map(lens => () =>
    tryAgent(`${PRE}\nCROSS-TEAM ADVERSARIAL RE-VERIFICATION (lens: ${lens}). This Critical/High was flagged but not fully verified at team level. Default to REFUTING; confirm only if it holds against the actual files. Re-read the cited code/data.\nFINDING ${f.id}: ${f.title}\nFile: ${f.file}\nClaim: ${f.claim}\nEvidence: ${(f.evidence||'').slice(0,1500)}\nFrom the ${lens} angle: real? Verdict + reasoning citing file:line.`,
      {schema:VERDICT_SCHEMA, label:`reverify:${f.id}:${lens}`, phase:'Reverify'})
  )).then(vs => {
    const v = vs.filter(Boolean)
    const ref = v.filter(x=>x.verdict==='refuted').length
    const con = v.filter(x=>x.verdict==='confirmed').length
    return { id:f.id, title:f.title, file:f.file, severity:f.severity, verdict: ref>=2?'refuted':(con>=2?'confirmed':'uncertain'), details:v.map(x=>x.reasoning) }
  })
))).filter(Boolean)
const reverifyConfirmed = reverified.filter(r=>r.verdict==='confirmed')

// ---- Phase 3: Deliverables ----
phase('Deliverables')
const themesJson = JSON.stringify(themeBriefs).slice(0, 60000)
const reverJson = JSON.stringify(reverified).slice(0, 20000)

const [execR, roadmapR, dossierR, coverageR] = await parallel([
  () => tryAgent(`${PRE}\nWrite the EXECUTIVE SUMMARY for the entire mega-review (10 teams, ~1,900 findings, ${THEMES.length} systemic themes).\nTHEME BRIEFS:\n${themesJson}\nRE-VERIFICATION (confirmed=${reverifyConfirmed.length}/${reverified.length}):\n${reverJson}\nProduce markdown: (1) one-paragraph bottom-line verdict on submission readiness; (2) the systemic-theme table (theme | severity | #findings | manuscript impact); (3) the count of confirmed Criticals and the 8-12 most submission-blocking issues with finding ids; (4) what is SOUND and survives (don't only list problems). Be direct and honest.`, {schema:MD_SCHEMA, label:'deliver:exec', phase:'Deliverables'}),
  () => tryAgent(`${PRE}\nWrite the REMEDIATION ROADMAP.\nTHEME BRIEFS:\n${themesJson}\nProduce a prioritized markdown roadmap ranked by (severity x blast-radius) / effort. Group into: (A) BLOCKERS — must fix before any resubmission (e.g. C2-rename crashes, sign-inversion, leakage, pseudoreplication, invalid convergence FDR); (B) NUMBER/DOC reconciliation (drifted numbers, retired-framing purge); (C) figure/portal regeneration; (D) defensible-but-document items. For each item: the canonical file(s), the concrete fix, the rebuild chain it triggers (e.g. 27a->75->217), and a rough effort tag (S/M/L). Reference finding ids.`, {schema:MD_SCHEMA, label:'deliver:roadmap', phase:'Deliverables'}),
  () => tryAgent(`${PRE}\nWrite the REVIEWER-DEFENSE DOSSIER. Read ${DIR}/team_H9_rigor_redteam.md and the theme briefs.\nTHEME BRIEFS:\n${themesJson}\nFor each of the strongest hostile-reviewer attacks (thesis F1->F2-smallest inconsistency, pseudoreplication across CCC/chromVAR/DIALOGUE, prediction+F-stage leakage, CCC sign-inversion, cell-type-heritability direction error, convergence-FDR invalidity, COLOC PIP myths, stale/irreproducible numbers, no wet-lab): give the attack phrasing, the underlying finding id(s), the current defense, the gap, and the recommended pre-submission hardening. Markdown, ranked by danger.`, {schema:MD_SCHEMA, label:'deliver:dossier', phase:'Deliverables'}),
  () => tryAgent(`${PRE}\nWrite the COVERAGE MATRIX confirming the review's completeness. Verify every major pipeline (bulk RNA-seq, genetics/causal, atlas/convergence, single-cell core, sc-programs/CCC, spatial, epigenomic, drug, proteomics, ncRNA, NMF, network, cross-species, prediction, figures, portal), every headline-number family, every method, and every one of the 20 known vulnerabilities was covered by >=1 team. List any GAP not covered by any team (candidates for a follow-up run). Markdown table + explicit gaps list.`, {schema:COVERAGE_SCHEMA, label:'deliver:coverage', phase:'Deliverables'}),
]).then(r=>r)

return {
  team:'CAPSTONE',
  themeBriefs,
  reverified,
  exec: execR && execR.markdown,
  roadmap: roadmapR && roadmapR.markdown,
  dossier: dossierR && dossierR.markdown,
  coverage: coverageR && { markdown: coverageR.markdown, gaps: coverageR.gaps },
}
