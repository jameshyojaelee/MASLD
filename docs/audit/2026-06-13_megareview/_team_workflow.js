// Generic mega-review team workflow (vertical + H8/H9 horizontal teams).
// Parameterized by `args`: { team, title, cartograph:[{label,prompt}], clusters:[{k,label,paths,focus}],
//   seams:[{k,label,focus}], repro:[string], verifyCap:int }
// Skeleton: Cartograph -> Find(cluster+seam) -> Critic(+gap finders) -> Reproduce -> Verify(3-lens) -> Synthesize.
export const meta = {
  name: 'megareview-team',
  description: 'Generic mega-review team: cartograph -> find -> critic -> reproduce -> adversarial-verify -> synthesize (parameterized by args)',
  phases: [
    { title: 'Cartograph', detail: 'map the territory' },
    { title: 'Find', detail: 'cluster + seam finders' },
    { title: 'Critic', detail: 'completeness gaps + gap-driven finders' },
    { title: 'Reproduce', detail: 'recompute flagship numbers from disk' },
    { title: 'Verify', detail: '3-lens adversarial refutation of High/Critical' },
    { title: 'Synthesize', detail: 'team report + tally' },
  ],
}

const A = typeof args === 'string' ? JSON.parse(args) : (args || {})
const ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'
const CTX = ROOT + '/docs/audit/2026-06-13_megareview/REVIEW_CONTEXT.md'
const PRE = `You are an auditor on the MASLD mega-review, TEAM ${A.team} (${A.title}). FIRST Read ${CTX} in full and obey its rubric, severity definitions, evidence standard, and the READ-ONLY / login-node-safe mandate. Work under project root ${ROOT}. Paths given are approximate — confirm with Grep/Glob/Read.`
const CAP = A.verifyCap || 18

const STR = { type:'string' }
const NUM = { type:'number' }
const ENUM = (e) => ({ type:'string', enum:e })
const ARR = (items) => ({ type:'array', items })
const S = (props, req) => ({ type:'object', additionalProperties:false, required:req, properties:props })
const FINDING_ITEM = S({ title:STR, severity:ENUM(['Critical','High','Medium','Low']), type:ENUM(['bug','rigor','repro','literature']), file:STR, line:STR, claim:STR, evidence:STR, recommendation:STR, lit_question:STR }, ['title','severity','type','file','claim','evidence','recommendation'])
const FINDINGS_SCHEMA = S({ findings: ARR(FINDING_ITEM) }, ['findings'])
const REPRO_ITEM = S({ number_name:STR, claimed:STR, reproduced:STR, source_file:STR, matches:ENUM(['yes','no','partial','could_not']), notes:STR }, ['number_name','claimed','reproduced','source_file','matches'])
const REPRO_SCHEMA = S({ checks: ARR(REPRO_ITEM) }, ['checks'])
const VERDICT_SCHEMA = S({ verdict:ENUM(['confirmed','refuted','uncertain']), corrected_severity:ENUM(['Critical','High','Medium','Low','drop']), reasoning:STR }, ['verdict','reasoning'])
const GAP_ITEM = S({ area:STR, why:STR, probe:STR }, ['area','probe'])
const GAPS_SCHEMA = S({ gaps: ARR(GAP_ITEM) }, ['gaps'])
const MAP_SCHEMA = S({ summary:STR, risks:STR }, ['summary'])
const TALLY = S({ critical:NUM, high:NUM, medium:NUM, low:NUM }, ['critical','high','medium','low'])
const REPORT_SCHEMA = S({ executive_summary:STR, tally:TALLY, markdown_report:STR }, ['executive_summary','tally','markdown_report'])

// transient-failure retry: a dead agent (terminal API error -> null) is retried once before becoming a gap
const _agent = agent
const tryAgent = async (p, o) => { let r = await _agent(p, o); if (r == null) r = await _agent(p, o); return r }

log(`Team ${A.team}: ${A.title} — ${A.clusters.length} clusters, ${A.seams.length} seams, ${A.repro.length} repro tasks`)

phase('Cartograph')
const maps = (await parallel(A.cartograph.map(c => () =>
  tryAgent(`${PRE}\n${c.prompt}`, {schema:MAP_SCHEMA, label:c.label, phase:'Cartograph'})
))).filter(Boolean)
const mapSummary = maps.map(m=>m.summary).join('\n\n')

const round1 = (await parallel(A.clusters.map(c => () =>
  tryAgent(`${PRE}\nYour cluster: "${c.label}".\nLikely files: ${c.paths}\nFocus: ${c.focus}\nDo a DEEP review across code-correctness/logic, scientific rigor, and reproducibility. Read the actual scripts; trace data flow, thresholds, formulas, joins, indexing, NA/allele handling, seeds, and consistency with the canonical numbers in REVIEW_CONTEXT.md. For any method whose appropriateness needs external-literature checking, add a type='literature' finding with a precise lit_question. Cite file:line + exact evidence. Be skeptical and specific.\nMap context:\n${mapSummary}`,
    {schema:FINDINGS_SCHEMA, label:`find:${c.k}`, phase:'Find'})
))).filter(Boolean).flatMap(r=>r.findings)

const round2 = (await parallel(A.seams.map(s => () =>
  tryAgent(`${PRE}\nYou audit a CROSS-CUTTING SEAM: "${s.label}".\nFocus: ${s.focus}\nRead across all relevant scripts AND on-disk outputs. Cite file:line and CSV+column. Return findings per schema.\nMap context:\n${mapSummary}`,
    {schema:FINDINGS_SCHEMA, label:`seam:${s.k}`, phase:'Find'})
))).filter(Boolean).flatMap(r=>r.findings)

let allFinds = [...round1, ...round2]

phase('Critic')
const titles = allFinds.map(f=>`[${f.severity}/${f.type}] ${f.title} (${f.file})`).join('\n')
const gapsRes = await tryAgent(`${PRE}\nYou are the COMPLETENESS CRITIC for team ${A.team}.\nMAP:\n${mapSummary}\nFINDINGS SO FAR (${allFinds.length}):\n${titles}\nIdentify concrete GAPS: scripts/claims/methods NOT examined, or dimensions (correctness/rigor/repro/literature) under-covered. Specific enough for a follow-up agent to act. Up to 8 gaps.`, {schema:GAPS_SCHEMA, label:'critic', phase:'Critic'})
const gaps = (gapsRes && gapsRes.gaps ? gapsRes.gaps : []).slice(0,8)
if (gaps.length) {
  const round3 = (await parallel(gaps.map((g,i) => () =>
    tryAgent(`${PRE}\nGap-driven deep probe: "${g.area}".\nWhy: ${g.why||''}\nProbe: ${g.probe}\nRead the relevant scripts/outputs; return findings per schema with file:line evidence.`,
      {schema:FINDINGS_SCHEMA, label:`gap:${i}`, phase:'Find'})
  ))).filter(Boolean).flatMap(r=>r.findings)
  allFinds = [...allFinds, ...round3]
}

const seen = new Set(); const finds = []
for (const f of allFinds) { const key=(f.file||'')+'|'+(f.title||''); if(!seen.has(key)){ seen.add(key); finds.push(f) } }

phase('Reproduce')
const repro = (await parallel(A.repro.map((t,i) => () =>
  tryAgent(`${PRE}\nNUMBER-REPRODUCTION task (lightweight, login-node-safe; read CSV/RDS, count/filter/correlate only — NO heavy re-runs).\n${t}\nReport claimed vs reproduced with the exact source file. If a file is missing or not cheaply reproducible, mark could_not and explain.`,
    {schema:REPRO_SCHEMA, label:`repro:${i}`, phase:'Reproduce'})
))).filter(Boolean).flatMap(r=>r.checks)

phase('Verify')
const order = {Critical:0, High:1, Medium:2, Low:3}
const highCrit = finds.filter(f=>f.severity==='Critical'||f.severity==='High').sort((a,b)=>order[a.severity]-order[b.severity]).slice(0,CAP)
const LENSES = ['correctness','statistics','reproduce-it']
const verified = (await parallel(highCrit.map(f => () =>
  parallel(LENSES.map(lens => () =>
    tryAgent(`${PRE}\nADVERSARIAL VERIFIER (lens: ${lens}). Default to REFUTING this finding; only confirm if the evidence truly holds against the actual files. Re-read the cited code/data yourself.\nFINDING: ${JSON.stringify(f)}\nFrom the ${lens} angle, is this real? Return verdict, corrected_severity, reasoning citing file:line.`,
      {schema:VERDICT_SCHEMA, label:`verify:${f.severity}`, phase:'Verify'})
  )).then(vs => {
    const v = vs.filter(Boolean)
    const ref = v.filter(x=>x.verdict==='refuted').length
    const con = v.filter(x=>x.verdict==='confirmed').length
    const verdict = ref>=2 ? 'refuted' : (con>=2 ? 'confirmed' : 'uncertain')
    return { key:(f.file||'')+'|'+(f.title||''), verdict }
  })
))).filter(Boolean)
const verdictByKey = new Map(verified.map(v=>[v.key, v.verdict]))
for (const f of finds) { const k=(f.file||'')+'|'+(f.title||''); f.adversarial_verdict = verdictByKey.get(k) || (f.severity==='Critical'||f.severity==='High' ? 'uncertain' : 'unverified') }

phase('Synthesize')
const survivors = finds.filter(f=>f.adversarial_verdict!=='refuted')
const report = await tryAgent(`${PRE}\nYou are the SYNTHESIZER for team ${A.team} (${A.title}). Produce the team report.\nFINDINGS (post-verification; ${finds.length} total, ${survivors.length} not-refuted):\n${JSON.stringify(finds).slice(0,90000)}\nNUMBER-REPRODUCTION CHECKS:\n${JSON.stringify(repro)}\nWrite markdown_report per repo audit rubric: (1) Executive summary + bottom-line verdict on this domain's correctness/rigor/reproducibility; (2) Tally table by severity (count NOT-refuted in headline; note how many refuted); (3) Number-reproduction table; (4) Findings grouped by severity with file:line, evidence, adversarial_verdict, recommendation; (5) Literature questions to route to team H10. Fill the structured tally with NOT-refuted counts.`,
  {schema:REPORT_SCHEMA, label:'synthesize', phase:'Synthesize'})

return { team:A.team, findings:finds, repro, report }
