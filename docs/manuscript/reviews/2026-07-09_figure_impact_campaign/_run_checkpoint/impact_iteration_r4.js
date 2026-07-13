export const meta = {
  name: 'masld-review-iteration-completion',
  description: 'Complete the loop-until-dry iteration of the MASLD figure/writing impact review + re-check the just-applied edits',
  phases: [
    { title: 'Critic', detail: 'completeness critic over the 123 round-1 findings + current edited files' },
    { title: 'Augment', detail: 'per flagged unit: new findings + drafted fixes' },
    { title: 'Synthesize', detail: 'strategy addendum + convergence verdict' },
  ],
}

const ROOT = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'
const OUT = ROOT + '/docs/manuscript/reviews/2026-07-09_figure_impact_campaign'
const MAX_EXTRA_ROUNDS = 1  // FINAL verification round (round 4); rounds 1-3 already ran

const CTX_FILES =
  OUT + '/00_STRATEGY.md  (round-1 strategy: chosen headline, Fig1->Fig5 arc, top-10 moves, cross-cutting fixes, reviewer rebuttal)\n' +
  OUT + '/spine.md  (chosen headline + arc + venue verdict + competitor lines + buried gems + rubric)\n' +
  OUT + '/findings_digest.txt  (ALL 136 findings from rounds 1-3, one line each — this is what is ALREADY COVERED; do NOT re-report these)'

const EDITED_FILES =
  ROOT + '/docs/manuscript/working/abstract.md  (NEW + reconciled through round 3: drops 1,162, anchors ">90%" to the 1,031 effector set, names the 1,244 two-way transcriptomic-genetic intersection for gsMap, scopes THRB-above-NR1H4, "35 GWAS strata")\n' +
  ROOT + '/docs/manuscript/working/01_intro.md  (round-1 numbers 473/1,031, 35 GWAS, five ancestries, one Broadaway, pooled-not-mega, single-transition removed, Figure 1 callout; round-3 cash-out now MATCHES the abstract: both approved targets, THRB-above-NR1H4, RORA/HKDC1 receipts, gsMap 1.6-1.8 tissue-localization)\n' +
  ROOT + '/docs/manuscript/working/fig3.md  (round 1: deleted the orthogonality "expected" hedge)\n' +
  ROOT + '/docs/manuscript/working/fig5_discussion.md  (round 1: "expected"/TREAT removed; round 3: +Nelson 2015/Open-Targets prior-art paragraph after the "agreement is informative" lesson, +competitor-positioning paragraph (Feng/Li/Tzouanas, none named except Hmgcs2 conceded) before Limitations)\n' +
  ROOT + '/docs/manuscript/working/METHODS.md  (round 3: "mega-analysis"/"mega cohorts" prose -> "pooled"; the include_in_mega code flag preserved)\n' +
  ROOT + '/docs/paper_outline.md  (round 2: Working Title -> convergence-with-receipts spine; :25 breadth-is-the-value -> breadth-is-required; headline-reframe blockquote atop §1 demoting the cascade to a Fig3 supporting finding)'

const GROUNDING = [
  'TARGET VENUE: Nature Genetics / Nature Metabolism. The chosen headline = genetics/expression ORTHOGONALITY (92%, 1,072/1,162) -> convergence-as-necessity -> label-blind drug-gradient recovery (2 approved..1,504 preclinical) + prospective RORA/HKDC1 receipts -> inherited risk tissue-real ONLY at the convergent intersection (gsMap OR 1.55-1.78) while genetic-only is depleted (0.45-0.58).',
  'STYLE/RULES the edits and any new fix must obey: NO "TREAT" in the manuscript (say "adjusted P<0.05" or move interval-null detail to Methods); NO "mega" (say "pooled, cohort-adjusted"); human-only (no cross-species); active voice, finding-first; cite REAL files; FLAG any critique that targets a stale/superseded number (numbers churned across >=4 DEG defs, >=5 convergence builds); method<->modality (no cross-assay tool transplants).',
  'KNOWN OPEN ITEMS the critic should USE but not re-litigate as new: move #4 (Fig5C label-blind recovery is currently a COUNT, not a proven rank-enrichment; on-file null suggests P~1 -> the abstract deliberately says "recovered both approved targets + RORA/HKDC1 receipts", NOT "recovers the entire gradient"); move #5 (drug pair atlas {THRB,SLC5A2} vs text {THRB,GLP1R} unreconciled; SLC5A2 renal); move #9 ("96% non-coding of 31,036 variants" has a shaky source; gene-level 3.3% coding is sound); coloc cardinality bridge 473 subset 677 subset 1,031 subset 3,038, and the 1,162 orthogonality denominator has no clean bridge yet.',
  'COMPETITORS to test positioning against: Feng 2026 (closest scoop), Tzouanas 2026 (Cell, in-vivo HMGCS2 KO this paper cannot match), Li 2025, Broadaway 2024 (reused eQTL panel), Govaere (benchmark/circularity), Open Targets/Nelson 2015 (the "human genetics doubles success" prior art the drug-recovery headline invokes).',
].join('\n\n')

// unit -> the current file(s) an augment agent must read
const UNIT_FILE = {
  fig1: ROOT + '/docs/paper_outline.md (Fig 1 spec in §3) + note figures/main/fig1_atlas_overview/figure1_JL.png is the schematic composite',
  fig2: ROOT + '/docs/manuscript/working/fig2.md',
  fig3: ROOT + '/docs/manuscript/working/fig3.md',
  fig4: ROOT + '/docs/manuscript/working/fig4.md',
  fig5: ROOT + '/docs/manuscript/working/fig5_discussion.md (Results portion)',
  thesis_arc: ROOT + '/docs/paper_outline.md §1',
  abstract: ROOT + '/docs/manuscript/working/abstract.md',
  intro: ROOT + '/docs/manuscript/working/01_intro.md',
  discussion: ROOT + '/docs/manuscript/working/fig5_discussion.md (Discussion + Limitations portion)',
}

const CRITIC_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    gaps: { type: 'array', items: {
      type: 'object', additionalProperties: false,
      properties: {
        unit: { type: 'string', description: 'one of: fig1 fig2 fig3 fig4 fig5 thesis_arc abstract intro discussion (or "cross-cutting")' },
        gap: { type: 'string', description: 'a genuinely NEW high/medium gap NOT already in findings_digest.txt' },
        priority: { type: 'string', enum: ['high', 'medium', 'low'] },
        rationale: { type: 'string', description: 'why it is new + why it matters at the NG/NatMetab bar; cite the file' },
      }, required: ['unit', 'gap', 'priority', 'rationale'],
    } },
  }, required: ['gaps'],
}

const AUG_FINDING = {
  type: 'object', additionalProperties: false,
  properties: {
    target: { type: 'string' }, axis: { type: 'string' },
    severity: { type: 'string', enum: ['Critical', 'High', 'Medium', 'Low'] },
    leverage: { type: 'string', enum: ['high', 'medium', 'low'] },
    problem: { type: 'string' },
    concrete_fix: { type: 'string', description: 'FOR WRITING: the DRAFTED ready-to-paste sentence/subheading/legend (obey the no-TREAT / no-mega / active-voice rules). FOR PANELS: exact panel id + change. NO placeholders.' },
    effort: { type: 'string', enum: ['low', 'medium', 'high'] },
  }, required: ['target', 'axis', 'severity', 'leverage', 'problem', 'concrete_fix', 'effort'],
}
const AUG_SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    unit: { type: 'string' },
    extra_findings: { type: 'array', items: AUG_FINDING },
    addendum_md: { type: 'string', description: 'markdown block to append to this unit verdict' },
  }, required: ['unit', 'extra_findings', 'addendum_md'],
}

async function tryAgent(prompt, opts, attempts) {
  const n = attempts || 2
  for (let i = 0; i < n; i++) {
    const r = await agent(prompt, opts)
    if (r !== null && r !== undefined) return r
    log('retry ' + (opts.label || '') + ' attempt ' + (i + 1))
  }
  return null
}

function criticPrompt(roundNo, priorNewText) {
  return 'You are the COMPLETENESS CRITIC (round ' + roundNo + ' — FINAL VERIFICATION) for an adversarial IMPACT+ARGUMENT review of a MASLD atlas paper targeting Nature Genetics / Nature Metabolism. Rounds 1-3 produced 136 findings; the highest-leverage fixes have been applied (abstract + intro reconciled and mutually consistent; fig3/fig5 de-hedged; fig5 Discussion gained a Nelson/Open-Targets prior-art paragraph and a competitor-positioning paragraph; METHODS "mega"->"pooled"; outline retitled and reframed). This is the FINAL pass: CONFIRM the review has converged (return an EMPTY gaps array if so), or surface any genuinely remaining high/medium gap OR any NEW problem the round-2/3 edits introduced (a broken number, an internal contradiction between the just-edited abstract/intro/discussion, an over-cut hedge, a dangling reference). Be strict about "new" — anything already in findings_digest.txt does not count.\n\n' +
    'GROUNDING:\n' + GROUNDING + '\n\n' +
    'STEP 1 - read the round-1 context (know what is already covered):\n' + CTX_FILES + '\n\n' +
    'STEP 2 - read the CURRENT (post-edit) files and scrutinize them:\n' + EDITED_FILES + '\n\n' +
    (priorNewText ? ('Also already surfaced EARLIER THIS ITERATION (do not repeat):\n' + priorNewText + '\n\n') : '') +
    'TASK: Return ONLY genuinely NEW high/medium gaps that are NOT in findings_digest.txt and NOT listed above. Hunt specifically for: (a) a load-bearing claim never attacked; (b) a competitor comparison never made (Feng/Tzouanas/Li/Broadaway/Govaere/OpenTargets); (c) a CROSS-FIGURE inconsistency — including any the edits just created (e.g. the abstract or intro now disagreeing with a figure, the retitle/reframe not propagated into fig2/fig4/discussion, the 1,162 vs 1,031 vs 473 coloc denominators, modality 6-vs-7); (d) a place the paper still UNDER-claims a strong true result; (e) an arc-level gap (does Fig1->Fig5 now actually build to the headline?); (f) any NEW defect the edits introduced (a broken number, an over-reach, a dangling reference, an over-cut hedge). If you find nothing new of high/medium leverage, return an EMPTY gaps array — that is the correct answer if the review has converged. Do not manufacture low-value gaps to look busy.'
}

function augPrompt(unit, gaps) {
  const file = UNIT_FILE[unit] || '(cross-cutting — read whichever of the manuscript files the gaps reference)'
  return 'You are addressing completeness gaps a critic flagged for unit "' + unit + '" in an adversarial impact review of a MASLD atlas paper (Nature Genetics / Nature Metabolism target).\n\n' +
    'GROUNDING:\n' + GROUNDING + '\n\n' +
    'Read the current file(s) for this unit: ' + file + '\nAlso available for context: ' + OUT + '/spine.md and ' + OUT + '/00_STRATEGY.md\n\n' +
    'FLAGGED GAPS:\n- ' + gaps.join('\n- ') + '\n\n' +
    'TASK: For each gap, produce a surviving finding with a CONCRETE, ready-to-use fix: for writing = the DRAFTED replacement/added sentence, subheading, or legend (obey no-TREAT, no-mega, active-voice, human-only); for a panel/figure = the exact panel id + precise change. NO placeholders. Verify the target claim exists in the CURRENT text before asserting it (flag if a gap targets a stale number). Also return addendum_md: a short markdown block summarizing these additions for this unit verdict.'
}

function synthPrompt(newFindings, converged, roundsRun) {
  const digest = newFindings.length
    ? newFindings.map(f => '- [' + f.unit + '] (' + f.leverage + '/' + f.severity + '/' + f.axis + ') ' + f.target + ' -> ' + f.concrete_fix).join('\n')
    : '(no new findings — the review converged)'
  return 'You are the CAPSTONE STRATEGIST writing the ITERATION ADDENDUM to a MASLD atlas impact-review strategy (Nature Genetics / Nature Metabolism target).\n\n' +
    'Round 1 produced 123 findings + `00_STRATEGY.md`. This iteration ran ' + roundsRun + ' further completeness round(s) and ' + (converged ? 'CONVERGED (a round returned no new high-leverage findings).' : 'hit the round cap without a fully dry round.') + '\n\n' +
    'NEW findings surfaced this iteration:\n' + digest + '\n\n' +
    'For grounding, read ' + OUT + '/00_STRATEGY.md (the existing strategy) and ' + OUT + '/spine.md.\n\n' +
    'TASK: Write a markdown ADDENDUM with: (1) a one-paragraph CONVERGENCE VERDICT (did the review reach dry? what does that mean for confidence in the round-1 strategy?); (2) any REVISIONS to the round-1 Top-10 that the new findings force (promote/demote/add — be explicit, reference the round-1 move numbers); (3) a short ranked list of the NEW actionable items with their drafted fixes; (4) an explicit note on whether the recent edits (abstract.md, intro, fig3, fig5, outline) held up or introduced anything that must be corrected. Be decisive and concrete.'
}

// -------------------------------------------------------------------- run
const accumulated = []
const accumulatedAddenda = []
let priorNewText = ''
let converged = false
let roundsRun = 0

for (let r = 1; r <= MAX_EXTRA_ROUNDS; r++) {
  const roundNo = r + 3  // FINAL verification = round 4 (rounds 1-3 done)
  phase('Critic')
  const critic = await tryAgent(criticPrompt(roundNo, priorNewText),
    { label: 'critic:round' + roundNo, phase: 'Critic', schema: CRITIC_SCHEMA, effort: 'high' })
  const gaps = (critic && critic.gaps) ? critic.gaps.filter(g => g.priority !== 'low') : []
  log('round ' + roundNo + ' critic: ' + gaps.length + ' new high/med gaps')
  if (gaps.length === 0) { converged = true; roundsRun = roundNo; break }
  roundsRun = roundNo

  phase('Augment')
  const byUnit = {}
  for (const g of gaps) { (byUnit[g.unit] = byUnit[g.unit] || []).push(g.gap + ' [' + g.priority + ']' + (g.rationale ? (' -- ' + g.rationale) : '')) }
  const flagged = Object.keys(byUnit)
  const augs = await parallel(flagged.map(uk => () =>
    tryAgent(augPrompt(uk, byUnit[uk]), { label: 'aug:r' + roundNo + ':' + uk, phase: 'Augment', schema: AUG_SCHEMA, effort: 'high' })))
  const roundNew = []
  augs.forEach((a, i) => {
    if (!a) return
    const uk = flagged[i]
    for (const f of (a.extra_findings || [])) roundNew.push(Object.assign({}, f, { unit: uk, round: roundNo }))
    if (a.addendum_md) accumulatedAddenda.push('### ' + uk + ' (round ' + roundNo + ')\n' + a.addendum_md)
  })
  accumulated.push(...roundNew)
  priorNewText += '\n[ROUND ' + roundNo + ' additions]\n' + roundNew.map(f => '[' + f.unit + '] (' + f.leverage + '/' + f.severity + '/' + f.axis + ') ' + f.target).join('\n')
  log('round ' + roundNo + ' augment: +' + roundNew.length + ' findings across ' + flagged.length + ' units')
}

phase('Synthesize')
const addendum = await tryAgent(synthPrompt(accumulated, converged, roundsRun),
  { label: 'iteration-addendum', phase: 'Synthesize', effort: 'high' })

return {
  converged,
  roundsRun,
  new_findings_count: accumulated.length,
  new_findings: accumulated,
  unit_addenda: accumulatedAddenda,
  addendum_md: addendum || '(addendum synthesis failed)',
}
