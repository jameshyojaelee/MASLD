# Codex to Claude: C7 read-only documentation check, 2026-09-29

**Wrong or unsupported numbers in the assigned passages: none found.** No
authority was edited.

I checked `docs/STATUS.md:366-368`, the A1 identity paragraph, Model B test 1
and closure text in `docs/technical/NEW_MODELS_PROGRAM_2026-09-23.md`, and
`Analysis/MASLD_Model_Benchmark/MODELS.md:504` against their named inputs.

- A1 summary: 1,466 libraries, 1,417 with enough identity sites, 1,349
  individuals; within-cohort same-person pair counts 25/36/5/3; 10 first-degree,
  4 second-degree, 23 sex mismatches, 57 ambiguous expression-sex calls.
- The 844-fit identity roster gives at most 813 people; the Figure 3D source
  table has 416 libraries from 404 called people (303 advanced from 291 and
  113 controls from 113). The stage-contrast manifest reproduces every
  library→person count in the assigned text, including four same-person
  cross-arm placements in NAS 3–4 versus 5–8 and one in F0 versus F1.
- Local GSE213621 metadata contain 367 libraries from 342 identity-called
  people, including 299 staged libraries from 274 people. The related-pair
  table has 25 same-person pairs; 11 differ in stage labels with the stated
  6/4/1 category split. The 361-record upper bound of 342 people is valid.
- Model B's corrected input has 1,007 records, 898 with an identity call;
  860 distinct called people plus 109 records without a call give the stated
  at-most-969 bound if provisional GSE193066 calls hold. Treating six such
  calls as unconfirmed gives at most 975. The test-1 source gives gain
  0.0279663, p=1/201, and positive gain in five of seven cohorts. The
  901-record sensitivity excluding GSE193066 has weighted gain 0.0233987.
  Test 3's four p-values and BH q=0.981 agree with its summary. The NUTS
  source reports zero divergences and max R-hat 4.5506 with `converged=false`.
- The stated GSE193066 132/164 possibly-mixed count matches the seal's A1
  identity rule; its provisional qualification is retained. The Model B
  component attribution wording follows PRESPEC_B §7's GSE281797
  untestability declaration and does not claim characterized axes.

The source paper's 368-human-subject statement and the GEO record's 368
samples were not re-retrieved in this C7 pass; the local 367-library and
342-person counts above were checked directly. This limits any independent
claim about the missing GEO sample, but does not contradict the assigned text.
