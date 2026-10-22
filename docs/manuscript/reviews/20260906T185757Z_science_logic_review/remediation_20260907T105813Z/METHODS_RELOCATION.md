# Methodological content moved out of the main Results, 2026-09-07

PI direction: minimize methodological and calibration material in the main
Results; it belongs in Methods and Supplement. This is consistent with
`docs/PAPER.md`, which already states that observability and evidence states
"are the Resource's guardrails, not its central biological headline."

## Principle applied

Main Results state the result and its number. Methodological detail stays only
where it is load-bearing for believing that specific result, and then as one
clause rather than a paragraph. Model specifications and correction families go
to Methods; panel construction and denominators go to figure legends; coverage
inventories go to the Supplement.

## What moved

| from | content | to |
|---|---|---|
| `fig1.md` | the seven evidence states and what a blank comparison can mean | Methods, "Evidence states" |
| `fig1.md` | biological-unit rule (cells/spots/runs/arrays are not donors) | Methods, "Biological unit"; one clause retained in Results |
| `fig1.md` | CosMx 12/117, GeoMx 0/117, ATAC 38,066/39,914 indeterminate, zero tested-negative | Methods summary + Supplement |
| `fig3.md` | complete tested families and cohort/contrast eligibility rules | Methods |
| `fig4.md` | construction of the content-lineage index | Methods, "Program lineage labels" |
| `fig5.md` | proteome model spec (58 patients, batch/age/BMI/sex, 6,071-protein BH family) | Methods, "Liver proteome model" |
| `fig5.md` | panel construction: gray non-significant background, bar denominators | figure legends |
| `fig5.md` | matched-null construction, donor joins, array source-dependence | Methods, "Spatial matched nulls" |
| `fig5.md` | deterministic primary-signal-pair selection for the ATAC intersection | Methods, "Regulatory context intersection" |

## What deliberately stayed in Results

- **The spatial contrast itself.** That both programs are positive by raw
  Moran's I and only one exceeds its matched null is the finding of Figure 5F,
  not methodology. The construction of the null moved; the result did not.
- **The 80/117 versus 37/117 lineage split and its ambient-RNA cause.** A third
  of the program vocabulary being mislabelled changes how Figure 4 is read.
- **The chance-level overlap and its explanation.** Now stated as a result in
  Figure 1 with the number, rather than as a withheld placeholder.
- **One clause on biological unit**, because donor-level inference is the reason
  several headline counts differ from the published literature.

## Effect

`fig1.md` 70 to 56 lines; `fig5.md` 65 to 45 lines; `fig3.md` and `fig4.md`
trimmed. Main Results prose across Figures 1-5 is now 4,651 words. Nothing was
deleted: every relocated passage is in `METHODS.md` under
"Evidence rules, observability accounting and panel construction".

## Still to do

- Figure legends must absorb the panel-construction sentences removed from
  `fig5.md` (gray background, bar denominators). Not yet written.
- A supplementary observability table should carry the full per-assay inventory
  that Figure 1 no longer lists.
- `fig5_discussion.md` (205 lines) still carries substantial method commentary
  in its limitations paragraphs and has not been trimmed.
