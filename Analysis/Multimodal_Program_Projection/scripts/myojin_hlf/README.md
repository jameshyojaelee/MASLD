# Myojin HLF outcome-isolation firewall

This directory implements the sealed phases A/B firewall and a separately
frozen, one-pass phase-C executor for Plan 40. Phase C cannot run unless both
`SEALED` and `FIREWALL_READY` authorize the exact core specification and
`PHASE_C_CODE_READY` authenticates every executor, statistical, validation,
wrapper, fixture, environment, and frozen-input hash.

The trust boundary is explicit:

1. `00_mask_myojin_workbook.py` is the source custodian. It authenticates the
   raw workbook and exports only gene identifiers, the two neutral `num`
   fields, and missingness indicators. MAGeCK effects, scores, p values, FDRs,
   ranks, and `goodsgrna` values never enter the blind bundle.
2. `01_extract_depmap_hlf.py` extracts HLF expression and generic Chronos from
   one version-matched public DepMap 24Q4 triplet. Chronos is not relabeled as
   Myojin vehicle fitness.
3. `02_freeze_blind_spec.py` freezes the target universe, four evidence
   classes, Plan-20 programs and weights, expected directions, covariates,
   models, permutation families, seeds, known-hit exclusions, and the Figure 5
   magnitude gate. It writes `SEALED` only when every source gate passes;
   otherwise it writes `BLOCKED` and leaves phase C unauthorized.
4. `03_validate_blind_spec.py` independently re-derives hashes, membership,
   eligibility, leak checks, and seal logic. `FIREWALL_READY` is possible only
   for a fully sealed bundle.
5. `09_freeze_phase_c_code.py` authenticates the exact execution environment
   and freezes the reviewed code bundle while the outcome remains unread.
6. `10_execute_one_pass.py` reauthenticates all twelve source gates, the
   official DepMap triplet, HLF identity, raw-workbook hash, design rank, and
   code hashes before it writes the append-only unseal record. It then reads
   the screen once and executes every prespecified primary and sensitivity
   family without tuning.
7. `11_validate_phase_c.py` independently re-derives the tested universe,
   observed OLS/HC3 statistics, protective-hit calls, BH families, program
   scores, Figure-5 verdict, and release hashes. Biological positivity is not
   a software acceptance condition.

Pinned DepMap input is the public CC-BY-4.0 release `DepMap Public 24Q4`, DOI
`10.25452/figshare.plus.27993248.v1`. The complete matching triplet is used to
avoid combining the unversioned local Chronos snapshot with expression from a
different release. The Figshare API metadata fixes file IDs, byte sizes, and
published MD5 values in `01_extract_depmap_hlf.py`; the adapter independently
checks all three before writing `DEPMAP_READY`. The local Model/Chronos pair is
recorded only as a rejected preflight because the repository does not carry a
release manifest or a matching HLF expression file.

The pinned 24Q4 matrices contain two real cross-file annotation conflicts:
`MEF2B` is Entrez 100271849 in Chronos versus 4207 in expression, and `RLN2`
is 6019 versus 6013. The adapter never merges these by symbol. It records both
IDs in `depmap_gene_mapping_audit.tsv`, blanks both HLF covariates, and makes
the affected screen genes untestable. This is a gene-level fail-closed mapping
rule; it does not weaken the release-level source gate.

The 80% class gate measures finite numeric TPM plus Chronos availability before
the biological TPM>=1 floor. `covariate_coverage_audit.tsv` reports both values
side by side. Low-but-observed TPM is complete for source-power auditing but is
still excluded from the primary tested universe, preserving the detectability
rule without mislabeling low expression as missing data.

All outputs remain under:

```text
Analysis/Multimodal_Program_Projection/candidates/
  program-context-v2-candidate-2026-08-07/myojin_hlf/
```

The official DepMap gate passed and the blind specification is sealed at
`1f99b7a7a95d9aebbbb911bc2317860f4a3fada0d1da2a4932d5c1ab190c4cdd`.
The phase-C code bundle was frozen at
`51650e5d21f4a16f1e67ab30ee27bcd5d2ed3a697912b1a0518f43195a629d0b`
with outcomes unread. SLURM job `19639096` (`permutation`, 72-hour request)
completed in 10 minutes 23 seconds with exit code 0, then wrote
`PHASE_C_VALIDATED`; all 24 independent release checks passed. The mechanical
Figure-5 verdict is `FALSE`, so this release is a complete supplementary
assay-specific non-support result rather than a main-figure rescue. The
predeclared failure behavior remains: a pre-unseal gate failure creates
`phase_c_terminal_blocker.tsv` without opening the outcome workbook, while a
post-unseal error creates `phase_c_terminal_failure.tsv` and forbids rerunning
the executor under the same release ID.
