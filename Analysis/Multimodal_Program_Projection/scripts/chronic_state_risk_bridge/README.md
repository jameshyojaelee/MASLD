# Plan 43 execution

This directory implements the isolated `chronic-state-risk-bridge-2026-08-09`
candidate. Outcome-bearing source acquisition is prohibited until
`SEALED.json` and `SEAL_VALIDATED.json` exist.

Execution begins with:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/chronic_state_risk_bridge/run_freeze.sbatch
```

The freeze reconstructs Tier-1/2 genetic drivers from the 35 per-GWAS COLOC
tables, writes the corrected 326-locus registry, and freezes the signed
117-program chronic-state axis. It never overwrites Plans 41 or 42. Once its
independent validation passes, `run_acquire.sbatch` opens three isolated public
source lanes (`human`, `crop`, and `protein`); download failures become explicit
source-gate records rather than favorable-result-dependent substitutions.

## Terminal status

Plan 43 completed `complete_no_promotion` on 2026-08-09 EDT. The terminal validator
passed 51/51 checks and sealed 141 payloads under release-manifest SHA256
`81a800fe8b23498217e7e2a0c67309331e00e01a47a311a59ce2a861a97b2757`.
No canonical write, Figure 5 replacement, or journal escalation is authorized.

Execution order after the seal:

```text
run_acquire.sbatch
  -> run_human_gates.sbatch / run_crop_gate.sbatch
  -> run_human_reversal.sbatch
  -> run_protein.sbatch / run_context.sbatch / run_clcc1_correction.sbatch
  -> run_spatial_identity.sbatch
  -> run_adjudicate_validate.sbatch
```

The most important terminal boundaries are:

- GSE106737 and GSE83452 are one overlapping cohort, not independent
  replications; the independent GSE48452 reversal test was direction-discordant
  and null.
- GSE281160 failed before state-axis outcome access because only five corrected
  Tier-1/2 loci were fully oriented (four activity-increasing, one decreasing).
- The outcome-unseen PXD052787 secretome comparison was null and failed the
  prespecified balance and continuous-bulk-effect safeguards.
- The M8/M20 donor-level identity audit supports the labels
  `hepatocyte-cluster stromal ECM program` and `ductular injury program`, but
  the prespecified gene-class spatial comparison was unavailable.

Authoritative terminal files are `VALIDATED.json`, `promotion_verdict.tsv`,
`source_gate_status.tsv`, and `release_manifest.tsv` under the candidate root.
