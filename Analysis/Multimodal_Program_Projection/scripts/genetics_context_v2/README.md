# Genetics/context v2: source and power preflight

This directory implements the Wave-2, claim-independent portion of
`docs/archive/plans/2026-08-07_paper_program/30_GENETICS_CONTEXT_AND_FIG3.md`.
Current paper authority is `docs/PAPER.md`; the archived plan is retained for
method-level provenance only.
It freezes and re-derives the static genetics/state contract, audits the 35
primary GWAS phenotypes, and records what the deposited Broadaway liver-eQTL
lead table can and cannot establish. It also contains an isolated,
public-source-only GEN-03 audit of Hong et al. (*Nature Genetics*, 2025;
doi:10.1038/s41588-025-02237-8).

Neither the preflight nor GEN-03 runs context rescue, matching,
phenotype-stratified outcome tests, or figure assembly. In particular:

- the deposited Broadaway lead table establishes 9,013 significant signals
  for 6,564 source-significant eGenes;
- it does not provide the complete 18,322-gene source expression/testing
  universe described in the paper;
- the local 19,072-gene OTTERS annotation and 18,975-row downstream COLOC
  aggregation are derived analysis universes and are not substituted for that
  missing source denominator;
- absence from that lead table is therefore `indeterminate`, never a
  non-eGene or powered negative; and
- downstream COLOC `n_snps` is retained only as assay-coverage metadata, not
  as eGene strength.

The GEN-03 public-source gate is terminally `coverage_limited`:

- official Nature supplementary files, Zenodo record 14586466, 25 bounded
  archive entries, and official code commit
  `650e9fe46f57a0889b0b6ee0f7fd35f6fc7f4730` are pinned and verified;
- 48 source donors and the 44-donor eQTL subset are reproduced (23 control,
  4 MASL, 7 eMASH, and 10 aMASH among included donors);
- all 23 deposited interaction contexts and 1,299,254 state-specific tested
  gene–variant pairs are present;
- Zenodo and Supplementary Table 11 contain the same 601 unique regulatory
  quartet tuples;
- none of 16 prespecified flat/hierarchical BH reconstructions using either
  the specified likelihood-ratio statistic or the strongest rival
  interaction-coefficient statistic reproduces the reported 13,648 calls and
  deposited cell-type union flags; and
- complete base liver-eQTL backgrounds and required matching covariates were
  therefore not assembled in this bounded run.

Consequently, the source tables can be cited for cohort/method provenance and
the deposited 601 examples only. GEN-03 authorizes no source-significant gene
call, rescued fraction, context bridge, enrichment, or genetic/context
negative claim.

Identifier units are explicit: 6,564 is the source-wide unique Ensembl-eGene
count, 6,583 is the mapped symbol-annotation row count, and 6,562 is the count
of unique source Ensembl identifiers represented by those rows. These values
are not interchangeable in tables or captions.

All runtime products go beneath the isolated candidate root
`Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/genetics_context/`.
Canonical GWAS, atlas, and figure paths are rejected by a path guard.

## Entry point

Submit only after reviewing the pinned input hashes:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2/run_genetics_preflight.sbatch
```

The wrapper uses a 48-hour CPU allocation and a single-word `coloc` job name.
It was created for submission but is not submitted by this implementation.

The Hong source audit is a separate, fail-closed chain:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2/run_hong_source_acquisition.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2/run_hong_source_audit.sbatch
sbatch Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2/run_hong_source_validation.sbatch
```

These jobs use the single-word names `zenodo` and `eqtl`, request at least 48
hours, and write only beneath the isolated candidate root. The acquisition
uses HTTP byte ranges plus ZIP CRC32 verification; it does not download the
6.6 GB archive wholesale.

The terminal closure is a separate, non-mutating final step:

```bash
sbatch Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2/run_genetics_terminal_closure.sbatch
```

It pins the terminal Hong gate and its independent validation, preserves the
historical preflight gate byte-for-byte, and supersedes the stale GEN-03 rows
only for downstream gate selection. The 48-hour `eqtl` wrapper is prepared
for review and has not been submitted.

For lightweight fixtures only:

```bash
bash Analysis/Multimodal_Program_Projection/scripts/genetics_context_v2/tests/run_tests.sh
```

## Stage contract

1. `01_freeze_and_rederive.py` verifies pinned inputs and freezes the primary
   evidence classes.
2. `02_build_phenotype_registry.py` applies the prespecified phenotype/source
   registry and cross-ancestry limitation labels.
3. `03_build_eqtl_observability.py` reproduces the source-positive eGene
   deposit and assembles gene-level observability metadata.
4. `04_power_interface.R` emits the all-joint interface and explicit skip rows
   for source universes that cannot be reconstructed.
5. `05_validate_preflight.py` independently re-derives load-bearing counts,
   verifies the fail-closed states, and seals a checksummed preflight manifest.
6. `06_fetch_hong_public_sources.py` retrieves and CRC-verifies only the 23
   complete interaction tables, source annotation, and quartet table.
7. `07_audit_hong_public_contract.py` reproduces donor/context/quartet
   contracts and attacks the unspecified interaction multiplicity family.
8. `08_validate_hong_source_gate.py` independently re-parses the complete
   public interaction universe, re-derives the primary and strongest-rival
   families, checks quarantine exclusion, and seals the terminal verdict.
9. `09_close_genetics_terminal_gate.py` creates the hashed, non-mutating
   downstream closure and identifier-unit contract.
10. `10_validate_genetics_terminal_closure.py` independently rederives the
    authority hashes, stale-row supersession, and 6,564/6,583/6,562 boundary.

No favorable biological result is an acceptance criterion.
