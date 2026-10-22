# Receiving H3 read recipe: bounded independent review

Reviewed 2026-09-30. Source and public roster metadata only; no receiving read,
BAM, count, target or prediction values were opened. No jobs, numerical fits,
installations or changes to the owner's source were performed.

## Finding

The inspected recipe preserves the declared receiving assay and frozen prediction
comparison. No source-level reason to change the model, target units, fixed
regions or biological donor set was found. Successful acquisition/counting would
establish available measurements, not predictive accuracy or independent
participant identity. The owner's paired-profile report already documents the
recipe and queued job; this review adds the precise limits of its storage and
alignment guarantees.

Inspected snapshots:

- `scripts/analysis/in_biopsy_allelic/rec_codex/codex_receiving_h3_reads.py`:
  SHA256 `f641ea6485515b7f986f3ce279f449c5e6b85e99f676c3a5ac959fc93bf8b509`.
- `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_h3_reads.sbatch`:
  SHA256 `57603cd755e2948ed90d0ec78a76d8e9f9b283a44888b7fe0e2dc0c99f2d4eb8`.

## Contract preserved in source

The helper requires the completed reference index and B1 plus four exact RNA
groups before download (lines70–115,213). Each of the 17 donors must have native
RNA quantification and both frozen predictions, global `rrr` and local
`rrr_offset_cis`, using `raw` transport. Receipt, protocol and prediction bytes
are hash-bound; prediction arrays are not deserialized here. The prediction
producer, `codex_native_rna_pilot.py:172–194`, checks complete gene coverage and
the exact 96,460-region prediction axis. Hash checking does not itself rederive
that prediction calculation.

The H3 recipe joins all 202 runs to their frozen public donor/GSM/experiment/
BioSample, URL, MD5 and size metadata, explicitly requiring SINGLE ChIP-Seq
(lines116–140). It selects H3 rows, excluding Input. The declared source bytes
are 56,783,483,922. No alternate runs, assay pooling or donor substitutions are
available in this recipe.

The unchanged SAF is hash-bound. Bowtie2 sensitive/end-to-end single-read
alignment is followed by `samtools view -F2308 -q10`: unmapped, secondary and
supplementary records are excluded (lines243–249). FeatureCounts uses SAF,
unstranded primary reads, MAPQ10 and one-region assignment, with no `-O`, `-M`,
paired-fragment counting, Input subtraction or PCR deduplication (lines252–254).
Technical runs remain separate count columns with donor bindings. The helper
does not assemble the final axis, transform B99 targets, score predictions or
fit models (lines255–264). The declared count unit is selected mapped single
reads, not original molecular fragments; CUT&RUN-to-ChIP transport remains
assay-qualified (lines228–231). Biological n remains 17.

## Two limits to retain when reporting completion

1. **Storage is periodically monitored, not strictly quota-enforced.**
   `bound()` checks saved files against 230GB and donor caps; subprocess checks
   occur every five seconds (lines148–150,183–206). Downloads check every 64MiB
   against 250GB (lines168–173). File growth between checks is not bounded by
   these checks. Consequently `hard_cap_bytes=250GB` at line227 is a declared
   allowance, not proof that saved bytes can never exceed it. The 20GB buffer
   is precautionary. Metadata arithmetic gives aggregate donor caps
   317,917,419,610 bytes, above the global sampled stop. [Unverified] Completion
   of all 17 donors within the global allowance is not established by the cap
   formula. Report actual saved bytes and any stop without treating missing
   donors as successful or changing the assay recipe to finish.

2. **A selected primary alignment is not proof of unique locus origin.**
   Default one-best Bowtie2 output, MAPQ10 and absence of featureCounts `-M`
   preserve this prespecified recipe; they do not establish that every retained
   read has only one possible genomic origin. Keep the explicit primary-read
   terminology. Do not introduce unique-molecule or original-fragment claims.
   Changing reporting language needs no assay rerun; changing the filter would
   be a different count recipe.

The launcher requests one CPU-partition nslab job, eight CPUs, 64GB and 90h
(lines2–10), without an array. Alignment uses seven workers plus streamed BAM
conversion; featureCounts uses eight. Peak RAM, runtime, disk sufficiency and
actual receiving accuracy were not tested by this source review. Existing
primary full-MSE evaluation and assay/source identity qualifications remain
unchanged. No additional executable check was run.
