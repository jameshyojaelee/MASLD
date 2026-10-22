# GSE267145 H3K27ac region coordinates: frozen as 1-based inclusive by structural inference

Date: 2026-09-06

Status: PROPOSED — written by the analysis session, awaiting owner adoption. The frozen evidence
tree already exists and is immutable; this record documents it and the one config sentence that
describes it.

Decision: the deposited H3K27ac region labels of GSE267119 (`chr1:96466-96736`, 96,460 regions)
are read as **1-based inclusive**. The decision rests on two prespecified, falsifiable structural
tests that agree, not on a primary-source statement, and it is voided by any later primary-source
statement that contradicts it.

Evidence, from `Analysis/MASLD_Model_Benchmark/executions/gse267145-coordinate-convention-21473307/`
(ARTIFACTS sha256 `23931bcfd5d4116efbcc95bd8b0af337006a99aadb6950f77c5cab3d34a9f754`):

- Test A, chromHMM grid and offset, on the official Supplementary Table 3 annotation columns
  against the independent 0-based Roadmap E066 segmentation: 14,320 rows tested; inclusive widths
  are multiples of 200 in **14,140** rows against **58** under a half-open reading; `start − 1`
  reproduces an E066 breakpoint **7,596** times, `start` **0** times, `start + 1` **0** times.
  Under a 0-based file these four counts swap wholesale, so the test can fail.
- Test B, `GenomicRanges::reduce` postcondition, on the deposited axis itself: 96,436
  same-contig adjacent pairs, minimum gap **2**, **zero** pairs at gap 0 or 1, five at gap 2.
  `bedtools merge` on 0-based intervals leaves gap-1 pairs (about seven expected here); observing
  none is p of order 1e-3.
- GEO `!Sample_data_processing` records MACS2 2.2.5 → ENCODE IDR v2.2.2 → merged in R 4.2.1 →
  Rsubread 2.10.5, whose SAF annotation is 1-based inclusive and whose GeneID is the deposited key.

What changed, as one coordinated edit:

- `scripts/audit_gse267145_coordinate_activation.py` gains `chromhmm_grid_and_offset_test`,
  `reduce_postcondition_gap_test` and `coordinate_structural_evidence`. `coordinate_claims()` is
  byte-identical.
- `tests/unit/test_audit_gse267145_coordinate_activation.py` gains five tests, including one that
  proves the structural tests return `zero_based_half_open` on a 0-based fixture and one that
  proves a contradicting primary-source sentence voids the freeze. The five existing tests are
  untouched.
- `scripts/freeze_gse267145_coordinate_convention.py` and
  `slurm/freeze_gse267145_coordinate_convention.sbatch` produce the frozen tree.
- `config/datasets/gse267145_znf469_human_liver.toml`: `blockers[0]` now states the resolution
  and cites this record. No other authority file changes.
- `config/evaluation/gse267145_gene_anchored_h3k27ac_state_task.toml` and its promotion gate are
  registered as a new sibling TaskSpec. The two existing GSE267145 TaskSpecs, whose hashes are
  pinned in six places each, are not edited.

What did **not** change, deliberately:

- `coordinate_semantics_resolved` stays `False` everywhere. That flag has always meant "a primary
  source stated the convention," and three guards (`activate_gse267145_histology_state.py:390`,
  `build_gse267145_paired_molecular_fixture.py:507`, the coordinate contract test) require it to
  be `False`; flipping it would make the paired fixture unbuildable. The freeze is recorded under
  its own keys: `coordinate_convention_frozen`, `resolution_class`,
  `primary_source_statement_found = false`.
- `bounds_compatibility_used_as_semantic_evidence` and
  `tool_default_or_file_extension_used_as_semantic_evidence` are asserted `false` in the frozen
  tree even though the freeze succeeded. Neither was used.
- The activation ARTIFACTS `e5042eb4…` and TaskSpec `8bf05fe6…` hashes are unchanged.

Residual uncertainty, stated plainly. Test B proves the merge and count steps ran under 1-based
arithmetic on the integers as given. It cannot exclude that the IDR narrowPeak entered R without
the `+1` that `rtracklayer::import` applies; in that case the labels are BED numbers
mis-provenanced but counted 1-based. At gene scale this is invisible — 7 of 96,460 regions change
block at 25 kb. For sequence, `fetch(chrom, start − 1, end)` returns the interval the reads were
counted over under either scenario. Test A speaks to the supplement's annotation columns and
reaches the deposited keys only through the inference that one R session wrote both; Test B is the
one that touches the deposited axis directly. The freeze rests on their conjunction.

What this authorises:

- Tier 1, now: interval overlap, gene-anchored aggregation, genomic context, distance to TSS, with
  a declared ±1 bp tolerance and no threshold applied at an exact boundary without recording the
  tie count.
- Tier 2, with a declaration: base-resolution sequence extraction where the artifact records
  `sequence_semantics = "counted_interval"` and the ±1 bp residual.
- Still blocked: any assertion of base-resolution correspondence to the *called* peak, and motif
  calls at region edges without a ±1 bp sensitivity check.

A query has been drafted for the depositor asking the base-count question directly. A reply
confirming 1-based upgrades `resolution_class` to `primary_source_statement`. A reply contradicting
it voids this record and every number derived under it; that consequence is accepted in advance.

Revoking this: delete this record, restore `blockers[0]` in the dataset TOML to its 2026-08-24
text, and remove the sibling TaskSpec and gate. The frozen tree stays on disk as a record of what
was measured; nothing in the manuscript authorities, figures or validator depends on it yet.
