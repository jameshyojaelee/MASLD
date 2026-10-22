# Independent Codex review of the actual four RNA batches

This review used Python source, small public metadata and actual protocol
receipts only. No FASTQ, estimated-expression table, prediction array,
receiving H3 measurement or predictive outcome was read. No live script,
job, or other root's document was changed.

## Population and recipe checks passed

The actual four receipts are under
`codex_native_rna_batch_{21998814,21998815,21998816,21998817}/protocol.json`
in the REC execution root. Their combined donor population contains all
sixteen remaining fixed donors exactly once and excludes completed B1.
Every donor has eight technical runs and sixteen FASTQs in the unchanged
frozen public roster. The totals independently reproduce:

| Actual job | Fixed donors | Compressed raw bytes | Declared saved cap |
|---|---|---:|---:|
| 21998814 | B7,B8,B9,B10 | 17,592,037,059 | 25,612,037,059 |
| 21998815 | B12,B15,B19,B21 | 19,745,789,824 | 27,765,789,824 |
| 21998816 | B22,B24,B26,B36 | 17,426,600,755 | 25,446,600,755 |
| 21998817 | B38,B41,B46,B47 | 18,774,126,559 | 26,794,126,559 |
| **Total** | **16 unique donors** | **73,538,554,197** | **105,618,554,197** |

Actual wrapper source SHA256 is
`a4980388792007611dac5cae3c5c094ebdb120835cab61dea1863dd199d64247`.
All four archived wrapper copies and their copied helper sources match the
receipt hashes. The input-hash dictionaries, successful index-summary hash,
and native index population/length checks match B1's saved protocol exactly.

Source inspection confirms the same per-donor paired-file pooling, Salmon
1.10.3, library A before mates, eight threads, no new bias flags, unchanged
native transcript/gene admission and exact tximport expected-count aggregation.
The same frozen scorer requires all 42,163 modeled genes and generates raw
global `rrr` and local `rrr_offset_cis` profiles on 96,460 regions. No new
calibration, fitting, target-dependent selection, row padding or H3 read is
present. Failures remain recorded, and full-cohort completion stays false
in these batch summaries.

Each launcher requests eight CPUs, 64 GB and 90 hours, with donors processed
serially inside its allocation. This matches the successful pilot's resource
request; all four together request 32 CPUs/256 GB. Peak memory was not measured
by this review. No claim of predictive validation follows from recipe agreement.

## Bounded limitations and completion checks

Line references below concern the actual `codex_native_rna_batch.py` source.

1. **Keep the original recipe files fixed while jobs run.** Lines 51–61
   execute the copied pilot bytes but its ordinary imports resolve the index
   and fixed-artifact helpers from their original paths. The R adapter and
   scorer also retain original paths. The source-hash checks before/after
   each donor detect persistent changes, but copied files alone do not make
   these jobs independent of later edits to the original recipe files.
   No source drift was observed in the inspected receipts/current copies.

2. **Shared index checks establish metadata identity, not complete payload
   hashes.** Lines 70–79 verify successful build metadata, recipe, summary,
   native names and lengths. They do not hash every search-index payload
   file or reread the full reference FASTAs. These checks match B1's successful
   qualification, but should not be described as checking all binary index
   bytes before/after. Preserve the completed index unchanged while the four
   batches use it. No index alteration was observed or inferred here.

3. **Measure final saved sizes after writing summaries.** Lines 209–217 check
   donor/group bytes before writing the donor summary; lines 225–229 write
   final state/summary without a final inclusive cap assertion. The definite
   omission is that donor-summary bytes are absent from that donor's recorded
   completion size. The 20 MB group receipt buffer makes whole-group excess
   from ordinary JSON unlikely; no overrun is claimed. After completion, a
   metadata-only check of each whole donor directory and group directory,
   including summaries, should verify their declared caps before accepting
   completion. It requires no expression or prediction read and no rerun.

The source/metadata review found no current frozen-recipe, donor-coverage,
feature-unit or prediction-form discrepancy. Numerical outputs and full
cohort completion were not evaluated. The independently prepared single-job
wrapper is operationally superseded and must remain unsubmitted.

## Queued evaluation and RNA comparability source review

The inspected evaluation source is
`scripts/analysis/in_biopsy_allelic/rec_codex/codex_receiving_profile_evaluation.py`,
SHA256 `e12c2e617faf1c35e3093ec715bb5f1af2c018eb8cba92869ce393c330efff33`.
Its launcher SHA256 is
`5cc5f246750abde4f105d8f3c211d391617a6a8c2a872a7d05cee7fc2e0002b5`.
SLURM inspection confirms job 21998897 is pending Dependency; no evaluation
result is claimed. Only source, public identity metadata, and small protocol/
measurement receipts were inspected. No molecular array, count table,
prediction value, H3 outcome, or restricted participant value was opened.

**The RNA transformation and modeled denominator agree with the released
global/local models.** This is stronger than merely accepting the API schema:

- Native adapter `codex_native_tximport_gene_counts.R:70–91` imports Salmon
  NumReads with `countsFromAbundance="no"`, checks transcript estimates against
  native NumReads, sums them to genes, and checks direct sums/conservation.
  It does not substitute TPM, length-scaled abundance, or rounded read counts.
- `codex_native_rna_pilot.py:155–179` normalizes identifiers, sums duplicate
  normalized IDs, requires every one of the 42,163 frozen modeled genes,
  writes that exact axis, and supplies the same matrix to `rrr/raw` and
  `rrr_offset_cis/raw`, with coverage 1.0. The larger native gene universe does
  not enter the model's denominator.
- Released `score.py:167–190` joins the frozen axis and computes
  `log2(1 + 1e6 * modeled_count / sum(modeled_counts))` per donor.
  Lines 191–202 apply the fixed training transformations; lines 225–234 apply
  fixed RRR and the cis residual correction. Raw transport does not estimate
  receiving-batch means/scales. Because the caller admits the complete axis,
  its unused partial-input mean filling at lines 192–199 does not fabricate
  missing receiving genes.
- The training source `chromatin-beyond-histology-20260907T161521Z/common.py:49–51`
  defines the same log2(CPM+1) formula, and lines 70–73 require the modeled
  99-person by 42,163-gene RNA matrix. The stable producer
  `chromatin-stable-rrr-forms-20260908T192024Z/transfer_functions.py:113–114`
  constructs its training X with that function. The release builder
  `chromatin-state-v12-build-20260908T172230Z/build_v12.py:34,45–59` preserves
  the gene axis and deposits the stable RRR/cis factors. The small molecular
  receipt separately records no prior RNA normalization.

**This does not establish equivalent training/receiving RNA measurement.**
The native protocol declares Salmon 1.10.3 default quantification, paired-run
pooling and no added bias flags. Its `private_PISCES_equivalence` is false.
The saved training measurement audit
`gse267145-rna-measurement-21065097/measurement_evidence.json` likewise sets
`public_code_is_exact_private_generator=false`. Gene-summed native Salmon
expected counts have the admitted continuous-count units, but matching units,
annotation IDs, denominator and model transformations do not prove equivalent
transcript filtering, assignment or private preprocessing. No measurement
equivalence or biological accuracy conclusion follows from this source review;
no extra input requirement is introduced.

The complete **42,163 genuine native gene-ID/estimated-count axis** remains
required. Native structural-zero estimates are retained as algorithmic
unobservability under the declared quantifier, not biological absence. Missing
native rows cannot be fabricated, inferred from training zeros, or removed
from the modeled denominator.

The queued evaluation driver requires all 17 completed H3 donors and RNA
prediction hashes bound before H3 processing (lines 80–120). It enforces the
frozen donor/region/form/raw-transport axes and equality of the training-mean
baseline (lines 157–176). Target construction uses the fixed B99 helper
(lines 141–155), not a receiving fit. Its protocol explicitly describes
approximate RNA and ChIP-versus-CUT&RUN transport and excludes assay-shift
uncertainty from the conditional donor bootstrap (lines 146–151).
No concrete normalization, denominator, model-arm, or target-transform defect
was found in these inspected paths. Numerical and scientific transport checks
were not run by this review.

One bounded storage omission remains in the evaluation source: lines 183–192
check the 200 MB saved-output allowance before writing `summary.json` and
record the size as `saved_output_bytes_before_summary`. A final inclusive
directory-size check is still needed after completion, as with the earlier
batch summaries. No storage excess is claimed and no source was changed.
