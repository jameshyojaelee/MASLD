# Codex: bounded RNA/frozen-prediction completion check

Root ran the metadata check **21998896**, completed in two seconds,exit0;
both B1 arms and all six invalid-archive refusals passed. Root queued full
completion check **21998908**, afterok on21998814–21998817. It has not run.
This task does not invoke H3,
acquisition, fitting, predictive evaluation or publication. It adds no
prerequisite to the separately queued H3 owner workflow.

Both individual jobs use name1581,cpu/nslab,one CPU,2GB,one hour with unique
actual-ID logs. Executed checker SHA256 is
`2e7a815ab64e305c76d1feb475984daefe18585b1607e2e88e628c90c947b4c8`.
The metadata check did not deserialize prediction values or evaluate all17
readiness. SLURM's zero RSS reading is unusable. Root verified the completion
dependency matches all four remaining producer jobs; B1 is already terminal
successful. No producer was rerun and no new H3 prerequisite was introduced.

Files:

- `scripts/analysis/in_biopsy_allelic/rec_codex/check_native_rna_readiness.py`
- `scripts/analysis/in_biopsy_allelic/rec_codex/run_check_native_rna_readiness.sbatch`

## Actual schema evidence

The checker follows the successful B1 pilot 21998663 and actual four batch
jobs 21998814–21998817, with fixed donor source order:
`B1,B7,B8,B9,B10,B12,B15,B19,B21,B22,B24,B26,B36,B38,B41,B46,B47`.
The batches cover the remaining sixteen once, in the already launched groups.
No donor is substituted or excluded.

Observed source and receipts establish these fields and storage layouts:

| Artifact | Actual evidence used |
|---|---|
| B1 `summary.json` | Completed native quantification, modeled/native gene counts, both prediction paths/hashes/forms, protocol hash; no predictive validation |
| Batch `batch_state.json` and `summary.json` | Same state dictionary, `group_complete`, four `status` entries and each successful donor summary hash; `cohort_complete` remains false |
| Donor `summary.json` | Native/model axis metadata, both frozen arms, saved protocol hash, completion/qualification flags |
| `geneCounts.report.json` | 227,368 native transcript rows, 60,623 native genes, exact expected-count aggregation, zero-estimate qualifications, output hashes and source/environment |
| Prediction NPZ | Producer line185–187 stores `profile`, `training_mean`, `region_key`, `sample_id`, `form`, `transport` |
| Prediction dimensions | Producer line180–183 requires profile `(1,96460)`, complete 42,163 modeled genes, frozen full region identity and finite prediction values |

Source line references concern `codex_native_rna_pilot.py`, hash
`413b9617160e0f472682aecc21a629002fae0f54556e1244ffc3aac044e232c5`.
Completion mode does not invent additional producer fields or recompute
prediction values. Its new readiness flag belongs only to its own report.

## Completion mode

With no CLI flag, compute invocation requires all five actual producers to
be `COMPLETED` with `0:0` in `sacct`. It then requires both actual batch
terminal dictionaries to match, every fixed donor entry to be succeeded,
all seventeen full donor summaries, and exact native/model metadata.

It verifies recorded original/copied source, B1 recipe, mapping, scorer,
weights, binary and existing index metadata hashes; then quantifier,
aggregated counts, native-zero diagnostics and prediction byte hashes.
No expression table or numerical prediction payload is deserialized.
The producer's finite-value admission remains the evidence for numerical
finiteness; this check does not repeat it.

The handwritten NPY reader uses bounded literal headers and rejects object
or structured dtypes. Numeric profile and training-mean members are checked
only for float64 headers, dimensions and byte lengths. Only the producer's
encoded Unicode identity metadata is decoded: all 96,460 region keys must
match the frozen public feature-axis order; sample ID, raw transport and
global/local form must match the corresponding donor/arm. No NumPy import,
`np.load`, pickle, or fitted array calculation is present.

Final donor and group inventories include all files and summaries and must
fit the actual declared caps. Inventories inspect sizes only. Raw FASTQs
are not rehashed: their original acquisition receipts and file-size
metadata are checked. No whole-genome or full search-index payload rehash
is added. Existing frozen reference/index metadata qualification is retained,
without claiming a complete pre-pilot search-payload hash record.

A fresh `codex_native_rna_readiness_<jobid>/` holds the executed checker,
launcher and result. Any missing, failed, changed, malformed, differently
ordered or over-cap artifact produces an explicit failure with phase,
donor and traceback, and the job raises. Existing output is refused.
No producer file is modified. Success records only all-seventeen RNA and
frozen-prediction readiness, not predictive usefulness, biological
independence, private PISCES equivalence or independent validation.

Native zero estimates remain valid rows. There is no per-gene positivity
rule, padding, renormalization, donor omission or lowered coverage threshold.

## Bounded metadata-decoder check mode

`--check-metadata` runs on compute before the all-seventeen completion job.
It verifies the completed B1 prediction byte hashes, then checks both real
NPZ headers and sample/form/transport/full-region identity metadata without
reading numeric profile or baseline values.

It also constructs one tiny valid three-region metadata archive and six
explicit refusal cases: wrong profile shape, missing required member,
wrong region order, wrong donor ID, object dtype, and truncated identity
payload. Synthetic numeric entries contain placeholder bytes and are never
deserialized. Production calls retain the fixed 96,460-region default;
three-region calls are confined to this synthetic metadata mode.

This mode writes a fresh `codex_native_rna_metadata_check_<jobid>/` with
actual B1 metadata and refusal results. It does not evaluate all-seventeen
readiness or invoke downstream work. A failure is retained with phase and
traceback.

## Static checks and root submission commands

Executed on login: `check_native_rna_readiness.py --inspect` passed actual
small receipt/source setup, without scheduler calls, inventories, output
hashing, archive reads or numerical deserialization. Python parsed the full
helper. `bash -n` passed for the launcher. Focused whitespace checking passed.
The installed `sacct --help` confirms the used accounting flags.

Root submitted the metadata-decoder check as job **21998896 / 1581** with
`--check-metadata`. Root verified its successful result described above.
The all-seventeen completion mode is queued as21998908 and **not run**. The shared launcher
passes optional CLI flags:

```bash
sbatch /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_check_native_rna_readiness.sbatch --check-metadata
```

Header: one CPU, 2 GiB, one hour, cpu partition, nslab account/QOS; no
array. Reusable name1581 has unique `1581_%j.out/err` logs. Root may override
the name/log prefix when submitting. After the bounded check passed, root
submitted completion mode with its already authorized four-batch dependencies:

```bash
sbatch --dependency=afterok:21998814:21998815:21998816:21998817 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_check_native_rna_readiness.sbatch
```

The helper itself also checks B1 terminal status and all five accounting
states. Submission/dependency choice and job-index recording belong to root.
No live producer/shared source, existing H3 queue, other root's report or
job index was changed by this preparation.

## Frozen submission sources

Root confirmed the helper SHA256 immediately before submitting 21998896:
`2e7a815ab64e305c76d1feb475984daefe18585b1607e2e88e628c90c947b4c8`.
The launcher SHA256 is
`b6d999a6e228c827744521063215b8022fa40fac35dfca7fa8d1f3b009cd598b`.
Both sources are frozen after handoff. Root independently repeated
`--inspect` and `bash -n`; both passed. This submission does not establish
all-seventeen readiness or add a prerequisite to the independently queued H3
work. Root owns job-index recording and review of the compute result.

## Independent queued-check review, 2026-10-01 02:20 UTC

Root submitted completion job **21998908 / 1581**. Live `scontrol` showed
`PENDING / Dependency`, with precisely
`afterok:21998814:21998815:21998816:21998817`, cpu/nslab, one CPU, 2 GiB and
one hour. Its unique logs are `1581_21998908.out/err`. The checker also
requires successful terminal B1 accounting; omitting B1 from the dependency
does not omit that check (`check_native_rna_readiness.py:107`).

The bounded source/JSON review found no likely false refusal or false
admission against the actual saved producer schemas. This is source and
metadata evidence, not an all-seventeen completion result.

| Producer/check | Saved donor/check completion | Live state |
| --- | --- | --- |
| B1 21998663 | 1/1 ████████ | COMPLETED, 0:0 |
| 21998814 | 1/4 ██░░░░░░ | RUNNING |
| 21998815 | 1/4 ██░░░░░░ | RUNNING |
| 21998816 | 2/4 ████░░░░ | RUNNING |
| 21998817 | 2/4 ████░░░░ | RUNNING |
| Decoder 21998896 | 6/6 refusals ████████ | COMPLETED, 0:0 |
| Completion 21998908 | 0/17 verified ░░░░░░░░ | PENDING, four dependencies |

Snapshot: four running, one pending, two completed jobs; no saved donor
failure in the inspected batch states. Job-level MaxRSS was blank in the
bounded accounting query, so actual memory use is unverified here. The
independently queued H3 work is unchanged.

The four groups at checker lines 30–33 match every actual protocol exactly:
B7/B8/B9/B10, B12/B15/B19/B21, B22/B24/B26/B36 and B38/B41/B46/B47. With
B1 first, they give the fixed seventeen donor source order. The six saved
remaining-donor summaries present at inspection have the actual required
field names and native 227,368-transcript, 60,623-gene and 42,163-modeled-gene
metadata. No expression or prediction values were loaded.

Checker lines 385–397 require the final batch summary to equal the final
batch-state dictionary, successful status for every group donor and the saved
donor-summary hashes. Producer lines 225–228 write those two identical final
state dictionaries. Final summaries are correctly absent while these batches
remain active; the checker refuses incomplete states after dependency release.

Checker lines 471–483 count all donor and group files, including summaries,
against their saved caps. Every actual protocol matches the producer formula
of raw-file bytes plus 2,000,000,000 bytes per donor, then a 20,000,000-byte
group buffer. The four inclusive group caps are respectively 25,612,037,059,
27,765,789,824, 25,446,600,755 and 26,794,126,559 bytes. B1's saved cap is
50,000,000,000 bytes. Actual final inventories and compliance are **not run**
by this source review; they remain the queued compute check's responsibility.

The current frozen checker/launcher hashes match the submission hashes above.
The current four original recipe sources, batch source/launcher, and all four
archived small-source copies match their immutable expected hashes. Actual
batch input hashes, index-summary hash and native-index metadata agree with
B1. No new reference-payload or raw-FASTQ hashing was introduced.

The saved `codex_native_rna_metadata_check_21998896/summary.json` reports
`metadata_decoder_checks_passed=true`, both real B1 arms checked, all six
invalid cases refused, `numerical_values_deserialized=false` and
`all17_readiness_evaluated=false`. Live accounting confirms that check exited
0:0. No NPZ was opened in this independent review. Frozen sources and live
producers were not edited; no job was submitted or rerun.
