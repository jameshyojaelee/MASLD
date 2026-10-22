# Codex: remaining fixed-donor native RNA calculation

**Operationally superseded; do not submit this single-job wrapper.** Another
concurrent root launched the same sixteen donors as four separate jobs
21998814–21998817 using `codex_native_rna_batch.py` and its launcher.
The prepared single-job wrapper was never submitted and must not duplicate
those acquisitions. Its preparation evidence below is retained as a record
of the independent recipe/population checks.

The B1 pilot
21998663 completed native RNA quantification and frozen global/local
predictions, with no receiving H3 measurements or predictive accuracy test.
That is one-donor feasibility, not independently validated performance.

## Frozen population and metadata

The remaining workflow retains all sixteen donors in the existing seventeen-
donor order, excluding only the already completed B1. Every donor has eight
technical runs and sixteen paired FASTQs. The frozen public roster is
`paired_h3k27ac_raw_roster_21997742/raw_files.tsv`, SHA256
`f91c00fb37387020a624a02a388ff00d440efc7fa5bdac5ef41da2a8b1ae7e78`.
Mate file numbers must be 1/2, and source biosample and experiment must match
within each run. RNA URLs and local basenames are unique; all seventeen
donors' 272 URLs are also distinct. This is file/sample-label identity,
not proof of independent biological identity.

| Remaining donor | Runs | FASTQs | Frozen compressed bytes |
|---|---:|---:|---:|
| B7 | 8 | 16 | 5,920,709,635 |
| B8 | 8 | 16 | 5,324,813,796 |
| B9 | 8 | 16 | 5,120,570,029 |
| B10 | 8 | 16 | 1,225,943,599 |
| B12 | 8 | 16 | 5,399,183,119 |
| B15 | 8 | 16 | 4,672,408,560 |
| B19 | 8 | 16 | 5,318,269,771 |
| B21 | 8 | 16 | 4,355,928,374 |
| B22 | 8 | 16 | 4,506,412,013 |
| B24 | 8 | 16 | 4,665,473,611 |
| B26 | 8 | 16 | 4,080,098,258 |
| B36 | 8 | 16 | 4,174,616,873 |
| B38 | 8 | 16 | 4,377,630,736 |
| B41 | 8 | 16 | 4,976,458,108 |
| B46 | 8 | 16 | 4,607,620,831 |
| B47 | 8 | 16 | 4,812,416,884 |
| **Total remaining** | **128** | **256** | **73,538,554,197** |

All seventeen donors total 79,800,518,482 compressed RNA bytes; B1 contributes
6,261,964,285 bytes and will not be downloaded or overwritten again.

## Executable implementation

New files confined to the Codex code directory:

- `run_remaining_native_rna.py`
- `run_remaining_native_rna.sbatch`

The wrapper first verifies the completed B1 summary, protocol and aggregation
receipt, inherited pilot/source/adapter hashes, and successful native index
summary/protocol. It does **not** invoke the B1-specific `main()` on other
donors. Each donor's run list, mates and bytes come directly from the frozen
roster.

Before imports on compute it archives exact bytes of the pilot, index helper,
fixed-artifact helper, scorer and R adapter. Python module execution preserves
the original `__file__` and repository-root conventions. Only explicitly
recorded globals change: donor per fixed iteration, shared sampled disk stop,
and paths to identical archived scorer/adapter source. The unchanged functions
perform acquisition, native index identity checks, Salmon quantification,
tximport aggregation and frozen prediction. No model coefficient, feature
threshold, target, transform, loss or statistical margin changes.

The recipe remains Salmon 1.10.3, native k31 decoy-aware index with duplicate
sequences kept and poly-A clipping disabled; paired runs pooled per donor;
library type A before mate arguments; no added bias flags. The unchanged R
adapter requires all 227,368 exact native transcript rows, retains native
expected counts with `countsFromAbundance=no`, and aggregates them by exact
versioned transcript-to-gene identity. The frozen scorer requires all 42,163
modeled genes and returns both raw `rrr` and raw `rrr_offset_cis` profiles
on all 96,460 regions. No absent rows are manufactured.

Native structural zero estimates and identical-sequence ambiguity retain
their existing qualifications. This approximate RNA recipe is not established
as the private PISCES measurement used in training; zero is not automatically
biological absence. Multithread quantification determinism is not asserted.

The wrapper hash-checks original and archived sources, fixed reference,
binary, mapping, scorer, adapter, weights and index metadata before/after
execution. It also records full current index payload hashes and checks
that payload and its file inventory at completion. These full payload hashes
describe the current qualified index; a complete payload hash record from
before the B1 pilot was not available. B1 gene-count/quantifier outputs and
both frozen predictions are checked against their saved output hashes and
reused without modifying B1.

## Resource and failure handling

One individual CPU job, no array, using the existing index:

```text
--job-name=1580 --account=nslab --qos=nslab --partition=cpu
--cpus-per-task=8 --mem=64G --time=90:00:00
```

The reusable name 1580 has distinct `1580_%j.out` and `1580_%j.err` logs.
Outputs are fresh `remaining_native_rna_<jobid>/` directories, containing
fresh per-donor directories and the executed sources/protocol. No deletion,
index rebuilding or B1 acquisition is performed.

The remaining workflow has a 250,000,000,000-byte hard saved-output cap and
230,000,000,000-byte sampled stop. Before each acquisition the complete
expected next file must fit under the sampled stop; file acquisition is
strictly size/MD5/SHA256/CRC checked. The inherited subprocess monitor
checks the entire remaining workflow every 30 seconds and terminates its
own subprocess on excess, preserving partial outputs. A subprocess can grow
between polls; completion and JSON writes therefore also require the sampled
stop, with summary bytes included. No success can be reported above the
hard cap. This is a local allocation within the approved 1.5 TB campaign.

Per-donor summaries and the progress summary are written only after native
axis admission, both predictions, source checks and disk checks pass. Final
success requires all sixteen remaining donors in fixed order plus unchanged,
hash-qualified B1 outputs. A failure preserves files, records the actual
phase/traceback when space permits, and re-raises; it does not substitute a
failed donor, missing row or guessed result.

## Checks and remaining scientific work

Executed: small-file/source `run_remaining_native_rna.py --inspect`, Python
syntax inspection and `bash -n` on its launcher. The metadata checks passed
for the complete fixed donor/run/mate/byte population. No raw receiving
reads, numerical arrays, quantification, prediction, installation or job
submission was performed by this preparation agent.

Even complete all-seventeen RNA predictions would establish only the RNA
side of the prospective comparison. Receiving H3 counting, the unchanged
B99 residual target on the full region axis, donor/assay qualifications,
and the prespecified paired-donor performance comparison remain required.
No predictive or independent-validation verdict is established here.

Historical prepared command, **not to be submitted while the four actual
batches cover these donors**:

```bash
sbatch /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_remaining_native_rna.sbatch
```

Root records the actual job ID in the code directory's job index. This agent
did not edit another worker's pilot/adapter, any live P3 diagnostic/arithmetic
replay, or the job index. No overlapping remaining-driver file was present
when these uniquely named files were created.

## Actual execution replaces this prepared wrapper

The actual batches cover B7/B8/B9/B10 (21998814), B12/B15/B19/B21
(21998815), B22/B24/B26/B36 (21998816), and B38/B41/B46/B47 (21998817).
Independent small-file/source checks confirm all sixteen remaining donors
are covered once and all four actual receipts match the completed B1 recipe's
frozen input hashes, native index metadata and installed versions. Total
remaining raw bytes are still 73,538,554,197; the four actual saved-output
caps total 105,618,554,197 bytes. No actual raw reads, estimated counts,
prediction arrays or receiving H3 measurements were opened by this review.

The actual batch implementation has its own tighter per-donor/group caps
and checks. The 230/250 GB policy described above belongs only to the
unsubmitted single-job preparation. It must not be represented as the live
batch cap policy. Details and bounded limitations are in
`2026-09-30_codex_native_rna_batch_independent_review.md`.
