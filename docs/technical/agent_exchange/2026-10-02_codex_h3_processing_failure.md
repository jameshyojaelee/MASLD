# Fixed H3 receiving processing: download failure and prepared recovery

2026-10-02. Codex inspected processing source, public run ownership metadata,
SLURM accounting, stderr/stdout, and processing receipts. No receiving count
values, predictions, target arrays, accuracy, or private genotype/allelic files
were opened. This agent submitted no jobs and changed no original outputs.
After reviewing the recovery, the root agent submitted **22023804** and recorded
it in `scripts/analysis/in_biopsy_allelic/rec_codex/job_index.tsv`.

## Finding

Job **21998882** stopped during HTTPS acquisition of **B36 / SRR8702094**.
`response.read(1<<20)` raised `BrokenPipeError: [Errno 32] Broken pipe` in the
original processor's `acquire()` function, before this run's alignment or
counting. [Inference] This is a transport interruption; the receipt does not
identify the remote or local cause of the broken connection.

The frozen file length is **688,988,160 bytes**. Its retained incomplete file is
**536,207,727 bytes**. No successful acquisition receipt exists for SRR8702094.
B36 has nine preceding run acquisition/filter receipts (SRR8702085–2093);
SRR8702095 was not acquired. B36 has no completed donor receipt.

Twelve donors have completed processing receipts: B1, B7, B8, B9, B10, B12,
B15, B19, B21, B22, B24, B26. The remaining fixed donors are B36, B38, B41,
B46, B47. There is no whole-run completion receipt. Missing processing does
not constitute failed or successful prediction, and does not authorize donor
omission or a different endpoint.

### Scheduler evidence

Read from `sacct` on 2026-10-02; this table is a historical diagnostic, not live
scheduler authority.

| Job | Recorded status | Processing completion | Errors / memory | Downstream |
| --- | --- | --- | --- | --- |
| 21998882 / 1597 | FAILED, exit 1:0; elapsed 13:13:23 | 12/17 ███████░░░ | One acquisition exception; batch MaxRSS 3,477,468 KiB | Whole-run completion absent |
| 21998897 / 1598 | CANCELLED by 91825; elapsed 0 | 0/1 ░░░░░░░░░░ | No execution or MaxRSS recorded | Primary fixed evaluation not run |
| 21999365 / 1584 | CANCELLED by 91825; elapsed 0 | 0/1 ░░░░░░░░░░ | No execution or MaxRSS recorded | Descriptive PCA comparison not run |

These three historical jobs have running/pending/done counts **0/0/0**, with
one failed and two cancelled. Accounting does not establish why the two jobs
were cancelled.

## Prepared recovery

New source:

- `scripts/analysis/in_biopsy_allelic/rec_codex/codex_receiving_h3_reads_recovery.py`
  SHA256 `733f2c793206fcb122250d3f39aa2c19a9705d926858ed32bd3aa5f36180f01d`.
- `scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_h3_reads_recovery.sbatch`
  SHA256 `5a87be11e043090f30ebf63d639ed7198006d8e93b2939e180a4530bea5a7857`.

The original processor remains unchanged, SHA256
`f641ea6485515b7f986f3ce279f449c5e6b85e99f676c3a5ac959fc93bf8b509`.
The recovery imports that fixed helper. Original protocol, failure, archived
source, and each successful donor receipt have explicit SHA256 pins. On compute,
it rechecks original RNA/reference prerequisites, executable hashes/versions,
run ownership, every alignment/filter/count command, and all twelve successful
count-table and count-summary hashes. Count payloads are read only as byte streams
for hashing; they are **not parsed or copied**. The original run's raw files,
BAMs, count tables, metadata, and partial download remain in place.

The new job creates a separate output directory and reuses the twelve completed
donor measurements at their original paths. It repeats all B36 runs and processes
the four remaining donors. It preserves the exact 17-person roster, 202 public
read files, fixed 96,460-region SAF, reference, alignment seed/settings,
primary-read/MAPQ filter, unstranded counting, and frozen RNA predictions.
No model fit, target assembly, accuracy calculation, calibration, or selection
occurs in this recovery.

Acquisition uses up to three attempts at the same full URL, with 10/20-second
delays after transient connection failures. Each failed attempt's new partial is
retained under a distinct suffix. There is no byte-range continuation, alternate
URL, redirect admission, checksum relaxation, or automatic retry of an HTTP
status, identity, checksum, CRC, or storage refusal. A file must still pass its
original length, MD5 and complete gzip CRC checks before alignment.

Original retained files currently occupy **123,141,051,877 bytes**, measured from
file metadata. Periodic global cap checks dynamically add all original retained
bytes to all new bytes; per-donor checks likewise combine original and new files
for B36. The original 230GB sampled stop and 250GB download allowance are
unchanged. As in the original helper, periodic checks are not filesystem quotas;
growth between checks remains possible. Failure preserves all partials and does
not create a whole-run successful receipt.

### Downstream metadata contract

The recovery's successful `summary.json` contains the exact ordered 17 donors,
`all17_reads_processed=true`, `biological_n=17`, and:

- `donor_receipt_paths`: twelve original `...21998882/<donor>/summary.json`
  paths and five new `...recovery_<jobid>/<donor>/summary.json` paths.
- `donor_summary_sha256`, `protocol_sha256`, original/recovery source hashes,
  original/new/combined saved bytes, and `reused_donors`.
- `count_payloads_parsed_or_copied=false`,
  `predictive_accuracy_evaluated=false`, `predictions_refitted=false`.

The protocol binds current prerequisites and the recovery source/launcher under
`input_sha256`; it separately binds original receipts/count bytes under
`reused_input_sha256` and records the original directory and retry policy.

The evaluator and descriptive PCA supplement now have **metadata-only input
binding changes**, prepared for root review. No scientific calculations changed.
The primary evaluator requires the exact recovery job **22023804** and reviewed
source SHA via CLI. It loads the pinned recovery definitions, requires a complete
17-person receipt with no recovery failure receipt, validates the exact original
twelve/new five receipt origins, verifies the exact 39 original metadata/count
artifact bindings, and includes both hash maps in checks before and after the
existing calculation. Count-command, column ownership and 17-person/96,460-region
checks remain. Primary source/launcher pins:

- `codex_receiving_profile_evaluation.py`:
  `97c7150fe7e5ede6cd2ec1bd0f14a20025826e62a1ea6d0f9505eafc9865bf1e`.
- `run_codex_receiving_profile_evaluation.sbatch`:
  `ff5bb7d55de7aff3508b3e5f37799bf1ff2993b3bdb61fdbdbd0b7e1f8e21142`.

The PCA supplement requires the replacement primary job and reviewed primary
source SHA via CLI. It preserves and checks the original descriptive declaration
on disk, SHA256 `bf305fdfb60c8383787b7041193b41ff950ba351332d0f7093c25efb7f458b46`.
Only its cancelled primary job/path/source hash are superseded in memory, recorded
in a separate `primary_input_binding.json`. The supplement requires the exact
recovery-origin metadata in the successful primary protocol and verifies all
primary output hashes, then calculates the same two descriptive contrasts. Its
obsolete dependency on 21998897 was removed from the launcher; the root must
supply an `afterok` dependency on the actual replacement primary job. Neither
downstream job was submitted by this agent. After reviewing the patches, the
root submitted primary **22024055** and supplement **22024060**. Final supplement
source/launcher pins:

- `residual_pca_supplement.py`:
  `d67148b6bebff84d7b77a088858977dac5328fc1a42a200dd7692e8afa2dbf14`.
- `run_residual_pca_supplement.sbatch`:
  `8a99e088a7a42282a9ff687c2c2e2717ff18606049894059acc8f6ab8fcbdb91`.

## Submission commands

One job; account/QOS nslab, partition cpu, **8 CPUs, 64GB, 90h**. No array.
Header fields are in the new launcher, with separate `1597_h3_recovery_%j` logs.

```bash
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
sbatch scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_h3_reads_recovery.sbatch 733f2c793206fcb122250d3f39aa2c19a9705d926858ed32bd3aa5f36180f01d
```

The root submitted the command above as **22023804** after review. After any
source change, recompute/review its expected SHA; that submitted source is fixed.

Root-submitted primary **22024055**: **one cpu/nslab job, 1 CPU, 4GB, 1h**, after successful
completion of recovery. The fixed completion and hash checks run before numeric
deserialization:

```bash
sbatch --dependency=afterok:22023804 scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_profile_evaluation.sbatch 22023804 733f2c793206fcb122250d3f39aa2c19a9705d926858ed32bd3aa5f36180f01d
```

Root-submitted PCA supplement **22024060**: **one cpu/nslab job, 1 CPU, 2GB, 1h**,
after successful replacement primary **22024055**:

```bash
sbatch --dependency=afterok:22024055 scripts/analysis/in_biopsy_allelic/rec_codex/run_residual_pca_supplement.sbatch 22024055 97c7150fe7e5ede6cd2ec1bd0f14a20025826e62a1ea6d0f9505eafc9865bf1e
```

## Checks executed

- `sacct -j 21998882,21998897,21999365 --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem,AllocCPUS,NodeList -P`: the states above.
- Python AST parsing of the new source: passed. Source/launcher and original
  small metadata SHA256 checks executed; original source remained unchanged.
- `bash -n scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_h3_reads_recovery.sbatch`: passed.
- AST parsing of the patched primary/PCA sources and `bash -n` for their launchers:
  passed. `git diff --check` for all affected scripts: passed; explicit new-file
  whitespace inspection also passed. The original processing source and original
  PCA declaration hashes remained unchanged.

Numerical tests, live network retry, compute-side reused count byte verification,
recovery execution, and receiving evaluation: **not run by this agent**. The root
started the recovery/evaluation chain; this diagnosis does not claim its completion
or accuracy. Static checks establish source syntax and the recorded failure, not
successful recovery, predictive utility, or independent participant identity.

## Independent verification of completed primary loss arithmetic: prepared

New Codex-owned source and launcher, prepared for root review without changing
the running recovery or queued primary/PCA sources:

- scripts/analysis/in_biopsy_allelic/rec_codex/codex_verify_receiving_profile_losses.py
  (20,622 bytes), SHA256
  57a2361192d1de84d08da02e074f98d0915f9d7831966878e88f7ce67b5b56a0.
- scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_verify_receiving_profile_losses.sbatch
  (1,025 bytes), SHA256
  863b276f4b13f26496d77efc6301598f4aa18028e188a9998c2a9781d7091dc3.

This verifier is restricted to successful primary job **22024055**, executed
source SHA256 97c7150fe7e5ede6cd2ec1bd0f14a20025826e62a1ea6d0f9505eafc9865bf1e,
the exact ordered 17 donors, and all 96,460 fixed regions. It requires a
root-supplied SHA256 of the completed primary summary.json, rather than an
invented future digest. Failure receipts in either primary or recovery prevent
admission. The completion scope, executed source, primary protocol, saved target
and result hashes, all 34 frozen prediction hashes, native fixed baseline, region
axis, and original/recovery pre-H3 prediction bindings are checked before array
values are opened. It verifies only hashes of consumed permitted files:
**raw receiving counts and descriptor outputs are neither opened nor hashed**.
Their presence in the complete primary receipt is not independent verification
of their values or construction.

## Second failure: B36 retained-file allowance, 2026-10-02

SLURM reports recovery **22023804 FAILED, exit1:0, elapsed03:42:04**;
batch MaxRSS is3,475,944KiB. The failure receipt records B36 alignment,
12 completed donors, and preserved partials. Its SHA256 is
`8773a9f5a1e41d7708883fc72e8d302e716e108d65e37eeb365fd79dee5e1ab4`;
the recovery protocol SHA256 is
`34e6c6cd524986d9965f1beb870b0fd8a54404837329dfa988f5fabafdd95e56`.
SRR8702094 and2095 acquisition receipts passed their original size/MD5/CRC
requirements. The final2095 alignment is partial; B36 has no donor completion
receipt. This is a storage stop, not another documented download failure.

The metadata-only census, independently checked by a parallel agent, is:

| Scope | Retained bytes | Sampled limit | Difference |
| --- | ---: | ---: | ---: |
| Original21998882 plus recovery22023804 |135,722,050,761|230,000,000,000|94,277,949,239 available|
| B36 across those two directories |22,047,696,522|22,047,488,435|208,087 over|

B36's eleven H3-only compressed files total4,009,497,687 bytes; the original
operational allowance is five times that total plus2GB. Repeating its earlier
completed runs while retaining both attempts exhausted this donor allowance.
No unchanged-cap continuation can pass its initial check with these files
retained. Global storage remains below its sampled stop.

Primary22024055 is pending with DependencyNeverSatisfied; descriptive PCA22024060
and independent arithmetic22024813 remain pending downstream. None has executed.
Processing remains **12/17 ███████░░░**; receiving accuracy remains unread.

Root is preparing a separate recovery with a **40,000,000,000-byte B36 retention
allowance**, preserving both failed directories and all twelve completed donor
measurements. The combined230GB sampled stop,250GB acquisition allowance,
other donor limits and broader user-authorized1.5TB campaign budget are unchanged.
This is a prospective operational storage amendment before receiving results,
not a change to the fixed cohort, target, models or statistical comparison.
[Estimate] Repeating B36 adds about13GB, putting its combined retention near35GB;
this estimate does not guarantee that all remaining processing fits.
The new helper and downstream bindings require root review before submission.
No new recovery has yet been submitted and no old job cancelled.

### Reviewed recovery submitted

Root reviewed the entire separate helper/launcher and their exact final hashes.
New source `codex_receiving_h3_reads_storage_recovery.py` SHA256 is
`b71d1f1bdd080bc157d24b2d121d558048bea6d5b835c212f99c22e1cf8811c3`;
launcher `run_codex_receiving_h3_reads_storage_recovery.sbatch` SHA256 is
`7b17e1645208cba780bb1524b0cd1f5f90289af3fdf8680ac6bbdc208717bdec`.
It reuses the pinned prior processing functions and original39 receipt/count
bindings. Failed22023804 source/protocol/failure metadata are independently
bound under input_sha256; all files in both retained roots count toward each
applicable cap. No H3 data values were parsed during preparation.

Root independently checked AST, exact donor/guard dictionary agreement and
launcher syntax, then submitted **22026732**,1597,cpu/nslab,8CPU/64G/90h:

    sbatch scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_h3_reads_storage_recovery.sbatch b71d1f1bdd080bc157d24b2d121d558048bea6d5b835c212f99c22e1cf8811c3

The actual job was appended to the owned job index under an exclusive lock.
The obsolete pending primary22024055/PCA22024060/check22024813 were first
confirmed PENDING, then cancelled with --state=PENDING; all show elapsed0.
Original live synthetic21985999 was untouched and now has11/50 P3 receipts.

[Estimate] The twelve completed donors retain113,674,251,183 bytes for
33,055,944,465 compressed raw bytes, a3.43884 saved/raw ratio. Applying that
ratio to the remaining four donors'19,718,041,770 raw bytes and adding13GB
for repeated B36 gives216,529,329,279 bytes combined, below230GB. This is an
empirical storage estimate, not a guarantee or a relaxation of either bound.
The latest B36 recovery actually retains12.58GB; the original incomplete copy
retains9.47GB. There is no observed20GB-per-attempt figure.

Claude performed a tools-disabled, read-only terminal review of supplied storage
facts and completed-track qualifications. Its saved response SHA256 is
`bcd9dd3024c4b7fbaf3f9c8891c5ba266c6dae33fc34bd462e4f31508ea8d0a7`,
is_error=false, one turn. Root retained cap-semantics, empirical storage,
exact new B36 origin and environment-binding checks. Root rejected its
restricted-genotype check suggestion, new decision thresholds/multiplicity rules,
and a dependency between this fixed transport bootstrap and the unrelated
synthetic REC calibration. No new hypothesis, endpoint or statistical rule was
added. The reviewers' operational advice is not scientific validation.

Primary binding was subsequently reviewed and queued as**22026740**,1598,
cpu/nslab,1CPU/4G/1h, afterok22026732. Its source SHA256 is
`7ae9d401789880bd8ca40cbdbe75b54d7e73f43d1008048f6f390ee6308c9895`.
The guarded source table retains the old733 processing definitions and then
archives the newb71 processor. Before numeric access it requires exact new17-donor
completion, original12/new5 receipt origins, both retained directories, failed
receipt pins and the recorded40GB amendment with unchanged global bounds.
The entire evaluate_completed_transport function AST is unchanged; only
admitted_tables and main plus source constants changed. Root read the complete
admission changes, checked AST/hash/bash syntax, submitted and indexed the job.

    sbatch --dependency=afterok:22026732 scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_receiving_profile_evaluation.sbatch 22026732 b71d1f1bdd080bc157d24b2d121d558048bea6d5b835c212f99c22e1cf8811c3

No endpoint, prediction, target transformation, inference unit, bootstrap draw,
interval or loss calculation changed. New PCA/check bindings are still under
review; receiving accuracy remains unopened.

The parallel agent completed the exact PCA and independent-loss metadata
bindings. Root inspected those changes, verified source hashes/ASTs, shell
syntax and embedded launcher admission Python, then queued both afterok22026740:

- PCA**22026751**,1584,cpu/nslab,1CPU/2G/1h. Source SHA256
  `03ea44956d0f1519ecf5f9d01d76f15fe9995a5ea0557c7da8ed4ea22b361d9b`.
- Independent arithmetic**22026752**,1596,cpu/nslab,1CPU/2G/1h. Source SHA256
  `913d74e8de358232877ed9253d23dfeb7b1c2eae5004759370a017c7e10732d1`;
  launcher SHA256
  `2927c3272c06ab8fb9693b1f417028e1213736d7ba963a21f1e9fe5934bdf8eb`.

Both were appended to the owned job index. They require the exact newly reviewed
primary source and storage processor bindings before arrays. PCA's numerical
suffix and every independent-verifier numerical function AST remain identical.
PCA's original descriptive declaration, fixed17 donors,96,460 regions,10,000
paired donor draws, seed and1e−12 arithmetic tolerances remain unchanged.

    sbatch --dependency=afterok:22026740 scripts/analysis/in_biopsy_allelic/rec_codex/run_residual_pca_supplement.sbatch 22026740 7ae9d401789880bd8ca40cbdbe75b54d7e73f43d1008048f6f390ee6308c9895
    sbatch --dependency=afterok:22026740 scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_verify_receiving_profile_losses.sbatch 913d74e8de358232877ed9253d23dfeb7b1c2eae5004759370a017c7e10732d1

Current bounded chain has2 running (synthetic21985999 and recovery22026732),
3 pending,0 newly completed, and one historical failed recovery22023804.
Recovery completion is12/17 ███████░░░, original synthetic P3 is11/50 ██░░░░░░░░,
and downstream execution is0/3 ░░░░░░░░░░. Latest live MaxRSS is3,469,844KiB
for recovery and55,224KiB for synthetic. No new error is reported; actual
recovery success and all three downstream numerical checks remain **not run**
to completion. No paper/figure/release authority changed, so the Resource scope
check was not run for these technical-only changes.

### B36 completed in the separate storage recovery

Live job22026732 produced B36/summary.json and started B38. Root inspected only
the small donor receipt: B36 identity, all11 run acquisitions and count columns,
gzip-CRC flags, exact new count path, and no-refit/no-target-evaluation flags
passed. Receipt SHA256 is
`ccb5915903afdd0e41f92573dc9090c207e849f5e0706e4ea18d2310a84f83fb`.
The new B36 outputs retain13,643,434,380 bytes; combined with both older B36
directories, retained donor output is35,691,130,902 bytes, below40GB.
No count payload, receiving target, prediction or accuracy value was opened.
This metadata check does not independently reproduce counting or validate utility.
The first inspection command compared an absolute count path to a relative
inspection path; resolving the latter corrected that inspection assertion.

Fixed cohort processing is now**13/17 ████████░░**; B38/B41/B46/B47 remain.
Live bounded jobs remain2 running and3 pending, with no new failure. Serial
synthetic P3 has13/50 ███░░░░░░░ saved receipts; three downstream evaluations
remain0/3 ░░░░░░░░░░. Latest MaxRSS is3,475,984KiB for recovery and55,228KiB
for synthetic. The whole17-donor processing summary is still absent.

### B38 completed, 2026-10-03

Live recovery22026732 subsequently produced B38/summary.json. Root checked
donor identity, all15 acquisition CRC flags and count-column donor bindings,
the exact new count path, and no-refit/no-target-evaluation flags: passed.
Receipt SHA256 is
`89a05350e7b03ed4b90ee659fec7c0c62fa39f6e016af3cb338cfecc4f4719aa`.
The receipt records12,148,698,247 saved output bytes. No count payload,
receiving target, prediction or accuracy was read; metadata does not reproduce
counting or establish model utility. Fixed processing is now14/17 ████████░░;
B41/B46/B47 remain. Synthetic21985999 has15/50 saved P3 receipts. The two
jobs remain running, with three downstream jobs pending and no new failure.

### B41 completed, 2026-10-03

Live recovery22026732 subsequently produced B41/summary.json. Root checked
donor identity, all12 acquisition CRC flags and count-column donor bindings,
the exact new count path, and no-refit/no-target-evaluation flags: passed.
Receipt SHA256 is
`6d0509c409f5d9de9fc22600d8c3053c4040f312af1f40a958dea2d0e18ee695`;
it records4,983,266,688 saved output bytes. No count payload, receiving target,
prediction or accuracy was read. Fixed processing is now15/17 █████████░;
B46/B47 remain. Synthetic21985999 has16/50 saved P3 receipts. Both jobs remain
running, with three downstream jobs pending and no new failure. Completion
metadata does not reproduce counting or establish independent model utility.

### B46 completed, 2026-10-03

Live recovery22026732 subsequently produced B46/summary.json. Root checked
donor identity, all11 acquisition CRC flags and count-column donor bindings,
the exact new count path, and no-refit/no-target-evaluation flags: passed.
Receipt SHA256 is
`9c44799d9579773c613cb01d0b92a0fac1fd81a2de3669704cfd758ae7e31651`;
it records8,887,004,352 saved output bytes. No count payload, receiving target,
prediction or accuracy was read. Fixed processing is now16/17 █████████░;
B47 remains. Synthetic21985999 has17/50 saved P3 receipts. Both jobs remain
running, with three downstream jobs pending and no new failure. Completion
metadata does not reproduce counting or establish independent model utility.

After admission it uses the saved fixed B99 residual target, the two already
frozen profiles for each donor, and the unchanged one-dimensional source
training mean. It independently sums squared signed residuals using scalar
math.fsum within 8,192-window chunks, then sums chunks and averages all windows.
There is no clipping, nonnegativity requirement, centering, renormalization,
calibration, model fitting, new input selection, or donor/region omission.
It reconstructs the four existing donor and cohort MSEs (global, local, fixed
training mean, zero), the global-minus-local donor/cohort contrast, and the
existing descriptive skills versus the fixed training mean.

The same prespecified 10,000 paired donor draws use PCG64, seed 20260930, and
int64 indices. Instead of the producer's array-indexing reduction, the verifier
counts each repeated donor's multiplicity explicitly and uses scalar summation.
It sorts the bootstrap values and manually interpolates positions
(10000-1)*0.025 and (10000-1)*0.975 for the same linear percentile interval.
No new interval, hypothesis, p-value, selection rule, or effect margin is added.
Comparison tolerances (absolute 1e-12, relative 1e-12) are prospectively declared
floating-point reconstruction tolerances. Root review tightened the relative
tolerance before any external result was available. These are not scientific
success thresholds.

All consumed SHA256, file size/time/inode/device identities and array/NPZ member
headers are recorded and rechecked after arithmetic. Object arrays and changed
or incomplete axes are refused. The verifier imports no model, decomposition,
target-transform or producer functions. Its outputs consist of a small loss
reconstruction, consumed-input receipt and completion/failure receipt; saved
outputs plus logs are checked against **4 MiB**.

Resource header: one job named **1596**, account/QOS nslab, partition cpu,
**1 CPU, 2GB, 1h**; no array. Logs use 1596_h3_verify_%j.out/.err.
[Estimate] Sequential prediction reads and a target memory map should fit within
the requested memory; peak memory and runtime have not been measured.

Root submitted job22024813 after source review, with dependency on successful
primary22024055. Its launcher now admits the exact complete17-donor primary
receipt and executed-source hash on compute, captures the summary SHA256,
then supplies that hash to the unchanged independent verifier before arrays
are read. This permits the check to run automatically when the primary finishes.
No model, comparison, donor, region, bootstrap or tolerance changed.

Executed submission:

    sbatch --dependency=afterok:22024055 scripts/analysis/in_biopsy_allelic/rec_codex/run_codex_verify_receiving_profile_losses.sbatch 57a2361192d1de84d08da02e074f98d0915f9d7831966878e88f7ce67b5b56a0

The actual job is indexed. Updated launcher SHA256 is
`80e5c73d7660ad614aaf0ca94370da5643e18142d8c6c560c0165857034f8366`.
Root checked its shell syntax and embedded admission Python AST. Numerical
execution remains pending; no receiving values were inspected during this change.

Python AST parsing, launcher bash -n, explicit new-file whitespace inspection
and focused git diff --check: passed. No NumPy import, numerical self-test,
receiving-array read, executable verification, submission, or job index append
was performed by this agent. Root reviewed the complete helper and launcher;
the recorded final source pin includes its requested relative-tolerance adjustment.

A passing result would verify **saved primary loss arithmetic only**. It would
not establish correct B99 target construction, assay equivalence/transport,
source-person disjointness, independent validation, or practical predictive
utility. Those qualifications remain explicit in the prepared completion receipt.

## Recovery completed (2026-10-03)

Job 22026732 completed all 17 donors with exit 0:0 in 06:34:09,
maximum RSS 3,479,912KiB. All previous partial files remain preserved.
Combined retained output is 216,343,027,693 bytes, below the unchanged sampled
230GB global stop. The B36-only 40GB amendment was declared before outcomes;
other caps and counting choices were unchanged. Root verified all 17 donor
completion-metadata hashes. Recovery summary SHA256 is
`55aa1a2ba242f34eda623fdcc6e8c07d6c60241646cf34230f2a82778bb395d4`.

Fixed primary 22026740, descriptive PCA 22026751 and independent saved-loss
reconstruction 22026752 then completed with exit 0:0. Reconstruction agrees
within 4.440892098500626e−16. Local predictions were worse than global in all
17 donors; both had higher donor-average error than the fixed training mean.
No independent useful transport model is established. Exact numbers and
scientific qualifications are in
[the completed fixed H3 comparison](2026-09-30_codex_paired_profile_evaluation.md#completed-fixed-h3-comparison-2026-10-03).
No receiving refit, recalibration or model selection followed these outcomes.
