# GSE223843: proposed feature/barcode identity acquisition

Initial proposal prepared 2026-10-01, metadata and primary methods only.
Root subsequently ran acquisition job 21999459; see the completed acquisition
and offline qualification below. Job 21999912 completed with 17/24 files
qualified; five fully read ATAC feature axes cannot supply the unchanged full
306,706-interval target. The seven unsupported files remain unqualified, not
biologically absent. This is a
six-person pediatric, cross-etiology source proposal, not adult MASLD validation
or model/evaluation admission. Person independence and historical exposure remain
unverified. No clinical fields, count matrices, raw reads or scores were read.

## Exact source and allowed members

The saved identity record is
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_gse223843_identity_metadata/identity_metadata.json`,
SHA256 `172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018`.
The [official file inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/filelist.txt)
was rechecked: 2,783 bytes, SHA256
`d0e3fe3012e57aa3d0364378d90357e400dbdfb6346a2ebd892d6d13a74c295a`.
Its header is `#Archive/File, Name, Time, Size, Type` with tabs.

The [official supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/)
lists only `GSE223843_RAW.tar` and `filelist.txt`. The verified downloadable URL is
[GSE223843_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/GSE223843_RAW.tar),
declared 1,066,270,720 bytes. The `File` rows name members of that archive;
individual file download URLs are **not verified** and must not be fabricated.
The exact archive URL, member names, roles and sizes are also saved in
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_gse223843_geometry_acquisition_proposal.json`,
11,465 bytes, SHA256
`35107ad4b17a7463bf28831191d7327ea895e28776b652f85d58ec20360cd7e8`.

| Author person | RNA accession / filename stem | RNA features bytes | RNA barcodes bytes | ATAC accession / filename stem | ATAC features bytes | ATAC barcodes bytes |
|---|---|---:|---:|---|---:|---:|
| Control1 | GSM6997741_GEX_Ctrl | 1,034,803 | 53,078 | GSM6997754_ATAC_Ctrl | 1,394,870 | 6,640 |
| Control3 | GSM6997744_GEX_Ctrl3 | 2,317,012 | 92,718 | GSM6997756_ATAC_Ctrl3 | 1,423,504 | 59,432 |
| Fontan1 | GSM6997746_GEX_Fontan1 | 1,453,275 | 35,835 | GSM6997758_ATAC_Fontan1 | 885,740 | 11,479 |
| Fontan2 | GSM6997748_GEX_Fontan2 | 1,565,546 | 35,054 | GSM6997760_ATAC_Fontan2 | 1,187,479 | 20,047 |
| Fontan3 | GSM6997750_GEX_Fontan3 | 1,354,139 | 60,744 | GSM6997762_ATAC_Fontan3 | 1,461,271 | 19,639 |
| Fontan4 | GSM6997752_GEX_Fontan4 | 1,097,809 | 32,729 | GSM6997764_ATAC_Fontan4 | 1,059,794 | 6,120 |

Each RNA stem expands to `_features.tsv.gz` and `_barcodes.tsv.gz`; each ATAC
stem expands to `_features.txt.gz` and `_barcodes.txt.gz`. These exact 24 names
are in the official inventory. Compressed totals are RNA features 8,822,584,
RNA barcodes 310,158, ATAC features 7,412,658, ATAC barcodes 123,357 bytes;
aggregate **16,668,757 bytes**. There is no separately named peaks BED or
fragments file in this inventory. Availability elsewhere remains unverified.

## Assay and coordinate meaning

The [primary methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11103255/),
sections “snRNA-ATAC-seq and initial data processing,” “snATAC-seq data
integration” and “snRNA-ATAC-seq data integration” (rendered lines 340–342,
349–350 and 353), describe Cell Ranger ARC 2.0.0 count processing against
GRCh38, followed by scATAC-pro 1.4.2 using FASTQ/fragments. Peaks/cells were
called per sample with defaults, yielding a nucleus-by-peak count matrix.
Integration uses unions of called peaks separately within the Control and
Fontan groups. Per-sample RNA and ATAC were linked through shared cell barcodes.
The methods do not establish the precise assignment/deduplication unit of
these deposited matrix entries or a common six-person peak axis. No source
count-scale claim follows from the inventory's MTX label.

The fixed development axis is 306,706 custom merged intervals, not the native
Cell Ranger peak axis. Existing source
`scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse296875_metadata_census.py:235–266`
declares GRanges one-based closed coordinates and the corresponding BED
`[start-1,end)` representation. External feature strings must be inspected
before choosing a coordinate parser; hg38/GRCh38 assembly agreement alone
does not fix off-by-one conventions. Do not infer those conventions from the
unopened TXT or TSV filenames.

If the external features contain every unchanged fixed interval in an exact
one-to-one coordinate join for each person, the identity phase can establish
coordinate coverage only. Native integer count compatibility remains a later
question. If exact fixed intervals are absent or only overlap different
called peaks, those source peak counts cannot establish the full unchanged
target counts: absent rows are not observed zeros, and overlapping peak counts
cannot generally reconstruct fragment assignment to new intervals. A fixed-axis
fragment recount would then be required; fragment acquisition is not proposed
in this phase and no deposited fragment URL is currently verified.

## Bounded proposed acquisition, pending root review

Installed lightweight parsing support is `/usr/bin/python3` 3.10.12, stdlib
`gzip` and zlib compile/runtime 1.2.11. Source compression is declared by the
`.gz` names; actual gzip magic, CRC, decoded field layout and coordinate
convention checks are not run. No package installation is needed.

1. Compute-only preparation would verify the saved identity hash, recheck the
   filelist hash and declared archive size, then request only TAR headers by
   exact HTTP byte ranges from the one verified archive URL. Require status
   206 and exact `Content-Range`/body length before reading; if the server
   ignores ranges, stop and close the response without consuming the archive.
   Range support is **unverified**. Validate TAR checksums, member types, names,
   sizes, unique names, padding, bounded PAX/GNU extensions and the archive end.
   Seek past all matrix bodies by arithmetic, never reading their payload.
2. Acquire only the 24 exactly named members after header sizes match the
   unchanged official inventory. No archive fallback, whole-archive stream,
   matrix, clinical file or raw-read route is authorized by this proposal.
   Set aggregate network cap 20,000,000 response-body bytes including at most
   2,000,000 header/extension bytes; expected allowed payload is 16,668,757.
   Bounds are proposed caps, not measured archive-header costs.
3. Stream gzip verification and identity-only parsing on compute. Proposed
   caps: 64 MiB decoded per member, 512 MiB decoded aggregate, 2,000,000 rows
   per member and 4,096 bytes per line. Stop on a bound, malformed fields,
   duplicate identities or ambiguous coordinates; preserve failure receipts.
   Do not silently truncate or fill rows. Treat RNA feature types explicitly
   if present, preserve source order, and namespace barcodes by author person.
   Compare full barcode sets per person, reporting both unmatched directions;
   matching barcode suffixes across people do not identify shared nuclei.
4. Report exact full-axis coordinate matches and missing/extra intervals by
   person against the frozen 306,706 geometry, with declared/source-convention
   uncertainty retained. Do not open source matrices or evaluate predictions.

[Estimate] One compute CPU, 2 GB memory, one hour is sufficient for this
identity-only task; no launcher or job is created here. Compressed retained
files should stay below 20 MB; a future helper should enforce a 64 MiB total
saved-output cap, emit only bounded diagnostic samples plus summary counts,
and avoid saving all decompressed text. No geometry or barcode acquisition has
occurred. This proposal does not admit a model, target measurement or cohort.

## Prepared implementation, not submitted

The bounded executable helper is
`scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_feature_barcode_ranges.py`;
launcher is `run_codex_gse223843_feature_barcode_ranges.sbatch` in the same
directory. Job label 1592 requests CPU/nslab, one CPU, 2 GB and one hour, with
actual `%j` stdout/stderr names. Root reviews and submits; acquisition and
runtime checks remain **not run**.

Reviewed preparation hashes: helper
`7f9845fd6ba4e756502b246fafd1500a55aa4b3987720e0569fc82c174224691`;
launcher `ad8c9b62aa7c984d94f142861a26732841052ee208ddc0b4efc9897da3e043bb`.
`ast.parse`, `bash -n` and exact launcher/helper hash-pin checks passed.
No functional acquisition, gzip or identity comparison was executed on login.

The launcher supplies the reviewed helper SHA256, which the helper checks
before proceeding. Frozen proposal, identity and 20,635,809-byte full fixture
geometry TSV are hash-guarded on compute; the official 2,783-byte filelist is
refetched and checked before TAR access. Fixture geometry SHA256 is
`50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3`.
The helper never opens the fixture count H5 or any molecular dataset.

TAR access requires exact 206/Content-Range and consistent Content-Length
before each body read, with redirects and content encoding refused. Regular
member names/sizes must match the whole official inventory; unselected regular
members must be its declared matrices. All regular bodies are first skipped
arithmetically during the header scan. Only the 24 allowlisted payload ranges
are subsequently requested. Local UTF-8 PAX path/size metadata is supported
with bounded safe parsing; links, sparse files, global/GNU extensions and
unsupported PAX keys stop closed. Two end blocks are verified; remaining
padding is not read. There is no whole-archive or guessed individual-URL route.

Schema checking precedes saving any row examples. Barcodes require one
canonical 10x barcode field; RNA supports narrow one-to-three identity fields
with optional symbol and explicit `Gene Expression` type; ATAC supports one
interval string or three chromosome/start/end fields. Unsupported or numeric
expression-like layouts stop and record field count without rejected field
contents. The original compressed allowlisted members remain available.
Decoded row contents are not printed to logs.

All exact barcode-set joins are namespaced by author person and report both
unmatched directions. RNA axes have original-byte and ENSG-row ordered
digests, feature-type counts and at most three verified identity examples;
no gene-version stripping or artificial rows occur. ATAC reports exact
original-string matches and separately labeled zero-based-half-open and
one-based-closed coordinate hypotheses. Neither coordinate hypothesis is
selected by overlap, and coordinate origin remains unverified. Matching
geometry is not count-unit equivalence. No all-axis zero filling or count
intersection is produced.

Root's parser review corrected single-base intervals before submission:
`start == end >= 1` remains valid under the one-based closed hypothesis;
zero-based half-open comparisons require `end > start`. The parser rejects
`start == end == 0` and reversed intervals. This preserves origin uncertainty
without discarding valid one-based geometry. AST, launcher syntax and updated
source hash-pin checks passed after this correction; acquisition remains unrun.

## Actual acquisition 21999459 and first-line classification

Root submitted job 21999459. It failed after 11 seconds at the first RNA
feature schema, whose first row has six fields. The source failure record
and TAR offsets independently establish successful exact-range acquisition
of all 24 compressed allowlisted members, with 16,690,996 total body bytes:
16,668,757 member bytes, 19,456 TAR metadata bytes and 2,783 filelist bytes.
HTTP range support is now observed, and no matrix body was read. This failure
does not mean RNA features or external ATAC coordinates are absent.

The original job and source helper remain unchanged. Original output is
`REC/codex_gse223843_feature_barcode_ranges_21999459/`;
`failure.json` SHA256
`ad5218429c3bd239a0cc46ae3e81b69b71480a901c1aa68fa87ee561fba2799a`
and `tar_identity_offsets.json` SHA256
`8e0467bcb07978e943cbbe69046fc359e4860677b61511a39bba1569aaf4f878`.

Root authorized a lightweight first-line lexical classification on login.
Only the first 61-byte decoded line of
`GSM6997741_GEX_Ctrl_features.tsv.gz` was examined; no row contents or numeric
values were printed/saved and no other feature rows were inspected. Its
1,034,803 compressed bytes and SHA256
`96c0f53ba3d3bb979ca9580f269503128c27163387b2e24d94b3f0295dcb0c38`
match the successful source range receipt. The six field patterns are:
Ensembl gene ID, gene-symbol syntax, recognized gene-expression feature-type
annotation, chromosome identity, unsigned-integer coordinate syntax and
unsigned-integer coordinate syntax. Thus the first row matches gene identity
plus coordinate annotation. The full-file schema, gzip CRC and coordinate
origin remain **unverified**; this is not a count-scale check.

Pattern-only evidence is saved in
`REC/codex_gse223843_first_feature_schema_20261001.json`. The original compressed
gene-feature files have been acquired, with their full schemas unproven.
No count matrix, clinical values, lineage assignments or model scores were read.

Proposed next step, awaiting root review: a separate offline helper should
guard original failed-job receipts and every saved member size/SHA against its
successful range receipt, use no network, and validate all rows of the strict
six-field RNA annotation schema before exposing identity examples. Coordinate
origin must remain qualified. Per-file schema failures should retain a field
count/reason and continue independent ATAC feature and barcode identity
assessment, so unsupported RNA annotation does not suppress all ATAC
observability reporting. Original source/job files and all count boundaries
remain unchanged. No further implementation or compute job is run here.

## Approved offline implementation prepared, not submitted

Root approved a distinct offline helper after the first-line classification.
`scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_offline_identity.py`
SHA256 is `2307979d63a864b080acd6275e514f18e648b2afc044b34ca500ab088e5b6906`;
`run_codex_gse223843_offline_identity.sbatch` SHA256 is
`d178fc54e2d3f1eb8cea3fca8252e7cbf13d058cafcf69bfc620a972565e9c4d`.
Original acquisition source, failed job, files and receipts are preserved.
This helper contains no network, subprocess, HDF5 or molecular-array imports.

Before opening gzip members, it pins the original failure/offset/filelist
receipts and executed acquisition source, the approved proposal and identity
metadata, the first-line classification receipt and the frozen full geometry.
All 24 source member sizes and hashes must match their exact successful range
receipts and approved per-person assay/role mapping. Source hashes are checked
again at completion. No matrix file, raw read, clinical record or model output
is an input.

RNA accepts the earlier narrow one-to-three identity-field layouts and the
strict six-field annotation layout. Six-field rows require an Ensembl gene ID,
gene-symbol syntax, `Gene Expression`, chromosome identity and two unsigned
integer coordinates in nondecreasing order. Unknown feature types are rejected
without filtering rows. Every row is checked before identity examples can be
included in the completed per-file result. Coordinate origin remains unverified.

Per-file schema, gzip-integrity and individual decoding-limit failures record
only reason/field count and mark CRC/full schema unverified, then continue
other files. No rejected values, valid-prefix examples or partial feature axis
are saved as if complete. A global aggregate-limit violation stops the job;
already written per-file progress receipts remain. Successfully read files
verify gzip CRC at EOF. Limits remain 64 MiB decoded per file, 512 MiB total,
2,000,000 rows and 4,096 bytes per line. No further rows of an unsupported
file are parsed after its failure.

The summary reports each author's exact sample-accession barcode join and both
unmatched directions. It does not infer well identities or shared nuclei across
people from barcode suffixes. ATAC comparisons report full fixed-axis matches,
missing fixed coordinates, unmatched source coordinates and source intervals
invalid under each of the separate one-based-closed and zero-based-half-open
hypotheses. One-based single-base intervals remain valid. Neither hypothesis
is selected by its match count. A missing coordinate is not a biological zero;
no count target or filtering is produced.

Launcher label 1592 requests one CPU, 2 GB, one hour on CPU/nslab with `%j`
logs and fresh `codex_gse223843_offline_identity_JOBID` output. Compressed
members are reused at their original locations; only source/metadata/results
are saved, under a 64 MiB inclusive cap. Completion means all files were
attempted; `full_identity_qualification_completed` remains false if any file
fails. Count payload, clinical/genotype/model values, lineage assignment,
person independence, model selection and evaluation admission remain closed
or unverified in explicit receipt fields.

`ast.parse`, `bash -n`, exact source/launcher hash pin and static import
inspection passed. Full offline gzip parsing, joins and geometry coverage are
**not run**. Root reviews and submits the single compute job.

## Completed offline qualification 21999912

The current result supersedes the prepared/not-run statements above. Root
submitted job 21999912; independent `sacct` inspection confirms COMPLETED,
11 seconds, exit 0:0. The original acquisition job 21999459 remains FAILED,
11 seconds, exit 1:0. Both sources and all original files are preserved.

| Job / phase | Terminal jobs | Running | Pending | Done | Errors | Memory | Downstream |
|---|---|---:|---:|---:|---|---|---|
|21999459 acquisition|1/1 ████████|0|0|0|Initial six-field RNA schema failure|Peak RSS unavailable in accounting|24 allowlisted compressed members retained; no matrix bodies read|
|21999912 offline identity|1/1 ████████|0|0|1|Seven per-file unsupported schemas retained|Peak RSS unavailable in accounting|17/24 ██████░░ files qualified; full target incompatible for five fully read ATAC axes|

Accounting prints zero batch MaxRSS and provides no usable peak-memory
measurement; it is not evidence of zero actual memory use. Offline stderr
is empty. The job requested one CPU, 2 GB and one hour.

Actual output is
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/codex_gse223843_offline_identity_21999912/`.
`summary.json` is 45,423 bytes, SHA256
`0ef0b29278da9e38512204fa3ba05a9c5311fd6d71b7aed425e6692cb046268d`.
Inclusive saved output is **127,964 bytes**, below 64 MiB. Total decoded
identity bytes are **34,131,578**, below 512 MiB. The executed helper and
launcher hashes independently match the frozen preparation hashes above.

All 24 per-file progress receipts reconcile with the summary. The source
member size/SHA receipts reconcile with the original successful range
receipts; person/assay/role names match the frozen 24-member proposal. Decoded
byte totals, qualification totals, barcode shared/unmatched arithmetic and
fixed-coordinate missing/unmatched arithmetic reconcile. This review parsed
only bounded receipt JSON and accounting. No compressed members, geometry
payload, matrix or rejected-line contents were reopened; identity row examples
were excluded from interpretation and not printed.

### Five complete ATAC feature axes do not cover the fixed full target

All five successfully parsed ATAC feature files are one-field interval-string
axes, with full gzip CRC verified. The sixth, Control3, stopped at unsupported
interval-string syntax and is not assigned a complete axis count or coverage.

| Author person | Complete source intervals | Exact original-string matches | Zero-based half-open coordinate matches / 306,706 | One-based closed coordinate matches / 306,706 |
|---|---:|---:|---:|---:|
|Control1|173,326|0|1|0|
|Control3|Unqualified|Unqualified|Unqualified|Unqualified|
|Fontan1|108,788|0|1|0|
|Fontan2|146,467|3|1|3|
|Fontan3|181,993|0|4|0|
|Fontan4|131,006|1|1|1|

Coordinate origin remains unresolved; neither hypothesis was selected by the
number of matches. Each of the five complete axes has fewer than 306,706
intervals, and neither qualified interpretation contains the unchanged full
fixed axis. Thus these deposited per-person peak axes cannot validate the
existing full-target task. Resolving the remaining Control3 or RNA schemas
cannot repair this incompatibility for the five completely inspected sources.
This is a coordinate-observability limitation, not evidence that missing
regions are biologically inaccessible. Matrix values remain closed; zero
filling, intersecting the target or converting overlapping called-peak counts
into a claimed full-axis measurement is not justified. A fragment recount
would require a separately verified and authorized fragment route, currently
not established in this inventory. No such acquisition or recount is proposed
or run here.

### Barcode linkage and retained unsupported schemas

All 12 barcode files passed the canonical one-field identity schema and full
gzip CRC. Exact joins use the explicit author sample/person namespace and
paired GEO accessions; no cross-person suffix matching or invented well
identity was used. These are barcode identities, not independent biological
replicates or RNA expression measurements.

| Author person | RNA barcodes | ATAC barcodes | Shared exact barcodes | RNA only | ATAC only |
|---|---:|---:|---:|---:|---:|
|Control1|10,686|1,346|1,319|9,367|27|
|Control3|19,506|15,275|15,186|4,320|89|
|Fontan1|6,977|2,482|2,462|4,515|20|
|Fontan2|6,844|4,674|4,605|2,239|69|
|Fontan3|12,360|4,569|4,538|7,822|31|
|Fontan4|6,359|1,218|1,217|5,142|1|

All six RNA feature files stopped on unsupported symbol syntax after
1,888,312 decoded prefix bytes **per file**, with six fields. Control3 ATAC
features stopped on unsupported interval-string syntax after 3,334,677
decoded prefix bytes, with one field. No rejected field contents are saved in
the reports, and no partial axes are reported as complete. Full-file schema
and gzip CRC remain unverified for these seven files. Unsupported annotation
syntax does not establish absent genes, inaccessible peaks or source-file
biological invalidity. No parser widening, rejected-row classification or
retry is needed to establish the five-source full-target incompatibility and
none was performed.

The summary correctly distinguishes completed inspection of all 24 files from
full identity qualification: the former is true, the latter false. Source
compressed gene-feature files were acquired but their full schemas remain
unproven. Count matrix/payload, clinical/genotype/model values, lineage
assignment, fits, model selection, network requests in the offline job,
zero-fill/intersection targets and evaluation admission are all explicitly
false. Person independence, historical exposure and cross-source disjointness
remain unverified. This is a pediatric source-observability result, not a
successful external model validation. No further jobs or source changes were
made in this reconciliation.

Completion receipts explicitly state count payload, clinical values, lineage
assignment, model scores, person independence and evaluation admission were
not inspected or established. Fresh output uses actual SLURM ID. Compressed
partial/error phase files are preserved; failure never writes a preexisting
directory. The 64 MiB saved-output cap includes sources, compressed members,
receipts and final summary, with reserved failure-receipt space during member
writes. Static AST and launcher syntax checks were performed on login before
submission. HTTP ranges ran in job 21999459; offline job 21999912 attempted
all 24 files, with 17 fully CRC/schema qualified and seven still unverified.
