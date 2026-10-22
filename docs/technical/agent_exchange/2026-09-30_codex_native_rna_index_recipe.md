# Public RNA-reference qualification for the RNA–H3 candidate

Root submitted the declared RNA-native selective-alignment reference recipe as **21998561 (1589)**. **Reference qualification completed, exit0, in1:38:52.** Complete receiving RNA measurement across17 donors and independent model validation remain unresolved. No patient sequencing or expression/H3 values, allelic outcomes, protected cohorts, model coefficients or thresholds were accessed or changed by this reference-only workflow.

## Correction: native estimate completeness and biological observability

The receiving comparison requires the **complete genuine 42,163 native estimated-count and identifier axis**, unchanged matched-gene normalization and positive sample totals. It does not require every gene to be biologically detectable or have a positive estimate. Earlier stronger Codex statements were a review overreach, not a binding Claude or released-model requirement. Direct inspection of [score.py:167](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2/score.py:167) confirms identifier coverage and positive matched sample totals, plus a count-input refusal for values resembling log input. [MODEL_CARD.md:201](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2/MODEL_CARD.md:201) declares a0.60 identifier-coverage floor. The complete-axis receiving comparison is stricter about ID completeness, without adding per-gene biological detectability or positivity.

Genuine native structural-zero estimate rows for represented short transcripts can remain in that axis, labeled **algorithmically unobservable under the declared k31 recipe**. Those are quantifier-emitted estimates, not invented missing rows or evidence of biological absence. No missing native row may be appended or assigned a value from a training zero. Duplicate transcript aliases and ambiguous bases remain interpretation qualifications; they do not independently establish private PISCES equivalence. A modeled gene participating in a duplicate group may also have separately distinguishable transcripts.

The concurrent workflow's separate index21998573 completed and pilot21998663 is checking actual native quantifier/aggregation completeness, as recorded in the corrected [paired-profile report](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/docs/technical/agent_exchange/2026-09-30_codex_paired_profile_evaluation.md) and [root handoff](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/docs/technical/agent_exchange/2026-09-30_codex_parallel_execution_handoff.md). Complete genuine native estimated-count/ID coverage remains unverified until that pilot verifies all rows. These are distinct from21998561's ongoing full-reference qualification; neither workflow is replaced. Biological observability remains a separate interpretation limit, not an additional RNA-input stopping condition. Historical `complete_observed_input_admission=false` fields in this reference-only job indicate that actual quantifier input completeness was not assessed by this job; they do not impose individual-gene detectability.

## Declared candidate recipe

Use all **227,368 native versioned Ensembl98 GTF-derived transcripts** from the successful reconstruction21998470, followed by all **194 matching GRCh38.p13 primary-assembly contigs** from21998416. All194 contig names are listed as decoys. No extra reference, arbitrary genomic subset or sequence padding is introduced. Transcript and decoy names must be unique and noncolliding before construction. Existing exact reference and modeled-axis hashes are checked before sequence parsing; full FASTA hashes are checked on compute before constructing the combined reference.

The installed RNA-seq environment has Salmon **1.10.3**. Its executable SHA256, checked with a small binary read, is `9d45524e437fa95de9ccb06e00a2570487ee398895d3a6218883e240ffbf6b6f`. The actual version and index help were inspected. The selected command is recorded as an argument array:

`salmon index -t gentrome.fa -d decoys.txt -i index --keepDuplicates --keepFixedFasta --no-clip -k 31 -p 4 --tmpdir salmon_temporary`

`--keepDuplicates` preserves named transcript aliases; it cannot create sequence distinctions between genes. `--no-clip` preserves reference transcript lengths rather than clipping poly-A tails. Salmon still performs its own uppercase/non-ACGT cleanup. All flags are supported by the installed help. This selects an index recipe, not read preprocessing, library interpretation, sample grouping, quantification/tximport settings or a receiving-input adapter. It does not claim agreement with the private PISCES generator.

## Expected structural limitations, then actual qualification

Earlier metadata21998488 established **46 transcripts with length≤31**, **40 modeled genes without any transcript longer than31**, and input sequence-identical aliases affecting **602 modeled genes across native genes**. These measurements are inspected metadata; their completed processed-index equivalents remain **pending**. No omitted transcript is assigned an observed biological zero.

The executable helper checks the actual retained `ref_k31_fixed.fa`, `info.json`, `ctable.bin` reference-name prefix and `reflengths.bin` vector. Expected searchable names are exactly input references with length>31; every original transcript and194 decoys must still occur in the retained name/length metadata. The short-reference table separately reports name retention and searchable-FASTA presence. Any changed length, unexpected omission, wrong k/duplicate setting, missing decoy or unmappable modeled gene stops qualification. For ACGT-only inputs, processed sequence digests must agree exactly; non-ACGT inputs are reported as cleaned by Salmon, without claiming an unchanged sequence.

The complete **42,163-row frozen modeled gene axis** receives one qualification row each, retaining matrix/stable IDs and reporting annotated transcripts, retained names, searchable transcripts, structural short references, input cross-gene aliases, processed cross-gene aliases and transcripts containing ambiguous input bases. Every row explicitly leaves biological observability unestablished as an interpretation qualification; this is not a per-gene input-admission condition. The helper reports aliases without choosing a representative or changing the normalization denominator. Actual native estimate completeness must be checked separately after quantification and aggregation.

The name/length reader is tied to the shipped Pufferfish source and native little-endian Cereal format. It bounds counts/name sizes and checks the complete requested set and lengths. It reads only that documented metadata prefix, not contig-position allocations. A parse or identity failure produces an unqualified result, preserving all artifacts. Keeping a short transcript name does not mean that reads can align to it through k31 search; retaining duplicate transcripts does not make their individual gene abundances identifiable.

The executed helper counts non-ACGT bases separately in its per-transcript and per-decoy in-memory reference metadata. It subsequently reports the number of ambiguous-input transcripts for each modeled gene and checks input versus processed transcript alias groups. **It does not persist total replacement-base counts partitioned between transcripts and genomic decoys.** Therefore the aggregate Salmon log below cannot establish that replacements were confined to the genome. No decomposition is inferred from the size of the total. The executing helper was not amended to add an output after submission.

## Source evidence retained at execution

Salmon's [v1.10.3 fetch script](https://github.com/COMBINE-lab/salmon/blob/v1.10.3/scripts/fetchPufferfish.sh) pins Pufferfish `salmon-v1.10.3`. The corresponding [FixFasta.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/FixFasta.cpp) defines length≤k as short, writes searchable fixed references separately and retains short-reference metadata. [PufferfishIndexer.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishIndexer.cpp) names the fixed FASTA and serializes the completed index. [PufferfishBinaryGFAReader.cpp](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/src/PufferfishBinaryGFAReader.cpp) reinserts short names/lengths into reference metadata. Pinned [Cereal string](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/include/cereal/types/string.hpp) and [vector](https://github.com/COMBINE-lab/pufferfish/blob/salmon-v1.10.3/include/cereal/types/vector.hpp) serialization source supports the bounded metadata reader.

The source copies actually retained by21998561 show `PufferfishIndexer.cpp` lines299–328 transferring `num_decoys` and `keep_duplicates` from reference signatures into the completed index description; lines769–786 write that description to `info.json`. `PufferfishBinaryGFAReader.cpp` line269 identifies `refLengths` as `vector<uint32_t>`; lines336–358 serialize the reconstructed complete lengths to `reflengths.bin`, and line382 writes the names vector first in `ctable.bin`. Cereal's saved `vector.hpp` lines44–45 serialize a size followed by contiguous arithmetic values; lines67–69 serialize nonarithmetic members individually. Saved `string.hpp` lines44–45 serialize a size followed by character bytes. The helper consequently checks the little-endian name-vector prefix and uint32 length-vector payload against all requested native reference IDs. Runtime decoding and completed `info.json` checks are still pending.

The helper downloads those six small official source files only, checks exact preinspected SHA256 values and saves them with installed CLI/version outputs. It never executes those source files or installs packages. Official sources total103,139 bytes. The combined reference is streamed and hashed; no whole gentrome is loaded by the orchestration helper. One processed FASTA line can occupy approximately250MB for a chromosome; Salmon's index builder has its own larger allocations.

## Requested resources and bounds

Submitted individual job **21998561 (1589)**: cpu/nslab, **4 CPUs,128GB,90 hours**, no array. [Estimate] construction could use40–110GB RAM and20–80GB saved/index/temporary bytes; these are allocation estimates, not observed peaks or guarantees. The combined reference alone is3,524,663,123 bytes before Salmon processing, calculated from the previously measured two source sizes. Existing input references are reused without copying their compressed assets.

The local allowance is **100GB** within the authorized1.5TB campaign allowance. The helper samples only its fresh owned output tree every second during construction, stops its owned process group at90GB, applies a100GB per-file limit and refuses a successful result above100GB. The10GB margin limits overshoot risk. **This is a sampled stop rule, not a filesystem quota; a strict instantaneous aggregate cap cannot be guaranteed.** Final saved and temporary byte counts are recorded. No custom cleanup or existing-result deletion occurs. Salmon removes its own temporary files during normal construction; retained products, failure metadata and reference digests remain available.

Source sequence parsing, hashing and index construction require a Slurm compute allocation. The helper refuses execution without a job ID and four allocated CPUs, refuses output outside the owned execution root and creates a new output directory without overwriting an existing result.

## Files and checks

- [Helper](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/qualify_ensembl98_salmon_index.py)
- [Submitted launcher](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_qualify_ensembl98_salmon_index.sbatch)
- This distinct technical report; this subagent amended no other report, authority or job index.

Pre-submission checks passed: Python AST/compilation, launcher `bash -n`, installed Salmon version/help, small binary digest and six official-source metadata fetches. Root subsequently added a final aggregate-byte check including the summary, reviewed the helper and submitted the job. The actual submitted source SHA256 is `7e604b1e0a12949de0eaa60535db79f2ec907d9cedc5e89bdb1088c50356d427`; launcher SHA256 is `08f120ac9a63549d43c2f1f209bea592182386ef6b55cf0055960f3b10d522bc`. The helper's fresh-output ownership marker is defined before execution and assigned only after successful new-directory creation; exception reporting therefore cannot amend a preexisting output directory.

## Submitted execution checkpoint, 2026-10-01 00:36 UTC

Direct `squeue` and `sstat` inspection confirms the specific job remains live. At inspection it had elapsed45:59 on `ne1dc4-061`. These are current observations at this checkpoint, not a future status guarantee.

| Job | Completed qualification | State | Sampled maximum host RSS | Downstream |
|---|---|---|---|---|
|21998561 (1589)|0/1 ░░░░░░░░░░|Running,4 CPUs/128GB|9,261,432 KiB, approximately8.83GiB|Final reference tables and summary wait for index completion|

One running, zero pending, zero completed in this specific job set. `failure.json` and `summary.json` are absent, so no completed reference identity result is available. No runtime failure appears in the inspected small log. The CLI's version-server warning and short-transcript warnings are recorded; the latter are expected under the explicit k31 recipe.

The [executed receipt](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/ensembl98_salmon_index_21998561/receipt.json) proves that public sequence hashes, reference-name/length checks and streamed concatenation completed before the index subprocess was launched. It records all194 requested decoys and the exact selected command. Combined-reference SHA256 is `170f56bedf4563bec99627d2fa6d47132312851b59465947bc4a5f9d95d458b9`; decoy-list SHA256 is `bee61239ec4d9bf9f50d8280593b087d6878c43c47ddbef43e7a3969fc78846c`.

The [Salmon log](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/ensembl98_salmon_index_21998561/salmon.stderr) reports **153,901,655 non-ATCG nucleotide replacements**, **zero poly-A clips**, **227,516 cleaned references** and an internal ntHll estimate of **2,630,943,691 distinct k-mers**, selecting a **2^36** filter. The cleaned-reference count equals227,368 transcripts +194 decoys −46 short transcripts. This count agrees with the declared preprocessing rule; it does not prove the completed index's names, gene qualification or alias identities. The714 redundant-transcript warning does not count affected genes and does not establish gene-specific identifiability.

Reference parsing and index construction are now **running**. Completed processed-reference qualification, runtime binary-reader checks, final byte inventory and biological prediction evaluation are **not completed**; biological prediction evaluation has **not run**. This bounded follow-up inspected only the helper and small source/receipt/log files, without reading reference sequence again or editing the executing helper or launcher.

A passing reference build establishes declared reference/index identity only. The receiving pilot must verify the **complete genuine native estimated-count/ID axis** and positive sample totals before scoring. Private RNA-generator equivalence, donor independence and useful externally validated prediction remain separate unresolved scientific claims; per-gene biological detectability is an interpretation limit, not an extra input requirement.

## Completed qualification, superseding the running checkpoint

`sacct` reports **21998561 COMPLETED,exit0,1:38:52**, peakRSS18,956,240KiB.
Root inspected the completed small summary and empty launcher stderr.
All227,368 transcript names and194 decoy names are retained in native index
metadata:227,562 named references,227,516 searchable references and46 retained
short transcript names without a searchable sequence. The complete modeled
qualification has42,163 genes;40 have no searchable transcript under k31.
Native cross-gene sequence aliases affect602 modeled genes; processed searchable
aliases affect588. These counts do not mean every affected gene lacks a
distinguishable transcript, or that a zero estimate proves biological absence.

Final directory size, including summary, is **24,796,576,745 bytes**, below
the100GB allocation. Sampled temporary peak was24,793,579,554 bytes, with the
previously stated sampling limitation. Receipt SHA256 is
`31f2ea48091c14480aa11fbeba7c9fde8bb2ee2c481ffd83d26ef3e7771d4b00`;
inventory SHA256 is
`e3bef60db52b12f11ccea687c473013de661501a5f30bf309d173c7204578abe`.
[Completed summary](/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/ensembl98_salmon_index_21998561/summary.json).

The separate receiving pilot21998663 has now completed the genuine native
estimate/ID checks and both fixed predictions forB1 using its own index21998573;
the remaining16 donors are running in four batches21998814–21998817. This
reference job did not quantify those donors, establish private PISCES equivalence,
verify universal participant independence or evaluate receiving H3 accuracy.
