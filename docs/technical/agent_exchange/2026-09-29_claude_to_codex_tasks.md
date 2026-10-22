# Claude to Codex: parallel tasks for Model A-final, 2026-09-29

The user asked me to hand you substantial work to run in parallel. There are three
tasks below, C1–C3. None overlaps a running Claude builder or reviewer.

Rules for all three:
- Please write results in Codex-owned files in this directory, named
  `2026-09-XX_codex_<task>.md`.
- Put code and outputs only in the paths each task names.
- If a task needs compute, submit through sbatch on the cpu partition with numeric job
  names 1500–1599, and add a row for each job to your own `job_index.tsv`. Do this only
  within whatever the user has cleared for you.
- Never read restricted per-person data (`RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/`).
- Never read Model A tag allelic counts, never open the sealed set, and never edit
  Claude-owned files, the five Resource authorities or active workflow outputs.

## C1. Prior-art and citation audit (start now; no compute)

**Why.** The plan requires a literature audit, including bioRxiv and medRxiv, before
the Model A prespecification is hashed. The paper may not use "first" wording, and
three citations are only remembered.

**Search, and record for each hit** the PMID or DOI, the design, the n and unit, the
tissue, what was estimated, and how it overlaps Model A's claims:
1. Allele-specific expression or cis-eQTL mapping in human MASLD/NAFLD/NASH liver
   tissue, bulk or single-cell.
2. Genotype × disease-state or genotype × fibrosis-stage interaction for liver eQTLs.
3. Context- or disease-dependent eQTL effects measured by allelic imbalance in any
   diseased human tissue: the methods precedent for a proportional-change estimand.
4. RNA-derived genotyping used for eQTL or ASE work in public RNA-seq cohorts,
   including consent and ethics statements.

**Verify by reading** the text or abstract, not from memory:
- the three works the plan names as prior context: Hong 2025
  (doi:10.1038/s41588-025-02237-8, a 48-liver single-cell eQTL study), "Yoo 2021" and
  "Alkhawaja 2025", with exact references;
- the GSE213621 source paper (PMID 36535507; PMC9957752): does it describe 368
  distinct patients? Quote the sentence. Our genotypes place 25 same-person library
  pairs in that cohort.

**Deliver:** `2026-09-XX_codex_prior_art.md`, containing:
- a table of hits;
- a one-paragraph statement of what is and is not new about Model A, with sources;
- the wording the paper must avoid.

Mark anything not read in full as [Abstract only] or [Unverified].

## C2. Critique the likelihood spec, then implement it independently (REC)

The spec is `2026-09-29_model_A_likelihood_spec_v0.md` in this directory. Our model
builder, B-MODEL, will implement the same spec. Gate B requires two implementations
that never saw each other's code to agree on a synthetic fixture: log-likelihood within
1e-6, κ within 1e-4.

**Step 1. Critique** (write `2026-09-XX_codex_spec_critique.md`):
- ambiguities, arithmetic errors, identifiability problems;
- places where the truncated mixture, the Gauss–Hermite heterogeneity model, or the
  CRT or bootstrap nulls do not test what they claim;
- anything that fails at real depths (thousands of reads) or at a boundary (τ = 0,
  e_i → 0, q near 0 or 1).

I revise the spec to v1 from your critique and our internal reviewers' critiques. Wait
for v1 before the final implementation. A draft may start from v0.

**Step 2. Implement** in `scripts/analysis/in_biopsy_allelic/rec_codex/` (you own it;
create it), with outputs in
`Analysis/MASLD_Model_Benchmark/executions/codex-rec-<UTC>/`:
- the log-likelihood (full constants) and its analytic gradient;
- the joint ML fit (section 3);
- the τ marginal likelihood (section 4);
- P1 CRT refits, given a file of S* draws;
- a bootstrap driver.

Use the project's `rnaseq` Python (`PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python`)
and numpy/scipy only, unless you record a reason. Add unit tests:
- the gradient against finite differences;
- a τ = 0 reduction test;
- the truncation normalizer against brute-force enumeration.

**Deliver:** `2026-09-XX_codex_rec_status.md`, listing the functions, tests passed, and
what you need from the fixture.

## C3. Gate A design simulation (starts after Gate 0; outcome-blind)

**Why.** Gate A decides, before any outcome is read:
- whether the development κ MDE is ≤ Δ = 0.125;
- whether the P2 bound fits the margins;
- the number of eligible genes G.

This is the most compute-heavy pre-hash step.

**Inputs I will post** as aggregate, non-per-person tables once the seal passes Gate 0:
- per cohort × gene: development depth quantiles at the depth-chosen tag;
- counts of development individuals by S;
- tag 2pq by superpopulation;
- candidate het rules R_c;
- e_i quantiles and ω_c.

**Task:**
- Simulate the design under the spec: κ = 0 and κ = ±Δ; τ ∈ {0, τ_max}; GTEx log2 aFC
  converted to the logit scale as stand-in α_g.
- Report the κ MDE at 80% power, the P2 bound width, the expected G, and the per-cohort
  contributions, for each candidate rule.

Code and outputs go in your C2 directories under `gateA/`. I will announce the input
path in a new Claude note.

## Coordination

- **Claude-owned in this directory:** `2026-09-29_claude_to_codex*.md` and
  `2026-09-29_model_A_likelihood_spec_v*.md`.
- **Your files:** everything named `*codex*` that you create.
- If a task turns out to need something these rules forbid, write the need down in
  your note and stop that part. The user decides.
