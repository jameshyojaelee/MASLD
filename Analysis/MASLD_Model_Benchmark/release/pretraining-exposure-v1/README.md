# pretraining-exposure-v1

Answers **"was my evaluation cohort in this model's pretraining data?"** when the corpus was
published **de-identified** and the metadata cannot answer it.

The corpus this was built against ships 581,503 x 20,010 expression values with an `obs` of
exactly two columns — a row index and a train/test split — and an empty `uns`. No GSE, GSM, SRR,
PRJNA or E-MTAB anywhere. So "our cohort is not in the corpus metadata" is not a finding: there is
no namespace for it to appear in. Today essentially every benchmark of a pretrained encoder simply
assumes the answer.

This tool measures it, in expression space. Applied to the producing project it found that
**GSE268273 — the external cohort every headline number there rests on — sits at corpus rows
141384..141492**, and that **7 of 8 cohorts are in the corpus, all in its TRAIN split**.

## Run it

```bash
python3 exposure_scan.py \
    --query   query.tsv \          # rows = genes, cols = samples (--query-orient samples for the transpose)
    --cohorts cohorts.tsv \        # sample_id <TAB> cohort
    --corpus  PreBULK.h5ad \       # .h5ad/.h5, .npy, .npz or .tsv(.gz), streamed in chunks
    --out     exposure.json
```

Needs **numpy**. `h5py` only for `.h5ad` corpora, `torch` only for a GPU. The corpus is read in
chunks and never held whole — the reference corpus was 46 GB.

Nothing about the tool is organism-, tissue- or assay-specific. The only requirement is that the
query matrix and the corpus **share a gene vocabulary**.

## What comes back

A per-cohort table and a state from the shared exposure vocabulary, plus the reasoning behind every
call: injectivity, the match window, the aligned run and its exceedance count, the same-sample band,
and mutual-argmax block confirmation for anything called.

```
GSE268273   encoder_seen  injectivity 1.0, 90% of matches inside a window 0.908x the cohort
                          size, median r 0.9866, band lower bound 0.9851, aligned run 109/109
                          with 0/2000 exceedances
GSE276114   unknown       injectivity 0.7797, ... median r 0.9363
```

⛔ **A negative reads "no near-duplicate detected", never "not in the corpus."** Read
`MODEL_CARD.md` before quoting anything.

## Test it

```bash
./test_release.sh
```

Cold start, **no external data**. It builds three seeded synthetic corpora with planted ground
truth — including a merely-similar cohort that must be refused — and runs 47 checks. Three of them
flip a single threshold and require the verdict to change, so a criterion that was never consulted
gets caught. Verified to fail: disabling leg 2, disabling the window test, and truncating the
corpus stream each turn it red (4, 4 and 11 failing checks).

## Dependency you must have

The per-study state vocabulary is **imported, not restated**: `../encoder-benchmark-v1/exposure.py`,
or wherever `$EXPOSURE_VOCABULARY` points. A convention copied instead of imported is a defect that
no reading of either copy can find, so the tool refuses to start without it rather than defining a
second, silently divergent copy.

## Files

| file | what |
|---|---|
| `exposure_scan.py` | the method and the CLI |
| `corpus_io.py` | streaming readers: h5ad (dense + CSR), npy, npz, tsv |
| `make_fixture.py` | the synthetic fixture with planted ground truth |
| `test_release.sh` | 47 checks, cold start, no external data |
| `expected.json` | the reference band, the reference application's numbers, the fixture's planted truth |
| `reference/reproduce_reference_call.py` | re-derives the reference application's call from its stored matches |
| `reference/reference_call_reproduction.json` | that result: 0 disagreements over 83 comparisons |
| `env/seeds.json`, `env/pip_freeze.txt` | seeds and the environment |
| `MODEL_CARD.md` | what it detects, what it cannot, and the limits that matter |
| `MANIFEST.sha256.json` | hashes of everything above |
