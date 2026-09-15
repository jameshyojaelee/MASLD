# masld-panel-v1

A 20-gene-pair log-ratio score that **orders liver biopsies by fibrosis severity**, runnable on
**20 qPCR wells**.

## Run it

```bash
# RNA-seq counts: rows = genes (Ensembl stable IDs), columns = samples
python3 score.py --counts my_counts.tsv.gz --out scores.tsv

# qPCR: rows = the 20 panel genes, columns = samples, values = Cq
python3 score.py --cq my_cq.tsv --out scores.tsv
```

Needs numpy only. Refuses rather than guesses: below 0.70 pair coverage it exits 2 with a message
naming the missing genes.

## What you carry to deploy

A pair list, one coefficient per pair, one intercept. That is all — no quantile reference, no
per-gene mean or standard deviation, no gene-axis-length array. Above the detection floor the score
is exactly a weighted sum of **delta-Cq** values, verified to 6.7e-16 against the counts path on 109
real samples and invariant to an arbitrary global cycle offset to 1.1e-15.

## What it is

**An ordering, not a stage.** It emits no interval and no probability, it is not calibrated to the
Kleiner scale, and it is not a diagnostic. It retains **97%** of a 26,629-gene model out of cohort
(0.973 / 1.007 / 0.974 across three external cohorts).

It does **not** beat SomaSignal (0.90) or FibroScan (0.85), both of which need no biopsy. Read
`MODEL_CARD.md` section 3 before quoting any number.

## Files

| file | what |
|---|---|
| `score.py` | the scorer, counts or Cq |
| `weights/masld_panel_v1.npz` | the frozen artifact: pairs, coefficients, intercept, floor |
| `weights/panel_pairs.tsv` | the same, human-readable, with gene symbols |
| `weights/POOL_V2_tsp_pairs_k20_rank_within_panel.npz` | a POST-HOC alternative, not promoted |
| `MODEL_CARD.md` | performance, limits, refusal contract, provenance |
| `expected.json` | the number `test_release.sh` asserts |
| `test_release.sh` | end-to-end check a stranger can run |
| `MANIFEST.sha256.json` | hashes of everything above |
