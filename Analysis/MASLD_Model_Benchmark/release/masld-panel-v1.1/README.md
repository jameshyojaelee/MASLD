# masld-panel-v1.1

A 20-gene-pair log-ratio score that **orders liver biopsies by fibrosis severity**, runnable on
**targeted liver RNA-seq or qPCR wells** with abundant protein-coding transcripts and
sample-independent missing-value handling.

## What changed from v1 (N07 fixes)

1. **Label-free abundant protein-coding candidate pool (`POOL_V2` only)**:
   - Restricted candidate genes strictly on `POOL_V2` (`557` participants, five cohorts) by
     **`median CPM >= 5.0`**, **`frac_zero <= 0.05`**, **`gene_biotype == "protein_coding"`**,
     non-`chrY`, and exclusion of olfactory receptor / testis-specific gene families (such as `OR2I1`).
   - Zero external test cohorts (including `GSE193066`, `GSE268273`, `GSE213621`, `GSE276114`)
     touched candidate selection or fitting.
2. **Sample-independent missing-pair handling and direct Cq path**:
   - Evaluates pair finiteness **per sample** (`ok_ps = isfinite(lg[ra, s]) & isfinite(lg[rb, s])`),
     never coupling samples through batch-wide `np.all` drops or `nansum(2^(-Cq))` totals.
   - `is_cq=True` directly evaluates `Cq_b - Cq_a` after averaging technical replicates per gene/sample.
   - Missing pairs (above the 0.70 pair coverage contract) are centered at their frozen training
     pair mean `pair_mean` (`m_p`) via `b0_raw + (P / used_s) * sum_{ok_s} coef_p * (diff_{p, s} - m_p)`,
     guaranteeing **0.0 score drift** on co-batched samples when any other patient has `NaN`s.
3. Preserves [`../masld-panel-v1/`](../masld-panel-v1/) completely untouched on disk.

## Run

```bash
python3 score.py --counts my_counts.tsv.gz --out scores.tsv
python3 score.py --cq my_cq.tsv --out scores.tsv
bash test_release.sh
```
