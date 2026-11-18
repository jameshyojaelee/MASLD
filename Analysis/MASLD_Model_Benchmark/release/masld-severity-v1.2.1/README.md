# masld-severity-v1.2.1

A latent fibrosis severity score for bulk liver RNA-seq, plus the activity ordering
head and label-free shift diagnostic introduced in v1.2, packaged with
[`kleiner_sensitivity.py`](kleiner_sensitivity.py) in place of the retired `kleiner.py`.

## What changed from v1.2

- **`score.py` and `weights/` are bitwise identical to v1.2** (`latent_severity`,
  `latent_activity`, shift diagnostic, input-contract bounds, and all refusal contracts
  are unchanged; out-of-cohort GSE268273 Spearman remains **0.663119**).
- **Replaced `kleiner.py` with `kleiner_sensitivity.py`**:
  - Requires explicit `--sigma-e` and `--assumption` arguments rather than treating
    serial-biopsy variability (`0.4655--0.6433`) as a qualified pure-noise range.
  - Adds a **feasibility guard (`|r_latent| <= 1.0`)**: refuses by default when
    disattenuation produces an impossible latent correlation `|r_latent| > 1.0`
    unless `--allow-infeasible-assumption` is explicitly supplied.
  - Adds a **predictor-vs-response guard**: refuses `--stage-role response` and
    `--stage-response-beta`, because dividing an OLS slope by `rho` applies only
    when stage is an error-prone **predictor**, not when stage is the **response**.
- Preserves [`../masld-severity-v1.2/`](../masld-severity-v1.2/) untouched on disk
  except for its `SUPERSEDED BY v1.2.1` README banner; does not alter candidate `v1.3`.

## Run

```bash
python score.py --counts my_counts.tsv --out scores.tsv --json report.json
python kleiner_sensitivity.py --stage-var 1.2 --sigma-e 0.4655 \
  --assumption "hypothetical pure reader noise sensitivity check" --r 0.62
bash test_release.sh
```

## Files

- `score.py` — bulk liver RNA-seq scorer (bitwise identical to v1.2)
- `kleiner_sensitivity.py` — guarded conditional sensitivity calculator
- `weights/` — frozen severity, activity, and shift-diagnostic weights (bitwise v1.2)
- `MODEL_CARD.md` — scope, performance, input contract, and measurement limits
- `expected_ood.json` — GSE268273 smoke-test target (`0.663119`)
- `test_release.sh` — end-to-end verification including `kleiner_sensitivity.py` guards
- `MANIFEST.sha256.json` — SHA-256 verification manifest for all release files
