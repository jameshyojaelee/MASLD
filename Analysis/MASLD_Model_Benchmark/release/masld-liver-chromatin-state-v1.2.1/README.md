# masld-liver-chromatin-state-v1.2.1

Predicts, from a bulk liver RNA counts matrix, ten directions of liver H3K27ac covariance learned in
one paired cohort (GSE267145, 99 participants) plus one **fixed steatosis-associated chromatin axis**,
and (via `--profile`) a per-region residualised H3K27ac profile across 96,460 regions.

## What changed from v1.2 (N07 fix)

- **Form-specific `prof_reliable_<form>` reliability masks (`weights/chromatin_state_v1_2_1.npz`)**:
  - In v1.2, `score.py:238` always returned `W["prof_reliable"]` (the 4,847-region `states_cis` mask)
    regardless of `--profile-form` (`rrr`, `rrr_cis`, `rrr_offset_cis`, `rrr_offset_states_cis`, or `states_cis`).
  - v1.2.1 stores `prof_reliable_<form>` boolean masks computed from out-of-fold per-fold skill
    (`skill_f > 0.2` across all 5 outer folds in `chromatin-stable-rrr-forms-20260908T192024Z/out/stable_oof.npz`
    and `chromatin-v12-exact-shipped-recipe-20260908T163444Z/out/fixed_config_oof.npz`) for ALL 5 forms,
    and `score.py` returns `np.asarray(W[f"prof_reliable_{form}"], bool)`.
  - In particular, the default `--profile-form rrr_offset_cis` returns its own **4,364**-region reliable mask
    while `--profile-form states_cis` returns its **4,847**-region mask.
- **Every v1.1 state head, steatosis head, and profile prediction matrix is bitwise identical to v1.2.**
- Preserves [`../masld-liver-chromatin-state-v1.2/`](../masld-liver-chromatin-state-v1.2/) untouched on disk.

```bash
python score.py --counts my_counts.tsv --out states.tsv --profile profile.npz --profile-form rrr_offset_cis
bash test_release.sh
```
