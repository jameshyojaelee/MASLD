#!/usr/bin/env python3
"""Render the four-arm aspect-axis result. Every number is recomputed from the
producing artifacts, never transcribed from a report."""
from __future__ import annotations
import json, pathlib, sys
import pandas as pd

res = pathlib.Path(sys.argv[1])
payload = json.loads((res / "aspect_axis_result.json").read_text())
frame = pd.read_csv(res / "directions.tsv", sep="\t")
pd.set_option("display.width", 200)

print("== per-arm ==")
for a, s in payload["per_arm"].items():
    print(f"{a} {s['cohort']:<11} n={s['n']:<4} family={s['family_size']:<6} "
          f"{s['three_aspect_outcome']:<15} axes={s['named_aspect_axes_supported']} "
          f"{s['activity_versus_fibrosis']:<27} "
          + "  ".join(f"{k.split('_')[0]}={v}" for k, v in s["aspect_status"].items()))
    c = s["controls"]
    print(f"    controls: sex={None if not c['positive_sex'] else c['positive_sex']['verdict']}"
          f"  N1={None if not c['negative_N1_duplicate'] else c['negative_N1_duplicate']['verdict']}"
          f"  N2={None if not c['negative_N2_shuffled'] else (c['negative_N2_shuffled']['verdict'], c['negative_N2_shuffled']['count'])}")

print("\n== conditional (family F) and harmonised (H) directions ==")
cols = ["arm", "family", "tag", "n", "family_size", "count_bh05",
        "null_count_p95", "null_count_max", "perm_p",
        "residual_rank_variance_fraction", "max_abs_rho", "null_maxrho_p95",
        "verdict"]
sub = frame[frame.family.isin(["F", "H1", "H2", "H3", "H4", "H5", "MARGINAL",
                               "CONTROL_POSITIVE", "CONTROL_NEGATIVE_N1",
                               "CONTROL_NEGATIVE_N2"])]
print(sub[cols].to_string(index=False))

print("\n== harmonised meta ==")
for tag, m in payload["harmonised_meta"].items():
    if "combined_z" not in m:
        print(f"{tag:<44} {m.get('reason')}")
        continue
    print(f"{tag:<44} arms={'+'.join(m['arms'])} n={m['pooled_n']} "
          f"z={m['combined_z']} p={m['combined_p']:.3e} "
          f"max_arm_share={m['concentration']['largest_single_arm_z_share']}")
    for pa in m["per_arm"]:
        print(f"      {pa['arm']} {pa['cohort']:<11} n={pa['n']:<4} "
              f"count={pa['count']:<6} p={pa['perm_p']:.4f} "
              f"z_share={pa['z_share']:<7} {pa['verdict']}")
    print(f"      leave-one-out: " + "  ".join(
        f"{k}={v['combined_z']}" for k, v in m["leave_one_arm_out"].items()))

print("\n== secondary (raw-exposure) null where computed ==")
for _, r in frame.iterrows():
    s = r.get("secondary_axis_permutation")
    if isinstance(s, str) and s.strip() not in ("", "{}"):
        d = json.loads(s.replace("'", '"')) if s.startswith("{'") else None
        print(f"  {r.arm} {r.tag}: {s[:200]}")
