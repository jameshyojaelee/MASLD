"""Family-size sensitivity for the 7B MPRA rows, from this review's own bootstrap draws (seed 20260825, B=10,000)."""
import sys, numpy as np, pandas as pd
d = np.load(sys.argv[1] + "/mpra_recomputed_draws.npz")
ref = d["available_simple_controls::allele_identity_ridge"]
seq = [k for k in d.files if not k.startswith("available_simple_controls")]
rows = []
for k in seq:
    g = d[k] - ref
    r = {"candidate": k}
    for m in (1, 7, 11):
        r[f"lower_two_sided_family{m}"] = np.quantile(g, 0.05 / m / 2)
    # Holm on bootstrap p (two-sided, proportion of draws on the far side of zero, doubled)
    r["p_boot_two_sided"] = min(1.0, 2 * min((g <= 0).mean(), (g >= 0).mean()))
    rows.append(r)
t = pd.DataFrame(rows).sort_values("p_boot_two_sided")
for fam, keys in {"7_displayed": [k for k in seq if k.endswith("delta_ridge") or k.startswith(("sei", "enformer"))],
                  "11_all_heads": seq}.items():
    sub = t[t.candidate.isin(keys)].sort_values("p_boot_two_sided").reset_index(drop=True)
    m = len(sub)
    sub[f"holm_adj_{fam}"] = np.minimum(1, np.maximum.accumulate(sub.p_boot_two_sided * (m - np.arange(m))))
    t = t.merge(sub[["candidate", f"holm_adj_{fam}"]], on="candidate", how="left")
t.to_csv(sys.argv[1] + "/mpra_family_sensitivity.tsv", sep="\t", index=False, float_format="%.5f")
print(t.to_string())
