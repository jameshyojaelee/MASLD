"""Check 5F spatial numbers: rerun (job 9091) vs job 1731 vs release; donor-level and program-level contrasts."""
import itertools
from pathlib import Path
import numpy as np, pandas as pd
P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
C = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F5ef"
J = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/1731/out"
REL = P / "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
mine_p = pd.read_csv(C / "spatial_rerun/out/program_level.tsv", sep="\t")
j_p = pd.read_csv(J / "program_level.tsv", sep="\t")
mine_u = pd.read_csv(C / "spatial_rerun/out/unit_level.tsv", sep="\t")
j_u = pd.read_csv(J / "unit_level.tsv", sep="\t")
num = ["residual_moran_i", "null_mean", "null_sd", "excess", "k_null_ge_obs", "p", "q_bh_two_program_within_dataset"]
print("program-level max |rerun - 1731|:", {c: float(np.max(np.abs(mine_p[c] - j_p[c]))) for c in num})
print("unit-level max |rerun - 1731|:", {c: float(np.max(np.abs(mine_u[c] - j_u[c]))) for c in ["residual_moran_i", "null_mean", "excess", "p"]})
print(mine_p[["dataset", "program_name", "design", "residual_moran_i", "null_mean", "null_sd", "excess", "k_null_ge_obs", "p", "q_bh_two_program_within_dataset"]].to_string())

# release program-level from release null draws (independent count)
nd = pd.read_csv(REL / "full_composition_null_draws.tsv.gz", sep="\t")
rp = pd.read_csv(REL / "full_composition_program_results.tsv", sep="\t")
out = []
for _, r in rp.iterrows():
    d = nd[(nd.dataset == r.dataset) & (nd.program_id == r.program_id)].residual_moran_i.to_numpy()
    k = int((d >= r.residual_moran_i).sum())
    out.append(dict(dataset=r.dataset, program=r.program_name, n=len(d), obs=r.residual_moran_i, null_mean=d.mean(), null_sd=d.std(ddof=1),
                    excess=r.residual_moran_i - d.mean(), k=k, p=(k + 1) / (len(d) + 1)))
o = pd.DataFrame(out)
for ds, g in o.groupby("dataset"):
    p = g.p.to_numpy(); m = len(p); idx = np.argsort(p)
    q = np.empty(m); q[idx] = np.minimum.accumulate((p[idx] * m / np.arange(1, m + 1))[::-1])[::-1]
    o.loc[g.index, "q"] = np.minimum(q, 1)
print("\nRelease null draws, recount:\n", o.to_string())

# program-level between-program difference test using release draws (draws independent across programs; pair by draw_id)
print("\nProgram-level SECM minus DI difference against paired independent matched-gene nulls (release draws):")
for ds in ["GSE192741", "Vu_et_al_2025"]:
    a = nd[(nd.dataset == ds) & (nd.program_id == "hotspot_hepatocytes_f05c535ae5bbc0b9")].sort_values("draw_id").residual_moran_i.to_numpy()
    b = nd[(nd.dataset == ds) & (nd.program_id == "hotspot_hepatocytes_48f39dd4d817a10e")].sort_values("draw_id").residual_moran_i.to_numpy()
    oa = rp[(rp.dataset == ds) & (rp.program_id == "hotspot_hepatocytes_f05c535ae5bbc0b9")].residual_moran_i.iloc[0]
    ob = rp[(rp.dataset == ds) & (rp.program_id == "hotspot_hepatocytes_48f39dd4d817a10e")].residual_moran_i.iloc[0]
    obs_diff = (oa - a.mean()) - (ob - b.mean())
    null_diff = (a - a.mean()) - (b - b.mean())
    k = int((null_diff >= obs_diff).sum())
    print(f"  {ds}: excess difference {obs_diff:.4f}; null SD {null_diff.std(ddof=1):.4f}; k={k}; one-sided P={(k+1)/(len(null_diff)+1):.4f}")

# release per-section observed residual Moran's I (all-16) -> donor-level sign counts
ru = pd.read_csv(REL / "full_composition_unit_results.tsv", sep="\t")
ru["prog"] = ru.program_id.map({"hotspot_hepatocytes_f05c535ae5bbc0b9": "SECM", "hotspot_hepatocytes_48f39dd4d817a10e": "DI"})
print("\nRelease per-section all-16 residual Moran's I:")
print(ru[ru.dataset == "GSE192741"][["prog", "sample_id", "individual", "residual_moran_i"]].to_string())
don = ru.groupby(["dataset", "prog", "individual"]).residual_moran_i.mean().reset_index()
print("donor/array residual I > 0 counts:", don.groupby(["dataset", "prog"]).residual_moran_i.apply(lambda x: f"{(x>0).sum()}/{len(x)}").to_dict())
print(don[don.dataset == "GSE192741"].to_string())


def signflip(x):
    x = np.asarray(x, float); obs = x.sum()
    sums = [np.dot(s, np.abs(x)) for s in itertools.product([-1, 1], repeat=len(x))]
    return np.mean(np.array(sums) >= obs - 1e-12)


# donor-level excess from rerun (all16)
u = mine_u[(mine_u.design == "all16") & (mine_u.unit_type == "donor_or_array")].copy()
u["prog"] = u.program_name.str.split(" \\(").str[0]
for ds in ["GSE192741", "Vu_et_al_2025"]:
    s = u[(u.dataset == ds) & (u.prog == "Stromal ECM")].set_index("unit_id")
    d = u[(u.dataset == ds) & (u.prog == "Ductular injury")].set_index("unit_id")
    print(f"\n{ds} SECM unit excess:", s.excess.round(4).to_dict(), "pos", (s.excess > 0).sum(), "/", len(s), "signflip P", signflip(s.excess))
    print(f"{ds} SECM unit residual I:", s.residual_moran_i.round(4).to_dict(), "pos", (s.residual_moran_i > 0).sum())
    print(f"{ds} DI unit excess:", d.excess.round(4).to_dict(), "pos", (d.excess > 0).sum(), "signflip P", signflip(d.excess))
    print(f"{ds} DI unit residual I pos", (d.residual_moran_i > 0).sum(), "/", len(d))
    diff = (s.excess - d.excess.reindex(s.index))
    print(f"{ds} SECM-DI unit diff:", diff.round(4).to_dict(), "pos", (diff > 0).sum(), "signflip P", signflip(diff))
sec = mine_u[(mine_u.design == "all16") & (mine_u.unit_type == "section") & (mine_u.dataset == "GSE192741") & mine_u.program_name.str.startswith("Stromal")]
print("\nGSE192741 Stromal ECM sections all16:\n", sec[["unit_id", "individual", "residual_moran_i", "null_mean", "excess", "k_null_ge_obs", "p"]].to_string())
secl = mine_u[(mine_u.design == "lineage_only") & (mine_u.unit_type == "section") & (mine_u.dataset == "GSE192741") & mine_u.program_name.str.startswith("Stromal")]
print("GSE192741 Stromal ECM sections lineage_only:\n", secl[["unit_id", "individual", "residual_moran_i", "excess"]].to_string())
