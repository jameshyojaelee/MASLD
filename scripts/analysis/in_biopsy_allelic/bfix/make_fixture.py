#!/usr/bin/env python3
"""Model A-final B-FIX: the shared synthetic fixture for the Gate B agreement check (spec v1 section 11).

Written by the lead, independently of both likelihood implementations (REC by Codex, B-MODEL by a
later Claude agent); it computes no likelihood, estimate or p-value. Synthetic data only.

Generative model (spec v1 sections 1-4), seed 20260930:
  3 cohorts with their own het rule and omega_c; 60 individuals (20 per cohort) with S in {0,1,2};
  z columns b, d, ffpe (partly missing in C2), sex, age10 (all missing in C3), c; one relative
  cluster inside C1 and one cluster crossing C1 and C2; 8 genes with 1-2 tags each (both
  orientations, one tag with a third allele, one rare tag). Per individual x gene: depths for every
  tag, chosen tag = largest n (ties to the smallest position), genotype class from the merged prior,
  reads from the oriented class distribution, then the integer het predicate. Only included rows are
  written; the individual file carries the stage-law depth predictors over all chosen rows with
  n >= d_c. Planted: kappa -0.08, omega_S 0.02, rho_S 0.10, rho_b 0.05, rho_d -0.05, delta small,
  u_g ~ N(0, 0.08^2).
Draw files for the agreement of P1 (S* draws), P2 (cluster resamples and lambda draws, one invalid)
and P3 (uniforms for the genotype class and the count) are written beside the data.
"""
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import betabinom

SEED = 20260930
PHI = 0.0439
COHORTS = pd.DataFrame({"cohort": ["C1", "C2", "C3"], "min_dp": [10, 10, 10], "min_minor_reads": [1, 2, 1],
                        "min_minor_pct": [10, 15, 5], "omega_c": [0.24, 0.095, 0.08]})
GENES = pd.DataFrame({"gene_id": [f"G{i}" for i in range(1, 9)],
                      "alpha": [0.8, -0.6, 0.4, 1.2, -0.3, 0.0, 0.9, -1.0],
                      "r": [logit(0.03)] * 8,
                      "w_fib": [0.05, 0.10, 0.30, 0.02, 0.20, 0.10, 0.60, 0.05],
                      "w_mye": [0.05, 0.05, 0.10, 0.03, 0.30, 0.10, 0.10, 0.05]})
# tags: gene, pos, orientation, p_ref, p_alt, p_x
TAGS = [("G1", 1000, 1, 0.55, 0.45, 0.0), ("G1", 1500, -1, 0.40, 0.60, 0.0),
        ("G2", 2000, -1, 0.70, 0.30, 0.0), ("G3", 3000, 1, 0.50, 0.45, 0.05),
        ("G4", 4000, 1, 0.60, 0.40, 0.0), ("G4", 4200, 1, 0.65, 0.35, 0.0),
        ("G5", 5000, -1, 0.97, 0.03, 0.0), ("G6", 6000, 1, 0.50, 0.50, 0.0),
        ("G7", 7000, -1, 0.35, 0.65, 0.0), ("G8", 8000, 1, 0.52, 0.48, 0.0),
        ("G8", 8300, -1, 0.50, 0.50, 0.0)]
Z_COLS = ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10", "z_c"]
TRUTH = {"kappa": -0.08, "omega_S": 0.02, "rho_S": 0.10, "rho_b": 0.05, "rho_d": -0.05, "tau": 0.08,
         "delta": {"z_b_log_e": 0.02, "z_d_dup": -0.10, "z_ffpe": 0.0, "z_sex_female": 0.03, "z_age10": 0.01, "z_c": -0.5}}


def het_lo(n, k, F):
    return max(k, -(-F * n // 100))


def main(out):
    out = Path(out)
    (out / "draws").mkdir(parents=True, exist_ok=True)
    (out / "params").mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)

    ind = []
    for c in COHORTS["cohort"]:
        for j in range(20):
            ind.append({"individual_id": f"{c}_I{j:02d}", "cohort": c, "S": int(j % 3)})
    ind = pd.DataFrame(ind)
    n_i = len(ind)
    ind["e_i"] = np.clip(np.exp(rng.normal(np.log(0.002), 0.5, n_i)), 1e-4, None)
    ind["z_b_log_e"] = np.log(ind["e_i"])
    ind["z_d_dup"] = rng.uniform(0.05, 0.4, n_i)
    ind["z_ffpe"] = rng.normal(0.0, 0.003, n_i)
    ind.loc[(ind["cohort"] == "C2") & (ind.index % 4 == 0), "z_ffpe"] = np.nan
    ind["z_sex_female"] = rng.integers(0, 2, n_i).astype(float)
    ind["z_age10"] = rng.normal(5.2, 1.1, n_i)
    ind.loc[ind["cohort"] == "C3", "z_age10"] = np.nan
    ind["z_c"] = rng.normal(0.003, 0.002, n_i)
    ind["cluster_id"] = ind["individual_id"]
    ind.loc[ind["individual_id"].isin(["C1_I01", "C1_I02"]), "cluster_id"] = "C1_I01"
    ind.loc[ind["individual_id"].isin(["C1_I10", "C2_I10"]), "cluster_id"] = "C1_I10"

    # centred covariates for generation (spec 1.1: within cohort over all individuals; mean-filled)
    def centred(col):
        v = ind[col].copy()
        v = v.fillna(ind.groupby("cohort")[col].transform("mean"))
        return (v - v.groupby(ind["cohort"]).transform("mean")).fillna(0.0)
    St = centred("S")
    Zt = {c: centred(c) for c in Z_COLS}
    dz = sum(TRUTH["delta"][c] * Zt[c] for c in Z_COLS)
    u = rng.normal(0, TRUTH["tau"], len(GENES))

    tags = pd.DataFrame(TAGS, columns=["gene_id", "pos", "orientation", "p_ref", "p_alt", "p_x"])
    tags["tag_id"] = [f"{g}_T{p}" for g, p in zip(tags["gene_id"], tags["pos"])]
    tags["s_t"] = np.where(tags["orientation"] == -1, 1, -1)
    rules = COHORTS.set_index("cohort")

    rows, cand = [], {i: [] for i in ind["individual_id"]}
    for ii, r_i in ind.iterrows():
        rule = rules.loc[r_i["cohort"]]
        for gi, g in GENES.iterrows():
            gt = tags[tags["gene_id"] == g["gene_id"]]
            depth = np.maximum(0, np.round(np.exp(rng.normal(np.log(35), 1.1, len(gt))))).astype(int)
            if rng.random() < 0.05:
                depth[0] = int(rng.integers(1500, 3001))
            order = np.lexsort((gt["pos"].to_numpy(), -depth))
            t, n = gt.iloc[order[0]], int(depth[order[0]])
            if n >= rule["min_dp"]:
                cand[r_i["individual_id"]].append(n)
            if n == 0:
                continue
            w = np.array([2 * t["p_ref"] * t["p_alt"], t["p_ref"] ** 2 + 2 * t["p_ref"] * t["p_x"],
                          t["p_alt"] ** 2 + 2 * t["p_alt"] * t["p_x"]])
            cls = rng.choice(3, p=w / w.sum())  # 0 H, 1 R (hom tag-REF), 2 A (hom tag-ALT)
            if cls == 0:
                eta = g["alpha"] + (TRUTH["kappa"] * g["alpha"] + u[gi]) * St[ii] + g["alpha"] * dz[ii]
                p = expit(eta + t["s_t"] * (rules.loc[r_i["cohort"], "omega_c"] + TRUTH["omega_S"] * St[ii]))
                rho = 1e-6 + (1 - 2e-6) * expit(g["r"] + TRUTH["rho_S"] * St[ii] + TRUTH["rho_b"] * Zt["z_b_log_e"][ii]
                                                  + TRUTH["rho_d"] * Zt["z_d_dup"][ii])
                a_or = betabinom.rvs(n, p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho, random_state=rng)
                alt = a_or if t["orientation"] == 1 else n - a_or
            else:
                e = r_i["e_i"]
                other = betabinom.rvs(n, e * (1 - PHI) / PHI, (1 - e) * (1 - PHI) / PHI, random_state=rng)
                alt = other if cls == 1 else n - other
            ref = n - alt
            m = min(ref, alt)
            if n < rule["min_dp"] or m < het_lo(n, rule["min_minor_reads"], rule["min_minor_pct"]):
                continue
            a = alt if t["orientation"] == 1 else ref
            rows.append({"individual_id": r_i["individual_id"], "cohort": r_i["cohort"], "gene_id": g["gene_id"],
                         "tag_id": t["tag_id"], "pos": int(t["pos"]), "orientation": int(t["orientation"]), "s_t": int(t["s_t"]),
                         "n": n, "ref_reads": int(ref), "alt_reads": int(alt), "a": int(a),
                         "p_ref": t["p_ref"], "p_alt": t["p_alt"], "p_x": t["p_x"], "e_i": r_i["e_i"], "true_class": "HRA"[cls]})
    rows = pd.DataFrame(rows).sort_values(["individual_id", "gene_id"]).reset_index(drop=True)
    rows.insert(0, "row_id", np.arange(len(rows)))
    ind["stage_law_log1p_candidate_rows"] = [np.log1p(len(cand[i])) for i in ind["individual_id"]]
    ind["stage_law_log1p_median_n"] = [np.log1p(np.median(cand[i])) if cand[i] else 0.0 for i in ind["individual_id"]]

    truth_class = rows.pop("true_class")
    rows.to_csv(out / "rows.tsv", sep="\t", index=False)
    pd.DataFrame({"row_id": rows["row_id"], "true_class": truth_class}).to_csv(out / "rows_true_class.tsv", sep="\t", index=False)
    ind.to_csv(out / "individuals.tsv", sep="\t", index=False, na_rep="NA")
    COHORTS.to_csv(out / "cohorts.tsv", sep="\t", index=False)
    GENES[["gene_id", "w_fib", "w_mye"]].to_csv(out / "genes.tsv", sep="\t", index=False)
    tags.to_csv(out / "tags.tsv", sep="\t", index=False)

    # design columns after spec rules: z_age10 is all-missing in C3 (0 after centring there), z_ffpe is
    # partly missing in C2 (mean fill + indicator); indicators follow the base columns in base order
    design = Z_COLS + ["z_ffpe_missing"]
    G = len(GENES)
    theta = {
        "theta1_start_like": {"alpha": [0.5] * G, "r": [float(logit(0.03))] * G, "kappa": 0.0, "delta": [0.0] * len(design),
                              "omega_S": 0.0, "rho_S": 0.0, "rho_b": 0.0, "rho_d": 0.0, "v": 0.0},
        "theta2_truth_v0": {"alpha": GENES["alpha"].tolist(), "r": GENES["r"].tolist(), "kappa": TRUTH["kappa"],
                            "delta": [TRUTH["delta"][c] for c in Z_COLS] + [0.0], "omega_S": TRUTH["omega_S"],
                            "rho_S": TRUTH["rho_S"], "rho_b": TRUTH["rho_b"], "rho_d": TRUTH["rho_d"], "v": 0.0},
        "theta3_truth_vpos": None,
        "theta4_perturbed": {"alpha": [0.7, -0.5, 0.3, 1.0, -0.2, 0.05, 0.8, -0.9], "r": [float(logit(0.05))] * G,
                             "kappa": 0.1, "delta": [0.05, -0.2, 3.0, 0.02, -0.03, -1.0, 0.04], "omega_S": -0.03,
                             "rho_S": -0.2, "rho_b": 0.1, "rho_d": 0.2, "v": 0.02 ** 2},
    }
    theta["theta3_truth_vpos"] = dict(theta["theta2_truth_v0"], v=TRUTH["tau"] ** 2)
    for k, v in theta.items():
        (out / "params" / f"{k}.json").write_text(json.dumps(v, indent=1))

    # draw files
    B = 50
    sstar = []  # within-cohort permutations; both implementations read the draws from this file
    for b in range(B):
        g1 = np.random.default_rng([SEED, 1, b])
        v = ind["S"].to_numpy().copy()
        for c in COHORTS["cohort"]:
            m = (ind["cohort"] == c).to_numpy()
            v[m] = g1.permutation(v[m])
        sstar.append(v)
    sstar = np.array(sstar)
    pd.DataFrame(sstar, columns=ind["individual_id"]).rename_axis("draw").to_csv(out / "draws" / "p1_Sstar.tsv", sep="\t")
    rs = []
    cl = ind.groupby("cluster_id").agg(cohort=("cohort", "first"), first=("individual_id", "min")).reset_index()
    cl["stratum"] = [ind.loc[ind["individual_id"] == f, "cohort"].iloc[0] for f in cl["first"]]
    for b in range(B):
        g2 = np.random.default_rng([SEED, 2, b])
        for c, sub in cl.groupby("stratum"):
            pick = g2.integers(0, len(sub), len(sub))
            rs += [{"draw": b, "stratum": c, "slot": j, "cluster_id": sub["cluster_id"].iloc[p]} for j, p in enumerate(pick)]
    pd.DataFrame(rs).to_csv(out / "draws" / "p2_resample.tsv", sep="\t", index=False)
    lam = []
    for b in range(B):
        g3 = np.random.default_rng([SEED, 4, b])
        lam.append({"draw": b, "source": "GSE193066", "lambda": 0.03 if b == 7 else float(np.clip(g3.normal(0.60, 0.08), 0.06, 1))})
        lam.append({"draw": b, "source": "PXD051911_A1", "lambda": float(np.clip(g3.normal(0.75, 0.09), 0.06, 1))})
    pd.DataFrame(lam).to_csv(out / "draws" / "lambda_draws.tsv", sep="\t", index=False)
    np.save(out / "draws" / "p3_uniforms.npy", np.random.default_rng([SEED, 3]).random((B, len(rows), 2)))

    glob = {"phi": PHI, "Delta": 0.125, "tau_max_fixture": 0.05, "B": B, "design_columns": design,
            "gene_order": GENES["gene_id"].tolist(), "seed": SEED, "truth": TRUTH, "u_planted": u.tolist(),
            "p3_uniform_use": "column 0 picks the genotype class by inverse CDF over (H, R, A) in that order with weights P(G)C_G(n); column 1 picks a by inverse CDF over K_c(n) in increasing a",
            "p2_resample_rule": "clusters are resampled within the stratum of the cohort of their member with the smallest individual_id; S and z are centred within each member's own cohort using multiplicities",
            "lambda_rule": "draw 7 of GSE193066 is invalid (<= 0.05) on purpose",
            "rows": int(len(rows)), "individuals": int(n_i)}
    (out / "globals.json").write_text(json.dumps(glob, indent=1))
    print(json.dumps({k: glob[k] for k in ["rows", "individuals", "design_columns"]}), rows.groupby("cohort").size().to_dict(),
          pd.Series(truth_class).value_counts().to_dict())


if __name__ == "__main__":
    main(sys.argv[1])
