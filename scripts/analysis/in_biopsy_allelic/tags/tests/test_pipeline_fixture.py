"""End-to-end run of t5a -> t5b -> t6 on a synthetic fixture with known answers.

Fixture: 3 cohorts, 240 individuals, 40 genes with one tag each. Each gene has a true
oriented log2 allelic ratio; genes 0-7 also change it by 0.4 log2 per stage step.
Checks: main effects track the truth, GTEx sign agreement is high, the stage genes
rank at the top of the stage test, and the pipeline writes every expected file.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
PY = sys.executable
CT = ["Endothelial cells", "Hepatocytes", "Plasma cells", "T cells", "Cholangiocytes", "Fibroblasts", "Macrophages",
      "Circulating NK/NKT", "Resident NK", "Mono+mono derived cells", "Basophils", "B cells", "cDC1s", "cDC2s", "pDCs",
      "Neutrophils"]


def build(d, rng):
    n_ind, n_gene = 240, 40
    cohorts = np.array(["GSE135251", "GSE213621", "GSE126848"])[rng.integers(0, 3, n_ind)]
    ind = [f"IND{i:05d}" for i in range(n_ind)]
    run = [f"SRR{i:06d}" for i in range(n_ind)]
    S = rng.integers(0, 3, n_ind)
    cond = np.where(cohorts == "GSE213621", np.array(["Fibrosis_F0F1", "Fibrosis_F2", "Fibrosis_F3F4"])[S], "NASH")
    fib = np.where(S == 0, rng.integers(0, 2, n_ind), np.where(S == 1, 2, rng.integers(3, 5, n_ind))).astype(float)
    fib[cohorts == "GSE126848"] = np.nan
    fib[cohorts == "GSE213621"] = S[cohorts == "GSE213621"] + 1
    pd.DataFrame({"sample_id": run, "dataset": cohorts, "condition": cond, "fibrosis_stage": fib}).to_csv(d / "meta.csv", index=False)
    idt = d / "identity"; idt.mkdir()
    pd.DataFrame({"run": run, "cohort": cohorts, "individual_id": ind, "autosomal_called": 5000, "expr_sex": "F",
                  "chrX_called": 100, "chrX_het_rate": 0.2}).to_csv(idt / "sample_to_individual.tsv", sep="\t", index=False)
    pd.DataFrame({"#IID1": [run[0]], "IID2": [run[1]], "relation": ["first_degree"]}).to_csv(idt / "kinship_related_pairs.tsv", sep="\t", index=False)
    seal = d / "seal"; seal.mkdir()
    pd.DataFrame({"individual_id": ind, "home_cohort": cohorts, "sealed": False}).to_csv(seal / "sealed_individuals.tsv", sep="\t", index=False)
    pd.DataFrame({"run": ["none"]}).to_csv(seal / "excluded_possibly_mixed_libraries.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": ["x"], "is_first_biopsy": [True]}).to_csv(d / "p193.tsv", sep="\t", index=False)
    pd.DataFrame({"run": ["y"], "source_control_status": ["Control"]}).to_csv(d / "c130.tsv", sep="\t", index=False)
    hep = rng.uniform(0.6, 0.95, n_ind)
    comp = pd.DataFrame({c: (1 - hep) / (len(CT) - 1) for c in CT}); comp["Hepatocytes"] = hep
    comp.insert(0, "run", run); comp.insert(0, "cohort", cohorts)
    comp.to_csv(d / "comp.tsv", sep="\t", index=False)
    json.dump({"chosen_rule": {"min_dp": 20, "min_minor_reads": 2, "min_minor_frac": 0.02, "het_sensitivity": 0.95}},
              open(d / "t2.json", "w"))
    truth = rng.normal(0, 0.8, n_gene)
    stage_eff = np.where(np.arange(n_gene) < 8, 0.4, 0.0)
    orient = rng.choice([-1, 1], n_gene)
    tags, ad, gtf, eg = [], [], [], []
    for gi in range(n_gene):
        gid, pos = f"ENSG{gi:011d}.1", 1_000_000 + gi * 10_000
        tags.append((gid, f"G{gi}", f"chr1_{pos}_A_G_b38", "chr1", pos, "A", "G", 0.95 * orient[gi], 0.95 * orient[gi],
                     0.3, 0.3, orient[gi], False, gid.split(".")[0]))
        gtf.append(f"chr1\tX\texon\t{pos - 500}\t{pos + 500}\t.\t+\t.\tgene_id \"{gid}\"; gene_name \"G{gi}\";\n")
        eg.append((gid, 0.5 * truth[gi], truth[gi], 0.001))
        het = rng.random(n_ind) < 0.45
        depth = rng.poisson(60, n_ind) + 5
        mu = truth[gi] + stage_eff[gi] * S + rng.normal(0, 0.2, n_ind)
        p_leadalt = 2 ** mu / (1 + 2 ** mu)
        a_lead = rng.binomial(depth, p_leadalt)
        alt_reads = np.where(orient[gi] == 1, a_lead, depth - a_lead)
        alt_reads = np.where(het, alt_reads, rng.binomial(depth, 0.002))
        for i in range(n_ind):
            ad.append((cohorts[i], run[i], "chr1", pos, "A", "G", int(depth[i] - alt_reads[i]), int(alt_reads[i]), 0))
    pd.DataFrame(tags, columns=["gene_id", "gene_name", "lead_variant_id", "chrom", "pos", "ref", "alt", "r_eur", "r_eas",
                                "maf_eur", "maf_eas", "orientation", "tag_is_lead", "gene_key"]).to_csv(d / "tags.tsv", sep="\t", index=False)
    pd.DataFrame(ad, columns=["cohort", "run", "chrom", "pos", "ref", "alt", "ref_reads", "alt_reads", "other_reads"]).to_csv(
        d / "ad.tsv.gz", sep="\t", index=False, compression="gzip")
    (d / "genes.gtf").write_text("".join(gtf))
    (d / "red.tsv").write_text("chr1\t1390000\n")  # gene 39's tag: must be dropped
    pd.DataFrame(eg, columns=["gene_id", "slope", "log2_aFC", "qval"]).to_csv(d / "egenes.tsv", sep="\t", index=False)
    pb = pd.DataFrame(rng.uniform(1, 100, (n_gene, len(CT))), index=[f"G{g}" for g in range(n_gene)], columns=CT)
    pb.loc["G1", "Macrophages"] = 5000
    pb.to_csv(d / "pb.tsv", sep="\t")
    return truth


def test_pipeline():
    rng = np.random.default_rng(11)
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        truth = build(d, rng)
        subprocess.run([PY, str(HERE / "t5a_allelic_table.py"), "--tags", str(d / "tags.tsv"), "--ad", str(d / "ad.tsv.gz"),
                        "--t2", str(d / "t2.json"), "--identity-dir", str(d / "identity"), "--seal-dir", str(d / "seal"),
                        "--metadata", str(d / "meta.csv"), "--gse130970-controls", str(d / "c130.tsv"),
                        "--gse193066-placement", str(d / "p193.tsv"), "--composition", str(d / "comp.tsv"),
                        "--gtex-egenes", str(d / "egenes.tsv"), "--gtf", str(d / "genes.gtf"), "--rediportal-hits", str(d / "red.tsv"),
                        "--restricted-out", str(d / "r"), "--out", str(d / "t5a")], check=True, capture_output=True)
        me = pd.read_csv(d / "t5a" / "main_allelic_effects.tsv", sep="\t")
        s = json.loads((d / "t5a" / "t5a_summary.json").read_text())
        est = me.set_index("gene_id")["log2_allelic_ratio"]
        tr = pd.Series(truth, index=[f"ENSG{g:011d}.1" for g in range(40)])
        r = np.corrcoef(est, tr[est.index])[0, 1]
        assert r > 0.9, r
        assert s["g2"]["sign_agreement_all"] > 0.85, s["g2"]
        assert s["tag_rows_at_rediportal_sites_dropped"] == 1 and "ENSG00000000039.1" not in set(me["gene_id"])
        subprocess.run([PY, str(HERE / "t5b_stage_lineage.py"), "--participant-gene", str(d / "r" / "participant_gene_allelic.tsv.gz"),
                        "--t5a-dir", str(d / "t5a"), "--identity-dir", str(d / "identity"), "--out", str(d / "t5b")],
                       check=True, capture_output=True)
        st = pd.read_csv(d / "t5b" / "stage_tests.tsv", sep="\t")
        top = set(st.sort_values("p_perm").head(8)["gene_id"])
        hits = len(top & {f"ENSG{g:011d}.1" for g in range(8)})
        assert hits >= 6, (hits, st.head(10).to_string())
        sb = json.loads((d / "t5b" / "t5b_summary.json").read_text())
        assert "g3" in sb and "lineage" in sb
        subprocess.run([PY, str(HERE / "t6_predict.py"), "--participant-gene", str(d / "r" / "participant_gene_allelic.tsv.gz"),
                        "--t5a-dir", str(d / "t5a"), "--identity-dir", str(d / "identity"), "--pseudobulk", str(d / "pb.tsv"),
                        "--out", str(d / "t6")], check=True, capture_output=True)
        s6 = json.loads((d / "t6" / "t6_summary.json").read_text())
        assert s6["pairs"] > 20, s6
        print(json.dumps({"main_r": r, "g2": s["g2"]["sign_agreement_all"], "stage_top8_hits": hits,
                          "stage_summary": sb["stage"], "g3": sb["g3"]["stage"], "t6": {k: s6[k] for k in ("pairs", "l1_genes", "gain_L0_over_GTEx_relative_deviance", "gain_L1_over_L0_relative_deviance", "shuffle_null_95th")}}, indent=1, default=str))


if __name__ == "__main__":
    try:
        test_pipeline(); print("ok test_pipeline")
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode()[-3000:], e.stderr.decode()[-3000:]); raise
