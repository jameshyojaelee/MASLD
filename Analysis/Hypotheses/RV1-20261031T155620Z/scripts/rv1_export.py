"""RV1 step 1: export BF3 donor pseudobulk for three lineages and the Hotspot LOO match table."""
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
BUILD = ROOT / "Analysis/Hypotheses/BF3-20261023T154841Z/build-20261023T155101Z"
HS = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes"
ROSTER = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata/current_source_diagnosis_donor_roster.tsv"
OUT = Path(sys.argv[1])
LINEAGES = ["Hepatocytes", "Fibroblasts", "Cholangiocytes"]

z = np.load(BUILD / "donor_type_counts.npz", allow_pickle=True)
counts, donors, types, genes = z["counts"], list(z["donors"]), list(z["types"]), list(z["genes"])

donor_rows = list(csv.DictReader(open(BUILD / "donors.tsv"), delimiter="\t"))
assert [r["donor"] for r in donor_rows] == donors
corrected = {r["donor"]: r["disease_stage_coarse"] for r in csv.DictReader(open(ROSTER), delimiter="\t")}
with open(OUT / "donors.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["donor", "dataset", "diagnosis_bf3", "diagnosis_corrected", "sex"] + [f"nuclei_{l}" for l in LINEAGES] + [f"umi_{l}" for l in LINEAGES])
    for i, r in enumerate(donor_rows):
        umi = [int(counts[i, types.index(l)].sum()) for l in LINEAGES]
        w.writerow([r["donor"], r["dataset"], r["diagnosis"], corrected[r["donor"]], r["sex"]] + [r[f"nuclei_{l}"] for l in LINEAGES] + umi)

for l in LINEAGES:
    m = counts[:, types.index(l), :]
    with gzip.open(OUT / f"counts_{l}.tsv.gz", "wt") as fh:
        fh.write("gene\t" + "\t".join(donors) + "\n")
        for j, g in enumerate(genes):
            fh.write(g + "\t" + "\t".join(str(int(v)) for v in m[:, j]) + "\n")


def load_modules(path):
    d = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        d[int(r["module"])].append((float(r["weight"]), r["gene"]))
    return {m: [g for _, g in sorted(v, reverse=True)] for m, v in d.items()}


full = load_modules(HS / "module_genes.tsv")
loo = load_modules(HS / "loo/GSE244832/module_genes.tsv")
with open(OUT / "B1_loo_match.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["full_module", "n_genes", "top_genes", "best_loo_module", "jaccard_top50", "n_shared_top50", "loo_top_genes"])
    for m in sorted(full):
        a = set(full[m][:50])
        scored = [(len(a & set(loo[k][:50])) / len(a | set(loo[k][:50])), k) for k in loo]
        jac, k = max(scored)
        w.writerow([m, len(full[m]), ",".join(full[m][:10]), k, f"{jac:.4f}",
                    len(a & set(loo[k][:50])), ",".join(loo[k][:10])])
print("export done", OUT)
