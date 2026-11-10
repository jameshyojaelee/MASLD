"""Check that every pseudobulk column maps to one biological donor before the 4D content calls.

Reports, per lineage and dataset, the number of columns (runs/samples), the number mapped through a
donor_pairing.csv, and the number of distinct donors after collapse. Unmapped runs from a dataset that
has a pairing file would stay separate units (run-level pseudoreplication).
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
SC = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
meta = pd.read_csv(SC / "sample_celltype_metadata.tsv", sep="\t", dtype=str)
run2ds = dict(zip(meta["sample"], meta["dataset"]))

pairing = {}
for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103"):
    p = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for d, runs in zip(p["donor_id"], p["rna_srrs"]):
        for r in str(runs).split(";"):
            r = r.strip()
            if r and r != "nan":
                if r in pairing and pairing[r] != f"{ds}:{d}":
                    print("RUN MAPS TO TWO DONORS", r, pairing[r], ds, d)
                pairing[r] = f"{ds}:{d}"

rows = []
for ln in ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]:
    cols = pd.read_csv(SC / f"pseudobulk_{ln}.tsv.gz", sep="\t", nrows=0).columns[1:]
    dup_cols = int(pd.Index(cols).duplicated().sum())
    df = pd.DataFrame({"run": cols})
    df["dataset"] = df.run.map(run2ds)
    df["donor"] = df.run.map(pairing)
    for ds, g in df.groupby(df.dataset.fillna("NA")):
        rows.append(dict(lineage=ln, dataset=ds, n_columns=len(g), n_mapped=int(g.donor.notna().sum()),
                         n_unmapped=int(g.donor.isna().sum()),
                         n_units_after_collapse=int(g.donor.fillna(g.run).nunique()), dup_columns=dup_cols))
out = pd.DataFrame(rows)
print(out.to_string(index=False))
print(out.groupby("lineage").n_units_after_collapse.sum())
out.to_csv(sys.argv[1], sep="\t", index=False)
