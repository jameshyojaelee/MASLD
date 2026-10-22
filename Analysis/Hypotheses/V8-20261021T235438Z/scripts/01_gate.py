"""V8: GSE344087 CITE-seq gate. QC, CellTypist labels, hepatocyte counts per donor,
and the observability row for all 117 programs. Follows PRESPEC_V8.md.
No program score is computed. Usage: python -I 01_gate.py
"""
import gzip
import hashlib
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/V8-20261021T235438Z"
DL = RUN / "downloads_untrusted"
OUT = RUN / "outputs"
HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
MEMBERSHIP, REGISTRY = HOT / "program_membership_v2.tsv", HOT / "program_registry_v2.tsv"
CT_MODEL = Path("/gpfs/commons/home/jameslee/.celltypist/data/models/Healthy_Human_Liver.pkl")
LINEAGE_LABEL = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "cholangiocytes": "Cholangiocytes",
                 "macrophages": "Macrophages", "tcells": "T cells"}
SEED, MIN_HEP, DET_MIN, MIN_GENES, MIN_WEIGHT, MIN_LINEAGE = 20260923, 30, 0.01, 8, 0.80, 100


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


assert sha256(RUN / "PRESPEC_V8.md") == (RUN / "PRESPEC_V8.md.sha256").read_text().split()[0], "prespec hash"
assert sha256(MEMBERSHIP).startswith("feec3fc9"), "membership hash"
assert sha256(CT_MODEL) == "1b9999ba474e9423c583802a3482a900537c488d0e5973dbe20fa702446b8431", "CellTypist model hash"
OUT.mkdir(exist_ok=False)

# ---- design from SOFT -------------------------------------------------------------------
rows, cur = [], {}
for line in gzip.open(DL / "GSE344087_family.soft.gz", "rt"):
    if line.startswith("^SAMPLE"):
        cur = {"gsm": line.split("= ")[1].strip()}
        rows.append(cur)
    m = re.match(r"!Sample_characteristics_ch1 = (sample id|disease|fibrosis): (.+)", line.strip())
    if m:
        cur[m.group(1).replace(" ", "_")] = m.group(2)
design = pd.DataFrame(rows)
assert len(design) == 20, len(design)

# ---- GEX, QC, labels ---------------------------------------------------------------------
adatas = []
for f in sorted((DL / "raw").glob("GSM*_filtered_feature_bc_matrix.h5")):
    gsm = f.name.split("_")[0]
    a = sc.read_10x_h5(f, gex_only=True)
    a.var_names_make_unique()
    a.obs_names = [f"{gsm}_{b}" for b in a.obs_names]
    a.obs["gsm"] = gsm
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
    a.obs["raw_cells"] = a.n_obs
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.pct_counts_mt < 25)].copy()
    adatas.append(a)
A = ad.concat(adatas, join="inner")
A.X = sparse.csr_matrix(A.X)
counts = A.X.copy()
sc.pp.normalize_total(A, target_sum=1e4)
sc.pp.log1p(A)
import celltypist  # noqa: E402

np.random.seed(SEED)
pred = celltypist.annotate(A, model=celltypist.models.Model.load(model=str(CT_MODEL)), majority_voting=True)
A.obs["label"] = pred.predicted_labels["majority_voting"].astype(str).values

lin = pd.crosstab(A.obs.gsm, A.obs.label)
lin.to_csv(OUT / "lineage_counts_per_donor.tsv", sep="\t")
dq = design.merge(A.obs.groupby("gsm").agg(raw_cells=("raw_cells", "first"), cells_qc=("gsm", "size")).reset_index(), on="gsm")
hep_counts = lin["Hepatocytes"] if "Hepatocytes" in lin.columns else pd.Series(dtype=int)
dq["n_hep"] = dq.gsm.map(hep_counts).fillna(0).astype(int)
dq.to_csv(OUT / "donor_hepatocytes.tsv", sep="\t", index=False)
nofib = dq[(dq.disease == "MASLD") & (dq.fibrosis == "no fibrosis")]
fib = dq[(dq.disease == "MASLD") & (dq.fibrosis == "fibrosis")]
gate = (nofib.n_hep >= MIN_HEP).all() and len(nofib) == 3 and (fib.n_hep >= MIN_HEP).sum() >= 8

# ---- observability row, 117 programs ----------------------------------------------------
mem = pd.read_csv(MEMBERSHIP, sep="\t")
mem = mem[mem.mapped_symbol_status == "gencode_v49_unique_symbol_confirmed"]
reg = pd.read_csv(REGISTRY, sep="\t", usecols=["cell_type", "module", "program_uid", "module_name"])
genes = pd.Index(A.var_names)
cov = []
for ct, label in LINEAGE_LABEL.items():
    m = (A.obs.label == label).to_numpy()
    det = np.asarray((counts[m] > 0).mean(axis=0)).ravel() if m.sum() else np.zeros(len(genes))
    det = pd.Series(det, index=genes)
    for _, r in reg[reg.cell_type == ct].iterrows():
        w = mem[mem.program_uid == r.program_uid][["mapped_symbol", "original_l1_weight"]]
        on = w[w.mapped_symbol.isin(genes)]
        dd = on[det.reindex(on.mapped_symbol).to_numpy() >= DET_MIN]
        state = "lineage not captured" if m.sum() < MIN_LINEAGE else (
            "testable" if len(dd) >= MIN_GENES and dd.original_l1_weight.sum() >= MIN_WEIGHT else "not measured")
        cov.append(dict(program_uid=r.program_uid, cell_type=ct, module=r.module, module_name=r.module_name,
                        lineage_cells=int(m.sum()), n_members=len(w), axis_genes=len(on),
                        axis_weight=on.original_l1_weight.sum(), detected_genes=len(dd),
                        detected_weight=dd.original_l1_weight.sum(), state=state))
cov = pd.DataFrame(cov)
cov.to_csv(OUT / "observability_117_programs.tsv", sep="\t", index=False)
pd.DataFrame({"value": {
    "n_donors": len(dq), "n_nofib": len(nofib), "n_fib": len(fib),
    "nofib_min_hep": int(nofib.n_hep.min()), "fib_donors_ge30_hep": int((fib.n_hep >= MIN_HEP).sum()),
    "gate_pass": bool(gate), "cells_qc": int(A.n_obs), "hepatocytes_total": int((A.obs.label == "Hepatocytes").sum()),
    "programs_testable": int((cov.state == "testable").sum())}}).to_csv(OUT / "gate_V8.tsv", sep="\t")
print(dq.to_string(), flush=True)
print(cov.groupby(["cell_type", "state"]).size().to_string(), flush=True)
print("gate_pass", gate)
