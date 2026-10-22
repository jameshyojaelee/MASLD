"""WP5: 117-program observability across assays (candidate, not adopted).

Adds six new datasets to the 7 dataset rows of spatial-resource-candidate-2026-08-11, using the
same rule (spatial_resource_lib.coverage_state: >= 8 genes and >= 0.20 of fixed L1 weight among
genes detected in >= 1% of all QC cells, spots or segments). A second column applies the 0.80
replication-scoring rule. Detection rates only; no program score or outcome is computed.
Usage: python -I 01_observability_panel.py
"""
import gzip
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = Path(__file__).resolve().parents[1]
OUT = RUN / "outputs"
H = ROOT / "Analysis/Hypotheses"
sys.path.insert(0, str(ROOT / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"))
from spatial_resource_lib import DETECTION_MIN_FRACTION, PROGRAM_RELEASE_ID, coverage_state, load_frozen_programs  # noqa: E402

HOT = ROOT / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
EXISTING = ROOT / "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/coverage/spatial_program_coverage.tsv"
ANNOT = ROOT / ("RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/"
                "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz")
CONTROL_PREFIXES = ("Negative", "NegPrb", "SystemControl", "FalseCode", "NegativeAdd")
REPLICATION_MIN_WEIGHT = 0.80


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def frac_detected(X):
    X = sparse.csr_matrix(X)
    return np.asarray((X > 0).sum(axis=0)).ravel() / X.shape[0]


def cosmx_qc(X):
    X = sparse.csr_matrix(X)
    total = np.asarray(X.sum(axis=1)).ravel()
    ngen = np.asarray((X > 0).sum(axis=1)).ravel()
    return X[(total >= 20) & (ngen >= 10)]


def targets_only(df):
    return df[[c for c in df.columns if not str(c).startswith(CONTROL_PREFIXES)]]


def det_flex():
    import anndata as ad
    a = ad.read_h5ad(H / "V5-20261021T122723Z/derived/gse298719_qc_cells.h5ad")
    return pd.Series(frac_detected(a.X), index=a.var_names), "QC cells, 13 donors"


def det_citeseq():
    import anndata as ad
    import scanpy as sc
    n_cells, hits, genes = 0, None, None
    for f in sorted((H / "V8-20261021T235438Z/downloads_untrusted/raw").glob("GSM*_filtered_feature_bc_matrix.h5")):
        a = sc.read_10x_h5(f, gex_only=True)
        a.var_names_make_unique()
        X = sparse.csr_matrix(a.X)
        ngen = np.asarray((X > 0).sum(axis=1)).ravel()
        mt = a.var_names.str.startswith("MT-")
        tot = np.asarray(X.sum(axis=1)).ravel()
        pmt = np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(tot, 1) * 100
        X = X[(ngen >= 200) & (pmt < 25)]
        h = np.asarray((X > 0).sum(axis=0)).ravel()
        if genes is None:
            genes, hits = a.var_names, h
        else:
            assert genes.equals(a.var_names)
            hits = hits + h
        n_cells += X.shape[0]
        del a, X
    return pd.Series(hits / n_cells, index=genes), "QC cells (>= 200 genes, < 25% MT), 20 samples"


def det_geomx():
    k4 = pd.read_csv(H / "V6-20261021T234432Z/outputs/K4_gene_detection_alb_vin.tsv", sep="\t")
    return pd.Series(k4.frac_above_loq.to_numpy(), index=k4.gene), "AOIs above LOQ, Albumin+ and Vinculin+ segments"


def det_cosmx_atomx(files):
    parts = []
    for f in files:
        e = pd.read_csv(f)
        e = e[e.cell_ID != 0].drop(columns=["fov", "cell_ID"])
        parts.append(targets_only(e))
    genes = sorted(set.intersection(*[set(p.columns) for p in parts]))
    X = sparse.vstack([sparse.csr_matrix(p[genes].to_numpy(np.float32)) for p in parts])
    return pd.Series(frac_detected(cosmx_qc(X)), index=genes)


def det_gse292268():
    d = H / "V7-20261021T235717Z/downloads_untrusted/gse292268"
    s = det_cosmx_atomx([d / "GSM8855042_Slide1_exprMat_file.csv.gz", d / "GSM8855043_Slide2_exprMat_file.csv.gz"])
    return s, "QC cells (>= 20 counts, >= 10 genes), 2 slides"


def det_fukumoto():
    d = H / "V7-20261021T235717Z/downloads_untrusted/fukumoto"
    parts = [targets_only(pd.read_csv(d / f, index_col=0)) for f in
             ("CosMx_expression_data_MASH_sample1_NT.csv", "CosMx_expression_data_MASH_sample_2_NT.csv")]
    genes = sorted(set.intersection(*[set(p.columns) for p in parts]))
    X = sparse.vstack([sparse.csr_matrix(p[genes].to_numpy(np.float32)) for p in parts])
    return pd.Series(frac_detected(cosmx_qc(X)), index=genes), "QC cells (>= 20 counts, >= 10 genes), 2 sections"


def det_xenium():
    d = H / "V7-20261021T235717Z/downloads_untrusted/scp3359/SCP3359/expression/68f4186f4e6da7a7e1b0edbe"
    genes = pd.read_csv(d / "counts_features.tsv", header=None)[0].astype(str).to_numpy()
    X = scipy.io.mmread(d / "counts_matrix.mtx").tocsr()
    if X.shape[0] == len(genes):
        X = X.T.tocsr()
    assert X.shape[1] == len(genes), X.shape
    return pd.Series(frac_detected(cosmx_qc(X)), index=genes), "QC cells (>= 20 counts, >= 10 genes), 18 donors"


NEW = [  # dataset_id, assay_id, loader, biological_unit, lineage_context
    ("GSE298719", "flex_probe_rna_dissociated", det_flex, "donor", "dissociated_whole_cell_all_lineages"),
    ("GSE344087", "citeseq_rna_dissociated", det_citeseq, "donor", "dissociated_sorted_cells"),
    ("GSE255312", "geomx_wta_probe_rna", det_geomx, "roi_no_donor_key", "albumin_and_vinculin_segments"),
    ("Fukumoto2025_CosMx", "cosmx_targeted_rna", det_fukumoto, "section", "segmented_cells"),
    ("GSE292268", "cosmx_targeted_rna", det_gse292268, "biopsy", "segmented_cells"),
    ("SCP3359", "xenium_targeted_rna", det_xenium, "donor", "segmented_cells"),
]

if __name__ == "__main__":
    assert sha256(HOT / "program_membership_v2.tsv").startswith("feec3fc9"), "membership hash"
    OUT.mkdir(exist_ok=False)
    programs, weights = load_frozen_programs(HOT)
    union = sorted({g for w in weights.values() for g in w})
    annot = pd.read_csv(ANNOT, sep="\t", usecols=["gene_name", "gene_biotype"]).drop_duplicates("gene_name")
    biotype = dict(zip(annot.gene_name, annot.gene_biotype))

    rows, gene_rows, audit = [], [], []
    for ds, assay, loader, unit, ctx in NEW:
        det, scope = loader()
        det = det.groupby(level=0).max()
        measured, detected = set(det.index), set(det.index[det >= DETECTION_MIN_FRACTION])
        audit.append(dict(dataset_id=ds, assay_id=assay, gene_axis=len(measured), detected_genes=len(detected),
                          detection_scope=scope))
        print(ds, len(measured), len(detected), flush=True)
        for g in union:
            gene_rows.append(dict(dataset_id=ds, gene=g, measured=g in measured,
                                  detection_fraction=float(det.get(g, np.nan)), detected=g in detected))
        for uid, prog in programs.items():
            w = weights[uid]
            on, dg = set(w) & measured, set(w) & detected
            ret = sum(w[g] for g in dg)
            status, reason = coverage_state(len(dg), ret)
            rows.append(dict(dataset_id=ds, assay_id=assay, program_uid=uid, program_label=prog["module_name"],
                             cell_type=prog["cell_type"], n_program_genes=len(w), n_genes_on_axis=len(on),
                             axis_l1_weight=sum(w[g] for g in on), n_genes_detected=len(dg), retained_l1_weight=ret,
                             coverage_status=status, testability_reason=reason,
                             replication_rule_0p80=bool(len(dg) >= 8 and ret >= REPLICATION_MIN_WEIGHT),
                             biological_unit=unit, lineage_context=ctx, source="WP5 " + scope))
    ex = pd.read_csv(EXISTING, sep="\t")
    ex = ex.assign(axis_l1_weight=np.nan, replication_rule_0p80=(ex.n_genes_detected >= 8) & (ex.retained_l1_weight >= 0.80),
                   source="spatial-resource-candidate-2026-08-11")
    cols = list(pd.DataFrame(rows).columns)
    panel = pd.concat([ex.reindex(columns=cols), pd.DataFrame(rows)], ignore_index=True)
    assert panel.groupby(["dataset_id", "assay_id"]).program_uid.nunique().eq(117).all()
    panel.to_csv(OUT / "observability_panel_117.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(OUT / "gene_detection_new_datasets.tsv.gz", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(OUT / "gene_axis_audit_new.tsv", sep="\t", index=False)

    summ = (panel.groupby(["dataset_id", "assay_id"])
            .agg(observable=("coverage_status", lambda s: int((s == "observable").sum())),
                 partial=("coverage_status", lambda s: int((s == "partial").sum())),
                 untestable=("coverage_status", lambda s: int((s == "untestable").sum())),
                 replication_0p80=("replication_rule_0p80", "sum")).reset_index())
    summ.to_csv(OUT / "observability_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))

    # lncRNA share of fixed weight that each new gene axis carries (gene presence only)
    lnc = []
    for ds, assay, *_ in NEW:
        g = pd.DataFrame(gene_rows)
        g = g[g.dataset_id == ds].set_index("gene")
        for uid, w in weights.items():
            lw = {k: v for k, v in w.items() if biotype.get(k) == "lncRNA"}
            if not lw:
                continue
            lnc.append(dict(dataset_id=ds, assay_id=assay, program_uid=uid, n_lncrna=len(lw), lncrna_weight=sum(lw.values()),
                            lncrna_weight_on_axis=sum(v for k, v in lw.items() if g.measured.get(k, False)),
                            other_weight=sum(v for k, v in w.items() if k not in lw),
                            other_weight_on_axis=sum(v for k, v in w.items() if k not in lw and g.measured.get(k, False))))
    pd.DataFrame(lnc).to_csv(OUT / "lncrna_axis_presence.tsv", sep="\t", index=False)
    lt = pd.DataFrame(lnc).groupby(["dataset_id", "assay_id"])[["lncrna_weight", "lncrna_weight_on_axis", "other_weight", "other_weight_on_axis"]].sum()
    lt["lncrna_on_axis_frac"] = lt.lncrna_weight_on_axis / lt.lncrna_weight
    lt["other_on_axis_frac"] = lt.other_weight_on_axis / lt.other_weight
    lt.to_csv(OUT / "lncrna_axis_presence_summary.tsv", sep="\t")
    print(lt.to_string())
