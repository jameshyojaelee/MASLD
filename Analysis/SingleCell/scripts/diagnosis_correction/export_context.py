"""Export corrected scRNA context to a local review copy, without publishing."""
import json
import sys
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
sys.path.insert(0, str(root / "scripts/portal"))
from _portal_io import dump_json, load_ensembl_symbol_map, write_parquet
from generate_singlecell_data import build_pseudobulk_parquet, build_hep_markers_parquet, CELL_TYPES

out = task / "portal"
out.mkdir(exist_ok=False)
props = pd.read_csv(task / "cell_type_proportions.csv")
meta = pd.read_csv(task / "metadata/library_metadata_source_F.tsv", sep="\t")
props = props.drop(columns=[c for c in ("donor", "binary_eligible", "source_histology", "qualification_reason") if c in props])
props = props.merge(meta[["sample", "donor", "binary_eligible", "source_histology", "qualification_reason"]], on="sample", validate="one_to_one")
present = [c for c in CELL_TYPES if c in props]
composition = []
for _, r in props.iterrows():
    composition.append(dict(sample=r["sample"], donor=r["donor"], dataset=r["dataset"], condition=r["condition"],
        condition_harmonized=r["condition_harmonized"], binary_eligible=bool(r["binary_eligible"]),
        source_histology=r["source_histology"] if pd.notna(r["source_histology"]) else "",
        qualification_reason=r["qualification_reason"], fractions={c: float(r[c]) for c in present}))
sub = pd.read_csv(task / "composition/hepatocyte_subtype_results.tsv", sep="\t")
summary = dict(n_cells=int(meta.n_cells.sum()), n_datasets=int(meta.dataset.nunique()), n_samples=len(meta),
    n_source_donor_identifiers=int(meta.donor.nunique()), cell_types=present, composition=composition,
    lineage_composition_test="Donor-level arcsine-square-root proportion model, adjusted for dataset; BH over 16 lineages",
    celltype_disease_activity=[], fgsea_transitions=[], hep_subtypes=[],
    historical_score_state="Disease scores and fibrosis-proxy enrichments from superseded labels are withdrawn; no current Resource claim uses them.")
for _, r in sub.iterrows():
    summary["hep_subtypes"].append(dict(subtype=str(r["feature"]), donor_n=int(r["donor_n"]),
        adjusted_arcsine_effect=float(r["logFC"]), donor_padj=float(r["adj.P.Val"]),
        difference_percentage_points=float(r["difference_percentage_points"])))
dump_json(summary, out / "singlecell_summary.json")
emap = load_ensembl_symbol_map(root)
build_pseudobulk_parquet(root / "Analysis/SingleCell/results_gpu_v2", out, emap)
build_hep_markers_parquet(root / "Analysis/SingleCell/results_gpu_v2", out)
write_parquet(pd.read_csv(task / "composition/lineages_results.tsv", sep="\t"), out / "sc_lineage_composition.parquet", sort_col="feature")
(out / "source.json").write_text(json.dumps(dict(diagnosis="source-diagnosis-correction-20261026", published=False,
    binary="pseudobulk_de_v2", stage="pseudobulk_de_coarse", allcell="allcell_v2"), indent=2))
