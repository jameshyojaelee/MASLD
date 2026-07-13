#!/usr/bin/env python
"""04b_merge_annotate.py -- merge the 12 per-donor step1 h5ads (from 04a) and run
02's step2-5 (AnnDataSet -> LSI -> harmony -> leiden -> marker annotation) ONCE,
SKIPPING step6 macs3 (we pool over GSE244832 peaks, not GSE281367's own peaks).
Writes snapatac2_fast/snapatac2_processed.h5ad (.X=tile matrix, obs cell_type+donor_id).
Env: snapatac2."""
import os, sys, importlib.util, logging
import snapatac2 as snap

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
S = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome/scripts/02_snapatac2_processing.py")
spec = importlib.util.spec_from_file_location("snap02", S)
mod = importlib.util.module_from_spec(spec); sys.modules["snap02"] = mod
spec.loader.exec_module(mod)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("04b")

class A: pass
a = A()
a.input_dir  = os.path.join(ROOT, "Analysis/ATAC/Human_External/frag_links")
a.output_dir = os.path.join(ROOT, "Analysis/ATAC/Human_External/snapatac2_fast")
a.sample_csv = os.path.join(ROOT, "data/GSE281367/metadata/donor_pairing.csv")
a.n_features = 50000; a.leiden_resolution = 1.0

df = mod.load_sample_metadata(a.sample_csv)
donors_info = mod.discover_donors(a.input_dir, df, logger)
dcm = {d: cond for d, cond, _ in donors_info}

# reopen the 12 per-donor backed h5ads (retain fragment refs for make_gene_matrix)
processed = []
for d, cond, frag in donors_info:
    p = os.path.join(a.output_dir, "per_donor", f"{d}.h5ad")
    if not os.path.exists(p):
        logger.warning("missing per-donor h5ad for %s (%s) -> skip", d, p); continue
    processed.append((d, snap.read(p)))
logger.info("Reopened %d per-donor h5ads", len(processed))
assert len(processed) >= 10, f"only {len(processed)} donors available"

dataset = mod.step2_create_dataset(processed, a.output_dir, logger)
dataset = mod.step3_feature_selection_lsi(dataset, a, logger)
adata, use_rep = mod.step4_harmony_and_clustering(dataset, dcm, a, logger)
adata = mod.step5_cell_type_annotation(adata, dataset, logger)
# SKIP step6 macs3 (unneeded; we pool over GSE244832 peaks)
out = os.path.join(a.output_dir, "snapatac2_processed.h5ad")
adata.write(out)
logger.info("WROTE %s (macs3 skipped). cell_type counts:\n%s",
            out, adata.obs["cell_type"].value_counts().to_string())
