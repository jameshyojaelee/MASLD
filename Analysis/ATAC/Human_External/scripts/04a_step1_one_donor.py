#!/usr/bin/env python
"""04a_step1_one_donor.py -- parallel step1 (import/QC/tile) for ONE GSE281367 donor.
Reuses 02_snapatac2_processing.step1_import_qc_tile verbatim, one donor per SLURM
array task, so the ~43min/donor import+tile runs 12-way concurrent instead of serial.
Writes snapatac2_fast/per_donor/{donor}.h5ad (backed, tile matrix + QC, fragment refs).
Env: snapatac2."""
import os, sys, importlib.util, logging

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
S = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome/scripts/02_snapatac2_processing.py")
spec = importlib.util.spec_from_file_location("snap02", S)
mod = importlib.util.module_from_spec(spec); sys.modules["snap02"] = mod
spec.loader.exec_module(mod)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("04a")

class A: pass
a = A()
a.input_dir  = os.path.join(ROOT, "Analysis/ATAC/Human_External/frag_links")
a.output_dir = os.path.join(ROOT, "Analysis/ATAC/Human_External/snapatac2_fast")
a.sample_csv = os.path.join(ROOT, "data/GSE281367/metadata/donor_pairing.csv")
a.n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
a.min_tsse = 4.0; a.min_frags = 1000; a.max_frags = 100000
a.tile_size = 500; a.n_features = 50000; a.leiden_resolution = 1.0; a.skip_da = True
os.makedirs(os.path.join(a.output_dir, "per_donor"), exist_ok=True)

df = mod.load_sample_metadata(a.sample_csv)
donors = mod.discover_donors(a.input_dir, df, logger)
i = int(os.environ["SLURM_ARRAY_TASK_ID"]) - 1
one = [donors[i]]
donor_id = one[0][0]
out_h5 = os.path.join(a.output_dir, "per_donor", f"{donor_id}.h5ad")
if os.path.exists(out_h5):
    logger.info("%s already exists -> skip", out_h5); sys.exit(0)
logger.info("step1 import/QC/tile for donor %s", donor_id)
mod.step1_import_qc_tile(one, a, logger)
logger.info("DONE donor %s -> %s", donor_id, out_h5)
