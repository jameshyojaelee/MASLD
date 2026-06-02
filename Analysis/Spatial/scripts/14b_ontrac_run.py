#!/usr/bin/env python3
"""
14b_ontrac_run.py — DPT diffusion-pseudotime niche trajectory
(ONTraC fallback; ONTraC CLI unavailable).

PROVENANCE NOTE (F090): the niche-trajectory (NT) scores and niche clusters
this script emits are produced by a DPT diffusion-pseudotime FALLBACK rooted at
the highest-hepatocyte spot, NOT by ONTraC. The ONTraC CLI is not available in
this environment, so the ONTraC branch below has never executed successfully —
every production run has dropped into ``fallback_dpt_trajectory`` and stamped
``method=diffusion_pseudotime`` into ``parameter_sweep.csv``. Filenames, the
``ontrac/`` results dir, and the ``ONTraC`` job names are retained only for
file-path/back-compat reasons; the method is DPT diffusion-pseudotime.

The script still ATTEMPTS ONTraC first (so it would use the real tool if the CLI
were ever installed) and falls back to diffusion pseudotime on the cell-type
composition space when the CLI is absent or all ONTraC runs fail. When it falls
back, it logs a clear WARNING that the output is DPT and writes the honest
``method``/``reason`` provenance to ``parameter_sweep.csv``.

SLURM: --partition=gpu --gres=gpu:1 --cpus=8 --mem=64G --time=8:00:00
"""

import pathlib
import subprocess
import shutil
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, init_spatial_gpu,
    save_csv, print_header, print_step,
)


def load_ontrac_input():
    """Load ONTraC input files from 14a."""
    input_dir = RESULTS_DIR / "ontrac" / "input"

    meta_path = input_dir / "ontrac_metadata.csv"
    comp_path = input_dir / "ontrac_cell_type_composition.csv"

    if not meta_path.exists() or not comp_path.exists():
        print("  ERROR: Run 14a_ontrac_prepare_input.py first")
        sys.exit(1)

    metadata = pd.read_csv(meta_path, index_col=0)
    composition = pd.read_csv(comp_path, index_col=0)
    print(f"  Loaded: {len(metadata)} spots, {composition.shape[1]} cell types")
    return metadata, composition


def run_ontrac_sweep(metadata, composition, config):
    """Run ONTraC via CLI with parameter sweep over niche cluster counts.

    ONTraC is invoked as a subprocess for each k value. Output CSVs
    are read back to extract NT scores and cluster assignments.

    Returns: (best_k, nt_scores, niche_clusters, sweep_results)
    """
    k_range = config["niche_cluster_range"]
    n_neighbors = config["n_neighbors"]
    gnn_epochs = config["gnn_epochs"]

    sweep_results = []
    best_k = config["default_k"]
    best_silhouette = -1
    best_nt_scores = None
    best_clusters = None

    input_dir = RESULTS_DIR / "ontrac" / "input"
    meta_path = input_dir / "ontrac_metadata.csv"
    comp_path = input_dir / "ontrac_cell_type_composition.csv"

    for i, k in enumerate(k_range):
        print_step(f"ONTraC k={k}", i + 1, len(k_range))

        # Per-k output directories
        run_dir = RESULTS_DIR / "ontrac" / f"run_k{k}"
        nn_dir = run_dir / "NN"
        gnn_dir = run_dir / "GNN"
        nt_dir = run_dir / "NT"
        for d in [nn_dir, gnn_dir, nt_dir]:
            d.mkdir(parents=True, exist_ok=True)

        try:
            cmd = [
                "ONTraC",
                "--meta-input", str(meta_path),
                "--deconvoluted-ct-composition", str(comp_path),
                "--NN-dir", str(nn_dir),
                "--GNN-dir", str(gnn_dir),
                "--NT-dir", str(nt_dir),
                "--n-neighbors", str(n_neighbors),
                "-k", str(k),
                "--epochs", str(gnn_epochs),
                "--device", "cuda",
                "--patience", "100",
                "--min-epochs", "50",
            ]
            print(f"    CMD: ONTraC ... -k {k} --epochs {gnn_epochs}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

            if result.returncode != 0:
                print(f"    WARNING: ONTraC failed for k={k}")
                if result.stderr:
                    print(f"    stderr: {result.stderr[-300:]}")
                sweep_results.append({
                    "k": k, "silhouette_score": np.nan,
                    "n_clusters_found": 0, "nt_score_range": np.nan,
                    "nt_score_std": np.nan,
                })
                continue

            # Read NT scores from output
            nt_file = None
            for candidate in sorted(nt_dir.glob("*NTScore*")):
                nt_file = candidate
            if nt_file is None:
                csvs = sorted(nt_dir.glob("*.csv.gz")) + sorted(nt_dir.glob("*.csv"))
                if csvs:
                    nt_file = csvs[0]

            if nt_file is None:
                print(f"    WARNING: No NTScore output for k={k}")
                sweep_results.append({
                    "k": k, "silhouette_score": np.nan,
                    "n_clusters_found": 0, "nt_score_range": np.nan,
                    "nt_score_std": np.nan,
                })
                continue

            nt_df = pd.read_csv(nt_file)
            nt_col = next(
                (c for c in nt_df.columns if "NTScore" in c or "nt_score" in c.lower()),
                nt_df.columns[-1],
            )
            nt_scores = nt_df[nt_col].values

            # Read cluster assignments
            cl_files = sorted(gnn_dir.glob("*cluster*")) + sorted(nt_dir.glob("*cluster*"))
            if cl_files:
                cl_df = pd.read_csv(cl_files[0])
                cl_col = next(
                    (c for c in cl_df.columns if "cluster" in c.lower()),
                    cl_df.columns[-1],
                )
                clusters = cl_df[cl_col].values
            else:
                clusters = np.zeros(len(nt_scores), dtype=int)

            # Silhouette score
            n_comp = min(len(clusters), len(composition))
            if len(np.unique(clusters)) > 1 and n_comp > 1:
                sil = silhouette_score(composition.values[:n_comp], clusters[:n_comp])
            else:
                sil = -1.0

            sweep_results.append({
                "k": k,
                "silhouette_score": sil,
                "n_clusters_found": len(np.unique(clusters)),
                "nt_score_range": float(nt_scores.max() - nt_scores.min()),
                "nt_score_std": float(nt_scores.std()),
            })

            print(f"    Silhouette: {sil:.3f}, NT range: "
                  f"[{nt_scores.min():.3f}, {nt_scores.max():.3f}]")

            if sil > best_silhouette:
                best_silhouette = sil
                best_k = k
                best_nt_scores = nt_scores
                best_clusters = clusters

        except subprocess.TimeoutExpired:
            print(f"    WARNING: ONTraC timed out for k={k}")
            sweep_results.append({
                "k": k, "silhouette_score": np.nan,
                "n_clusters_found": 0, "nt_score_range": np.nan,
                "nt_score_std": np.nan,
            })
        except Exception as e:
            print(f"    WARNING: ONTraC failed for k={k}: {e}")
            sweep_results.append({
                "k": k, "silhouette_score": np.nan,
                "n_clusters_found": 0, "nt_score_range": np.nan,
                "nt_score_std": np.nan,
            })

    print(f"\n  Best k: {best_k} (silhouette={best_silhouette:.3f})")
    return best_k, best_nt_scores, best_clusters, pd.DataFrame(sweep_results)


def fallback_dpt_trajectory(metadata, composition):
    """Fallback: compute trajectory via DPT diffusion-pseudotime on composition space.

    Used when the ONTraC CLI is unavailable or all ONTraC runs fail. The output
    NT scores and niche clusters are NOT ONTraC outputs — they are DPT
    diffusion-pseudotime (rooted at the highest-hepatocyte spot) plus Leiden
    clusters on the same composition neighbor graph. This is the path that has
    actually run in production (F090).
    """
    print("  *** WARNING (F090): OUTPUT IS DPT DIFFUSION-PSEUDOTIME, NOT ONTraC. ***")
    print("  ONTraC CLI unavailable / all ONTraC runs failed — falling back to "
          "DPT diffusion-pseudotime on cell-type composition space.")
    print("  NT scores + niche clusters below are the DPT fallback; "
          "parameter_sweep.csv records method=diffusion_pseudotime as provenance.")

    # Build AnnData from composition matrix
    import anndata as ad
    # Fill NaN and remove zero-variance columns
    comp_clean = np.nan_to_num(composition.values, nan=0.0)
    col_var = comp_clean.var(axis=0)
    keep_cols = col_var > 1e-10
    if keep_cols.sum() < comp_clean.shape[1]:
        print(f"  Removed {(~keep_cols).sum()} zero-variance cell type columns")
        comp_clean = comp_clean[:, keep_cols]
    adata = ad.AnnData(X=comp_clean, obs=metadata.copy())
    adata.obs_names = metadata.index.astype(str)

    # PCA + neighbors on composition space
    n_comps = min(10, comp_clean.shape[1] - 1)
    if n_comps < 2:
        n_comps = 2
    sc.pp.pca(adata, n_comps=n_comps, svd_solver='arpack', zero_center=True)
    sc.pp.neighbors(adata, n_neighbors=25, n_pcs=n_comps)

    # Diffusion map
    sc.tl.diffmap(adata, n_comps=10)

    # Root: spot with highest hepatocyte proportion (most "normal")
    hep_cols = [c for c in composition.columns
                if "hepat" in c.lower() or "hep" in c.lower()]
    if hep_cols:
        hep_score = composition[hep_cols].sum(axis=1)
        root_idx = hep_score.idxmax()
    else:
        root_idx = composition.index[0]

    adata.uns["iroot"] = list(adata.obs_names).index(str(root_idx))
    sc.tl.dpt(adata)

    nt_scores = adata.obs["dpt_pseudotime"].values

    # Cluster via Leiden on composition neighbors
    sc.tl.leiden(adata, resolution=0.5, key_added="niche_cluster",
                 flavor="igraph", n_iterations=2, directed=False)
    clusters = adata.obs["niche_cluster"].astype(int).values

    # Silhouette
    if len(np.unique(clusters)) > 1:
        sil = silhouette_score(composition.values, clusters)
    else:
        sil = -1.0
    print(f"  Fallback silhouette: {sil:.3f}")

    return nt_scores, clusters


def main():
    print_header("14b: DPT Niche Trajectory (ONTraC fallback; CLI unavailable)")

    config = load_config()
    ontrac_config = config["ontrac"]
    output_dir = RESULTS_DIR / "ontrac"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize GPU
    print_step("Initializing GPU")
    init_spatial_gpu()

    # Load input from 14a
    print_step("Loading ONTraC input")
    metadata, composition = load_ontrac_input()

    # Check if ONTraC CLI is available
    use_fallback = False
    ontrac_check = shutil.which("ONTraC")
    if ontrac_check:
        print(f"  ONTraC CLI found: {ontrac_check}")
    else:
        print("  WARNING: ONTraC CLI not found in PATH")
        print("  Falling back to diffusion pseudotime on composition space")
        use_fallback = True

    if not use_fallback:
        # Attempt real ONTraC with parameter sweep (only if the CLI is present).
        print_step("Attempting ONTraC parameter sweep")
        best_k, nt_scores, clusters, sweep_df = run_ontrac_sweep(
            metadata, composition, ontrac_config
        )

        # Save sweep results
        save_csv(sweep_df, "parameter_sweep.csv", subdir="ontrac")

        if nt_scores is None:
            print("  *** WARNING (F090): all ONTraC runs failed — OUTPUT WILL BE "
                  "DPT DIFFUSION-PSEUDOTIME, NOT ONTraC. Falling back to DPT. ***")
            use_fallback = True

    if use_fallback:
        nt_scores, clusters = fallback_dpt_trajectory(metadata, composition)
        best_k = len(np.unique(clusters))

        # Save fallback marker — honest method provenance preserved so any
        # downstream reader can see the NT axis is DPT, not ONTraC (F090).
        fallback_df = pd.DataFrame([{
            "method": "diffusion_pseudotime",
            "reason": "ONTraC CLI not available or all runs failed",
            "k": best_k,
        }])
        save_csv(fallback_df, "parameter_sweep.csv", subdir="ontrac")

    # Build output DataFrames
    nt_df = pd.DataFrame({
        "barcode": metadata.index,
        "nt_score": nt_scores,
        "niche_cluster": clusters,
    })
    nt_df.index = metadata.index
    save_csv(nt_df, "niche_trajectory_scores.csv", subdir="ontrac")

    # Separate cluster summary
    cluster_df = pd.DataFrame({
        "barcode": metadata.index,
        "niche_cluster": clusters,
        "Sample": metadata["Sample"].values,
        "Cell_Type": metadata["Cell_Type"].values,
    })
    cluster_df.index = metadata.index
    save_csv(cluster_df, "niche_clusters.csv", subdir="ontrac")

    # Summary statistics
    print(f"\n  Summary:")
    print(f"    Optimal k: {best_k}")
    print(f"    NT score range: [{nt_scores.min():.3f}, {nt_scores.max():.3f}]")
    print(f"    NT score mean: {nt_scores.mean():.3f} +/- {nt_scores.std():.3f}")
    print(f"    Cluster sizes:")
    for c in sorted(np.unique(clusters)):
        n = (clusters == c).sum()
        print(f"      Cluster {c}: {n} spots ({n / len(clusters) * 100:.1f}%)")

    print_header("14b: Complete (DPT diffusion-pseudotime; ONTraC fallback)")


if __name__ == "__main__":
    main()
