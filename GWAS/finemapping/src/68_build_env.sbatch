#!/usr/bin/env bash
#SBATCH --job-name=chrombpnet
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chrombpnet_envbuild_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/chrombpnet_envbuild_%j.err

# Build isolated variant-scorer (ChromBPNet) TF env. Keras-2 (TF 2.15) so the
# ENCODE ChromBPNet H5 models load with compile=False; tensorflow-probability
# supplies multinomial_nll (registered custom object, never actually called at
# inference). Isolated from borzoi/rapids/alphagenome (its own env dir).
set -uo pipefail
ENVDIR=/gpfs/commons/home/jameslee/micromamba/envs/chrombpnet_vs
mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs
echo "=== build start $(date) on $(hostname) ==="
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"

micromamba create -y -p "$ENVDIR" -c conda-forge -c bioconda \
  python=3.10 bedtools ucsc-bigbedtobed
echo "=== conda create rc=$? $(date) ==="

"$ENVDIR/bin/pip" install --no-cache-dir --no-input \
  "tensorflow[and-cuda]==2.15.1" \
  "tensorflow-probability==0.23.0" \
  "numpy<2" pandas scipy h5py tqdm pybedtools pyfaidx pysam
echo "=== pip rc=$? $(date) ==="

echo "=== import + GPU check $(date) ==="
"$ENVDIR/bin/python" - <<'PY'
import numpy, pandas, scipy, h5py, tqdm, pyfaidx, pysam, pybedtools
print("numpy", numpy.__version__, "| pandas", pandas.__version__)
import tensorflow as tf
print("tf", tf.__version__, "| keras", tf.keras.__version__)
import tensorflow_probability as tfp
print("tfp", tfp.__version__)
print("bedtools:", pybedtools.helpers.get_bedtools_path() if hasattr(pybedtools.helpers,'get_bedtools_path') else 'ok')
PY
echo "=== BUILD COMPLETE $(date) ==="
