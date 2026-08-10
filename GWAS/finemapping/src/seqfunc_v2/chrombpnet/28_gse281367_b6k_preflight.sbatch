#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=gpu
#SBATCH --qos=interactive
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/b6k_preflight_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/b6k_preflight_%j.err
set -euo pipefail

module load CUDA/11.7.0
module load cuDNN/8.4.1.50-CUDA-11.7.0
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON=${ROOT}/.mamba/seqfunc_chrombpnet/bin/python

nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader
"${PYTHON}" -u - <<'PY'
import tensorflow as tf

print(f"tensorflow={tf.__version__}")
devices = tf.config.list_physical_devices("GPU")
print(f"gpu_devices={devices}")
if not devices:
    raise SystemExit("TensorFlow did not detect the allocated B6K GPU")
x = tf.random.normal([256, 256])
y = tf.linalg.matmul(x, x)
print(f"matmul_sum={float(tf.reduce_sum(y).numpy())}")
PY
