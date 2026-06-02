#!/bin/bash
# 401_env_cnmf.sh — install cNMF (Kotliar 2019) into rapids_singlecell env.
# cNMF is CPU-bound (sklearn NMF backend) but we install alongside rapids
# so GPU-accelerated data prep (scanpy on RAPIDS) and cNMF can share one env.
set -euo pipefail

ENV_NAME="${ENV_NAME:-rapids_singlecell}"
echo "[401] Installing cNMF into: ${ENV_NAME}"

micromamba run -n "${ENV_NAME}" pip install --no-cache-dir \
    "cnmf==1.5.3" \
    "nimfa==1.4.0"

echo "[401] Verifying import..."
micromamba run -n "${ENV_NAME}" python -c "
import cnmf
import nimfa
print(f'cnmf version: {cnmf.__version__}')
print(f'nimfa version: {nimfa.__version__}')
from cnmf import cNMF
print('cNMF class importable.')
"

echo "[401] DONE."
