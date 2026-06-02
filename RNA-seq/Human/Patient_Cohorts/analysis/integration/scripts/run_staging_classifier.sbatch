#!/bin/bash
# Master sbatch for staging classifier pipeline
# Submits Phase 1 scripts with proper dependencies, then Phases 2-4 in parallel
#
# Usage: bash run_staging_classifier.sbatch (NOT sbatch — this is a submission script, not a job)
# Or source it to get job IDs: source run_staging_classifier.sbatch

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/../results/staging_classifier/logs"
mkdir -p "${LOG_DIR}"

echo "=== Staging Classifier Pipeline ==="
echo "Script dir: ${SCRIPT_DIR}"
echo "Log dir: ${LOG_DIR}"
echo ""

# ============================================================
# PHASE 1: Shared Foundation
# ============================================================

# Script 60: One-vs-rest dream DEGs (bigmem, ~4h)
JOB60=$(sbatch --parsable \
  --job-name=stg60_ovr_dream \
  --partition=bigmem \
  --cpus-per-task=16 \
  --mem=200G \
  --time=48:00:00 \
  --output="${LOG_DIR}/60_ovr_dream_%j.out" \
  --error="${LOG_DIR}/60_ovr_dream_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 60_one_vs_rest_stage_dream.R'")
echo "Submitted Script 60 (OVR dream): ${JOB60}"

# Script 61: Feature engineering (cpu, ~1h; depends on Script 60 for OVR t-stats)
JOB61=$(sbatch --parsable \
  --dependency=afterok:${JOB60} \
  --job-name=stg61_features \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/61_features_%j.out" \
  --error="${LOG_DIR}/61_features_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 61_feature_engineering.R'")
echo "Submitted Script 61 (features): ${JOB61} (depends on ${JOB60})"

# Script 62: Prepare Python data (cpu, ~1h; depends on Script 61)
JOB62=$(sbatch --parsable \
  --dependency=afterok:${JOB61} \
  --job-name=stg62_pydata \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/62_pydata_%j.out" \
  --error="${LOG_DIR}/62_pydata_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && cd ${SCRIPT_DIR} && python3 62_prepare_python_data.py'")
echo "Submitted Script 62 (Python data): ${JOB62} (depends on ${JOB61})"

# ============================================================
# PHASE 2: Plan 1 — Ordinal Evidence Convergence
# ============================================================

# Script 63: Ordinal models (cpu, ~6h; depends on Script 62 for modeling_metadata.csv)
JOB63=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg63_ordinal \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/63_ordinal_%j.out" \
  --error="${LOG_DIR}/63_ordinal_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 63_ordinal_models.R'")
echo "Submitted Script 63 (ordinal models): ${JOB63} (depends on ${JOB61})"

# Script 64: Plan 1 model sweep (cpu, ~4h; depends on Script 62)
# Note: MLP uses CPU-only torch; no GPU needed
JOB64=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg64_sweep1 \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/64_sweep1_%j.out" \
  --error="${LOG_DIR}/64_sweep1_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 64_plan1_model_sweep.py'")
echo "Submitted Script 64 (Plan 1 sweep): ${JOB64} (depends on ${JOB62})"

# ============================================================
# PHASE 3: Plan 2 — Hierarchical Multi-Resolution
# ============================================================

# Script 65: Tier 1 disease detection (cpu, ~2h; depends on Script 62 for modeling_metadata.csv)
JOB65=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg65_tier1 \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/65_tier1_%j.out" \
  --error="${LOG_DIR}/65_tier1_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 65_tier1_disease_detection.R'")
echo "Submitted Script 65 (Tier 1): ${JOB65} (depends on ${JOB61})"

# Script 66: Tier 2 severity (cpu, ~4h; depends on Script 62 for modeling_metadata.csv)
JOB66=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg66_tier2 \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/66_tier2_%j.out" \
  --error="${LOG_DIR}/66_tier2_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 66_tier2_severity.R'")
echo "Submitted Script 66 (Tier 2): ${JOB66} (depends on ${JOB61})"

# Script 67: Tier 3 staging (cpu, ~6h; depends on Script 62 for modeling_metadata.csv)
JOB67=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg67_tier3 \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=64G \
  --time=48:00:00 \
  --output="${LOG_DIR}/67_tier3_%j.out" \
  --error="${LOG_DIR}/67_tier3_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 67_tier3_staging.R'")
echo "Submitted Script 67 (Tier 3): ${JOB67} (depends on ${JOB61})"

# Script 68: Hierarchical neural sweep (gpu, ~4h; depends on Script 62)
JOB68=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg68_hier_nn \
  --partition=gpu \
  --cpus-per-task=8 \
  --mem=64G \
  --gres=gpu:1 \
  --time=48:00:00 \
  --output="${LOG_DIR}/68_hier_nn_%j.out" \
  --error="${LOG_DIR}/68_hier_nn_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 68_hierarchical_sweep.py'")
echo "Submitted Script 68 (hierarchical NN): ${JOB68} (depends on ${JOB62})"

# ============================================================
# PHASE 4: Plan 3 — Foundation-Ready Embedding
# ============================================================

# Script 69: VAE training (gpu, ~6h; depends on Script 62)
JOB69=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg69_vae \
  --partition=gpu \
  --cpus-per-task=8 \
  --mem=64G \
  --gres=gpu:1 \
  --time=48:00:00 \
  --output="${LOG_DIR}/69_vae_%j.out" \
  --error="${LOG_DIR}/69_vae_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 69_vae_training.py'")
echo "Submitted Script 69 (VAE): ${JOB69} (depends on ${JOB62})"

# Script 70: Classical on embeddings (cpu, ~2h; depends on Script 69)
JOB70=$(sbatch --parsable \
  --dependency=afterok:${JOB69} \
  --job-name=stg70_emb_cls \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/70_emb_cls_%j.out" \
  --error="${LOG_DIR}/70_emb_cls_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 70_classical_on_embeddings.py'")
echo "Submitted Script 70 (embedding classifiers): ${JOB70} (depends on ${JOB69})"

# Script 71: Attention panel (gpu, ~4h; depends on Script 62)
JOB71=$(sbatch --parsable \
  --dependency=afterok:${JOB62} \
  --job-name=stg71_attn \
  --partition=gpu \
  --cpus-per-task=8 \
  --mem=64G \
  --gres=gpu:1 \
  --time=48:00:00 \
  --output="${LOG_DIR}/71_attn_%j.out" \
  --error="${LOG_DIR}/71_attn_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 71_attention_panel.py'")
echo "Submitted Script 71 (attention panel): ${JOB71} (depends on ${JOB62})"

# Script 72: Multimodal hooks (gpu, ~1h; depends on Script 69)
JOB72=$(sbatch --parsable \
  --dependency=afterok:${JOB69} \
  --job-name=stg72_mmhook \
  --partition=gpu \
  --cpus-per-task=4 \
  --mem=32G \
  --gres=gpu:1 \
  --time=48:00:00 \
  --output="${LOG_DIR}/72_mmhook_%j.out" \
  --error="${LOG_DIR}/72_mmhook_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 72_multimodal_hooks.py'")
echo "Submitted Script 72 (multimodal hooks): ${JOB72} (depends on ${JOB69})"

# ============================================================
# PHASE 5: Shared Outputs
# ============================================================

# All Phase 2-4 jobs must complete first
PHASE234_DEPS="${JOB63}:${JOB64}:${JOB65}:${JOB66}:${JOB67}:${JOB68}:${JOB69}:${JOB70}:${JOB71}:${JOB72}"

# Script 73: Minimal panel extraction (cpu, ~3h)
JOB73=$(sbatch --parsable \
  --dependency=afterok:${PHASE234_DEPS} \
  --job-name=stg73_panels \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/73_panels_%j.out" \
  --error="${LOG_DIR}/73_panels_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 73_minimal_panel_extraction.R'")
echo "Submitted Script 73 (panels): ${JOB73}"

# Script 74: Label imputation (cpu, ~1h)
JOB74=$(sbatch --parsable \
  --dependency=afterok:${PHASE234_DEPS} \
  --job-name=stg74_impute \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/74_impute_%j.out" \
  --error="${LOG_DIR}/74_impute_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rapids_singlecell && export LD_LIBRARY_PATH=\"\$(python -c \"import nvidia; import os; print(os.path.dirname(nvidia.__file__))\" 2>/dev/null)/cublas/lib:\${LD_LIBRARY_PATH:-}\" 2>/dev/null || true && cd ${SCRIPT_DIR} && python3 74_label_imputation.py'")
echo "Submitted Script 74 (imputation): ${JOB74}"

# ============================================================
# PHASE 6: Comparison & Synthesis
# ============================================================

# Script 75: Comparison synthesis (cpu, ~1h; depends on Phase 5)
JOB75=$(sbatch --parsable \
  --dependency=afterok:${JOB73}:${JOB74} \
  --job-name=stg75_compare \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=32G \
  --time=48:00:00 \
  --output="${LOG_DIR}/75_compare_%j.out" \
  --error="${LOG_DIR}/75_compare_%j.err" \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && cd ${SCRIPT_DIR} && Rscript 75_comparison_synthesis.R'")
echo "Submitted Script 75 (comparison): ${JOB75}"

# ============================================================
echo ""
echo "=== All jobs submitted ==="
echo "Phase 1: ${JOB60} -> ${JOB61} -> ${JOB62}"
echo "Phase 2: ${JOB63}, ${JOB64}"
echo "Phase 3: ${JOB65}, ${JOB66}, ${JOB67}, ${JOB68}"
echo "Phase 4: ${JOB69} -> ${JOB70}, ${JOB71}, ${JOB72}"
echo "Phase 5: ${JOB73}, ${JOB74}"
echo "Phase 6: ${JOB75}"
echo ""
echo "Monitor: squeue -u \$USER | grep stg"
echo "Logs: ${LOG_DIR}/"
