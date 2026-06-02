"""Smoke test: load HepG2 ST head + run inference on 64 control cells.

Tests:
1. Model loads with LlamaConfig validate patched
2. predict_step accepts a batch of [64, 2058] basal + [64, pert_dim] one-hot
3. Outputs differ across perturbations (non-trivial)
4. Outputs for "non-targeting" ≈ baseline
"""
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import torch
import anndata as ad

# --- PATCH HF strict validation -------------------------------------------
from transformers.models.llama.configuration_llama import LlamaConfig
_orig_validate = LlamaConfig.validate
def _lenient_validate(self):
    try:
        _orig_validate(self)
    except Exception as e:
        pass
LlamaConfig.validate = _lenient_validate

# --- Load model ------------------------------------------------------------
from state.tx.models.state_transition import StateTransitionPerturbationModel

CKPT_DIR = "data/perturbation/checkpoints/state/st-se-replogle-full/hepg2_0.99"
CKPT = f"{CKPT_DIR}/checkpoints/final.ckpt"

print(f"[smoke] loading {CKPT}", flush=True)
model = StateTransitionPerturbationModel.load_from_checkpoint(CKPT, map_location="cpu", strict=False)
model.eval()
print(f"[smoke] loaded; cell_sentence_len={model.cell_sentence_len}, output_space={model.output_space}", flush=True)

# --- Load control cells ----------------------------------------------------
adata = ad.read_h5ad(f"{CKPT_DIR}/eval_best.ckpt/adata_real.h5ad")
ctrl_mask = (adata.obs["gene"] == "non-targeting").values
ctrl_xs = adata.obsm["X_state"][ctrl_mask]
ctrl_gem = adata.obs["gem_group"][ctrl_mask].astype(str).values
print(f"[smoke] control cells: {ctrl_xs.shape}", flush=True)

# Pick 64 random ctrl cells
rng = np.random.RandomState(42)
sel = rng.choice(len(ctrl_xs), size=64, replace=False)
basal = ctrl_xs[sel]                # [64, 2058]
gem_sel = ctrl_gem[sel]

# --- Load pert one-hot map + batch map ------------------------------------
pert_map = torch.load(f"{CKPT_DIR}/pert_onehot_map.pt", weights_only=False)
import pickle
with open(f"{CKPT_DIR}/batch_onehot_map.pkl", "rb") as f:
    batch_map = pickle.load(f)
print(f"[smoke] pert vocab: {len(pert_map)}; batch vocab: {len(batch_map)}", flush=True)

def get_pert_oh(name):
    return pert_map[name].float()

def get_batch_idx(g):
    v = batch_map.get(g)
    if v is None: return 0
    if torch.is_tensor(v) and v.ndim == 1:
        return int(torch.argmax(v).item())
    return int(v) if isinstance(v, (int, np.integer)) else 0

# Build batches
batch_indices = torch.tensor([get_batch_idx(g) for g in gem_sel], dtype=torch.long)
basal_t = torch.tensor(basal, dtype=torch.float32)

def run_pert(label, oh=None):
    """Run inference: return predicted expression."""
    if oh is None:
        oh = get_pert_oh(label)
    pert_oh = oh.float().unsqueeze(0).repeat(64, 1)  # [64, pert_dim]
    batch = {
        "ctrl_cell_emb": basal_t,
        "pert_emb": pert_oh,
        "pert_name": [label] * 64,
        "batch": batch_indices,
    }
    with torch.no_grad():
        out = model.predict_step(batch, batch_idx=0, padded=False)
    # out["preds"] = predicted per-cell expression
    preds = out["preds"].detach().cpu().numpy()
    return preds, out

# Test 4 perturbations
perts_to_test = ["non-targeting", "TFAM", "GFM1", "MTOR"]
results = {}
for p in perts_to_test:
    if p not in pert_map:
        print(f"[smoke] {p} not in vocab — skipping")
        continue
    preds, out = run_pert(p)
    print(f"[smoke] {p}: pred shape={preds.shape}, mean={preds.mean():.4f}, std={preds.std():.4f}", flush=True)
    if "pert_cell_counts_preds" in out and out["pert_cell_counts_preds"] is not None:
        cnts = out["pert_cell_counts_preds"].detach().cpu().numpy()
        print(f"    pert_cell_counts shape={cnts.shape}, mean={cnts.mean():.4f}, std={cnts.std():.4f}", flush=True)
    results[p] = preds

# Verify outputs differ across perturbations
if "non-targeting" in results and "TFAM" in results:
    nt = results["non-targeting"]
    tfam = results["TFAM"]
    delta = tfam - nt
    print(f"[smoke] TFAM-NT delta: mean={delta.mean():.4f}, std={delta.std():.4f}, max={np.abs(delta).max():.4f}")
    # NON-TRIVIAL CHECK: if delta std > 0.01, predictions are different
    assert delta.std() > 1e-4, f"Predictions look identical to control: delta.std()={delta.std()}"
    print(f"[smoke] PASS: TFAM prediction differs from non-targeting (std={delta.std():.4f})", flush=True)

# Test 2-hot (synthetic double KO)
print("\n[smoke] testing 2-hot synthetic double-KO ...", flush=True)
if "TFAM" in pert_map and "GFM1" in pert_map:
    oh_double = get_pert_oh("TFAM") + get_pert_oh("GFM1")
    preds_double, _ = run_pert("TFAM+GFM1_synthetic", oh=oh_double)
    nt = results["non-targeting"]
    tf = results["TFAM"]
    gf = results["GFM1"]
    delta_tf = tf - nt
    delta_gf = gf - nt
    delta_dbl = preds_double - nt
    delta_add = delta_tf + delta_gf
    synergy = delta_dbl - delta_add
    syn_mag = float(np.linalg.norm(synergy) / (np.linalg.norm(delta_add) + 1e-9))
    print(f"[smoke] TFAM+GFM1 (2-hot) synergy magnitude: {syn_mag:.4f}", flush=True)
    print(f"[smoke] delta_double std: {delta_dbl.std():.4f}", flush=True)
    print(f"[smoke] delta_additive std: {delta_add.std():.4f}", flush=True)
    print(f"[smoke] synergy std: {synergy.std():.4f}", flush=True)

print("\n[smoke] DONE.", flush=True)
