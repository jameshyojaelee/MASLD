# RETRACTED: "BulkFormer" embeddings are actually original VAE embeddings due to Script 100 fallback. All results invalid. See review notes.
#!/usr/bin/env python3
"""
167_bulkformer_transfer.py
BulkFormer Transfer Learning — tests whether frozen foundation model
embeddings + classification head can match task-specific elastic net.

Comparison:
  1. LogisticRegressionCV on BulkFormer 64-dim embeddings (directly comparable)
  2. Elastic net on BulkFormer embeddings (same features, classical model)
  3. Elastic net on NAS-VAE embeddings (task-specific representation, 64-dim)
  4. DL head on BulkFormer embeddings (Linear→ReLU→Dropout→Linear; or CORAL for ordinal)
  5. Baselines: random embeddings (permuted BulkFormer), clinical-only (sex+age)

Targets:
  - F>=3 binary (primary)
  - Fibrosis ordinal F0-F4 (secondary, CORAL if torch+coral available)

LOCO-CV: 6 fibrosis folds (leave-one-cohort-out).
DL option: 5-fold inner CV on training data for early stopping epoch.

Input (from results/staging_classifier/):
  - bulkformer_embeddings.csv      1,444 x 64 (bf0..bf63)
  - nas_embeddings_all_samples.csv 1,444 x 64 (z0..z63)
  - modeling_metadata.csv          labels + LOCO folds

Output (to results/prognosis_v2/):
  - bulkformer_transfer_results.csv       per-fold AUROC/QWK per model variant
  - bulkformer_vs_nasvae_comparison.csv   head-to-head on same folds

Env:   rapids_singlecell (torch + coral-pytorch). Falls back to sklearn-only.
SLURM: io partition, 4 CPUs, 8G, 48h

Usage:
  sbatch --job-name=stg167_bf_transfer \
         --partition=io --cpus-per-task=4 --mem=8G --time=48:00:00 \
         --output=logs/167_bf_transfer_%j.out \
         --error=logs/167_bf_transfer_%j.err \
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \
                 micromamba activate rapids_singlecell && \
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \
                 python 167_bulkformer_transfer.py'"
"""

import os
import sys
import time
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------------------
# Try importing torch + coral (DL option)
# ---------------------------------------------------------------------------
HAS_TORCH = False
HAS_CORAL = False
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import TensorDataset, DataLoader
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
    HAS_TORCH = True
    try:
        from coral_pytorch.layers import CoralLayer
        from coral_pytorch.losses import coral_loss as coral_loss_fn
        from coral_pytorch.dataset import levels_from_labelbatch
        HAS_CORAL = True
    except ImportError:
        pass
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier",
)
OUTDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/prognosis_v2",
)
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def quadratic_weighted_kappa(y_true, y_pred, n_classes):
    """Compute QWK from integer labels."""
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(y_true, y_pred, labels=list(range(n_classes)))
    n = cm.sum()
    if n == 0:
        return 0.0
    w = np.zeros((n_classes, n_classes))
    for i in range(n_classes):
        for j in range(n_classes):
            w[i, j] = (i - j) ** 2 / (n_classes - 1) ** 2
    expected = np.outer(cm.sum(axis=1), cm.sum(axis=0)) / n
    num = (w * cm).sum()
    den = (w * expected).sum()
    return 1.0 - num / den if den > 0 else 1.0


def safe_auroc(y_true, y_score):
    """AUROC that returns NaN if only one class present."""
    if len(np.unique(y_true)) < 2:
        return np.nan
    return roc_auc_score(y_true, y_score)


# ===================================================================
# 1. Load data
# ===================================================================
log.info("=== 167: BulkFormer Transfer Learning ===")
log.info(f"Torch available: {HAS_TORCH}")
log.info(f"CORAL available: {HAS_CORAL}")

# --- Metadata ---
meta = pd.read_csv(os.path.join(SDIR, "modeling_metadata.csv"))
meta["sample_id"] = meta["sample_id"].astype(str).str.strip()
log.info(f"Metadata: {meta.shape[0]} samples, {meta.shape[1]} columns")

# --- BulkFormer embeddings (64-dim) ---
bf_df = pd.read_csv(os.path.join(SDIR, "bulkformer_embeddings.csv"))
bf_df["sample_id"] = bf_df["sample_id"].astype(str).str.strip()
bf_cols = [c for c in bf_df.columns if c.startswith("bf")]
bf_dim = len(bf_cols)
log.info(f"BulkFormer embeddings: {bf_df.shape[0]} samples, {bf_dim}-dim")

# --- NAS-VAE embeddings (64-dim) ---
nasvae_df = pd.read_csv(os.path.join(SDIR, "nas_embeddings_all_samples.csv"))
nasvae_df["sample_id"] = nasvae_df["sample_id"].astype(str).str.strip()
nv_cols = [c for c in nasvae_df.columns if c.startswith("z")]
nv_dim = len(nv_cols)
log.info(f"NAS-VAE embeddings: {nasvae_df.shape[0]} samples, {nv_dim}-dim")

# --- Merge ---
merged = meta.merge(bf_df, on="sample_id", how="inner")
merged = merged.merge(nasvae_df, on="sample_id", how="inner")
log.info(f"Merged dataset: {merged.shape[0]} samples")

# --- Filter to fibrosis-eligible (6 LOCO folds) ---
fib_data = merged[merged["loco_fold_fibrosis"] != "excluded"].copy()
fib_data = fib_data[fib_data["fib_stage"] >= 0].copy()
fib_folds = sorted(fib_data["loco_fold_fibrosis"].unique())
log.info(f"Fibrosis-eligible: {len(fib_data)} samples, {len(fib_folds)} folds")
log.info(f"Folds: {fib_folds}")
log.info(f"fib_ge3 distribution: {fib_data['fib_ge3'].value_counts().to_dict()}")
log.info(f"fib_stage distribution: {fib_data['fib_stage'].value_counts().sort_index().to_dict()}")

# --- Prepare clinical features ---
# Encode sex: M=0, F=1 (NaN treated as 0; GSE135251 + GSE240729 lack sex/age)
# This is a sanity-check baseline — expected to perform poorly
fib_data["sex_enc"] = (fib_data["sex"] == "F").astype(int)
# Age: fill missing with median (278/718 samples missing age)
age_median = fib_data.loc[fib_data["age"].notna(), "age"].median()
fib_data["age_fill"] = fib_data["age"].fillna(age_median)
clinical_cols = ["sex_enc", "age_fill"]

# --- Random embeddings baseline (permute bf features per-column) ---
rng = np.random.RandomState(SEED)
bf_arr = fib_data[bf_cols].values.copy()
random_emb = np.zeros_like(bf_arr)
for j in range(bf_arr.shape[1]):
    random_emb[:, j] = rng.permutation(bf_arr[:, j])
rand_cols = [f"rand_{i}" for i in range(bf_dim)]
for i, rc in enumerate(rand_cols):
    fib_data[rc] = random_emb[:, i]


# ===================================================================
# 2. Define model runners
# ===================================================================
def run_logistic_cv(X_train, y_train, X_test):
    """LogisticRegressionCV with elastic net (l1_ratio grid)."""
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)
    model = LogisticRegressionCV(
        Cs=10,
        cv=5,
        penalty="elasticnet",
        solver="saga",
        l1_ratios=[0.1, 0.5, 0.9],
        max_iter=5000,
        random_state=SEED,
        scoring="roc_auc",
        n_jobs=-1,
    )
    model.fit(X_tr, y_train)
    proba = model.predict_proba(X_te)[:, 1]
    pred = model.predict(X_te)
    return proba, pred


def run_logistic_simple(X_train, y_train, X_test):
    """Simple LogisticRegressionCV (L2)."""
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)
    model = LogisticRegressionCV(
        Cs=10, cv=5, penalty="l2", solver="lbfgs",
        max_iter=5000, random_state=SEED, scoring="roc_auc", n_jobs=-1,
    )
    model.fit(X_tr, y_train)
    proba = model.predict_proba(X_te)[:, 1]
    pred = model.predict(X_te)
    return proba, pred


def run_ordinal_logistic(X_train, y_train, X_test, n_classes):
    """Ordinal prediction via chained binary classifiers (cumulative model)."""
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)
    # Fit K-1 binary classifiers: P(Y > k) for k=0,..,K-2
    models = []
    for k in range(n_classes - 1):
        y_bin = (y_train > k).astype(int)
        if len(np.unique(y_bin)) < 2:
            models.append(None)
            continue
        m = LogisticRegression(
            C=1.0, penalty="l2", solver="lbfgs", max_iter=5000, random_state=SEED,
        )
        m.fit(X_tr, y_bin)
        models.append(m)
    # Predict cumulative probabilities -> ordinal class
    cum_probs = np.zeros((X_te.shape[0], n_classes - 1))
    for k, m in enumerate(models):
        if m is not None:
            cum_probs[:, k] = m.predict_proba(X_te)[:, 1]
        else:
            cum_probs[:, k] = 0.5
    # Convert cumulative P(Y>k) to class probabilities
    class_probs = np.zeros((X_te.shape[0], n_classes))
    class_probs[:, 0] = 1.0 - cum_probs[:, 0]
    for k in range(1, n_classes - 1):
        class_probs[:, k] = cum_probs[:, k - 1] - cum_probs[:, k]
    class_probs[:, n_classes - 1] = cum_probs[:, n_classes - 2]
    # Clip negatives (numerical)
    class_probs = np.clip(class_probs, 0, 1)
    row_sums = class_probs.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    class_probs = class_probs / row_sums
    pred = class_probs.argmax(axis=1)
    return class_probs, pred


# --- DL models (if torch available) ---
if HAS_TORCH:
    class BinaryHead(nn.Module):
        """Two-layer classification head for binary target."""
        def __init__(self, input_dim, hidden_dim=128, dropout=0.3):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, 1),
            )

        def forward(self, x):
            return self.net(x).squeeze(-1)

    def train_binary_head(X_train, y_train, X_test, input_dim, n_inner_folds=5):
        """Train DL binary head with inner CV for early stopping epoch.

        Applies StandardScaler internally — pass raw (unscaled) features.
        """
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Inner CV to determine best epoch
        inner_skf = StratifiedKFold(n_splits=n_inner_folds, shuffle=True, random_state=SEED)
        best_epochs = []

        for inner_tr_idx, inner_val_idx in inner_skf.split(X_train, y_train):
            sc_inner = StandardScaler()
            X_itr_np = sc_inner.fit_transform(X_train[inner_tr_idx])
            X_ival_np = sc_inner.transform(X_train[inner_val_idx])
            X_itr = torch.tensor(X_itr_np, dtype=torch.float32).to(device)
            y_itr = torch.tensor(y_train[inner_tr_idx], dtype=torch.float32).to(device)
            X_ival = torch.tensor(X_ival_np, dtype=torch.float32).to(device)
            y_ival = torch.tensor(y_train[inner_val_idx], dtype=torch.float32).to(device)

            model = BinaryHead(input_dim).to(device)
            optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
            criterion = nn.BCEWithLogitsLoss()

            best_val_loss = float("inf")
            best_epoch = 0
            patience_count = 0
            patience = 20

            for epoch in range(300):
                model.train()
                optimizer.zero_grad()
                logits = model(X_itr)
                loss = criterion(logits, y_itr)
                loss.backward()
                optimizer.step()

                model.eval()
                with torch.no_grad():
                    val_logits = model(X_ival)
                    val_loss = criterion(val_logits, y_ival).item()
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_epoch = epoch
                    patience_count = 0
                else:
                    patience_count += 1
                    if patience_count >= patience:
                        break

            best_epochs.append(best_epoch)

        n_epochs = int(np.median(best_epochs)) + 1
        n_epochs = max(n_epochs, 10)  # minimum 10 epochs
        log.info(f"    DL inner CV best epochs: {best_epochs} -> training for {n_epochs}")

        # Retrain on full training set
        scaler = StandardScaler()
        X_tr_sc = scaler.fit_transform(X_train)
        X_te_sc = scaler.transform(X_test)

        X_tr_t = torch.tensor(X_tr_sc, dtype=torch.float32).to(device)
        y_tr_t = torch.tensor(y_train, dtype=torch.float32).to(device)
        X_te_t = torch.tensor(X_te_sc, dtype=torch.float32).to(device)

        model = BinaryHead(input_dim).to(device)
        optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
        criterion = nn.BCEWithLogitsLoss()

        model.train()
        for epoch in range(n_epochs):
            optimizer.zero_grad()
            logits = model(X_tr_t)
            loss = criterion(logits, y_tr_t)
            loss.backward()
            optimizer.step()

        model.eval()
        with torch.no_grad():
            te_logits = model(X_te_t).cpu().numpy()
        proba = 1.0 / (1.0 + np.exp(-te_logits))  # sigmoid
        pred = (proba >= 0.5).astype(int)
        return proba, pred

    if HAS_CORAL:
        class CoralOrdinalHead(nn.Module):
            """Two-layer head with CORAL ordinal output."""
            def __init__(self, input_dim, n_classes, hidden_dim=128, dropout=0.3):
                super().__init__()
                self.encoder = nn.Sequential(
                    nn.Linear(input_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(dropout),
                )
                self.coral = CoralLayer(hidden_dim, n_classes)

            def forward(self, x):
                h = self.encoder(x)
                return self.coral(h)

        def train_coral_head(X_train, y_train, X_test, input_dim, n_classes, n_inner_folds=5):
            """Train CORAL ordinal head with inner CV for early stopping.

            Applies StandardScaler internally — pass raw (unscaled) features.
            """
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

            inner_skf = StratifiedKFold(n_splits=n_inner_folds, shuffle=True, random_state=SEED)
            best_epochs = []

            for inner_tr_idx, inner_val_idx in inner_skf.split(X_train, y_train):
                sc_inner = StandardScaler()
                X_itr_np = sc_inner.fit_transform(X_train[inner_tr_idx])
                X_ival_np = sc_inner.transform(X_train[inner_val_idx])
                X_itr = torch.tensor(X_itr_np, dtype=torch.float32).to(device)
                y_itr = torch.tensor(y_train[inner_tr_idx], dtype=torch.long).to(device)
                X_ival = torch.tensor(X_ival_np, dtype=torch.float32).to(device)
                y_ival = torch.tensor(y_train[inner_val_idx], dtype=torch.long).to(device)

                model = CoralOrdinalHead(input_dim, n_classes).to(device)
                optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)

                best_val_loss = float("inf")
                best_epoch = 0
                patience_count = 0
                patience = 20

                for epoch in range(300):
                    model.train()
                    optimizer.zero_grad()
                    logits = model(X_itr)
                    levels = levels_from_labelbatch(y_itr, num_classes=n_classes)
                    loss = coral_loss_fn(logits, levels)
                    loss.backward()
                    optimizer.step()

                    model.eval()
                    with torch.no_grad():
                        val_logits = model(X_ival)
                        val_levels = levels_from_labelbatch(y_ival, num_classes=n_classes)
                        val_loss = coral_loss_fn(val_logits, val_levels).item()
                    if val_loss < best_val_loss:
                        best_val_loss = val_loss
                        best_epoch = epoch
                        patience_count = 0
                    else:
                        patience_count += 1
                        if patience_count >= patience:
                            break

                best_epochs.append(best_epoch)

            n_epochs = int(np.median(best_epochs)) + 1
            n_epochs = max(n_epochs, 10)
            log.info(f"    CORAL inner CV best epochs: {best_epochs} -> training for {n_epochs}")

            scaler = StandardScaler()
            X_tr_sc = scaler.fit_transform(X_train)
            X_te_sc = scaler.transform(X_test)

            X_tr_t = torch.tensor(X_tr_sc, dtype=torch.float32).to(device)
            y_tr_t = torch.tensor(y_train, dtype=torch.long).to(device)
            X_te_t = torch.tensor(X_te_sc, dtype=torch.float32).to(device)

            model = CoralOrdinalHead(input_dim, n_classes).to(device)
            optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)

            model.train()
            for epoch in range(n_epochs):
                optimizer.zero_grad()
                logits = model(X_tr_t)
                tr_levels = levels_from_labelbatch(y_tr_t, num_classes=n_classes)
                loss = coral_loss_fn(logits, tr_levels)
                loss.backward()
                optimizer.step()

            model.eval()
            with torch.no_grad():
                te_logits = model(X_te_t).cpu().numpy()
            # CORAL: cumulative logits -> ordinal predictions
            cum_probs = 1.0 / (1.0 + np.exp(-te_logits))
            pred = (cum_probs > 0.5).sum(axis=1)
            pred = np.clip(pred, 0, n_classes - 1)
            return cum_probs, pred


# ===================================================================
# 3. LOCO-CV: F>=3 binary (primary target)
# ===================================================================
log.info("\n" + "=" * 70)
log.info("SECTION A: F>=3 Binary Classification (LOCO-CV)")
log.info("=" * 70)

results_binary = []

for fold in fib_folds:
    log.info(f"\n--- Fold: {fold} ---")
    test_mask = fib_data["loco_fold_fibrosis"] == fold
    train_mask = ~test_mask

    y_train = fib_data.loc[train_mask, "fib_ge3"].values
    y_test = fib_data.loc[test_mask, "fib_ge3"].values
    log.info(f"  Train: {train_mask.sum()} (pos={y_train.sum()}), Test: {test_mask.sum()} (pos={y_test.sum()})")

    # -- (A1) LogisticRegressionCV on BulkFormer embeddings --
    X_tr_bf = fib_data.loc[train_mask, bf_cols].values
    X_te_bf = fib_data.loc[test_mask, bf_cols].values

    proba, pred = run_logistic_cv(X_tr_bf, y_train, X_te_bf)
    auroc = safe_auroc(y_test, proba)
    acc = accuracy_score(y_test, pred)
    results_binary.append({
        "fold": fold, "model": "LR_BulkFormer", "target": "fib_ge3",
        "auroc": auroc, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  LR_BulkFormer: AUROC={auroc:.3f}, Acc={acc:.3f}")

    # -- (A2) Elastic net on BulkFormer embeddings --
    proba, pred = run_logistic_simple(X_tr_bf, y_train, X_te_bf)
    auroc = safe_auroc(y_test, proba)
    acc = accuracy_score(y_test, pred)
    results_binary.append({
        "fold": fold, "model": "ElasticNet_BulkFormer", "target": "fib_ge3",
        "auroc": auroc, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  ElasticNet_BulkFormer: AUROC={auroc:.3f}, Acc={acc:.3f}")

    # -- (A3) Elastic net on NAS-VAE embeddings --
    X_tr_nv = fib_data.loc[train_mask, nv_cols].values
    X_te_nv = fib_data.loc[test_mask, nv_cols].values

    proba, pred = run_logistic_cv(X_tr_nv, y_train, X_te_nv)
    auroc = safe_auroc(y_test, proba)
    acc = accuracy_score(y_test, pred)
    results_binary.append({
        "fold": fold, "model": "ElasticNet_NASVAE", "target": "fib_ge3",
        "auroc": auroc, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  ElasticNet_NASVAE: AUROC={auroc:.3f}, Acc={acc:.3f}")

    # -- (A4) Random embeddings baseline --
    X_tr_rand = fib_data.loc[train_mask, rand_cols].values
    X_te_rand = fib_data.loc[test_mask, rand_cols].values

    proba, pred = run_logistic_simple(X_tr_rand, y_train, X_te_rand)
    auroc = safe_auroc(y_test, proba)
    acc = accuracy_score(y_test, pred)
    results_binary.append({
        "fold": fold, "model": "Random_Embeddings", "target": "fib_ge3",
        "auroc": auroc, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  Random_Embeddings: AUROC={auroc:.3f}, Acc={acc:.3f}")

    # -- (A5) Clinical only (sex + age) --
    X_tr_clin = fib_data.loc[train_mask, clinical_cols].values
    X_te_clin = fib_data.loc[test_mask, clinical_cols].values

    # Skip if both clinical features have zero variance (e.g., all imputed)
    has_variance = any(np.std(X_tr_clin[:, j]) > 0 for j in range(X_tr_clin.shape[1]))
    if has_variance:
        proba, pred = run_logistic_simple(X_tr_clin, y_train, X_te_clin)
        auroc = safe_auroc(y_test, proba)
        acc = accuracy_score(y_test, pred)
    else:
        auroc = np.nan
        acc = np.nan
    results_binary.append({
        "fold": fold, "model": "Clinical_Only", "target": "fib_ge3",
        "auroc": auroc, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  Clinical_Only: AUROC={auroc:.3f}, Acc={acc:.3f}")

    # -- (A6) DL binary head on BulkFormer (if torch) --
    if HAS_TORCH:
        proba, pred = train_binary_head(X_tr_bf, y_train, X_te_bf, bf_dim)
        auroc = safe_auroc(y_test, proba)
        acc = accuracy_score(y_test, pred)
        results_binary.append({
            "fold": fold, "model": "DL_BulkFormer", "target": "fib_ge3",
            "auroc": auroc, "accuracy": acc, "n_test": len(y_test),
        })
        log.info(f"  DL_BulkFormer: AUROC={auroc:.3f}, Acc={acc:.3f}")

binary_df = pd.DataFrame(results_binary)

# Summary
log.info("\n--- F>=3 Binary Summary (mean +/- std across folds) ---")
for model in binary_df["model"].unique():
    sub = binary_df[binary_df["model"] == model]
    auroc_m = sub["auroc"].mean()
    auroc_s = sub["auroc"].std()
    acc_m = sub["accuracy"].mean()
    log.info(f"  {model:25s}: AUROC={auroc_m:.3f} +/- {auroc_s:.3f}, Acc={acc_m:.3f}")


# ===================================================================
# 4. LOCO-CV: Fibrosis ordinal F0-F4 (secondary target)
# ===================================================================
log.info("\n" + "=" * 70)
log.info("SECTION B: Fibrosis Ordinal F0-F4 (LOCO-CV)")
log.info("=" * 70)

N_FIB_CLASSES = 5  # F0, F1, F2, F3, F4
results_ordinal = []

for fold in fib_folds:
    log.info(f"\n--- Fold: {fold} ---")
    test_mask = fib_data["loco_fold_fibrosis"] == fold
    train_mask = ~test_mask

    y_train = fib_data.loc[train_mask, "fib_stage"].values
    y_test = fib_data.loc[test_mask, "fib_stage"].values
    log.info(f"  Train: {train_mask.sum()}, Test: {test_mask.sum()}")
    log.info(f"  Test stage dist: {dict(zip(*np.unique(y_test, return_counts=True)))}")

    # -- (B1) Ordinal LR on BulkFormer --
    X_tr_bf = fib_data.loc[train_mask, bf_cols].values
    X_te_bf = fib_data.loc[test_mask, bf_cols].values

    class_probs, pred = run_ordinal_logistic(X_tr_bf, y_train, X_te_bf, N_FIB_CLASSES)
    qwk = quadratic_weighted_kappa(y_test, pred, N_FIB_CLASSES)
    acc = accuracy_score(y_test, pred)
    results_ordinal.append({
        "fold": fold, "model": "OrdinalLR_BulkFormer", "target": "fib_ordinal",
        "qwk": qwk, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  OrdinalLR_BulkFormer: QWK={qwk:.3f}, Acc={acc:.3f}")

    # -- (B2) Ordinal LR on NAS-VAE --
    X_tr_nv = fib_data.loc[train_mask, nv_cols].values
    X_te_nv = fib_data.loc[test_mask, nv_cols].values

    class_probs, pred = run_ordinal_logistic(X_tr_nv, y_train, X_te_nv, N_FIB_CLASSES)
    qwk = quadratic_weighted_kappa(y_test, pred, N_FIB_CLASSES)
    acc = accuracy_score(y_test, pred)
    results_ordinal.append({
        "fold": fold, "model": "OrdinalLR_NASVAE", "target": "fib_ordinal",
        "qwk": qwk, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  OrdinalLR_NASVAE: QWK={qwk:.3f}, Acc={acc:.3f}")

    # -- (B3) Random baseline --
    X_tr_rand = fib_data.loc[train_mask, rand_cols].values
    X_te_rand = fib_data.loc[test_mask, rand_cols].values

    class_probs, pred = run_ordinal_logistic(X_tr_rand, y_train, X_te_rand, N_FIB_CLASSES)
    qwk = quadratic_weighted_kappa(y_test, pred, N_FIB_CLASSES)
    acc = accuracy_score(y_test, pred)
    results_ordinal.append({
        "fold": fold, "model": "Random_Embeddings", "target": "fib_ordinal",
        "qwk": qwk, "accuracy": acc, "n_test": len(y_test),
    })
    log.info(f"  Random_Embeddings: QWK={qwk:.3f}, Acc={acc:.3f}")

    # -- (B4) CORAL on BulkFormer (if available) --
    if HAS_TORCH and HAS_CORAL:
        cum_probs, pred = train_coral_head(X_tr_bf, y_train, X_te_bf, bf_dim, N_FIB_CLASSES)
        qwk = quadratic_weighted_kappa(y_test, pred, N_FIB_CLASSES)
        acc = accuracy_score(y_test, pred)
        results_ordinal.append({
            "fold": fold, "model": "CORAL_BulkFormer", "target": "fib_ordinal",
            "qwk": qwk, "accuracy": acc, "n_test": len(y_test),
        })
        log.info(f"  CORAL_BulkFormer: QWK={qwk:.3f}, Acc={acc:.3f}")

    # -- (B5) CORAL on NAS-VAE (if available) --
    if HAS_TORCH and HAS_CORAL:
        cum_probs, pred = train_coral_head(X_tr_nv, y_train, X_te_nv, nv_dim, N_FIB_CLASSES)
        qwk = quadratic_weighted_kappa(y_test, pred, N_FIB_CLASSES)
        acc = accuracy_score(y_test, pred)
        results_ordinal.append({
            "fold": fold, "model": "CORAL_NASVAE", "target": "fib_ordinal",
            "qwk": qwk, "accuracy": acc, "n_test": len(y_test),
        })
        log.info(f"  CORAL_NASVAE: QWK={qwk:.3f}, Acc={acc:.3f}")

ordinal_df = pd.DataFrame(results_ordinal)

log.info("\n--- Fibrosis Ordinal Summary (mean +/- std across folds) ---")
for model in ordinal_df["model"].unique():
    sub = ordinal_df[ordinal_df["model"] == model]
    qwk_m = sub["qwk"].mean()
    qwk_s = sub["qwk"].std()
    acc_m = sub["accuracy"].mean()
    log.info(f"  {model:25s}: QWK={qwk_m:.3f} +/- {qwk_s:.3f}, Acc={acc_m:.3f}")


# ===================================================================
# 5. Save results
# ===================================================================
log.info("\n" + "=" * 70)
log.info("SECTION C: Save Results")
log.info("=" * 70)

# Combine binary and ordinal results
all_results = pd.concat([binary_df, ordinal_df], ignore_index=True)
out_path = os.path.join(OUTDIR, "bulkformer_transfer_results.csv")
all_results.to_csv(out_path, index=False)
log.info(f"Saved: {out_path} ({len(all_results)} rows)")

# --- Head-to-head comparison table ---
comparison_rows = []

# Binary: BulkFormer vs NAS-VAE
for fold in fib_folds:
    bf_row = binary_df[
        (binary_df["fold"] == fold) & (binary_df["model"] == "LR_BulkFormer")
    ]
    nv_row = binary_df[
        (binary_df["fold"] == fold) & (binary_df["model"] == "ElasticNet_NASVAE")
    ]
    if len(bf_row) > 0 and len(nv_row) > 0:
        comparison_rows.append({
            "fold": fold,
            "target": "fib_ge3",
            "metric": "auroc",
            "bulkformer_value": bf_row["auroc"].values[0],
            "nasvae_value": nv_row["auroc"].values[0],
            "delta": bf_row["auroc"].values[0] - nv_row["auroc"].values[0],
        })

# Ordinal: BulkFormer vs NAS-VAE
for fold in fib_folds:
    bf_row = ordinal_df[
        (ordinal_df["fold"] == fold) & (ordinal_df["model"] == "OrdinalLR_BulkFormer")
    ]
    nv_row = ordinal_df[
        (ordinal_df["fold"] == fold) & (ordinal_df["model"] == "OrdinalLR_NASVAE")
    ]
    if len(bf_row) > 0 and len(nv_row) > 0:
        comparison_rows.append({
            "fold": fold,
            "target": "fib_ordinal",
            "metric": "qwk",
            "bulkformer_value": bf_row["qwk"].values[0],
            "nasvae_value": nv_row["qwk"].values[0],
            "delta": bf_row["qwk"].values[0] - nv_row["qwk"].values[0],
        })

comparison_df = pd.DataFrame(comparison_rows)
comp_path = os.path.join(OUTDIR, "bulkformer_vs_nasvae_comparison.csv")
comparison_df.to_csv(comp_path, index=False)
log.info(f"Saved: {comp_path} ({len(comparison_df)} rows)")

# --- Print final head-to-head summary ---
log.info("\n" + "=" * 70)
log.info("HEAD-TO-HEAD: BulkFormer vs NAS-VAE")
log.info("=" * 70)

for target in comparison_df["target"].unique():
    sub = comparison_df[comparison_df["target"] == target]
    metric = sub["metric"].values[0]
    bf_mean = sub["bulkformer_value"].mean()
    nv_mean = sub["nasvae_value"].mean()
    delta_mean = sub["delta"].mean()
    winner = "BulkFormer" if delta_mean > 0 else "NAS-VAE"
    log.info(
        f"  {target:15s} ({metric:5s}): BulkFormer={bf_mean:.3f}, "
        f"NAS-VAE={nv_mean:.3f}, delta={delta_mean:+.3f} -> {winner}"
    )

elapsed = time.time()
log.info(f"\nDone. Script 167 complete.")
