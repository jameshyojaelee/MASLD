#!/usr/bin/env bash
# upload-data.sh — Sync the atlas data tree to a public Hugging Face Dataset.
#
# Idempotent: creates the dataset repo on first run, uploads changed files
# only. Uses the `hf` CLI from huggingface_hub>=1.11 (installed at
# ~/.local/bin/hf). Run `hf auth login` once before running this script.
#
# Usage:
#   ./deploy/upload-data.sh              # full sync
#   ./deploy/upload-data.sh --dry-run    # show what would upload

set -euo pipefail

HF_USER="jameshyojaelee"
DATA_NAME="masld-atlas-data"
DATA_REPO="${HF_USER}/${DATA_NAME}"

PORTAL_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_SOURCE="$(cd "${PORTAL_DIR}/../masld-atlas-v2/public/data" 2>/dev/null && pwd || true)"

HF_BIN="${HF_BIN:-$HOME/.local/bin/hf}"

if [[ ! -x "$HF_BIN" ]]; then
  echo "ERROR: hf CLI not found at $HF_BIN." >&2
  echo "       Install with: pip install --user 'huggingface_hub>=1.11' 'typer>=0.18'" >&2
  exit 1
fi

if [[ -z "$DATA_SOURCE" || ! -d "$DATA_SOURCE" ]]; then
  echo "ERROR: data source not found." >&2
  echo "       Expected: ${PORTAL_DIR}/../masld-atlas-v2/public/data" >&2
  exit 1
fi

if ! "$HF_BIN" auth whoami >/dev/null 2>&1; then
  echo "ERROR: not logged in to Hugging Face. Run: hf auth login --add-to-git-credential" >&2
  exit 1
fi

echo "→ Data source: $DATA_SOURCE"
echo "→ Size:        $(du -sh "$DATA_SOURCE" | awk '{print $1}')"
echo "→ Repo:        https://huggingface.co/datasets/${DATA_REPO}"

if [[ "${1:-}" == "--dry-run" ]]; then
  echo "→ DRY RUN — showing file counts by subdir:"
  ( cd "$DATA_SOURCE" && find . -maxdepth 2 -type d | head -20 )
  echo "   Files total: $(find "$DATA_SOURCE" -type f | wc -l)"
  exit 0
fi

# Create the dataset repo (no-op if it already exists).
"$HF_BIN" repos create "$DATA_NAME" --repo-type dataset --no-private 2>&1 \
  | grep -Ev "(already exists|already created|409 Conflict)" || true

# Compress the dataset into a single tarball to bypass HF's 10k files/dir limit
ARCHIVE_PATH="/tmp/masld-atlas-data.tar.gz"
echo "→ Compressing dataset (this may take a minute)..."
(cd "$DATA_SOURCE" && tar -czf "$ARCHIVE_PATH" .)

# Upload the tarball
"$HF_BIN" upload "$DATA_REPO" "$ARCHIVE_PATH" "data.tar.gz" \
  --repo-type dataset \
  --commit-message "Sync atlas data $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Done. Browse at: https://huggingface.co/datasets/${DATA_REPO}/tree/main"
