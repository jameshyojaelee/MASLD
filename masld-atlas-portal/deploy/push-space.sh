#!/usr/bin/env bash
# push-space.sh — Deploy masld-atlas-portal code to a public Hugging Face Space.
#
# Strategy: clone the Space to a throwaway build dir, rsync portal code into
# it (preserving HF's default .gitattributes for LFS), commit, and push.
# Does NOT bundle the data tree — the Space fetches it at runtime via
# data/bootstrap.py from the linked dataset repo.
#
# Idempotent: creates the Space on first run, updates it on subsequent runs.
# Run `hf auth login --add-to-git-credential` once before running this script.
#
# Usage:
#   ./deploy/push-space.sh

set -euo pipefail

HF_USER="jameshyojaelee"
SPACE_NAME="masld-atlas"
SPACE_REPO="${HF_USER}/${SPACE_NAME}"
DATA_REPO_ID="${HF_USER}/masld-atlas-data"

PORTAL_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-/tmp/masld-atlas-space-build}"

HF_BIN="${HF_BIN:-$HOME/.local/bin/hf}"

if [[ ! -x "$HF_BIN" ]]; then
  echo "ERROR: hf CLI not found at $HF_BIN." >&2
  exit 1
fi

if ! "$HF_BIN" auth whoami >/dev/null 2>&1; then
  echo "ERROR: not logged in. Run: hf auth login --add-to-git-credential" >&2
  exit 1
fi

echo "→ Portal src: $PORTAL_DIR"
echo "→ Build dir:  $BUILD_DIR"
echo "→ Space:      https://huggingface.co/spaces/${SPACE_REPO}"

# Resolve the HF token from either an env var or the CLI cache so we can push
# via HTTPS without an interactive credential prompt.
HF_TOKEN_VAL="${HF_TOKEN:-$(cat "$HOME/.cache/huggingface/token" 2>/dev/null || true)}"
if [[ -z "$HF_TOKEN_VAL" ]]; then
  echo "ERROR: no HF token found. Run: hf auth login" >&2
  exit 1
fi
SPACE_URL="https://${HF_USER}:${HF_TOKEN_VAL}@huggingface.co/spaces/${SPACE_REPO}"

# --- 1. Ensure the Space exists (idempotent) ---------------------------------
# HF deprecated "streamlit" as an API sdk value — Streamlit spaces now run
# under the "docker" SDK with a Dockerfile at the repo root.
"$HF_BIN" repos create "$SPACE_NAME" --type space --space-sdk docker --no-private 2>&1 \
  | grep -Ev "(already exists|already created|409 Conflict)" || true

# --- 2. Fresh clone of the Space ---------------------------------------------
rm -rf "$BUILD_DIR"
git clone "$SPACE_URL" "$BUILD_DIR"

# --- 3. Mirror portal code into the clone ------------------------------------
# Excludes:
#   - .git            : never overwrite the clone's git dir
#   - .gitattributes  : keep HF's default (includes *.parquet LFS rule)
#   - __pycache__     : Python bytecode
#   - data_cache      : local data fallback, not meant for deployment
#   - deploy          : deploy scripts stay in the source repo
cd "$BUILD_DIR"
rsync -a --delete \
  --exclude '.git/' \
  --exclude '.gitattributes' \
  --exclude '__pycache__/' \
  --exclude 'data_cache/' \
  --exclude 'deploy/' \
  --exclude '.streamlit/secrets.toml' \
  --exclude 'atlas.parquet' \
  "${PORTAL_DIR}/" ./

# --- 4. Commit and push ------------------------------------------------------
git add -A
if git diff --cached --quiet; then
  echo "→ No changes to deploy."
else
  git commit -m "deploy: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  git push origin main
fi

echo
echo "Space pushed. Next: set runtime env vars in the Space settings:"
echo "  MASLD_DATA_REPO_ID = ${DATA_REPO_ID}"
echo
echo "Settings URL: https://huggingface.co/spaces/${SPACE_REPO}/settings"
echo "App URL:      https://huggingface.co/spaces/${SPACE_REPO}"
