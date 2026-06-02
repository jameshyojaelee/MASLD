"""Runtime data bootstrap for deployed hosts.

Local dev: no-op — ``data/paths.py`` resolves the sibling Next.js data
directory.

Deployed host: if ``MASLD_DATA_REPO_ID`` is set, snapshot that Hugging Face
repo (dataset or model) to a cache directory and export ``MASLD_DATA_DIR``
so every downstream loader reads from the fetched copy.

Env vars
--------
``MASLD_DATA_REPO_ID``   e.g. ``your-user/masld-atlas-data`` — required to trigger a fetch.
``MASLD_DATA_REPO_TYPE`` ``dataset`` (default) or ``model``.
``MASLD_DATA_REVISION``  Optional git revision / branch / tag; defaults to the repo's default branch.
``MASLD_DATA_DIR``       If already set, we honour it and skip the fetch entirely.
``HF_TOKEN``             Required only for private repos.

The cache is placed at ``/data/masld-atlas-data`` when ``/data`` is
writable (Hugging Face Space with persistent storage), else
``~/.cache/masld-atlas-data``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def _cache_root() -> Path:
    persistent = Path("/data")
    if persistent.is_dir() and os.access(persistent, os.W_OK):
        return persistent / "masld-atlas-data"
    return Path.home() / ".cache" / "masld-atlas-data"


def ensure_data() -> Path | None:
    """Ensure the data directory is populated; return its path, or ``None``
    if no remote repo is configured (local dev)."""
    if os.environ.get("MASLD_DATA_DIR"):
        return Path(os.environ["MASLD_DATA_DIR"]).expanduser().resolve()

    repo_id = os.environ.get("MASLD_DATA_REPO_ID")
    if not repo_id:
        return None

    try:
        from huggingface_hub import snapshot_download
    except ImportError:
        print(
            "[bootstrap] huggingface_hub is not installed; cannot fetch remote data.",
            file=sys.stderr,
        )
        return None

    repo_type = os.environ.get("MASLD_DATA_REPO_TYPE", "dataset")
    revision = os.environ.get("MASLD_DATA_REVISION")
    token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")

    target = _cache_root()
    target.mkdir(parents=True, exist_ok=True)

    sentinel = target / ".snapshot_complete"
    if not sentinel.exists():
        print(f"[bootstrap] Fetching {repo_id} → {target}", file=sys.stderr)
        # Only fetch the tarball — the repo may contain legacy individual
        # files from earlier upload attempts, and pulling all 10K+ of them on
        # every cold boot would add 5-10 minutes to startup for zero benefit
        # (the tarball contains everything we need).
        snapshot_download(
            repo_id=repo_id,
            repo_type=repo_type,
            revision=revision,
            local_dir=str(target),
            token=token,
            allow_patterns=["data.tar.gz"],
        )

        tar_path = target / "data.tar.gz"
        if tar_path.exists():
            import tarfile
            print(f"[bootstrap] Extracting {tar_path}...", file=sys.stderr)
            with tarfile.open(tar_path, "r:gz") as tar:
                tar.extractall(path=target)

        sentinel.touch()

    os.environ["MASLD_DATA_DIR"] = str(target)
    return target
