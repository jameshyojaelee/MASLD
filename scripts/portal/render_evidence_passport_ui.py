#!/usr/bin/env python3
"""Render manifest-bound local review images for the Plan 50 candidate UI."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import shutil
import struct
import subprocess
import tempfile
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Mapping, Sequence
from urllib.parse import urlencode


SCRIPT_PATH = Path(__file__).resolve()
LOGICAL_PRODUCER_ID = "scripts/portal/render_evidence_passport_ui.py"
UI_RELATIVE_ROOT = Path("portal_candidate")
REVIEW_MANIFEST_COLUMNS = (
    "review_id",
    "relative_path",
    "view",
    "gene_symbol",
    "width_px",
    "height_px",
    "sha256",
    "bytes",
    "renderer",
    "renderer_sha256",
    "browser_label",
    "browser_version",
    "browser_sha256",
    "mesa_module",
    "source_index_sha256",
    "ui_contract_sha256",
    "render_query",
    "visible_fields",
    "viewport_check",
    "max_required_bottom_px",
)
DEFAULT_CHROME = Path(
    "/gpfs/commons/home/jameslee/.cache/puppeteer/chrome/"
    "linux-150.0.7871.24/chrome-linux64/chrome"
)
EXPECTED_MESA_MODULE = "Mesa/23.1.9-GCCcore-13.2.0"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


class PassportUIRenderError(RuntimeError):
    """Raised when the browser review artifact cannot be proven."""


class _ReviewContractParser(HTMLParser):
    """Extract the browser-computed viewport contract from serialized DOM."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self._inside = False
        self._chunks: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "pre" and dict(attrs).get("id") == "review-viewport-contract":
            self._inside = True

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "pre" and self._inside:
            self._inside = False

    def handle_data(self, data: str) -> None:
        if self._inside:
            self._chunks.append(data)

    @property
    def payload(self) -> str:
        return html.unescape("".join(self._chunks)).strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(REVIEW_MANIFEST_COLUMNS),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != PNG_SIGNATURE or header[12:16] != b"IHDR":
        raise PassportUIRenderError(f"not a valid PNG with an IHDR header: {path}")
    return struct.unpack(">II", header[16:24])


def _chrome_path(explicit: Path | None) -> Path:
    if explicit is not None:
        candidate = explicit
    elif os.environ.get("PASSPORT_CHROMIUM"):
        candidate = Path(os.environ["PASSPORT_CHROMIUM"])
    elif DEFAULT_CHROME.is_file():
        candidate = DEFAULT_CHROME
    else:
        located = shutil.which("google-chrome") or shutil.which("chromium")
        if not located:
            raise PassportUIRenderError(
                "no approved Chrome/Chromium executable was found; set PASSPORT_CHROMIUM"
            )
        candidate = Path(located)
    candidate = candidate.resolve()
    if not candidate.is_file() or candidate.is_symlink() or not os.access(candidate, os.X_OK):
        raise PassportUIRenderError(f"unsafe or non-executable browser: {candidate}")
    return candidate


def _browser_version(chrome: Path) -> str:
    result = subprocess.run(
        [str(chrome), "--version"],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    version = result.stdout.strip() or result.stderr.strip()
    if not version:
        raise PassportUIRenderError("browser returned an empty version string")
    return version


def _require_mesa(expected: str) -> None:
    loaded = os.environ.get("LOADEDMODULES", "")
    if expected not in loaded.split(":"):
        raise PassportUIRenderError(
            f"required renderer module is not loaded: {expected}; LOADEDMODULES={loaded!r}"
        )


def _browser_command(chrome: Path, width: int, height: int) -> list[str]:
    """Return the minimal Chrome command proven on the project compute node."""

    return [
        str(chrome),
        "--headless=new",
        "--no-sandbox",
        "--disable-gpu",
        "--hide-scrollbars",
        f"--window-size={width},{height}",
    ]


def _viewport_contract(
    chrome: Path,
    render_uri: str,
    *,
    width: int,
    height: int,
    required_fields: Sequence[str],
) -> tuple[str, str]:
    command = [*_browser_command(chrome, width, height), "--dump-dom", render_uri]
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        raise PassportUIRenderError(
            "browser DOM audit failed: " + str(getattr(exc, "stderr", ""))
        ) from exc
    parser = _ReviewContractParser()
    parser.feed(result.stdout)
    if not parser.payload:
        raise PassportUIRenderError("browser DOM lacks review-viewport-contract")
    try:
        payload = json.loads(parser.payload)
    except json.JSONDecodeError as exc:
        raise PassportUIRenderError("invalid browser viewport contract JSON") from exc
    expected = list(required_fields)
    if payload.get("required") != expected:
        raise PassportUIRenderError(
            f"viewport required fields drift: {payload.get('required')!r} != {expected!r}"
        )
    measurements = payload.get("measurements")
    if not isinstance(measurements, dict) or set(measurements) != set(expected):
        raise PassportUIRenderError(
            f"viewport fields missing: {sorted(set(expected) - set(measurements or {}))}"
        )
    bottoms: list[float] = []
    for field in expected:
        value = measurements[field]
        if not isinstance(value, dict) or set(value) != {"top", "bottom"}:
            raise PassportUIRenderError(f"malformed viewport measurement: {field}")
        try:
            top = float(value["top"])
            bottom = float(value["bottom"])
        except (TypeError, ValueError) as exc:
            raise PassportUIRenderError(
                f"nonnumeric viewport measurement: {field}={value!r}"
            ) from exc
        if top < 0 or bottom <= top or bottom > height:
            raise PassportUIRenderError(
                f"required field is outside {width}x{height} review viewport: "
                f"{field}=({top},{bottom})"
            )
        bottoms.append(bottom)
    return ";".join(expected), f"{max(bottoms):.2f}"


def render_ui(
    bundle: Path,
    *,
    chrome_path: Path | None = None,
    mesa_module: str = EXPECTED_MESA_MODULE,
) -> dict[str, Any]:
    bundle = bundle.resolve()
    if not bundle.is_dir() or bundle.is_symlink():
        raise PassportUIRenderError(f"unsafe candidate bundle root: {bundle}")
    ui_root = bundle / UI_RELATIVE_ROOT
    index_path = ui_root / "index.html"
    contract_path = ui_root / "ui_contract.json"
    for path in (index_path, contract_path):
        if not path.is_file() or path.is_symlink():
            raise PassportUIRenderError(f"UI input is missing or symlinked: {path}")
    try:
        ui_contract = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PassportUIRenderError(f"invalid UI contract: {contract_path}") from exc
    required_views = ui_contract.get("required_views")
    if required_views != ["overview", "THRB", "HKDC1", "GLP1R", "MTARC1"]:
        raise PassportUIRenderError(f"unexpected required review views: {required_views!r}")
    width = int(ui_contract.get("render_width_px", 0))
    height = int(ui_contract.get("render_height_px", 0))
    if (width, height) != (1440, 1200):
        raise PassportUIRenderError(f"unexpected render dimensions: {(width, height)}")

    _require_mesa(mesa_module)
    chrome = _chrome_path(chrome_path)
    browser_version = _browser_version(chrome)
    browser_hash = sha256_file(chrome)
    review_root = ui_root / "review"
    review_root.mkdir(parents=True, exist_ok=True)
    if review_root.is_symlink():
        raise PassportUIRenderError(f"review output root is symlinked: {review_root}")

    rows = []
    source_index_sha = sha256_file(index_path)
    ui_contract_sha = sha256_file(contract_path)
    for view in required_views:
        gene = "" if view == "overview" else view
        query_values = {"review": "1"} if not gene else {"gene": gene, "review": "1"}
        query = f"?{urlencode(query_values)}"
        render_uri = index_path.as_uri() + query
        destination = review_root / f"{view}.png"
        descriptor, temporary_text = tempfile.mkstemp(
            prefix=f".{view}.", suffix=".tmp.png", dir=review_root
        )
        os.close(descriptor)
        temporary = Path(temporary_text)
        temporary.unlink()
        command = [
            *_browser_command(chrome, width, height),
            f"--screenshot={temporary}",
            render_uri,
        ]
        try:
            subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            temporary.unlink(missing_ok=True)
            raise PassportUIRenderError(
                f"browser rendering failed for {view}: {getattr(exc, 'stderr', '')}"
            ) from exc
        if not temporary.is_file() or temporary.is_symlink():
            raise PassportUIRenderError(f"browser did not create {temporary}")
        observed_dimensions = png_dimensions(temporary)
        if observed_dimensions != (width, height):
            temporary.unlink(missing_ok=True)
            raise PassportUIRenderError(
                f"{view} dimensions {observed_dimensions} != {(width, height)}"
            )
        os.replace(temporary, destination)
        required_fields = (
            ("gene_grain_boundary", "native_grain_navigation")
            if not gene
            else (
                "claim", "limitation", "call", "testability",
                "testability_reason", "provenance", "source_dependence",
                "grain_boundary", "next_experiment", "falsifier",
            )
        )
        visible_fields, max_bottom = _viewport_contract(
            chrome,
            render_uri,
            width=width,
            height=height,
            required_fields=required_fields,
        )
        relative = destination.relative_to(bundle).as_posix()
        rows.append(
            {
                "review_id": f"passport_ui_review:{view}",
                "relative_path": relative,
                "view": "overview" if not gene else "gene",
                "gene_symbol": gene,
                "width_px": width,
                "height_px": height,
                "sha256": sha256_file(destination),
                "bytes": destination.stat().st_size,
                "renderer": LOGICAL_PRODUCER_ID,
                "renderer_sha256": sha256_file(SCRIPT_PATH),
                "browser_label": "pinned_puppeteer_chrome",
                "browser_version": browser_version,
                "browser_sha256": browser_hash,
                "mesa_module": mesa_module,
                "source_index_sha256": source_index_sha,
                "ui_contract_sha256": ui_contract_sha,
                "render_query": query,
                "visible_fields": visible_fields,
                "viewport_check": "pass",
                "max_required_bottom_px": max_bottom,
            }
        )
    manifest_path = review_root / "review_manifest.tsv"
    _atomic_write_tsv(manifest_path, rows)
    return {
        "review_manifest": manifest_path.as_posix(),
        "n_views": len(rows),
        "browser_version": browser_version,
        "dimensions": [width, height],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--chrome", type=Path)
    parser.add_argument("--mesa-module", default=EXPECTED_MESA_MODULE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = render_ui(
        args.bundle, chrome_path=args.chrome, mesa_module=args.mesa_module
    )
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
