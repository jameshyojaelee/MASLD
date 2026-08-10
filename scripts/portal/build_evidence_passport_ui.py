#!/usr/bin/env python3
"""Build the local-only evidence-passport review interface.

This is a presentation adapter over the frozen Plan 50 candidate tables.  It
does not calculate, merge, or reinterpret a scientific call.  The resulting
HTML is self-contained so the review copy can be opened from an immutable
release without a server or network connection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
RENDERER_PATH = SCRIPT_PATH.with_name("render_evidence_passport_ui.py")
LOGICAL_PRODUCER_ID = "scripts/portal/build_evidence_passport_ui.py"
RENDERER_LOGICAL_PRODUCER_ID = "scripts/portal/render_evidence_passport_ui.py"
UI_RELATIVE_ROOT = Path("portal_candidate")
UI_CONTRACT_VERSION = "evidence_passport_candidate_ui_v1"
EMBEDDED_PAYLOAD_VERSION = "evidence_passport_embedded_payload_v1"
HERO_SYMBOLS = ("THRB", "HKDC1", "GLP1R", "MTARC1")
SOURCE_MANIFEST_COLUMNS = ("source_role", "relative_path", "sha256", "bytes")
SOURCE_TABLES = (
    ("gene_identity_and_claim", "passport_gene_index.parquet"),
    ("gene_evidence", "passport_evidence_long.parquet"),
    ("gene_testability_coverage", "passport_coverage_long.parquet"),
    ("next_experiment_and_falsifier", "passport_next_experiment.tsv"),
    ("domain_and_grain_boundary", "passport_domain_coverage.tsv"),
    ("program_inventory", "passport_program_index.parquet"),
    ("program_context_inventory", "passport_program_context.parquet"),
    ("assay_grain_inventory", "passport_assay_status.tsv"),
)


class PassportUIError(RuntimeError):
    """Raised when the frozen input contract cannot produce a faithful UI."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(SOURCE_MANIFEST_COLUMNS),
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


def _read(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path, engine="pyarrow")
    return pd.read_csv(path, sep="\t", dtype="string", keep_default_na=False)


def _string(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value)


def _optional_number(value: Any) -> int | float | None:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return None
    number = float(value)
    return int(number) if number.is_integer() else number


def _boolean(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = _string(value).strip().lower()
    if text in {"true", "1", "yes", "t"}:
        return True
    if text in {"false", "0", "no", "f", ""}:
        return False
    raise PassportUIError(f"invalid frozen boolean value: {value!r}")


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise PassportUIError(f"{label} lacks required columns: {missing}")


def _release_id(frames: Sequence[pd.DataFrame]) -> str:
    observed: set[str] = set()
    for frame in frames:
        if "analysis_release_id" in frame:
            observed.update(_string(value) for value in frame["analysis_release_id"])
    observed.discard("")
    if len(observed) != 1:
        raise PassportUIError(f"expected one analysis release, observed {sorted(observed)}")
    return next(iter(observed))


def _source_rows(bundle: Path) -> list[dict[str, Any]]:
    rows = []
    for source_role, relative in SOURCE_TABLES:
        path = bundle / relative
        if not path.is_file() or path.is_symlink():
            raise PassportUIError(f"UI source is missing or symlinked: {path}")
        rows.append(
            {
                "source_role": source_role,
                "relative_path": relative,
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    return rows


def _compact_evidence(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "id": _string(row["evidence_result_id"]),
        "domain": _string(row["evidence_domain"]),
        "assay": _string(row["assay"]),
        "dataset": _string(row["dataset_id"]),
        "phenotype": _string(row["phenotype"]),
        "context": _string(row["context"]),
        "grain": _string(row["biological_unit"]),
        "contrast": _string(row["contrast_or_exposure"]),
        "unit": _string(row["effect_unit"]),
        "estimate": _optional_number(row["estimate"]),
        "q": _optional_number(row["q_value"]),
        "direction": _string(row["direction"]),
        "testability": _string(row["testability_state"]),
        "testabilityReason": _string(row["testability_reason"]),
        "call": _string(row["call_state"]),
        "provenance": _string(row["provenance_state"]),
        "sourceDependent": _boolean(row["source_dependent"]),
        "sourceNode": _string(row["source_node_id"]),
        "sourceRelease": _string(row["source_release_id"]),
        "sourceRow": _string(row["source_input_row_id"]),
        "wording": _string(row["allowed_wording"]),
        "limitation": _string(row["limitation"]),
    }


def _compact_program_context(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "id": _string(row["program_context_id"]),
        "programId": _string(row["program_uid"]),
        "program": _string(row["program_label"]),
        "cellType": _string(row["cell_type"]),
        "dataset": _string(row["dataset_id"]),
        "assay": _string(row["assay"]),
        "grain": _string(row["biological_unit"]),
        "contrast": _string(row["contrast_or_exposure"]),
        "unit": _string(row["effect_unit"]),
        "estimate": _optional_number(row["estimate"]),
        "q": _optional_number(row["q_value"]),
        "direction": _string(row["direction"]),
        "testability": _string(row["testability_state"]),
        "testabilityReason": _string(row["testability_reason"]),
        "call": _string(row["call_state"]),
        "provenance": _string(row["provenance_state"]),
        "sourceDependent": _boolean(row["source_dependent"]),
        "sourceRelease": _string(row["source_release_id"]),
        "wording": _string(row["allowed_wording"]),
        "limitation": _string(row["limitation"]),
        "geneExpansionAuthorized": _boolean(row["gene_call_expansion_authorized"]),
    }


def _compact_assay_record(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "id": _string(row["assay_status_id"]),
        "grain": _string(row["result_grain"]),
        "entity": _string(row["entity_id"]),
        "assay": _string(row["assay"]),
        "dataset": _string(row["dataset_id"]),
        "status": _string(row["status"]),
        "call": _string(row["call_state"]),
        "testability": _string(row["testability_state"]),
        "estimate": _optional_number(row["estimate"]),
        "q": _optional_number(row["q_value"]),
        "unit": _string(row["effect_unit"]),
        "provenance": _string(row["provenance_state"]),
        "sourceDependent": _boolean(row["source_dependent"]),
        "sourceRelease": _string(row["source_release_id"]),
        "wording": _string(row["allowed_wording"]),
        "limitation": _string(row["limitation"]),
        "geneExpansionAuthorized": _boolean(row["gene_expansion_authorized"]),
    }


def build_payload(bundle: Path, *, fixture_mode: bool = False) -> dict[str, Any]:
    genes = _read(bundle / "passport_gene_index.parquet")
    evidence = _read(bundle / "passport_evidence_long.parquet")
    coverage = _read(bundle / "passport_coverage_long.parquet")
    experiments = _read(bundle / "passport_next_experiment.tsv")
    domains = _read(bundle / "passport_domain_coverage.tsv")
    programs = _read(bundle / "passport_program_index.parquet")
    program_context = _read(bundle / "passport_program_context.parquet")
    assay_status = _read(bundle / "passport_assay_status.tsv")

    _require_columns(
        genes,
        {
            "analysis_release_id", "passport_id", "ensembl_id", "symbol",
            "primary_evidence_class", "role_hypothesis", "allowed_wording",
            "limitation",
        },
        "gene index",
    )
    _require_columns(
        evidence,
        {
            "evidence_result_id", "passport_id", "evidence_domain", "assay",
            "dataset_id", "phenotype", "context", "biological_unit",
            "contrast_or_exposure", "effect_unit", "estimate", "q_value",
            "direction", "testability_state", "testability_reason", "call_state",
            "provenance_state", "source_dependent", "source_node_id",
            "source_release_id", "source_input_row_id", "allowed_wording",
            "limitation",
        },
        "gene evidence",
    )
    _require_columns(
        coverage,
        {"passport_id", "evidence_domain", "testability_state", "call_state"},
        "gene coverage",
    )
    _require_columns(
        experiments,
        {
            "passport_id", "biological_model", "context", "perturbation",
            "primary_readout", "falsifying_outcome", "disclaimer",
        },
        "next experiments",
    )
    _require_columns(
        domains,
        {
            "analysis_release_id", "evidence_domain", "gene_evidence_status",
            "program_context_status", "assay_status", "coverage_boundary",
            "allowed_wording", "limitation",
        },
        "domain coverage",
    )
    _require_columns(
        program_context,
        {
            "program_context_id", "program_uid", "program_label", "cell_type",
            "dataset_id", "assay", "biological_unit", "contrast_or_exposure",
            "effect_unit", "estimate", "q_value", "direction",
            "testability_state", "testability_reason", "call_state",
            "provenance_state", "source_dependent", "source_release_id",
            "allowed_wording", "limitation", "gene_call_expansion_authorized",
        },
        "program context",
    )
    _require_columns(
        assay_status,
        {
            "assay_status_id", "result_grain", "entity_id", "assay",
            "dataset_id", "status", "call_state", "testability_state",
            "estimate", "q_value", "effect_unit", "provenance_state",
            "source_dependent", "source_release_id", "allowed_wording",
            "limitation", "gene_expansion_authorized",
        },
        "assay status",
    )

    release_id = _release_id((genes, evidence, coverage, experiments, domains))
    if genes["passport_id"].duplicated().any():
        raise PassportUIError("gene index contains duplicate passport IDs")
    if experiments["passport_id"].duplicated().any():
        raise PassportUIError("next-experiment table contains duplicate passport IDs")
    if set(genes["passport_id"].astype(str)) != set(experiments["passport_id"].astype(str)):
        raise PassportUIError("next-experiment table does not match the gene universe")
    if not set(evidence["passport_id"].astype(str)).issubset(
        set(genes["passport_id"].astype(str))
    ):
        raise PassportUIError("evidence table references an unknown passport")
    if not set(coverage["passport_id"].astype(str)).issubset(
        set(genes["passport_id"].astype(str))
    ):
        raise PassportUIError("coverage table references an unknown passport")

    evidence_by_passport: dict[str, list[dict[str, Any]]] = {}
    for row in evidence.sort_values(
        ["passport_id", "evidence_domain", "assay", "evidence_result_id"],
        kind="stable",
    ).to_dict("records"):
        evidence_by_passport.setdefault(_string(row["passport_id"]), []).append(
            _compact_evidence(row)
        )
    experiment_by_passport = {
        _string(row["passport_id"]): {
            "model": _string(row["biological_model"]),
            "context": _string(row["context"]),
            "perturbation": _string(row["perturbation"]),
            "readout": _string(row["primary_readout"]),
            "falsifier": _string(row["falsifying_outcome"]),
            "disclaimer": _string(row["disclaimer"]),
        }
        for row in experiments.to_dict("records")
    }

    gene_rows = []
    ordered_genes = genes.assign(
        _symbol_key=genes["symbol"].astype(str).str.casefold(),
        _ensembl_key=genes["ensembl_id"].astype(str),
    ).sort_values(["_symbol_key", "_ensembl_key"], kind="stable")
    for row in ordered_genes.to_dict("records"):
        passport_id = _string(row["passport_id"])
        gene_evidence = evidence_by_passport.get(passport_id, [])
        gene_rows.append(
            {
                "passportId": passport_id,
                "symbol": _string(row["symbol"]),
                "ensembl": _string(row["ensembl_id"]),
                "primaryClass": _string(row["primary_evidence_class"]),
                "role": _string(row["role_hypothesis"]),
                "claim": _string(row["allowed_wording"]),
                "limitation": _string(row["limitation"]),
                "calls": sorted({item["call"] for item in gene_evidence}),
                "testability": sorted(
                    {item["testability"] for item in gene_evidence}
                ),
                "provenance": sorted({item["provenance"] for item in gene_evidence}),
                "sourceDependence": sorted(
                    {
                        "source_dependent" if item["sourceDependent"] else "independent"
                        for item in gene_evidence
                    }
                ),
                "domains": sorted({item["domain"] for item in gene_evidence}),
                "evidence": gene_evidence,
                "experiment": experiment_by_passport[passport_id],
            }
        )

    symbol_counts = genes["symbol"].astype(str).value_counts()
    hero_rows: dict[str, str | None] = {}
    fixture_fallback_passports = [row["passportId"] for row in gene_rows[:4]]
    for hero_index, symbol in enumerate(HERO_SYMBOLS):
        matches = genes.loc[genes["symbol"].astype(str) == symbol, "passport_id"]
        if len(matches) == 1 and int(symbol_counts[symbol]) == 1:
            hero_rows[symbol] = str(matches.iloc[0])
        elif fixture_mode:
            hero_rows[symbol] = fixture_fallback_passports[hero_index]
        else:
            raise PassportUIError(
                f"hero gene {symbol} must map to exactly one frozen passport"
            )

    domain_rows = []
    for row in domains.sort_values("evidence_domain", kind="stable").to_dict("records"):
        domain_rows.append(
            {
                "domain": _string(row["evidence_domain"]),
                "geneStatus": _string(row["gene_evidence_status"]),
                "programStatus": _string(row["program_context_status"]),
                "assayStatus": _string(row["assay_status"]),
                "boundary": _string(row["coverage_boundary"]),
                "wording": _string(row["allowed_wording"]),
                "limitation": _string(row["limitation"]),
            }
        )

    context_rows = [
        _compact_program_context(row)
        for row in program_context.sort_values(
            ["program_label", "dataset_id", "assay", "program_context_id"],
            kind="stable",
        ).to_dict("records")
    ]
    assay_rows = [
        _compact_assay_record(row)
        for row in assay_status.sort_values(
            ["result_grain", "entity_id", "assay_status_id"], kind="stable"
        ).to_dict("records")
    ]

    return {
        "payloadVersion": EMBEDDED_PAYLOAD_VERSION,
        "analysisReleaseId": release_id,
        "alphabeticalDefault": True,
        "fixtureMode": fixture_mode,
        "heroPassports": hero_rows,
        "grainBoundary": {
            "gene": "Gene pages show only frozen gene-grain evidence and coverage.",
            "program": (
                "Cell programs and spatial contexts remain at program grain; membership "
                "does not create a gene call."
            ),
            "assay": (
                "Functional challenge records remain at assay, class, or program grain; "
                "they do not create a gene call."
            ),
            "programs": len(programs),
            "programContexts": len(program_context),
            "assayRecords": len(assay_status),
        },
        "domains": domain_rows,
        "programContextRecords": context_rows,
        "assayGrainRecords": assay_rows,
        "genes": gene_rows,
    }


def _page(payload_json: str, release_id: str) -> str:
    escaped_release = html.escape(release_id, quote=True)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<title>MASLD evidence passports — local candidate review</title>
<style>
:root{{--ink:#17212b;--muted:#586570;--paper:#f5f2ec;--card:#fffefa;--line:#d7d2c8;--blue:#146b8c;--magenta:#a33b68;--gray:#707982;--pale:#e9f1f3;--warn:#f5e8df}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}}
header{{background:#102c3a;color:white;padding:28px max(24px,calc((100vw - 1240px)/2));border-bottom:6px solid #d4a547}}
header p{{max-width:850px;color:#dbe8ed;margin:.45rem 0 0}} header .eyebrow{{letter-spacing:.12em;text-transform:uppercase;font-size:12px;color:#e8ca80}}
h1{{margin:.2rem 0;font-size:clamp(30px,4vw,54px);line-height:1.03}} h2{{font-size:22px;margin:0 0 14px}} h3{{font-size:16px;margin:0 0 8px}}
main{{max-width:1240px;margin:0 auto;padding:24px}} .notice{{background:#fff4d8;border:1px solid #d9b75c;padding:13px 16px;border-radius:10px;margin-bottom:18px}}
.toolbar{{display:grid;grid-template-columns:minmax(240px,2fr) repeat(3,minmax(150px,1fr));gap:12px;background:var(--card);border:1px solid var(--line);padding:16px;border-radius:14px;position:sticky;top:0;z-index:2;box-shadow:0 4px 16px #1c2d3520}}
label{{display:grid;gap:5px;font-size:12px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}} input,select{{width:100%;border:1px solid #aeb5b8;border-radius:8px;padding:10px 11px;background:white;color:var(--ink);font:inherit}}
.summary{{display:flex;justify-content:space-between;align-items:center;gap:14px;margin:18px 2px 10px;color:var(--muted)}} .grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(270px,1fr));gap:12px}}
.native-nav{{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-bottom:16px}} .nav-card{{text-align:left;border:1px solid #89a8b4;border-radius:12px;padding:14px;background:#e8f1f4;color:inherit;cursor:pointer}} .nav-card strong{{display:block;font-size:17px}} .nav-card span{{color:var(--muted)}}
button.card{{text-align:left;border:1px solid var(--line);border-radius:12px;padding:15px;background:var(--card);color:inherit;cursor:pointer;min-height:138px}} button.card:hover,button.card:focus-visible{{border-color:var(--blue);box-shadow:0 0 0 3px #146b8c25;outline:none}}
.symbol{{font-size:23px;font-weight:800}} .ensembl{{font:12px ui-monospace,SFMono-Regular,monospace;color:var(--muted)}} .chips{{display:flex;gap:6px;flex-wrap:wrap;margin:10px 0}}
.chip{{font-size:11px;border-radius:999px;padding:3px 8px;border:1px solid #b9c2c5;background:var(--pale)}} .chip.call{{border-color:#c88aa5;background:#f5e8ee}} .chip.test{{border-color:#77a4b6;background:#e4f0f4}}
.detail-head{{display:flex;justify-content:space-between;gap:20px;align-items:flex-start}} .back{{border:1px solid #9da9ad;background:white;border-radius:8px;padding:9px 12px;cursor:pointer}}
.identity{{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:20px;margin-bottom:14px}} .claim-grid{{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:14px}}
.claim{{padding:14px;border-radius:10px;background:#e8f1f4;border-left:4px solid var(--blue)}} .limit{{padding:14px;border-radius:10px;background:var(--warn);border-left:4px solid var(--magenta)}}
.section{{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px;margin:14px 0}} .evidence{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}}
.evidence article{{border:1px solid var(--line);border-radius:10px;padding:14px}} dl{{display:grid;grid-template-columns:minmax(115px,145px) 1fr;gap:5px 12px;margin:8px 0}} dt{{font-weight:700;color:var(--muted)}} dd{{margin:0;overflow-wrap:anywhere}} .grain{{border-top:1px solid var(--line);margin-top:12px;padding-top:12px;color:var(--muted)}}
.domain-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px}} .domain{{padding:12px;border:1px solid var(--line);border-radius:9px;background:#fbfaf6}} .domain strong{{display:block;text-transform:capitalize}}
.experiment{{border-left:5px solid #d4a547}} .empty{{padding:28px;text-align:center;border:1px dashed #aab1b3;border-radius:12px;color:var(--muted)}} footer{{max-width:1240px;margin:0 auto;padding:12px 24px 30px;color:var(--muted);font-size:12px}}
.review-summary{{display:none}}
.review-layout header{{padding:12px max(18px,calc((100vw - 1240px)/2));border-bottom-width:3px}} .review-layout header h1{{font-size:30px}} .review-layout header p{{margin:.2rem 0;font-size:13px}} .review-layout main{{padding:10px 24px}} .review-layout .detail-head h1{{font-size:32px}} .review-layout .detail-head{{margin-bottom:5px}} .review-layout .identity,.review-layout .section{{padding:10px 13px;margin:7px 0}} .review-layout .claim-grid{{margin-top:7px;gap:7px}} .review-layout .claim,.review-layout .limit{{padding:8px 10px}} .review-layout .chips{{margin:5px 0}} .review-layout .experiment dl{{margin:4px 0;gap:2px 10px}} .review-layout .evidence{{display:none}} .review-layout .domain-section{{display:none}} .review-layout .grain{{margin:5px 0;padding-top:5px}} .review-layout footer{{display:none}}
.review-layout .review-summary{{display:grid;grid-template-columns:1fr 1fr;gap:5px;margin-top:7px}} .review-layout .review-summary article{{padding:5px 8px;border:1px solid var(--line);border-radius:7px;font-size:12px}} .review-layout .review-summary strong{{display:block}} .review-layout .review-summary p{{margin:2px 0}}
@media(max-width:850px){{.toolbar{{grid-template-columns:1fr 1fr;position:static}}.evidence,.claim-grid{{grid-template-columns:1fr}}}}
@media print{{header{{padding:18px 24px}}.toolbar,.back,footer{{display:none!important}}main{{max-width:none;padding:16px}}.section,.identity{{break-inside:avoid}}}}
</style>
</head>
<body data-ui-contract="{UI_CONTRACT_VERSION}" data-ready="false">
<header><div class="eyebrow">Candidate review · local only</div><h1>Evidence passports</h1><p>Source-preserving gene views with explicit testability, provenance, biological grain, limitations, and a discriminating next experiment.</p></header>
<main id="app" aria-live="polite"></main>
<footer>Frozen analysis release: <span id="release-footer">{escaped_release}</span>. Candidate review artifact; no canonical promotion is authorized.</footer>
<script id="passport-data" type="application/json">{payload_json}</script>
<script>
"use strict";
const DATA=JSON.parse(document.getElementById("passport-data").textContent);
const APP=document.getElementById("app");
const esc=(v)=>String(v??"").replace(/[&<>\"']/g,c=>({{"&":"&amp;","<":"&lt;",">":"&gt;",'\"':"&quot;","'":"&#39;"}}[c]));
const label=(v)=>String(v||"not stated").replaceAll("_"," ");
const fmt=(v)=>v===null||v===undefined||v===""?"not reported":Number.isFinite(Number(v))?Number(v).toPrecision(4):String(v);
const chips=(values,kind="")=>(values||[]).map(v=>`<span class="chip ${{kind}}">${{esc(label(v))}}</span>`).join("");
function params(){{return new URLSearchParams(location.search)}}
function goGene(passportId){{const url=new URL(location.href);url.search="";url.searchParams.set("passport",passportId);history.pushState({{}},"",url);render()}}
function goView(view){{const url=new URL(location.href);url.search="";url.searchParams.set("view",view);history.pushState({{}},"",url);render()}}
function goHome(){{const url=new URL(location.href);url.search="";history.pushState({{}},"",url);render()}}
function renderOverview(){{
  const calls=[...new Set(DATA.genes.flatMap(g=>g.calls))].sort();
  const tests=[...new Set(DATA.genes.flatMap(g=>g.testability))].sort();
  const prov=[...new Set(DATA.genes.flatMap(g=>g.provenance))].sort();
  APP.innerHTML=`<div class="notice" data-review-field="gene_grain_boundary"><strong>Interpretation boundary.</strong> Gene evidence is shown only at gene grain. Program and assay records remain separate and never create a member-gene call.</div>
  <section class="native-nav" data-review-field="native_grain_navigation" aria-label="Native-grain evidence views"><button class="nav-card" id="program-context-view"><strong>Program Context · ${{DATA.programContextRecords.length}} records</strong><span>Spatial, protein, and chromatin contexts at their native program grain.</span></button><button class="nav-card" id="assay-record-view"><strong>Assay / Class · ${{DATA.assayGrainRecords.length}} records</strong><span>Functional challenge status at assay, class, or program grain.</span></button></section>
  <section class="toolbar" aria-label="Gene passport filters"><label>Search gene or Ensembl ID<input id="search" type="search" autocomplete="off" placeholder="e.g. HKDC1 or ENSG…"></label>
  <label>Call<select id="call"><option value="">All calls</option>${{calls.map(v=>`<option value="${{esc(v)}}">${{esc(label(v))}}</option>`).join("")}}</select></label>
  <label>Testability<select id="test"><option value="">All states</option>${{tests.map(v=>`<option value="${{esc(v)}}">${{esc(label(v))}}</option>`).join("")}}</select></label>
  <label>Provenance<select id="prov"><option value="">All provenance</option>${{prov.map(v=>`<option value="${{esc(v)}}">${{esc(label(v))}}</option>`).join("")}}</select></label></section>
  <div class="summary"><span id="match-count"></span><span>Alphabetical by gene symbol · first 200 matches shown</span></div><section id="gene-grid" class="grid"></section>
  <section class="section"><h2>Domain and grain boundary</h2><div class="domain-grid">${{DATA.domains.map(d=>`<article class="domain"><strong>${{esc(label(d.domain))}}</strong><span>${{esc(label(d.boundary))}}</span><p>${{esc(d.wording)}}</p></article>`).join("")}}</div><p class="grain">${{esc(DATA.grainBoundary.program)}} ${{esc(DATA.grainBoundary.assay)}}</p></section>`;
  const update=()=>{{const q=document.getElementById("search").value.trim().toLowerCase(),call=document.getElementById("call").value,test=document.getElementById("test").value,prov=document.getElementById("prov").value;
    const found=DATA.genes.filter(g=>(!q||g.symbol.toLowerCase().includes(q)||g.ensembl.toLowerCase().includes(q))&&(!call||g.calls.includes(call))&&(!test||g.testability.includes(test))&&(!prov||g.provenance.includes(prov)));
    document.getElementById("match-count").textContent=`${{found.length.toLocaleString()}} matching genes`;
    document.getElementById("gene-grid").innerHTML=found.slice(0,200).map(g=>`<button class="card" data-passport="${{esc(g.passportId)}}" aria-label="Open ${{esc(g.symbol)}} (${{esc(g.ensembl)}}) evidence passport"><span class="symbol">${{esc(g.symbol)}}</span><div class="ensembl">${{esc(g.ensembl)}}</div><div class="chips"><span class="chip">${{esc(label(g.primaryClass))}}</span>${{chips(g.calls,"call")}}${{chips(g.testability,"test")}}</div><div>${{esc(g.claim)}}</div></button>`).join("")||`<div class="empty">No gene matches these filters.</div>`;
    document.querySelectorAll("button.card").forEach(b=>b.addEventListener("click",()=>goGene(b.dataset.passport)));
  }};
  document.getElementById("program-context-view").addEventListener("click",()=>goView("program-context"));
  document.getElementById("assay-record-view").addEventListener("click",()=>goView("assay-records"));
  ["search","call","test","prov"].forEach(id=>document.getElementById(id).addEventListener("input",update));update();
}}
function renderProgramContexts(){{
  APP.innerHTML=`<div class="detail-head"><div><div class="eyebrow">Native program grain</div><h1>Program Context</h1><p>These records remain program-level observations. They are not projected onto member genes.</p></div><button class="back" onclick="goHome()">← Gene index</button></div><section class="section"><div class="evidence">${{DATA.programContextRecords.map(r=>`<article><div class="chips"><span class="chip">program grain</span><span class="chip">${{esc(label(r.call))}}</span><span class="chip test">${{esc(label(r.testability))}}</span></div><h3>${{esc(r.program)}}</h3><dl><dt>Cell type</dt><dd>${{esc(r.cellType)}}</dd><dt>Dataset</dt><dd>${{esc(r.dataset)}}</dd><dt>Assay</dt><dd>${{esc(r.assay)}}</dd><dt>Biological unit</dt><dd>${{esc(r.grain)}}</dd><dt>Contrast</dt><dd>${{esc(r.contrast)}}</dd><dt>Effect</dt><dd>${{esc(fmt(r.estimate))}} ${{esc(r.unit)}}</dd><dt>Adjusted p</dt><dd>${{esc(fmt(r.q))}}</dd><dt>Provenance</dt><dd>${{esc(label(r.provenance))}} · ${{r.sourceDependent?"source dependent":"independent"}}</dd></dl><div class="claim">${{esc(r.wording)}}</div><div class="limit">${{esc(r.limitation)}}</div><p class="grain">Member-gene expansion authorized: ${{r.geneExpansionAuthorized?"yes":"no"}}</p></article>`).join("")}}</div></section>`;
}}
function renderAssayRecords(){{
  APP.innerHTML=`<div class="detail-head"><div><div class="eyebrow">Native assay / class / program grain</div><h1>Assay and class status</h1><p>Functional records remain at their frozen result grain and do not create gene calls.</p></div><button class="back" onclick="goHome()">← Gene index</button></div><section class="section"><div class="evidence">${{DATA.assayGrainRecords.map(r=>`<article><div class="chips"><span class="chip">${{esc(label(r.grain))}} grain</span><span class="chip">${{esc(label(r.status))}}</span><span class="chip call">${{esc(label(r.call))}}</span><span class="chip test">${{esc(label(r.testability))}}</span></div><h3>${{esc(r.entity)}}</h3><dl><dt>Dataset</dt><dd>${{esc(r.dataset)}}</dd><dt>Assay</dt><dd>${{esc(r.assay)}}</dd><dt>Effect</dt><dd>${{esc(fmt(r.estimate))}} ${{esc(r.unit)}}</dd><dt>Adjusted p</dt><dd>${{esc(fmt(r.q))}}</dd><dt>Provenance</dt><dd>${{esc(label(r.provenance))}} · ${{r.sourceDependent?"source dependent":"independent"}}</dd></dl><div class="claim">${{esc(r.wording)}}</div><div class="limit">${{esc(r.limitation)}}</div><p class="grain">Gene expansion authorized: ${{r.geneExpansionAuthorized?"yes":"no"}}</p></article>`).join("")}}</div></section>`;
}}
function renderDetail(passportId,heroSymbol){{
  const heroPassport=heroSymbol?DATA.heroPassports[heroSymbol]:null;
  const requestedPassport=passportId||heroPassport;
  const gene=DATA.genes.find(g=>g.passportId===requestedPassport);
  const requestedLabel=heroSymbol||passportId||"unknown passport";
  if(!gene){{APP.innerHTML=`<section class="identity"><button class="back" onclick="goHome()">← All genes</button><h2>${{esc(requestedLabel)}}</h2><div class="empty">This identity is not present in the frozen candidate gene universe.</div></section>`;return}}
  const ev=gene.evidence.map(e=>`<article><div class="chips"><span class="chip">gene grain</span><span class="chip">${{esc(label(e.domain))}}</span><span class="chip call">${{esc(label(e.call))}}</span><span class="chip test">${{esc(label(e.testability))}}</span></div><h3>${{esc(e.assay)}}</h3><dl><dt>Dataset</dt><dd>${{esc(e.dataset)}}</dd><dt>Phenotype</dt><dd>${{esc(e.phenotype)}}</dd><dt>Context</dt><dd>${{esc(e.context)}}</dd><dt>Biological unit</dt><dd>${{esc(e.grain)}}</dd><dt>Contrast</dt><dd>${{esc(e.contrast)}}</dd><dt>Estimate</dt><dd>${{esc(fmt(e.estimate))}} ${{esc(e.unit)}}</dd><dt>Adjusted p</dt><dd>${{esc(fmt(e.q))}}</dd><dt>Direction</dt><dd>${{esc(label(e.direction))}}</dd><dt>Testability reason</dt><dd>${{esc(e.testabilityReason)}}</dd><dt>Provenance</dt><dd>${{esc(label(e.provenance))}} · ${{e.sourceDependent?"source dependent":"independent"}}</dd><dt>Source release</dt><dd>${{esc(e.sourceRelease)}}</dd></dl><div class="claim">${{esc(e.wording)}}</div><div class="limit">${{esc(e.limitation)}}</div></article>`).join("");
  const reviewEv=gene.evidence.map(e=>`<article><strong>${{esc(e.assay)}}</strong><p>${{esc(label(e.call))}} · ${{esc(label(e.testability))}}</p><p>${{esc(e.testabilityReason)}}</p><p>${{esc(label(e.provenance))}} · ${{e.sourceDependent?"source dependent":"independent"}}</p></article>`).join("");
  const x=gene.experiment;
  APP.innerHTML=`<div class="detail-head"><div><div class="eyebrow">Gene-grain passport</div><h1>${{esc(gene.symbol)}}</h1><div class="ensembl">${{esc(gene.ensembl)}}</div></div><button class="back" onclick="goHome()">← All genes</button></div>
  <section class="identity"><div class="chips"><span class="chip">${{esc(label(gene.primaryClass))}}</span><span class="chip">${{esc(label(gene.role))}}</span><span data-review-field="call">${{chips(gene.calls,"call")}}</span><span data-review-field="testability">${{chips(gene.testability,"test")}}</span><span data-review-field="provenance">${{chips(gene.provenance)}}</span><span data-review-field="source_dependence">${{chips(gene.sourceDependence)}}</span></div><div class="claim-grid"><div class="claim" data-review-field="claim"><h3>Source-authorized claim</h3>${{esc(gene.claim)}}</div><div class="limit" data-review-field="limitation"><h3>Limitation</h3>${{esc(gene.limitation)}}</div></div><div class="review-summary" data-review-field="testability_reason">${{reviewEv}}</div><p class="grain" data-review-field="grain_boundary"><strong>Grain boundary:</strong> ${{esc(DATA.grainBoundary.gene)}} ${{esc(DATA.grainBoundary.program)}} ${{esc(DATA.grainBoundary.assay)}}</p></section>
  <section class="section experiment" data-review-field="next_experiment"><h2>Next discriminating experiment</h2><dl><dt>Model</dt><dd>${{esc(x.model)}}</dd><dt>Context</dt><dd>${{esc(x.context)}}</dd><dt>Perturbation</dt><dd>${{esc(x.perturbation)}}</dd><dt>Primary readout</dt><dd>${{esc(x.readout)}}</dd><dt>Falsifying outcome</dt><dd data-review-field="falsifier">${{esc(x.falsifier)}}</dd></dl><p class="grain">${{esc(x.disclaimer)}}</p></section>
  <section class="section"><h2>Testability, call, and provenance</h2><div class="evidence">${{ev||`<div class="empty">No accepted gene-grain evidence rows.</div>`}}</div><p class="grain"><strong>Grain boundary:</strong> ${{esc(DATA.grainBoundary.gene)}} ${{esc(DATA.grainBoundary.program)}} ${{esc(DATA.grainBoundary.assay)}}</p></section>
  <section class="section domain-section"><h2>Domain boundary</h2><div class="domain-grid">${{DATA.domains.map(d=>`<article class="domain"><strong>${{esc(label(d.domain))}}</strong><span>${{esc(label(d.boundary))}}</span><p>${{esc(d.wording)}}</p></article>`).join("")}}</div></section>`;
}}
function addReviewViewportContract(){{
  const required=document.body.dataset.view==="gene"?["claim","limitation","call","testability","testability_reason","provenance","source_dependence","grain_boundary","next_experiment","falsifier"]:["gene_grain_boundary","native_grain_navigation"];
  const measurements={{}};
  for(const field of required){{const node=document.querySelector(`[data-review-field="${{field}}"]`);if(!node)continue;const rect=node.getBoundingClientRect();measurements[field]={{top:Number(rect.top.toFixed(2)),bottom:Number(rect.bottom.toFixed(2))}}}}
  const marker=document.createElement("pre");marker.id="review-viewport-contract";marker.hidden=true;marker.textContent=JSON.stringify({{required,measurements,viewportHeight:innerHeight}});document.body.appendChild(marker);
}}
function render(){{document.body.dataset.ready="false";document.body.classList.toggle("review-layout",params().get("review")==="1");const p=params(),passportId=p.get("passport"),heroSymbol=p.get("gene"),view=p.get("view");if(passportId||heroSymbol)renderDetail(passportId,heroSymbol);else if(view==="program-context")renderProgramContexts();else if(view==="assay-records")renderAssayRecords();else renderOverview();document.body.dataset.view=(passportId||heroSymbol)?"gene":(view||"overview");document.body.dataset.passport=passportId||DATA.heroPassports[heroSymbol]||"";document.body.dataset.gene=heroSymbol||"";if(params().get("review")==="1")addReviewViewportContract();document.body.dataset.ready="true"}}
addEventListener("popstate",render);render();
</script>
</body></html>"""


def build_ui(bundle: Path, *, fixture_mode: bool = False) -> dict[str, Any]:
    bundle = bundle.resolve()
    if not bundle.is_dir() or bundle.is_symlink():
        raise PassportUIError(f"unsafe candidate bundle root: {bundle}")
    if not RENDERER_PATH.is_file() or RENDERER_PATH.is_symlink():
        raise PassportUIError(f"renderer producer is missing or symlinked: {RENDERER_PATH}")
    ui_root = bundle / UI_RELATIVE_ROOT
    ui_root.mkdir(parents=True, exist_ok=True)
    if ui_root.is_symlink():
        raise PassportUIError(f"UI output root is symlinked: {ui_root}")

    source_rows = _source_rows(bundle)
    source_manifest = ui_root / "ui_source_manifest.tsv"
    _write_tsv(source_manifest, source_rows)
    payload = build_payload(bundle, fixture_mode=fixture_mode)
    payload_json = stable_json(payload).replace("<", "\\u003c")
    payload_sha256 = hashlib.sha256(payload_json.encode("utf-8")).hexdigest()
    index_path = ui_root / "index.html"
    _atomic_write(index_path, _page(payload_json, payload["analysisReleaseId"]))

    contract_path = ui_root / "ui_contract.json"
    ui_contract = {
        "contract_version": UI_CONTRACT_VERSION,
        "analysis_release_id": payload["analysisReleaseId"],
        "alphabetical_default": True,
        "embedded_payload_version": EMBEDDED_PAYLOAD_VERSION,
        "embedded_payload_sha256": payload_sha256,
        "fixture_mode": fixture_mode,
        "filters": ["search", "call", "testability", "provenance"],
        "forbidden_gene_ordering": ["effect_order", "composite_order"],
        "hero_symbols": list(HERO_SYMBOLS),
        "required_views": ["overview", *HERO_SYMBOLS],
        "source_manifest_sha256": sha256_file(source_manifest),
        "generator": LOGICAL_PRODUCER_ID,
        "generator_sha256": sha256_file(SCRIPT_PATH),
        "renderer": RENDERER_LOGICAL_PRODUCER_ID,
        "renderer_sha256": sha256_file(RENDERER_PATH),
        "render_width_px": 1440,
        "render_height_px": 1200,
        "canonical_promotion_authorized": False,
        "scientific_call_recomputed": False,
    }
    _atomic_write(contract_path, stable_json(ui_contract) + "\n")
    return {
        "index": index_path,
        "contract": contract_path,
        "source_manifest": source_manifest,
        "analysis_release_id": payload["analysisReleaseId"],
        "n_genes": len(payload["genes"]),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--fixture-mode", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(stable_json(build_ui(args.bundle, fixture_mode=args.fixture_mode)))


if __name__ == "__main__":
    main()
