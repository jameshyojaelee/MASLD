#!/usr/bin/env python3
"""Adjudicate EpiBERT native preprocessing from checkpoint-bound primary sources."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any
from xml.etree import ElementTree

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


class NativePreprocessingAuditError(RuntimeError):
    """Raised when frozen evidence or a preprocessing requirement differs."""


EXPECTED = {
    "admission": "26e12915b05c2727c45eca091f85ec89f96c205e817c6ec9b64f181a563e7577",
    "runtime": "83e42afe350081391204703e9178ac30b52b7dfc15b863ba1af2e10611342aa4",
    "primary_evidence": "4047b17e387eaf1be1401a7e80bb9470670ef9af11b13e2f4209a25133e1e6e0",
    "revision": "dae61b434f885991c29e2fde3fbf7f70d0ff9b80",
    "motif": "aaa5f79457e991d470c69efd0066b9db25e395d4c4e8bd07bc4bffa3c2d6768c",
    "motif_means": "94a2177a4a2a9bea937fec095d3e2086a0e32b33299a16563e7da1f44486181d",
    "motif_std": "8394c3d002e6bdf07b8c990d9958f18f8bb30e7f8aa7c732d3bce2056a7e7123",
}


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def digest_bytes(value: bytes) -> str:
    return sha256(value).hexdigest()


def require_artifact(root: Path, expected: str) -> dict[str, Any]:
    try:
        value = verify_frozen_tree(root)
    except ArtifactError as error:
        raise NativePreprocessingAuditError(f"frozen artifact differs: {root}") from error
    if digest(root / "ARTIFACTS.json") != expected:
        raise NativePreprocessingAuditError(f"frozen artifact identity differs: {root}")
    return dict(value)


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise NativePreprocessingAuditError(f"expected JSON object: {path}")
    return value


def normalized_text(element: ElementTree.Element) -> str:
    return " ".join(" ".join(element.itertext()).split())


def paper_evidence(xml_path: Path, acquisition: dict[str, Any]) -> dict[str, Any]:
    root = ElementTree.parse(xml_path).getroot()
    full_text = normalized_text(root)
    required = {
        "fpm_label": "fragments-per-million (FPM) normalized",
        "background_50000": "randomly selecting 50000 peaks",
        "foreground_50k": "selected 50k random peaks",
        "motif_count": "693 consensus motif models",
        "min_max": "min-max scaled",
        "tn5_forward": "+4bp from the 5′ end",
        "tn5_reverse": "-5bp from the 3′ end",
    }
    absent = [name for name, value in required.items() if value not in full_text]
    if absent:
        raise NativePreprocessingAuditError(
            "paper method evidence differs: " + ",".join(absent)
        )
    members = acquisition["paper"]["article_evidence"]["supplementary_members"]
    if [row["href"] for row in members] != [
        "mmc1.pdf",
        "mmc2.xlsx",
        "mmc3.xlsx",
        "mmc4.pdf",
        "mmc5.pdf",
    ]:
        raise NativePreprocessingAuditError("paper supplement roster differs")
    return {
        "efetch_sha256": digest(xml_path),
        "method_facts": {
            "atac_signal_label": "fragments per million, with f defined as total fragments divided by 1e6",
            "tn5_shifts_bp": {"forward": 4, "reverse": -5},
            "motif_background": "union across pretraining datasets; 50000 random peaks present in at least 50 percent of datasets",
            "motif_foreground": "50000 random peaks per dataset",
            "motif_count": 693,
            "motif_summary": "SEA enrichment q-value per motif",
            "motif_normalization": "per-vector min-max scaling",
        },
        "supplementary_members": members,
        "supplement_binary_content_audited": False,
        "supplement_binary_blocker": "NCBI article-member endpoints are protected download surfaces; the acquisition freezes their JATS titles and exact PMC bindings but does not bypass the protection.",
    }


def run_git(repository: Path, *arguments: str, binary: bool = False) -> str | bytes:
    result = subprocess.run(
        ["/usr/bin/git", *arguments],
        cwd=repository,
        check=True,
        capture_output=True,
        text=not binary,
    )
    return result.stdout


def notebook_source(repository: Path, revision: str, path: str) -> str:
    raw = run_git(repository, "show", f"{revision}:{path}")
    if not isinstance(raw, str):
        raise NativePreprocessingAuditError("notebook source type differs")
    value = json.loads(raw)
    return "\n".join(
        "".join(cell.get("source", [])) for cell in value.get("cells", [])
    )


def repository_evidence(bundle: Path) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="epibert-audit-") as directory:
        repository = Path(directory) / "repository"
        subprocess.run(
            ["/usr/bin/git", "clone", "--quiet", str(bundle), str(repository)],
            check=True,
        )
        tagged = str(run_git(repository, "rev-parse", "v1.0.0^{}")).strip()
        head = str(run_git(repository, "rev-parse", "HEAD")).strip()
        count = int(str(run_git(repository, "rev-list", "--all", "--count")).strip())
        if tagged != EXPECTED["revision"] or count != 65:
            raise NativePreprocessingAuditError("repository history identity differs")

        example_commit = "53b290601cfcb4986ebb534a0559777d859cc29e"
        correction_commit = "097a0185beb70b3b584720265d461605a06e2fe6"
        example = notebook_source(repository, example_commit, "example_usage/data_processing.ipynb")
        corrected = notebook_source(repository, correction_commit, "example_usage/data_processing.ipynb")
        example_required = [
            "head -n 50000",
            "$2+$10-64",
            "$2+$10+64",
            "--thresh 1.0",
        ]
        corrected_required = [
            "head -n 50000",
            "$2+$10-64",
            "$2+$10+64",
            "--thresh 50000.0",
            "reads per 20 million",
            "-scale 0.1797",
        ]
        if any(value not in example for value in example_required) or any(
            value not in corrected for value in corrected_required
        ):
            raise NativePreprocessingAuditError("post-release example evidence differs")
        if subprocess.run(
            ["/usr/bin/git", "merge-base", "--is-ancestor", EXPECTED["revision"], example_commit],
            cwd=repository,
        ).returncode != 0:
            raise NativePreprocessingAuditError("post-release example ancestry differs")

        background = run_git(
            repository,
            "show",
            f"{example_commit}:example_usage/all_peaks_merged.counts.shared.centered.bed",
            binary=True,
        )
        motif_example = run_git(
            repository,
            "show",
            f"{example_commit}:example_usage/consensus_pwms.meme",
            binary=True,
        )
        sea_example = run_git(
            repository,
            "show",
            f"{correction_commit}:example_usage/ENCFF135AEX.motifs.tsv",
            binary=True,
        )
        assert isinstance(background, bytes)
        assert isinstance(motif_example, bytes)
        assert isinstance(sea_example, bytes)
        pinned_paths = str(
            run_git(repository, "ls-tree", "-r", "--name-only", EXPECTED["revision"])
        ).splitlines()
        invocation_candidates = [
            path
            for path in pinned_paths
            if Path(path).suffix.lower() in {".json", ".yaml", ".yml", ".inputs"}
        ]
        return {
            "bundle_sha256": digest(bundle),
            "head": head,
            "tag_v1_0_0": tagged,
            "commit_count_all_refs": count,
            "pinned_invocation_candidates": invocation_candidates,
            "post_release_example": {
                "commit": example_commit,
                "author_date": str(
                    run_git(repository, "show", "-s", "--format=%aI", example_commit)
                ).strip(),
                "top_peaks": 50000,
                "peak_window_bp": 128,
                "sea_threshold": 1.0,
                "background_peak_rows": len(background.splitlines()),
                "background_peaks_sha256": digest_bytes(background),
                "motif_file_sha256": digest_bytes(motif_example),
            },
            "post_release_correction": {
                "commit": correction_commit,
                "author_date": str(
                    run_git(repository, "show", "-s", "--format=%aI", correction_commit)
                ).strip(),
                "sea_threshold": 50000.0,
                "signal_scale_example": "20 divided by 111.279 million fragments equals 0.1797",
                "sea_output_sha256": digest_bytes(sea_example),
                "sea_output_lines": len(sea_example.splitlines()),
            },
            "checkpoint_binding": "none: both example commits postdate v1.0.0 and the December 2024 checkpoint archive",
        }


def motif_roster(path: Path) -> list[tuple[str, str]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MOTIF "):
            fields = line.split(maxsplit=2)
            rows.append((fields[1], fields[2] if len(fields) > 2 else ""))
    return rows


def comma_values(path: Path) -> list[float]:
    return [float(value) for value in path.read_text(encoding="utf-8").split(",")]


def pinned_source_evidence(source: Path) -> tuple[dict[str, Any], list[tuple[str, str]]]:
    signal = source / "data_processing/create_signal_tracks/fragments_to_bed_scores.wdl"
    producer = source / "data_processing/create_signal_tracks/bam_to_bed_ATAC.wdl"
    sea = source / "data_processing/motif_enrichment/meme_run_sea.wdl"
    motif = source / "data_processing/motif_enrichment/consensus_pwms_vierstra.meme"
    means = source / "src/motif_means_norm.tsv"
    std = source / "src/motif_std_norm.tsv"
    training = source / "training_utils_atac_pretrain.py"
    analysis = source / "analysis/interval_and_plotting_utilities.py"
    signal_text = signal.read_text(encoding="utf-8")
    producer_text = producer.read_text(encoding="utf-8")
    sea_text = sea.read_text(encoding="utf-8")
    training_text = training.read_text(encoding="utf-8")
    analysis_text = analysis.read_text(encoding="utf-8")
    required = [
        ("20.0 / scale_factor", signal_text),
        ("$2-5,$3+5", signal_text),
        ("print $1 / 1000000.0", producer_text),
        ("Float thresh", sea_text),
        ("Int top_peaks", sea_text),
        ("Int half_peak_width", sea_text),
        ("File background_peaks", sea_text),
        ("memesuite/memesuite:5.4.1", sea_text),
        ("(motif_activity - min_val) / (max_val - min_val)", training_text),
        ("(motif_activity - motif_means) / (motif_std + 1.0e-06)", training_text),
        ("processed_data.sort(key=lambda x: x[0])", analysis_text),
    ]
    if any(value not in text for value, text in required):
        raise NativePreprocessingAuditError("pinned preprocessing source differs")
    if any(
        re.search(rf"(?:Float|Int|File)\s+{name}\s*=", sea_text)
        for name in ("thresh", "top_peaks", "half_peak_width", "background_peaks")
    ):
        raise NativePreprocessingAuditError("SEA required input unexpectedly has a default")
    rows = motif_roster(motif)
    expected_ids = [f"AC{index:04d}:{rows[index - 1][0].split(':', 1)[1]}" for index in range(1, 694)]
    observed_ids = [row[0] for row in rows]
    if (
        digest(motif) != EXPECTED["motif"]
        or len(rows) != 693
        or observed_ids != expected_ids
        or observed_ids != sorted(observed_ids)
        or digest(means) != EXPECTED["motif_means"]
        or digest(std) != EXPECTED["motif_std"]
        or len(comma_values(means)) != 693
        or len(comma_values(std)) != 693
    ):
        raise NativePreprocessingAuditError("motif roster or normalization arrays differ")
    return (
        {
            "signal_wdl_sha256": digest(signal),
            "fragment_producer_wdl_sha256": digest(producer),
            "sea_wdl_sha256": digest(sea),
            "motif_file_sha256": digest(motif),
            "motif_means_sha256": digest(means),
            "motif_std_sha256": digest(std),
            "motif_count": len(rows),
            "motif_first": observed_ids[0],
            "motif_last": observed_ids[-1],
            "motif_order": "lexicographic sort of SEA field 3; identical to AC0001 through AC0693 order",
            "atac_scale_formula": "20.0 / scale_factor",
            "adjacent_producer_scale_units": "accepted BEDPE fragment rows divided by 1e6",
            "signal_expansion": "each Tn5 insertion interval is expanded by 5 bp on each side",
            "sea_required_inputs_without_defaults": [
                "top_peaks",
                "half_peak_width",
                "background_peaks",
                "thresh",
            ],
            "sea_container": "memesuite/memesuite:5.4.1",
            "motif_training_normalization": "per-vector min-max",
            "motif_validation_normalization": "released validation deserializer uses frozen means and standard deviations",
            "motif_analysis_normalization": "per-vector min-max",
        },
        rows,
    )


def archive_evidence(admission: Path) -> dict[str, Any]:
    models = load(admission / "inventories/epibert_models.inventory.json")
    code = load(admission / "inventories/epibert_code.inventory.json")
    model_regular = [row["path"] for row in models["members"] if row["type"] == "regular_file"]
    configurations = [
        path
        for path in model_regular
        if Path(path).suffix.lower() in {".json", ".yaml", ".yml", ".wdl", ".ipynb", ".tsv"}
    ]
    if models["member_count"] != 18 or configurations:
        raise NativePreprocessingAuditError("checkpoint archive configuration inventory differs")
    code_paths = [row["path"] for row in code["members"] if row["type"] == "regular_file"]
    invocation_configs = [
        path
        for path in code_paths
        if Path(path).suffix.lower() in {".json", ".yaml", ".yml", ".inputs"}
    ]
    if invocation_configs:
        raise NativePreprocessingAuditError("pinned code archive contains an invocation config")
    return {
        "model_archive_member_count": models["member_count"],
        "model_archive_regular_files": model_regular,
        "model_archive_preprocessing_configs": configurations,
        "code_archive_member_count": code["member_count"],
        "code_archive_invocation_configs": invocation_configs,
    }


def audit(args: argparse.Namespace) -> tuple[dict[str, Any], list[tuple[str, str]]]:
    require_artifact(args.admission, EXPECTED["admission"])
    require_artifact(args.runtime, EXPECTED["runtime"])
    require_artifact(args.primary_evidence, EXPECTED["primary_evidence"])
    acquisition = load(args.primary_evidence / "sources/acquisition_receipt.json")
    if (
        acquisition["repository"]["pinned_revision"] != EXPECTED["revision"]
        or acquisition["paper"]["evidence_mode"]
        != "ncbi_efetch_jats_plus_article_metadata"
    ):
        raise NativePreprocessingAuditError("primary acquisition contract differs")
    source = next(args.runtime.glob(f"source/EpiBERT-{EXPECTED['revision']}"))
    pinned, roster = pinned_source_evidence(source)
    paper = paper_evidence(
        args.primary_evidence / "sources/paper/PMC11872434.efetch.xml",
        acquisition,
    )
    history = repository_evidence(
        args.primary_evidence / "sources/repository/EpiBERT-all-refs.bundle"
    )
    archives = archive_evidence(args.admission)

    fields = {
        "atac_tn5_shifts": {
            "status": "authoritative",
            "value": {"forward_bp": 4, "reverse_bp": -5},
            "evidence": "paper STAR Methods and official post-release example agree",
        },
        "atac_scale_factor_units": {
            "status": "authoritative_units_but_transform_conflicted",
            "value": "total accepted fragments divided by 1e6",
            "blocker": "The paper describes FPM as 1/f, while pinned WDL and the official post-release example use 20/f after insertion expansion and call it reads per 20 million. No checkpoint-bound invocation selects one interpretation.",
        },
        "atac_scale_numerator": {
            "status": "unresolved_checkpoint_native_transform",
            "candidates": [1.0, 20.0],
            "blocker": "paper prose and executable WDL disagree",
        },
        "motif_top_peaks": {
            "status": "unresolved",
            "candidate": 50000,
            "blocker": "The paper says 50000 random foreground peaks, while the later example takes the top 50000 by score; pinned WDL has no default or invocation.",
        },
        "motif_peak_window": {
            "status": "unresolved",
            "candidate_bp": 128,
            "blocker": "The later example centers 128 bp on the narrowPeak summit, while pinned WDL trims both peak boundaries by an unbound half_peak_width.",
        },
        "motif_background_peaks": {
            "status": "unresolved",
            "paper_recipe": paper["method_facts"]["motif_background"],
            "post_release_example_rows": history["post_release_example"]["background_peak_rows"],
            "post_release_example_sha256": history["post_release_example"]["background_peaks_sha256"],
            "blocker": "The paper samples 50000 peaks, but the post-release file has a different row count and is not checkpoint-bound; seed and exact training background remain absent.",
        },
        "motif_sea_threshold": {
            "status": "unresolved",
            "post_release_candidates": [1.0, 50000.0],
            "blocker": "The first official example used 1.0 and a later correction changed it to 50000.0; pinned WDL has no default or invocation.",
        },
        "motif_roster_and_order": {
            "status": "authoritative",
            "count": 693,
            "sha256": pinned["motif_file_sha256"],
            "first": pinned["motif_first"],
            "last": pinned["motif_last"],
            "order": pinned["motif_order"],
        },
        "motif_normalization": {
            "status": "training_and_inference_minmax_with_validation_conflict",
            "value": "per-vector min-max scaling for training and released analysis utilities",
            "blocker": "The released validation deserializer instead applies frozen mean/standard-deviation normalization, so validation parity must not be assumed.",
        },
        "sea_runtime": {
            "status": "authoritative_for_pinned_workflow",
            "value": "MEME Suite 5.4.1",
            "note": "The later notebook invokes 5.5.6 and is not checkpoint-bound.",
        },
        "peak_calling_and_reference": {
            "status": "unresolved_for_new_context",
            "blocker": "No checkpoint-bound invocation freezes the exact background file, peak-selection seed, native FASTA, or complete peak-calling inputs for a new MASLD context.",
        },
    }
    unresolved = [name for name, row in fields.items() if row["status"] != "authoritative" and not row["status"].startswith("authoritative_for")]
    result = {
        "schema_version": "masld-bench-epibert-native-preprocessing-audit-v1",
        "status": "pass_terminal_blocker_frozen",
        "model_id": "epibert",
        "checkpoint_revision": EXPECTED["revision"],
        "input_artifacts_sha256": {
            "admission": EXPECTED["admission"],
            "runtime": EXPECTED["runtime"],
            "primary_evidence": EXPECTED["primary_evidence"],
        },
        "archives": archives,
        "pinned_source": pinned,
        "paper": paper,
        "repository_history": history,
        "field_disposition": fields,
        "unresolved_fields": unresolved,
        "native_preprocessing_exact": False,
        "native_preprocessing_fixture_built": False,
        "biological_data_read": False,
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
        "metrics_calculated": False,
        "gpu_queue_item_created": False,
        "task_scope": "observed_ATAC_only",
        "rna_conditioned_atac_eligible": False,
        "sealed_champion_eligible": False,
        "terminal_disposition": "checkpoint_runtime_ready_native_preprocessing_terminally_blocked_no_fixture",
        "terminal_reason": "At least the ATAC scale transform, foreground selection, peak-window mapping, exact background, and SEA threshold remain non-authoritative or contradictory at the checkpoint boundary.",
    }
    return result, roster


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--primary-evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise NativePreprocessingAuditError("output already exists")
    result, roster = audit(args)
    args.output.mkdir(parents=True, mode=0o750)
    (args.output / "native_preprocessing_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (args.output / "motif_roster.tsv").open("x", encoding="utf-8") as handle:
        handle.write("index\tmotif_id\tmotif_description\n")
        for index, (motif_id, description) in enumerate(roster, start=1):
            handle.write(f"{index}\t{motif_id}\t{description}\n")
    print(json.dumps({"status": result["status"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
