"""Command-line control plane for the MASLD model benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Sequence

from .artifacts import ArtifactError, verify_frozen_tree
from .hashing import sha256_file
from .campaign import (
    CampaignError,
    execute_admission_run,
    execute_run,
    submission_commands,
    verify_run_execution_attempt,
)
from .planner import (
    PlanningError,
    freeze_campaign,
    load_frozen_plan,
    validate_registry_tree,
)


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_ROOT = PACKAGE_ROOT / "config"


def _json(value: Any) -> None:
    print(json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False))


def _registry_validate(arguments: argparse.Namespace) -> int:
    _json(validate_registry_tree(arguments.config_root))
    return 0


def _reference_validate(arguments: argparse.Namespace) -> int:
    from .reference import ReferenceBundle

    bundle = ReferenceBundle.load_toml(Path(arguments.config_root) / "resources.toml")
    result = bundle.validate()
    _json(
        {
            "validated_paths": {
                key: Path(value).as_posix() for key, value in result.items()
            },
            "sequence_extraction_ready": bundle.sequence_extraction_ready,
            "sequence_extraction_blockers": list(
                bundle.sequence_extraction_blockers
            ),
        }
        if isinstance(result, dict)
        else bundle.to_dict()
    )
    return 0


def _checkpoint_stage_geneformer(arguments: argparse.Namespace) -> int:
    from .checkpoint_preflight import stage_geneformer_bundle

    target = stage_geneformer_bundle(
        model_id=arguments.model_id,
        source_root=arguments.source_root,
        output_root=arguments.output_root,
        manifest_path=arguments.manifest,
    )
    manifest = verify_frozen_tree(target)
    _json(
        {
            "bundle": target.as_posix(),
            "model_id": arguments.model_id,
            "artifact_count": len(manifest["artifacts"]),
            "manifest_sha256": sha256_file(target / "ARTIFACTS.json"),
            "downloads_performed": False,
            "valid": True,
        }
    )
    return 0


def _campaign_freeze(arguments: argparse.Namespace) -> int:
    target = freeze_campaign(
        arguments.campaign,
        config_root=arguments.config_root,
        package_root=PACKAGE_ROOT,
        output_root=arguments.output_root,
    )
    plan = load_frozen_plan(target)
    _json(
        {
            "candidate": target.as_posix(),
            "plan_sha256": plan["plan_sha256"],
            "campaign_sha256": sha256_file(target / "ARTIFACTS.json"),
            "resource_totals": plan["resource_totals"],
            "submitted": False,
        }
    )
    return 0


def _campaign_verify(arguments: argparse.Namespace) -> int:
    manifest = verify_frozen_tree(arguments.candidate)
    plan = load_frozen_plan(arguments.candidate)
    _json(
        {
            "candidate": Path(arguments.candidate).resolve().as_posix(),
            "plan_sha256": plan["plan_sha256"],
            "campaign_sha256": sha256_file(
                Path(arguments.candidate).resolve() / "ARTIFACTS.json"
            ),
            "artifact_count": len(manifest["artifacts"]),
            "valid": True,
        }
    )
    return 0


def _campaign_submit(arguments: argparse.Namespace) -> int:
    commands = submission_commands(
        candidate=arguments.candidate,
        approved_campaign_sha256=arguments.approved_campaign_sha256,
        execute=arguments.execute,
        submission_ledger_root=arguments.submission_ledger_root,
    )
    _json({"execute": arguments.execute, "job_count": len(commands), "commands": commands})
    return 0


def _run_admission(arguments: argparse.Namespace) -> int:
    output = execute_admission_run(
        candidate=arguments.candidate,
        run_id=arguments.run_id,
        output_root=arguments.output_root,
    )
    _json({"run_id": arguments.run_id, "attempt": output.as_posix(), "status": "complete"})
    return 0


def _run_execute(arguments: argparse.Namespace) -> int:
    output = execute_run(
        candidate=arguments.candidate,
        run_id=arguments.run_id,
        output_root=arguments.output_root,
    )
    receipt = verify_run_execution_attempt(output)
    _json(
        {
            "run_id": arguments.run_id,
            "attempt": output.as_posix(),
            "status": receipt["status"],
        }
    )
    return 0


def _run_verify(arguments: argparse.Namespace) -> int:
    receipt = verify_run_execution_attempt(
        arguments.attempt,
        require_succeeded=not arguments.allow_failed,
    )
    _json(receipt)
    return 0


def _selection_freeze(arguments: argparse.Namespace) -> int:
    from .selection import freeze_selection_lock, verify_selection_lock

    path = freeze_selection_lock(
        candidate=arguments.candidate,
        decisions_path=arguments.decisions,
        output_root=arguments.output_root,
    )
    lock = verify_selection_lock(path)
    _json({"selection": path.as_posix(), "lock_id": lock.lock_id})
    return 0


def _selection_verify(arguments: argparse.Namespace) -> int:
    from .selection import verify_selection_lock

    lock = verify_selection_lock(arguments.selection)
    _json({"lock_id": lock.lock_id, "valid": True})
    return 0


def _release_stage(arguments: argparse.Namespace) -> int:
    from .release import stage_release, verify_release

    path = stage_release(spec_path=arguments.spec, output_root=arguments.output_root)
    manifest = verify_release(path)
    _json({"release": path.as_posix(), "release_id": manifest["release_id"], "valid": True})
    return 0


def _release_verify(arguments: argparse.Namespace) -> int:
    from .release import verify_release

    manifest = verify_release(arguments.release)
    _json({"release_id": manifest["release_id"], "valid": True})
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    groups = parser.add_subparsers(dest="group", required=True)

    registry = groups.add_parser("registry", help="validate frozen registry contracts")
    registry_actions = registry.add_subparsers(dest="action", required=True)
    registry_validate = registry_actions.add_parser("validate")
    registry_validate.add_argument("--config-root", type=Path, default=DEFAULT_CONFIG_ROOT)
    registry_validate.set_defaults(handler=_registry_validate)

    reference = groups.add_parser("reference", help="validate sequence reference admission")
    reference_actions = reference.add_subparsers(dest="action", required=True)
    reference_validate = reference_actions.add_parser("validate")
    reference_validate.add_argument("--config-root", type=Path, default=DEFAULT_CONFIG_ROOT)
    reference_validate.set_defaults(handler=_reference_validate)

    checkpoint = groups.add_parser(
        "checkpoint", help="stage or verify model checkpoint bundles"
    )
    checkpoint_actions = checkpoint.add_subparsers(dest="action", required=True)
    checkpoint_geneformer = checkpoint_actions.add_parser("stage-geneformer")
    checkpoint_geneformer.add_argument(
        "--model-id",
        required=True,
        choices=(
            "geneformer_v1_10m",
            "geneformer_v2_104m",
            "geneformer_v2_316m",
        ),
    )
    checkpoint_geneformer.add_argument("--source-root", type=Path, required=True)
    checkpoint_geneformer.add_argument("--output-root", type=Path, required=True)
    checkpoint_geneformer.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_CONFIG_ROOT
        / "artifacts"
        / "models"
        / "geneformer"
        / "checkpoints.json",
    )
    checkpoint_geneformer.set_defaults(handler=_checkpoint_stage_geneformer)

    campaign = groups.add_parser("campaign", help="freeze, verify, or submit campaign plans")
    campaign_actions = campaign.add_subparsers(dest="action", required=True)
    campaign_freeze = campaign_actions.add_parser("freeze")
    campaign_freeze.add_argument("--campaign", type=Path, required=True)
    campaign_freeze.add_argument("--config-root", type=Path, default=DEFAULT_CONFIG_ROOT)
    campaign_freeze.add_argument("--output-root", type=Path, default=PACKAGE_ROOT / "candidates")
    campaign_freeze.set_defaults(handler=_campaign_freeze)
    campaign_verify = campaign_actions.add_parser("verify")
    campaign_verify.add_argument("--candidate", type=Path, required=True)
    campaign_verify.set_defaults(handler=_campaign_verify)
    campaign_submit = campaign_actions.add_parser("submit")
    campaign_submit.add_argument("--candidate", type=Path, required=True)
    campaign_submit.add_argument(
        "--approved-campaign-sha256",
        required=True,
        help="SHA-256 of the frozen candidate ARTIFACTS.json, including every job script",
    )
    campaign_submit.add_argument(
        "--execute",
        action="store_true",
        help="actually call sbatch; omitted means render exact commands only",
    )
    campaign_submit.add_argument(
        "--submission-ledger-root",
        type=Path,
        default=PACKAGE_ROOT / "executions" / "submission-ledgers",
        help="stable root for the single-use immutable initial-submission ledger",
    )
    campaign_submit.set_defaults(handler=_campaign_submit)

    run = groups.add_parser("run", help="execute one reviewed immutable run")
    run_actions = run.add_subparsers(dest="action", required=True)
    run_admission = run_actions.add_parser("admission")
    run_admission.add_argument("--candidate", type=Path, required=True)
    run_admission.add_argument("--run-id", required=True)
    run_admission.add_argument("--output-root", type=Path, required=True)
    run_admission.set_defaults(handler=_run_admission)
    run_execute = run_actions.add_parser("execute")
    run_execute.add_argument("--candidate", type=Path, required=True)
    run_execute.add_argument("--run-id", required=True)
    run_execute.add_argument("--output-root", type=Path, required=True)
    run_execute.set_defaults(handler=_run_execute)
    run_verify = run_actions.add_parser("verify")
    run_verify.add_argument("--attempt", type=Path, required=True)
    run_verify.add_argument(
        "--allow-failed",
        action="store_true",
        help="verify a terminal failed receipt as well as successful attempts",
    )
    run_verify.set_defaults(handler=_run_verify)

    selection = groups.add_parser("selection", help="freeze all choices before sealed inference")
    selection_actions = selection.add_subparsers(dest="action", required=True)
    selection_freeze = selection_actions.add_parser("freeze")
    selection_freeze.add_argument("--candidate", type=Path, required=True)
    selection_freeze.add_argument("--decisions", type=Path, required=True)
    selection_freeze.add_argument("--output-root", type=Path, default=PACKAGE_ROOT / "selections")
    selection_freeze.set_defaults(handler=_selection_freeze)
    selection_verify = selection_actions.add_parser("verify")
    selection_verify.add_argument("--selection", type=Path, required=True)
    selection_verify.set_defaults(handler=_selection_verify)

    release = groups.add_parser("release", help="stage or verify an open research release")
    release_actions = release.add_subparsers(dest="action", required=True)
    release_stage = release_actions.add_parser("stage")
    release_stage.add_argument("--spec", type=Path, required=True)
    release_stage.add_argument("--output-root", type=Path, default=PACKAGE_ROOT / "releases")
    release_stage.set_defaults(handler=_release_stage)
    release_verify = release_actions.add_parser("verify")
    release_verify.add_argument("--release", type=Path, required=True)
    release_verify.set_defaults(handler=_release_verify)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    arguments = parser.parse_args(argv)
    expected_errors: list[type[Exception]] = [
        ArtifactError,
        CampaignError,
        PlanningError,
        ValueError,
    ]
    if arguments.group == "checkpoint":
        from .checkpoint_preflight import CheckpointPreflightError

        expected_errors.append(CheckpointPreflightError)
    elif arguments.group == "selection":
        from .selection import SelectionError

        expected_errors.append(SelectionError)
    elif arguments.group == "release":
        from .release import ReleaseError

        expected_errors.append(ReleaseError)
    try:
        return int(arguments.handler(arguments))
    except tuple(expected_errors) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
