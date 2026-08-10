#!/usr/bin/env python3
"""Token-aware regression gate for every reachable paired featureCounts producer."""

from __future__ import annotations

import argparse
import os
import re
import shlex
from pathlib import Path


HUMAN_REQUIRED = (
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE130970.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE135251.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE174478.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk0.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk1.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk2.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk3.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE240729.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/scripts/run_featurecounts.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/scripts/run_featurecounts_parallel.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/run_featurecounts.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/run_featurecounts_array.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/run_featurecounts_contract.sh",
)
DYNAMIC_REQUIRED = (
    "RNA-seq/Human/Patient_Cohorts/pipelines/shared/Snakefile",
)
GUARDED_CUSTOM_SNAKEFILES = (
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE130970/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE135251/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/workflow/Snakefile",
)
SHARED_AFFECTED_CONFIGS = (
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/workflow/config.yaml",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/workflow/config.yaml",
)
ACTIVE_COUNT_ROUTES = (
    (
        "GSE130970",
        "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE130970/scripts/submit_pipeline.sh",
        "custom",
    ),
    (
        "GSE135251",
        "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE135251/scripts/submit_pipeline.sh",
        "custom",
    ),
    (
        "GSE174478",
        "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/scripts/submit_pipeline.sh",
        "shared",
    ),
    (
        "GSE213621",
        "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/scripts/submit_pipeline.sh",
        "custom",
    ),
    (
        "GSE240729",
        "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/run_pipeline.sh",
        "shared",
    ),
)
MOUSE_REQUIRED = (
    "RNA-seq/Mouse/InHouse_MCD/workflow/Snakefile",
    "RNA-seq/Mouse/Public_MCD/GSE205974/workflow/Snakefile",
)
RETIRED_REQUIRED = (
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_PRJNA512027.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_PRJNA512027_chunk0.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_PRJNA512027_chunk1.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_PRJNA512027_chunk2.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_PRJNA512027_chunk3.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/merge_PRJNA512027_chunks.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/scripts/submit_pipeline.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/scripts/submit_megabulk.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/scripts/launch_pipeline_wrapper.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/scripts/fix_and_resubmit.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/scripts/fix_and_resubmit_r2.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/workflow/Snakefile",
)
DISABLED_LIVE_PRODUCERS = (
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE130970.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE135251.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE174478.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk0.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk1.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk2.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE213621_chunk3.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/run_GSE240729.sh",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/merge_GSE213621_chunks.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/scripts/run_featurecounts.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/scripts/run_featurecounts_parallel.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/scripts/merge_featurecounts_and_launch.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/run_featurecounts.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/run_featurecounts_array.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/merge_featurecounts.sh",
)
CLASSIFIED = set(HUMAN_REQUIRED + DYNAMIC_REQUIRED + GUARDED_CUSTOM_SNAKEFILES + MOUSE_REQUIRED)
SOURCE_SUFFIXES = {".sh", ".bash", ".sbatch", ".smk", ".py", ".R", ".r"}
PRUNE_DIRS = {"archive", "results", "logs", ".snakemake", "__pycache__", "qc", "counts", "tests"}


def command_regions(text: str) -> list[str]:
    """Extract executable-looking command lines after joining shell continuations."""
    text = re.sub(r"(?m)^\s*#.*$", "", text)
    text = text.replace("\\\n", " ")
    regions: list[str] = []
    for line in text.splitlines():
        for match in re.finditer(r"\bfeatureCounts\b", line):
            prefix = line[: match.start()].strip()
            # Exclude diagnostics/documentation/output-header writers. Quoted
            # Snakemake shell strings are retained.
            if re.search(r"(?:^|\s)(?:echo|printf|print|write|writeLines|cat)\b", prefix):
                continue
            if any(token in prefix for token in ("fh.write", "handle.write", "f.write")):
                continue
            regions.append(line[match.start() :])
    return regions


def tokens(region: str) -> list[str]:
    cleaned = region.replace("{", "X").replace("}", "X")
    # Strip common trailing quote/parenthesis punctuation without changing
    # embedded flag tokens such as --countReadPairs=false.
    try:
        return shlex.split(cleaned)
    except ValueError:
        return cleaned.split()


def has_token(region: str, token: str) -> bool:
    return token in tokens(region)


def has_strand_two(region: str) -> bool:
    parsed = tokens(region)
    values = [parsed[index + 1] for index, value in enumerate(parsed[:-1]) if value == "-s"]
    return values == ["2"]


def shell_refusal_is_first_executable(text: str, required_message: str | None = None) -> bool:
    """Require an ERROR + exit-64 guard before any executable shell statement.

    SLURM directives are comments to the shell and must remain above the guard
    so ``sbatch`` can parse them.  Checking executable lines also avoids the
    previous brittle 500-byte preamble limit.
    """
    executable = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        executable.append(stripped)
        if len(executable) == 2:
            break
    if len(executable) != 2:
        return False
    message, refusal = executable
    return (
        message.startswith('echo "ERROR:')
        and message.endswith('>&2')
        and (required_message is None or required_message in message)
        and refusal == "exit 64"
    )


def iter_source_files(root: Path):
    scan_roots = (
        root / "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2",
        root / "RNA-seq/Human/Patient_Cohorts/pipelines",
        root / "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation",
        root / "RNA-seq/Mouse/InHouse_MCD/workflow",
        root / "RNA-seq/Mouse/Public_MCD/GSE205974/workflow",
    )
    for scan_root in scan_roots:
        if not scan_root.is_dir():
            continue
        for current, dirnames, filenames in os.walk(scan_root):
            dirnames[:] = [name for name in dirnames if name not in PRUNE_DIRS]
            for filename in filenames:
                path = Path(current) / filename
                if filename == "Snakefile" or path.suffix in SOURCE_SUFFIXES:
                    yield path


def validate_paired_region(relative: str, region: str, human_contract: bool) -> list[str]:
    failures = []
    if not has_token(region, "--countReadPairs"):
        failures.append(f"paired invocation missing --countReadPairs: {relative}")
    if human_contract and not has_token(region, "-B"):
        failures.append(f"paired invocation missing -B: {relative}")
    if human_contract and not has_strand_two(region):
        failures.append(f"paired human invocation lacks the single exact '-s 2' contract: {relative}")
    return failures


def rule_block(text: str, rule_name: str) -> str:
    match = re.search(
        rf"(?ms)^[ \t]*rule[ \t]+{re.escape(rule_name)}[ \t]*:[ \t]*.*?"
        rf"(?=^[ \t]*rule[ \t]+\w+[ \t]*:|\Z)",
        text,
    )
    return match.group(0) if match else ""


def validate_actual_count_routes(root: Path) -> list[str]:
    """Verify the launcher -> workflow -> config composition for all five cohorts."""
    failures: list[str] = []
    for custom_relative in GUARDED_CUSTOM_SNAKEFILES:
        custom_path = root / custom_relative
        if not custom_path.is_file():
            failures.append(f"missing guarded custom Snakefile: {custom_relative}")
            continue
        count_block = rule_block(custom_path.read_text(errors="replace"), "featurecounts")
        if (
            not count_block
            or "exit 64" not in count_block
            or "bg001_remediation/recount_array.sbatch" not in count_block
            or any(has_token(region, "-p") for region in command_regions(count_block))
        ):
            failures.append(f"custom live count rule is not fail-closed: {custom_relative}")

    shared_relative = DYNAMIC_REQUIRED[0]
    shared_path = root / shared_relative
    shared_text = shared_path.read_text(errors="replace") if shared_path.is_file() else ""
    shared_count = rule_block(shared_text, "featurecounts_per_sample")
    shared_merge = rule_block(shared_text, "merge_featurecounts")
    affected_block = shared_text.split("BG001_AFFECTED_DATASETS", 1)[-1].split("DATASET_ID", 1)[0]
    for dataset, launcher_relative, route_kind in ACTIVE_COUNT_ROUTES:
        launcher = root / launcher_relative
        if not launcher.is_file():
            failures.append(f"missing active count launcher: {launcher_relative}")
            continue
        launcher_text = launcher.read_text(errors="replace")
        if route_kind == "custom":
            expected_snakefile = f"RNA-seq/Human/Patient_Cohorts/pipelines/custom/{dataset}/workflow/Snakefile"
            if not re.search(r"(?:-s|--snakefile)\s+(?:[\"']?)workflow/Snakefile\b", launcher_text):
                failures.append(f"active launcher does not select guarded custom Snakefile: {launcher_relative}")
            custom_path = root / expected_snakefile
            if not custom_path.is_file():
                failures.append(f"missing guarded custom Snakefile: {expected_snakefile}")
                continue
        else:
            if "pipelines/shared/Snakefile" not in launcher_text and "../../shared/Snakefile" not in launcher_text:
                failures.append(f"active launcher does not select shared Snakefile: {launcher_relative}")
            config_token = "--configfile workflow/config.yaml"
            if config_token not in launcher_text:
                failures.append(f"active shared launcher does not select workflow/config.yaml: {launcher_relative}")
            if f'"{dataset}"' not in affected_block:
                failures.append(f"shared workflow BG-001 guard omits {dataset}")

    if not shared_path.is_file():
        failures.append(f"missing shared workflow: {shared_relative}")
    elif (
        "if BG001_RECOUNT_REQUIRED:" not in shared_text
        or "else:" not in shared_text
        or not re.search(
            r"BG001_COUNTS_PRESENT\s*=.*\.is_file\(\).*\.is_file\(\)",
            shared_text,
        )
        or not re.search(
            r"BG001_PER_SAMPLE_GUARD_OUTPUTS\s*=\s*\(\)\s*if\s+BG001_COUNTS_PRESENT\s+else",
            shared_text,
        )
        or not re.search(
            r"BG001_MERGE_GUARD_OUTPUTS\s*=\s*\(\)\s*if\s+BG001_COUNTS_PRESENT\s+else",
            shared_text,
        )
        or "BG001_PER_SAMPLE_GUARD_OUTPUTS" not in shared_count
        or "BG001_MERGE_GUARD_OUTPUTS" not in shared_merge
        or not shared_count
        or "BG001_RECOUNT_REFUSAL" not in shared_count
        or "exit 64" not in shared_count
        or bool(command_regions(shared_count))
        or not shared_merge
        or "BG001_RECOUNT_REFUSAL" not in shared_merge
        or "exit 64" not in shared_merge
        or bool(command_regions(shared_merge))
    ):
        failures.append("shared affected-dataset count/merge guards are incomplete")

    for config_relative in SHARED_AFFECTED_CONFIGS:
        config_path = root / config_relative
        if not config_path.is_file():
            failures.append(f"missing affected shared-workflow config: {config_relative}")
            continue
        config_text = config_path.read_text(errors="replace")
        if not re.search(r"(?m)^\s*strandedness:\s*2(?:\s|#|$)", config_text):
            failures.append(f"affected shared-workflow config is not exact -s 2: {config_relative}")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve(strict=True)
    failures: list[str] = []

    for relative in HUMAN_REQUIRED + MOUSE_REQUIRED:
        path = root / relative
        if not path.is_file():
            failures.append(f"missing classified producer: {relative}")
            continue
        regions = command_regions(path.read_text(errors="replace"))
        paired_regions = [region for region in regions if has_token(region, "-p")]
        if not paired_regions:
            failures.append(f"no executable paired featureCounts invocation found: {relative}")
            continue
        require_s2 = relative in HUMAN_REQUIRED
        for region in paired_regions:
            failures.extend(validate_paired_region(relative, region, require_s2))

    # The shared workflow expresses pairedness through a parameter; validate
    # both its definition and its use rather than treating the placeholder as a
    # scientific flag.
    for relative in DYNAMIC_REQUIRED:
        path = root / relative
        if not path.is_file():
            failures.append(f"missing classified dynamic producer: {relative}")
            continue
        text = path.read_text(errors="replace")
        if not re.search(r"pe_flag\s*=.*[\"']-p --countReadPairs -B[\"']", text):
            failures.append(f"dynamic paired flag contract is invalid: {relative}")
        if "{params.pe_flag}" not in text:
            failures.append(f"dynamic paired flag is not used: {relative}")

    failures.extend(validate_actual_count_routes(root))

    for relative in RETIRED_REQUIRED:
        path = root / relative
        if not path.is_file():
            failures.append(f"missing retired launcher guard: {relative}")
            continue
        text = path.read_text(errors="replace")
        guarded_shell = shell_refusal_is_first_executable(text, "retired/noncanonical")
        guarded_snakefile = relative.endswith("Snakefile") and "raise ValueError" in text[:500]
        if not (guarded_shell or guarded_snakefile):
            failures.append(f"retired PRJNA512027 launcher is not fail-closed at its preamble: {relative}")

    for relative in DISABLED_LIVE_PRODUCERS:
        path = root / relative
        if not path.is_file() or not shell_refusal_is_first_executable(path.read_text(errors="replace")):
            failures.append(f"legacy live count/merge producer is not fail-closed at its preamble: {relative}")

    # Discovery is independent of the classified list: a newly reachable '-p'
    # command cannot silently escape the regression gate. Single-end commands
    # contain no exact -p token and are explicitly exempt.
    discovered_paired_files: set[str] = set()
    for path in iter_source_files(root):
        relative = str(path.relative_to(root))
        if "PRJNA512027" in relative:
            continue
        for region in command_regions(path.read_text(errors="replace")):
            if not has_token(region, "-p"):
                continue
            discovered_paired_files.add(relative)
            require_s2 = relative.startswith("RNA-seq/Human/") and "pipelines/shared/Snakefile" not in relative
            failures.extend(validate_paired_region(relative, region, require_s2))
    unclassified = sorted(discovered_paired_files - CLASSIFIED)
    if unclassified:
        failures.extend(f"unclassified reachable paired producer: {path}" for path in unclassified)

    if failures:
        raise SystemExit("\n".join(sorted(set(failures))))
    print(
        f"PASS: {len(HUMAN_REQUIRED)} fixed human, {len(GUARDED_CUSTOM_SNAKEFILES)} guarded custom, "
        f"{len(DYNAMIC_REQUIRED)} guarded dynamic shared, "
        f"and {len(MOUSE_REQUIRED)} mouse paired producers satisfy the fragment-counting contract"
    )


if __name__ == "__main__":
    main()
