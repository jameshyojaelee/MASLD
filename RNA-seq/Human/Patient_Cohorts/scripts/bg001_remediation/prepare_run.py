#!/usr/bin/env python3
"""Create an immutable, candidate-only BG-001 run skeleton and source snapshot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
from datetime import datetime, timezone
from pathlib import Path

RUN_RE = re.compile(r"^bg001-fragment-v211-gencode49-[0-9]{8}T[0-9]{6}Z$")
SOURCE_FILES = (
    "config/human_datasets.yaml",
    "data/gencode_v49_gene_metadata.tsv.gz",
    "data/published_gene_panels/govaere_2020_panel.tsv",
    "data/published_gene_panels/niddk_pipeline_2024.tsv",
    "data/published_gene_panels/opentargets_masld_2025.tsv",
    "docs/archive/documentation_consolidation_2026-08-11/originals/docs/manuscript/NUMBERS.md",
    "docs/archive/documentation_consolidation_2026-08-11/originals/docs/paper_outline.md",
    "figures/README.md",
    "figures/main/fig3_RNAseq/README.md",
    "figures/main/fig4_validation/README.md",
    "scripts/figures/run_pub_figures.sh",
    "Analysis/Multimodal_Program_Projection/README.md",
    "Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv",
    "Analysis/Multimodal_Program_Projection/scripts/run_figures.sbatch",
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv",
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/01_sample_qc.R",
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/03_integrate_counts.R",
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05h_limma_voom_qw_canonical.R",
    "RNA-seq/Human/Patient_Cohorts/pipelines/shared/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/shared/scripts/merge_featurecounts_exact.py",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE130970/scripts/submit_pipeline.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE135251/scripts/submit_pipeline.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/workflow/config.yaml",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/workflow/config.yaml",
    "RNA-seq/results/convergence/triple_convergence_targets.csv",
    "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv",
    "RNA-seq/results/multi_evidence/convergence_evidence_inhibitor_target_candidates.csv",
    "RNA-seq/results/validation/positive_control_validation.csv",
)
SOURCE_GLOBS = (
    "docs/archive/documentation_consolidation_2026-08-11/originals/docs/manuscript/working/**/*.md",
    "scripts/figures/**/*.R",
    "RNA-seq/46d*.R",
    "RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/*.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE130970/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE135251/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/scripts/*.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE213621/scripts/*.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts/*.sh",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/workflow/Snakefile",
    "RNA-seq/Human/Patient_Cohorts/pipelines/custom/PRJNA512027/scripts/*.sh",
    "RNA-seq/Mouse/InHouse_MCD/workflow/Snakefile",
    "RNA-seq/Mouse/Public_MCD/GSE205974/workflow/Snakefile",
)
DEPENDENCY_PATTERNS = {
    "raw_counts": ("gene_counts.txt",),
    "qc": ("sample_qc_report.csv",),
    "merged_raw": ("merged_counts_raw.rds",),
    "dge": ("merged_dge.rds",),
    "canonical_deg": ("canonical_deg_results.csv",),
    "consensus": ("consensus_degs.csv",),
    "multi_evidence_atlas": ("multi_evidence_atlas.csv",),
    "convergence_evidence": ("convergence_evidence.csv",),
}
DEPENDENCY_SUFFIXES = {".R", ".r", ".py", ".sh", ".bash", ".sbatch", ".smk", ".yaml", ".yml"}
DEPENDENCY_EXCLUDED_PREFIXES = (
    "Cas13_Library_Design/", "archive/", "results/",
    "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/",
)
WRAPPER_TOKEN = re.compile(r"(?P<path>[A-Za-z0-9_./-]+[.](?:R|r|py|sh|bash|sbatch))")
ASSIGNMENT_RE = re.compile(r"^\s*(?:export\s+)?(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*(?:<-|=)\s*(?P<value>.+?)\s*$")
READ_HINT_RE = re.compile(r"\b(?:fread|readRDS|read[.]csv|read_csv|read[.]table|load|scan)\s*\(")
WRITE_HINT_RE = re.compile(r"\b(?:fwrite|saveRDS|write[.]csv|write_csv|write[.]table|writeLines)\s*\(")
VARIABLE_REF_RE = re.compile(r"\$(?:\{(?P<braced>[A-Za-z_][A-Za-z0-9_]*)\}|(?P<plain>[A-Za-z_][A-Za-z0-9_]*))")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def require_secure_directory(path: Path, label: str, *, owner_private: bool) -> None:
    """Reject a publication ancestor another account could replace or enter."""
    try:
        status = os.lstat(path)
    except FileNotFoundError as exc:
        raise SystemExit(f"Missing {label}: {path}") from exc
    if stat.S_ISLNK(status.st_mode) or not stat.S_ISDIR(status.st_mode):
        raise SystemExit(f"{label} must be a non-symlink directory: {path}")
    permissions = stat.S_IMODE(status.st_mode)
    if owner_private:
        if status.st_uid != os.geteuid() or permissions & 0o700 != 0o700 or permissions & 0o077:
            raise SystemExit(
                f"{label} must be owned by effective UID {os.geteuid()} with owner-only rwx: "
                f"{path} mode={permissions:#05o}"
            )
    elif permissions & 0o022:
        raise SystemExit(
            f"{label} must not be writable by group/other accounts: "
            f"{path} mode={permissions:#05o}"
        )


def run_text(command: list[str], cwd: Path) -> str:
    return subprocess.run(command, cwd=cwd, check=True, text=True, capture_output=True).stdout


def noncomment_lines(path: Path) -> list[tuple[int, str]]:
    lines = []
    for number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "//")):
            continue
        lines.append((number, line))
    return lines


def expand_known_variables(value: str, variables: dict[str, str]) -> str:
    for _ in range(5):
        updated = VARIABLE_REF_RE.sub(
            lambda match: variables.get(match.group("braced") or match.group("plain"), match.group(0)),
            value,
        )
        if updated == value:
            break
        value = updated
    return value


def infer_variables(project: Path, relative: str, lines: list[tuple[int, str]]) -> dict[str, str]:
    parent = (project / relative).parent
    variables = {
        "SCRIPT_DIR": str(parent),
        "PROJECT_ROOT": str(project),
    }
    for _, line in lines:
        match = ASSIGNMENT_RE.match(line)
        if not match:
            continue
        name, raw = match.group("name"), match.group("value").strip()
        if "dirname" in raw and ("$0" in raw or "BASH_SOURCE" in raw):
            variables[name] = str(parent)
            continue
        if "$(" in raw or "`" in raw:
            continue
        raw = raw.split(" #", 1)[0].strip().strip("\"'")
        expanded = expand_known_variables(raw, variables)
        if "$" not in expanded and expanded:
            variables[name] = expanded
    return variables


def resolve_wrapper_targets(
    project: Path,
    relative: str,
    lines: list[tuple[int, str]],
    sources_by_basename: dict[str, set[str]],
) -> set[str]:
    targets: set[str] = set()
    parent = (project / relative).parent
    variables = infer_variables(project, relative, lines)
    for _, line in lines:
        expanded_line = expand_known_variables(line, variables)
        for match in WRAPPER_TOKEN.finditer(expanded_line):
            token = match.group("path")
            candidates = (parent / token, project / token)
            resolved_target = False
            for candidate in candidates:
                try:
                    resolved = candidate.resolve(strict=True)
                    target = str(resolved.relative_to(project))
                except (FileNotFoundError, ValueError):
                    continue
                targets.add(target)
                resolved_target = True
                break
            if not resolved_target:
                matches = sources_by_basename.get(Path(token).name, set())
                if len(matches) == 1:
                    targets.update(matches)
    return targets


def role_access_modes(
    lines: list[tuple[int, str]],
    role: str,
    needles: tuple[str, ...],
) -> tuple[set[str], list[str]]:
    modes: set[str] = set()
    evidence: list[str] = []
    variables: set[str] = set()
    for line_number, line in lines:
        if not any(needle in line for needle in needles):
            continue
        evidence.append(f"{line_number}:{line.strip()[:500]}")
        assignment = ASSIGNMENT_RE.match(line)
        if assignment:
            variables.add(assignment.group("name"))
        if WRITE_HINT_RE.search(line):
            modes.add("write")
        if READ_HINT_RE.search(line):
            modes.add("read")
    for _, line in lines:
        for variable in variables:
            if not re.search(rf"\b{re.escape(variable)}\b", line):
                continue
            if WRITE_HINT_RE.search(line):
                modes.add("write")
            if READ_HINT_RE.search(line):
                modes.add("read")
    if evidence and not modes:
        # An artifact path assigned to a variable and subsequently passed to a
        # project-specific loader is a read even when the loader is not one of
        # the standard library functions recognized above.
        modes.add("read" if variables else "reference")
    return modes, evidence


def main() -> None:
    # Every candidate byte is private by construction, independent of the
    # submitting shell's umask.  The process exits after preparing one run, so
    # there is no caller state to restore.
    os.umask(0o077)
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--run-id")
    args = parser.parse_args()
    project = args.project_root.resolve(strict=True)
    remediation_src = project / "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
    source_validator = remediation_src / "validate_featurecounts_sources.py"
    live_source_gate = subprocess.run(
        ["python3", str(source_validator), "--project-root", str(project)],
        check=True,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    figure_validator = remediation_src / "validate_rendered_figure_labels.py"
    figure_inventory = remediation_src / "active_figure_panel_inventory.tsv"
    live_figure_gate = subprocess.run(
        ["python3", str(figure_validator), "--project-root", str(project)],
        check=True,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    with figure_inventory.open(newline="") as handle:
        figure_reader = csv.DictReader(handle, delimiter="\t")
        if tuple(figure_reader.fieldnames or ()) != (
            "panel_id", "pdf_path", "label_mode", "classification_reason",
            "pdf_sha256", "pdftotext_version", "extracted_text_sha256",
        ):
            raise SystemExit(f"Malformed active figure inventory: {figure_inventory}")
        figure_rows = list(figure_reader)
    if len(figure_rows) != 77 or len({row["pdf_path"] for row in figure_rows}) != 77:
        raise SystemExit("Active figure inventory must bind exactly 77 distinct PDFs")
    run_id = args.run_id or datetime.now(timezone.utc).strftime("bg001-fragment-v211-gencode49-%Y%m%dT%H%M%SZ")
    if not RUN_RE.fullmatch(run_id):
        raise SystemExit(f"Invalid BG-001 run ID: {run_id}")
    candidate_parent_path = project / "results/remediation/bg001"
    candidate_parent_path.mkdir(parents=True, exist_ok=True)
    candidate_parent = candidate_parent_path.resolve(strict=True)
    if candidate_parent != candidate_parent_path:
        raise SystemExit(f"Candidate parent traverses a symlink: {candidate_parent_path}")
    require_secure_directory(project / "results", "results directory", owner_private=False)
    require_secure_directory(
        project / "results/remediation", "remediation directory", owner_private=True
    )
    require_secure_directory(candidate_parent, "BG-001 candidate parent", owner_private=True)
    run_root = candidate_parent / run_id
    if run_root.exists() or run_root.is_symlink():
        raise SystemExit(f"Refusing existing run root: {run_root}")
    run_root.mkdir(mode=0o700)
    os.chmod(run_root, 0o700)
    require_secure_directory(run_root, "run root", owner_private=True)
    (run_root / ".bg001_candidate_root").write_text(run_id + "\n")
    for relative in (
        "contract",
        "manifests",
        "manifests/bam_hash_shards",
        "frozen_sets",
        "source_snapshot",
        "counts",
        "comparisons",
        "logs",
        "promotion",
    ):
        directory = run_root / relative
        directory.mkdir(mode=0o700)
        os.chmod(directory, 0o700)
    for arm in ("R0", "F_locked", "F_legacy", "F_five"):
        arm_root = run_root / "arms" / arm
        for relative in ("qc", "results/integration", "logs", "provenance"):
            (arm_root / relative).mkdir(parents=True, exist_ok=True)
        (arm_root / ".bg001_candidate_root").write_text(f"{run_id}\t{arm}\n")

    remediation_dst = run_root / "source_snapshot/RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation"
    shutil.copytree(remediation_src, remediation_dst, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    source_rows = []
    explicit_sources = set(SOURCE_FILES)
    explicit_sources.update(row["pdf_path"] for row in figure_rows)
    for pattern in SOURCE_GLOBS:
        explicit_sources.update(
            str(path.relative_to(project))
            for path in project.glob(pattern)
            if path.is_file()
        )
    for relative in sorted(explicit_sources):
        src = project / relative
        if not src.is_file():
            raise SystemExit(f"Missing source file: {src}")
        dst = run_root / "source_snapshot" / relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    # Freeze every tracked code/config consumer that explicitly reads a
    # BG-001-sensitive artifact. The post-recount dependency report is built
    # from this inventory, not from a concurrently changing worktree.
    tracked_raw = subprocess.run(
        ["git", "ls-files", "-z"], cwd=project, check=True, capture_output=True
    ).stdout.split(b"\0")
    untracked_raw = subprocess.run(
        ["git", "ls-files", "--others", "--exclude-standard", "-z"],
        cwd=project,
        check=True,
        capture_output=True,
    ).stdout.split(b"\0")
    tracked = {raw.decode(errors="surrogateescape") for raw in tracked_raw if raw}
    untracked = {raw.decode(errors="surrogateescape") for raw in untracked_raw if raw}
    sources_by_basename: dict[str, set[str]] = {}
    for relative in tracked | untracked:
        path = project / relative
        if path.suffix in DEPENDENCY_SUFFIXES or path.name == "Snakefile":
            sources_by_basename.setdefault(path.name, set()).add(relative)
    dependency_records: dict[str, dict[str, object]] = {}
    for relative in sorted(tracked | untracked):
        if relative.startswith(DEPENDENCY_EXCLUDED_PREFIXES):
            continue
        src = project / relative
        if src.suffix not in DEPENDENCY_SUFFIXES or not src.is_file() or src.stat().st_size > 20_000_000:
            continue
        lines = noncomment_lines(src)
        direct_roles: set[str] = set()
        consumer_roles: set[str] = set()
        direct_access: dict[str, set[str]] = {}
        evidence: list[str] = []
        for role, needles in DEPENDENCY_PATTERNS.items():
            modes, role_evidence = role_access_modes(lines, role, needles)
            if not modes:
                continue
            direct_roles.add(role)
            direct_access[role] = modes
            evidence.extend(role_evidence)
            if modes != {"write"}:
                consumer_roles.add(role)
        dependency_records[relative] = {
            "direct_roles": direct_roles,
            "roles": consumer_roles,
            "direct_access": direct_access,
            "evidence": evidence,
            "wrapper_targets": resolve_wrapper_targets(project, relative, lines, sources_by_basename),
            "source_kind": "tracked" if relative in tracked else "untracked",
        }

    # Propagate exact artifact dependencies through literal wrapper invocations.
    changed = True
    while changed:
        changed = False
        for record in dependency_records.values():
            inherited: set[str] = set()
            for target in record["wrapper_targets"]:
                child = dependency_records.get(target)
                if child:
                    inherited.update(child["roles"])
            before = set(record["roles"])
            record["roles"].update(inherited)
            changed = changed or record["roles"] != before

    dependency_rows = []
    for relative, record in dependency_records.items():
        roles = sorted(record["roles"])
        if not roles:
            continue
        src = project / relative
        snapshot_relative = Path("dependency_sources") / relative
        dst = run_root / "source_snapshot" / snapshot_relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        dependency_rows.append(
            {
                "consumer_path": relative,
                "snapshot_relative_path": str(snapshot_relative),
                "consumed_roles": "|".join(roles),
                "direct_roles": "|".join(sorted(record["direct_roles"])),
                "direct_access_modes": "|".join(
                    f"{role}:{','.join(sorted(modes))}"
                    for role, modes in sorted(record["direct_access"].items())
                ),
                "effective_access_modes": "|".join(
                    f"{role}:{','.join(sorted(record['direct_access'].get(role, {'transitive_read'})))}"
                    for role in roles
                ),
                "wrapper_targets": "|".join(sorted(record["wrapper_targets"])),
                "evidence_lines": " || ".join(record["evidence"]),
                "source_kind": record["source_kind"],
                "source_sha256": sha256(dst),
            }
        )
    dependency_rows.sort(key=lambda row: row["consumer_path"])
    with (run_root / "contract/dependency_source_inventory.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "consumer_path", "snapshot_relative_path", "consumed_roles", "direct_roles",
                "direct_access_modes", "effective_access_modes", "wrapper_targets", "evidence_lines",
                "source_kind", "source_sha256",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(dependency_rows)
    with (run_root / "contract/untracked_dependency_sources.tsv").open("w", newline="") as handle:
        fields = ("consumer_path", "source_sha256")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {"consumer_path": row["consumer_path"], "source_sha256": row["source_sha256"]}
            for row in dependency_rows
            if row["source_kind"] == "untracked"
        )
    for path in sorted((run_root / "source_snapshot").rglob("*")):
        if path.is_file():
            source_rows.append(
                {
                    "relative_path": str(path.relative_to(run_root / "source_snapshot")),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )

    commit = run_text(["git", "rev-parse", "HEAD"], project).strip()
    status = run_text(["git", "status", "--porcelain=v2", "--untracked-files=normal"], project)
    diff = run_text(["git", "diff", "--binary", "HEAD", "--", "."], project)
    (run_root / "contract/git_status_porcelain_v2.txt").write_text(status)
    (run_root / "contract/dirty_worktree.patch").write_text(diff)
    with (run_root / "contract/source_manifest.tsv").open("w") as handle:
        handle.write("relative_path\tsize_bytes\tsha256\n")
        for row in source_rows:
            handle.write(f"{row['relative_path']}\t{row['size_bytes']}\t{row['sha256']}\n")

    source_gate = subprocess.run(
        [
            "python3",
            str(remediation_dst / "validate_featurecounts_sources.py"),
            "--project-root",
            str(run_root / "source_snapshot"),
        ],
        check=True,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    (run_root / "contract/source_regression_check.txt").write_text(
        "LIVE WORKTREE VALIDATION\n"
        + live_source_gate.stdout
        + live_source_gate.stderr
        + "FROZEN SNAPSHOT VALIDATION\n"
        + source_gate.stdout
        + source_gate.stderr
    )

    frozen_figure_gate = subprocess.run(
        [
            "python3",
            str(remediation_dst / "validate_rendered_figure_labels.py"),
            "--project-root",
            str(run_root / "source_snapshot"),
        ],
        check=True,
        text=True,
        capture_output=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    figure_gate_path = run_root / "contract/rendered_figure_label_check.txt"
    figure_gate_path.write_text(
        "LIVE WORKTREE VALIDATION\n"
        + live_figure_gate.stdout
        + live_figure_gate.stderr
        + "FROZEN SNAPSHOT VALIDATION\n"
        + frozen_figure_gate.stdout
        + frozen_figure_gate.stderr
    )

    # Discovery is allowed only for this draft. A separate SLURM prehash phase
    # fills every included BAM digest and atomically publishes bam_manifest.tsv
    # before featureCounts may start.
    builder = remediation_dst / "build_bam_manifest.py"
    draft_manifest = run_root / "manifests/bam_manifest.draft.tsv"
    subprocess.run(
        ["python3", str(builder), "--project-root", str(project), "--output", str(draft_manifest)],
        check=True,
    )
    os.chmod(draft_manifest, 0o600)

    contract = {
        "schema_version": "1.2",
        "run_id": run_id,
        "project_root": str(project),
        "candidate_root": str(run_root),
        "reference_root": str(run_root / "source_snapshot"),
        "git_commit": commit,
        "dirty_patch_sha256": sha256(run_root / "contract/dirty_worktree.patch"),
        "source_manifest_sha256": sha256(run_root / "contract/source_manifest.tsv"),
        "source_regression_check_sha256": sha256(run_root / "contract/source_regression_check.txt"),
        "rendered_figure_label_check_sha256": sha256(figure_gate_path),
        "dependency_source_inventory_sha256": sha256(run_root / "contract/dependency_source_inventory.tsv"),
        "untracked_dependency_sources_sha256": sha256(run_root / "contract/untracked_dependency_sources.tsv"),
        "bam_manifest_draft_sha256": sha256(draft_manifest),
        "bam_manifest_finalization_required": True,
        "active_affected_cohorts": ["GSE130970", "GSE135251", "GSE174478", "GSE213621", "GSE240729"],
        "canonical_affected_cohorts": ["GSE130970", "GSE135251", "GSE213621"],
        "featurecounts": {
            "version": "2.1.1",
            "arguments": ["-p", "--countReadPairs", "-B", "-s", "2"],
            "gtf": "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz",
            "gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
        },
        "arms": {arm: str(run_root / "arms" / arm) for arm in ("R0", "F_locked", "F_legacy", "F_five")},
        "arm_modes": {
            "R0": {"counts": "read", "qc": "recompute", "filter_scope": "legacy_all", "gene_mode": "native"},
            "F_locked": {"counts": "fragment", "qc": "locked", "filter_scope": "legacy_all", "gene_mode": "locked"},
            "F_legacy": {"counts": "fragment", "qc": "recompute", "filter_scope": "legacy_all", "gene_mode": "native"},
            "F_five": {"counts": "fragment", "qc": "reuse_F_legacy", "filter_scope": "canonical_five", "gene_mode": "native"},
        },
        "analysis_model": {
            "formula": "~ dataset + inferred_sex + group_binary",
            "engine": "voomWithQualityWeights -> lmFit -> treat",
            "disease_coefficient": "group_binaryDisease",
            "treat_lfc": 0.25,
            "qc_seed": 42,
        },
        "thresholds": {
            "qc_symmetric_difference_max": 4,
            "gene_jaccard_min": 0.99,
            "sample_weight_spearman_min": 0.99,
            "offset_median_abs_delta_max": 0.01,
            "offset_p95_abs_delta_max": 0.03,
            "offset_max_abs_delta_max": 0.10,
            "offset_group_shift_max": 0.02,
            "logfc_spearman_min": 0.995,
            "expressed_logfc_median_abs_delta_max": 0.01,
            "expressed_logfc_p95_abs_delta_max": 0.05,
            "treat_jaccard_min": 0.95,
            "treat_count_relative_change_max": 0.05,
        },
    }
    (run_root / "contract/run_contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(run_root)


if __name__ == "__main__":
    main()
