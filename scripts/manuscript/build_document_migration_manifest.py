#!/usr/bin/env python3
"""Build the auditable documentation-consolidation migration manifest."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "docs/archive/documentation_consolidation_2026-08-11"
INPUT = ARCHIVE / "ALL_MARKDOWN_BEFORE_DIRECTORY_MOVE.tsv"
OUTPUT = ROOT / "docs/archive/DOCUMENT_MIGRATION_MANIFEST.tsv"


HUMAN_OWNERS = {
    "docs/CODEBASE_CURRENT_STATE.md": "docs/STATUS.md",
    "docs/RESOURCE_PAPER_SCOPE.md": "docs/PAPER.md",
    "docs/paper_outline.md": "docs/PAPER.md",
    "docs/paper_narrative.md": "docs/PAPER.md",
    "docs/progress.md": "docs/STATUS.md",
    "docs/manuscript/NUMBERS.md": "docs/RESULTS.md",
    "docs/manuscript/working/README.md": "docs/manuscript/README.md",
    "docs/manuscript/working/F_FIVE_RESOURCE_ADOPTION.md": "docs/STATUS.md",
    "docs/manuscript/working/RECONCILIATION_STATUS.md": "docs/STATUS.md",
    "docs/manuscript/working/key_points_by_section.md": "docs/PAPER.md",
    "docs/manuscript/working/METHODS.md": "docs/manuscript/METHODS.md",
    "docs/competitor_paper_details.md": "docs/literature/REVIEW.md",
    "docs/hmsma_access_request_2026-08-10.md": "docs/ROADMAP.md",
    "docs/single_cell_analysis.md": "docs/technical/SINGLE_CELL_HISTORY.md",
}

PRESERVED_IN_PLACE = {
    "docs/FIGURE_GUIDELINES.md",
    "docs/manuscript/05_figure_legends.md",
    "docs/manuscript/release/README.md",
}

NEW_CURRENT_EXACT = {
    "docs/README.md",
    "docs/PAPER.md",
    "docs/STATUS.md",
    "docs/RESULTS.md",
    "docs/ROADMAP.md",
    "docs/manuscript/README.md",
    "docs/manuscript/METHODS.md",
    "docs/archive/INDEX.md",
}
NEW_CURRENT_PREFIXES = (
    "docs/manuscript/draft/",
    "docs/technical/",
    "docs/literature/",
)

for name in (
    "title.md",
    "abstract.md",
    "01_intro.md",
    "fig1.md",
    "fig2.md",
    "fig3.md",
    "fig4.md",
    "fig5_discussion.md",
):
    HUMAN_OWNERS[f"docs/manuscript/working/{name}"] = f"docs/manuscript/draft/{name}"


def classify(path: str) -> tuple[str, str, str]:
    """Return disposition, current owner/location, and preserved copy."""

    protected = (
        "docs/archive/documentation_consolidation_2026-08-11/originals/" + path
    )
    if path in HUMAN_OWNERS:
        return "retired_to_current_owner", HUMAN_OWNERS[path], protected

    if path in PRESERVED_IN_PLACE:
        return "updated_in_place_original_preserved", path, protected

    if path.startswith("docs/plans/"):
        suffix = path.removeprefix("docs/plans/")
        return "archived_plan", f"docs/archive/plans/{suffix}", f"docs/archive/plans/{suffix}"

    if path.startswith("docs/progress/"):
        suffix = path.removeprefix("docs/progress/")
        return (
            "archived_progress",
            f"docs/archive/progress/{suffix}",
            f"docs/archive/progress/{suffix}",
        )

    if path.startswith("docs/audits/"):
        suffix = path.removeprefix("docs/audits/")
        return "archived_audit", f"docs/archive/audits/{suffix}", f"docs/archive/audits/{suffix}"

    if path.startswith("docs/manuscript/reviews/"):
        suffix = path.removeprefix("docs/manuscript/reviews/")
        return (
            "archived_review",
            f"docs/archive/reviews/{suffix}",
            f"docs/archive/reviews/{suffix}",
        )

    if path.startswith("docs/archive/"):
        return "existing_archive", path, path

    if path.startswith(("docs/reference/", "docs/competitor_analysis/")):
        if (ROOT / protected).is_file():
            return "updated_in_place_original_preserved", path, protected
        return "source_reference_kept_in_place", path, path

    if path.startswith("docs/manuscript/"):
        return "manuscript_or_release_kept_in_place", path, path

    return "technical_or_reference_kept_in_place", path, path


def main() -> None:
    with INPUT.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    fieldnames = [
        "original_path",
        "bytes",
        "sha256",
        "disposition",
        "current_owner_or_location",
        "preserved_copy",
    ]
    source_rows = [
        row
        for row in rows
        if row["path"] not in NEW_CURRENT_EXACT
        and not row["path"].startswith(NEW_CURRENT_PREFIXES)
    ]

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in source_rows:
            disposition, current, preserved = classify(row["path"])
            writer.writerow(
                {
                    "original_path": row["path"],
                    "bytes": row["bytes"],
                    "sha256": row["sha256"],
                    "disposition": disposition,
                    "current_owner_or_location": current,
                    "preserved_copy": preserved,
                }
            )

    print(f"DOCUMENT_MIGRATION_MANIFEST\t{OUTPUT.relative_to(ROOT)}")
    print(f"ROWS\t{len(source_rows)}")


if __name__ == "__main__":
    main()
