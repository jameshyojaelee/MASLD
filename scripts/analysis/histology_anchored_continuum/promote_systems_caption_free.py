#!/usr/bin/env python3
"""Promote caption-free Figure 4G and S4Q-S4U molecular-systems panels."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", ""))
if not ROOT.is_dir():
    raise SystemExit("MASLD_PROJECT_ROOT is unset or invalid")

FAMILY_SOURCE = ROOT / (
    "figures/candidates/"
    "histology-continuum-molecular-systems-caption-free-20260825T130145Z"
)
S4U_SOURCE = ROOT / (
    "figures/candidates/"
    "histology-continuum-molecular-systems-s4u-caption-free-20260825T130145Z"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
SOURCE_DESTINATION = FIG4 / (
    "source_tables/current_candidate/molecular_systems_caption_free_20260825"
)
PROMOTION_MANIFEST = FIG4 / (
    "manifests/molecular_systems_caption_free_promotion_20260825.tsv"
)


@dataclass(frozen=True)
class Panel:
    callout: str
    figure_id: str
    source_root: Path
    source_name: str
    destination: Path
    old_sha256: str
    new_sha256: str


PANELS = (
    Panel(
        "4G", "continuum_effect_map", FAMILY_SOURCE,
        "s4_molecular_systems_continuum_effect_map.pdf",
        FIG4 / "panels/fig4g_molecular_systems_continuum_effect_map.pdf",
        "b44b1d8c295b540d4b3c012716782686c47322f476fe0105972ab4684c306706",
        "5dede3189a33366361e928b3b0771f519b19b5bc1f33e717cb636f83af985b42",
    ),
    Panel(
        "S4Q", "membership_atlas", FAMILY_SOURCE,
        "s4_molecular_systems_membership_atlas.pdf",
        FIG4 / "panels/supplementary/figs4q_molecular_systems_membership_atlas.pdf",
        "ed67c86b0e4ffe8d1fd8c8780b3b61e861cd20ec1a7207b2a5b81c9537f2d551",
        "c862c8dbd978c4c0d56c1b39f8bf579717cfdf7f8876aac2fa5e4f26917314f7",
    ),
    Panel(
        "S4R", "fibrosis_stage_map", FAMILY_SOURCE,
        "s4_molecular_systems_fibrosis_stage_effect_map.pdf",
        FIG4 / "panels/supplementary/figs4r_molecular_systems_fibrosis_stage_map.pdf",
        "f86f4ee6e76bebcd78481c57aeace1e65c2c67d03c5e64f5329afac2a5fec639",
        "2f4ceb5e8225f39778d6836ed935579b9abee37d11eb60d3c2c8ea827153638b",
    ),
    Panel(
        "S4S", "system_stage_continuum_summary", FAMILY_SOURCE,
        "s4_molecular_systems_stage_vs_continuum_summary.pdf",
        FIG4 / "panels/supplementary/figs4s_molecular_systems_stage_vs_continuum.pdf",
        "a45da244ced048a819ec7d69db41f17b5e5ce31ffa3f8e4444f465a873efad28",
        "524a65064d1e82ac30b22cbe241a5eeb1432e01e4f0c910fb4ccf1344802139d",
    ),
    Panel(
        "S4T", "five_cohort_sliding_windows", FAMILY_SOURCE,
        "s4_molecular_systems_five_cohort_sliding_windows.pdf",
        FIG4 / "panels/supplementary/figs4t_molecular_systems_continuum_windows.pdf",
        "e06153eeb7e3c87ba63cac5b67992af77fbc7d88339e1e8feba159deaad19382",
        "5661febc1df2618866f6ff0bbe55306af128a316566d06b6c86e42cec8860563",
    ),
    Panel(
        "S4U", "pooled_continuum_windows_with_inference", S4U_SOURCE,
        "s4_molecular_systems_pooled_continuum_windows_with_inference.pdf",
        FIG4 / "panels/supplementary/figs4u_molecular_systems_pooled_continuum_windows.pdf",
        "a615371846b451a455dc1aff59d34f27db336177da86163af457fefaa9298626",
        "07419c3b34b247716d60122d76cddfed4fc1ae32f31e0b6050d866a0fce06b8f",
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def staged_copy(source: Path, destination: Path) -> Path:
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    shutil.copy2(source, temporary)
    if sha256(temporary) != sha256(source):
        temporary.unlink(missing_ok=True)
        raise SystemExit(f"Staged panel hash mismatch: {source}")
    return temporary


def validate_candidate(root: Path) -> dict[str, dict[str, str]]:
    if not (root / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit(f"Candidate is not review-ready: {root}")
    audit = read_tsv(root / "provenance/audit.tsv")
    if not audit or any(row["passed"].upper() != "TRUE" for row in audit):
        raise SystemExit(f"Candidate audit failed: {root}")
    figures = {row["figure_id"]: row for row in read_tsv(root / "figure_manifest.tsv")}
    if not figures or any(
        row.get("embedded_explanatory_text", "TRUE").upper() != "FALSE"
        for row in figures.values()
    ):
        raise SystemExit(f"Caption-free figure contract drift: {root}")
    for row in read_tsv(root / "provenance/output_checksums.tsv"):
        path = Path(row["path"])
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise SystemExit(f"Candidate output drift: {path}")
    return figures


def main() -> None:
    family_figures = validate_candidate(FAMILY_SOURCE)
    s4u_figures = validate_candidate(S4U_SOURCE)
    if SOURCE_DESTINATION.exists() or PROMOTION_MANIFEST.exists():
        raise SystemExit("Caption-free promotion namespace already exists")

    staged: list[tuple[Path, Path]] = []
    try:
        for panel in PANELS:
            figures = family_figures if panel.source_root == FAMILY_SOURCE else s4u_figures
            row = figures.get(panel.figure_id)
            source = panel.source_root / "panels" / panel.source_name
            if row is None or Path(row["path"]).resolve() != source.resolve():
                raise SystemExit(f"Candidate manifest path drift: {panel.callout}")
            if sha256(source) != panel.new_sha256:
                raise SystemExit(f"Reviewed panel hash drift: {panel.callout}")
            if not panel.destination.is_file() or sha256(panel.destination) != panel.old_sha256:
                raise SystemExit(f"Active panel hash drift: {panel.callout}")
            staged.append((staged_copy(source, panel.destination), panel.destination))
        for temporary, destination in staged:
            os.replace(temporary, destination)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)

    SOURCE_DESTINATION.mkdir(parents=True)
    copied: list[Path] = []
    for label, source_root in (("family", FAMILY_SOURCE), ("s4u", S4U_SOURCE)):
        destination_root = SOURCE_DESTINATION / label
        destination_root.mkdir()
        for relative in (
            "figure_manifest.tsv", "READY_FOR_REVIEW.ok", "source_tables", "provenance"
        ):
            source = source_root / relative
            destination = destination_root / relative
            if source.is_dir():
                shutil.copytree(source, destination)
                copied.extend(path for path in destination.rglob("*") if path.is_file())
            else:
                shutil.copy2(source, destination)
                copied.append(destination)

    with PROMOTION_MANIFEST.open("x", newline="") as handle:
        fields = (
            "callout", "destination", "source", "old_sha256", "new_sha256",
            "embedded_explanatory_text", "data_values_changed",
            "geometry_changed", "display_order_changed", "revision",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for panel in PANELS:
            writer.writerow({
                "callout": panel.callout,
                "destination": str(panel.destination.relative_to(ROOT)),
                "source": str(
                    (panel.source_root / "panels" / panel.source_name).relative_to(ROOT)
                ),
                "old_sha256": panel.old_sha256,
                "new_sha256": panel.new_sha256,
                "embedded_explanatory_text": "FALSE",
                "data_values_changed": "FALSE",
                "geometry_changed": "FALSE",
                "display_order_changed": "FALSE",
                "revision": "external_caption_only_remove_embedded_bottom_prose",
            })

    checksum_path = SOURCE_DESTINATION / "promoted_source_checksums.tsv"
    with checksum_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "sha256"))
        for path in sorted(copied):
            writer.writerow((str(path.relative_to(ROOT)), sha256(path)))
    (SOURCE_DESTINATION / "PROMOTED.ok").write_text(
        "user_approved_caption_free_molecular_systems_family\n"
    )


if __name__ == "__main__":
    main()
