#!/usr/bin/env python3
"""Promote hash-pinned five-cohort descriptive continuum panels."""

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

SOURCE = ROOT / (
    "figures/candidates/"
    "histology-continuum-five-cohort-context-20260824T175647Z/panels"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
FIGS3 = ROOT / "figures/supplementary/figS03_bulk_transcriptomics"


@dataclass(frozen=True)
class Panel:
    callout: str
    source_name: str
    destination: Path
    old_sha256: str
    new_sha256: str


PANELS = (
    Panel("4F", "fig4f_continuum_program_trajectories.pdf",
          FIG4 / "panels/fig4f_continuum_program_trajectories.pdf",
          "eefddbd0c72ae1591a5eaf169c70a546d4e52e2c8846522e6381ba7595bae73b",
          "51f8d5ad4d82f112f36a224d93af7710e6914934bdb1e1248cff03efcf3eba30"),
    Panel("S3I", "s3_continuum_calibration_by_stage.pdf",
          FIGS3 / "panels/figs3i_continuum_calibration_by_stage.pdf",
          "a71bd46f29f5dc523189b000fef7c4762e5c985ad705eb27475e7a0436f0835a",
          "f485cbf14d2e3523e8ed7e169c6f9e164e84fee34bf7d1913b7e9c9e67b5596f"),
    Panel("S3K", "s3_all_continuum_scores_by_stage.pdf",
          FIGS3 / "panels/figs3k_all_continuum_scores_by_stage.pdf",
          "7e62e9535ea8d21b8149e500624f4238311bf6a401603ff66aa91efd281364b2",
          "db059d94465866f30cd1dd21ea5870e092a3af40a78b3b48a313e9a8e81a50a6"),
    Panel("S3P", "s3_fixed_gene_continuum_trajectories.pdf",
          FIGS3 / "panels/figs3p_fixed_gene_continuum_trajectories.pdf",
          "14b477605588ac8991bd07b78fbb72292d91eef0f320a46dd97208bf7191c909",
          "34dd64a20fbdd256b62ef999fb63378fd09021555360b85c286b4ca0521c9e53"),
    Panel("S3T", "s3_six_hallmarks_by_continuum_windows.pdf",
          FIGS3 / "panels/figs3t_six_hallmarks_by_continuum_windows.pdf",
          "c25d7577fefb196391cb8b547db1821d70e5f419c48c4b94052ee52d3ad288ff",
          "bd568999d9a801b577593c18f0536d4d0cc6c3fd5c1c4fcce6bfea93c30a8669"),
    Panel("S3U", "s3_six_hallmarks_by_fibrosis_stage.pdf",
          FIGS3 / "panels/figs3u_six_hallmarks_by_fibrosis_stage.pdf",
          "6ec47c1d6fc1b0206248818460d3865581db2606d4e58880dd9371f9ac1e7890",
          "c60fedf2370c97b0ae4583ca078da645e053d42c40e801d420b813c2f4d999f9"),
    Panel("S4F", "s4_focal_hotspot_by_continuum_windows.pdf",
          FIG4 / "panels/supplementary/figs4f_focal_hotspot_by_continuum_windows.pdf",
          "2a64e1bb3aec98369ab20a54ec701a0606cf7299bc7c12dcae7bca488bd5db3b",
          "e48f6229a0c7fc4f0a55396fc477ef02c6c6acdd03239e53ca19d973a9b32333"),
    Panel("S4G", "s4_focal_hotspot_by_fibrosis_stage.pdf",
          FIG4 / "panels/supplementary/figs4g_focal_hotspot_by_fibrosis_stage.pdf",
          "ca1b6682081600c2747b9d76d49c273a392ea87da388e58861ab8476b2c13874",
          "14631a108ed340513bff864db5c606aa9541247aee9da2047bea7db6dc381a5d"),
    Panel("S4H", "s4_all117_hotspot_continuum_windows_labeled.pdf",
          FIG4 / "panels/supplementary/figs4h_all117_hotspot_continuum_windows_labeled.pdf",
          "b6726b02a768ddeb9b9db0e0f83712c1bcb47630f1674872aaafbe4b7af8594d",
          "a1dd08aa970042e2eb0769fc58c4be6fc27bba09cc9062114e76f03fb1a59974"),
    Panel("S4I", "s4_continuum_associated_hotspot_windows_labeled.pdf",
          FIG4 / "panels/supplementary/figs4i_supported_hotspot_continuum_windows_labeled.pdf",
          "2e325899643cbe782577b5669c5cbd7964aca22aa9b4505652d672e6531ef77a",
          "406ed019b3505813482dbf4c03e034e10020ed93986d2f0b9e214f87a78b6862"),
    Panel("S4M", "s4_all10_nmf_by_continuum_windows.pdf",
          FIG4 / "panels/supplementary/figs4m_all10_nmf_by_continuum_windows.pdf",
          "03505dd112ac95f29c1728e0f2bf647ac6ab123451bec0514cf7aa31d13e4272",
          "bb3131142e37b004586c7b68ee7d0314f359a86b580f553a38f7924178103947"),
    Panel("S4N", "s4_all10_nmf_by_fibrosis_stage.pdf",
          FIG4 / "panels/supplementary/figs4n_all10_nmf_by_fibrosis_stage.pdf",
          "599301edfa21e5e7ea11c43dfc4e87c82da8d105766235fa7067dc8da68bd623",
          "3a9d704b0c4a5d686196cbf0a78e6eafecdda54bfc44ab0e895d3112541b09f7"),
    Panel("S4O", "s4_four_tf_regulon_continuum_windows.pdf",
          FIG4 / "panels/supplementary/figs4o_four_tf_regulon_continuum_windows.pdf",
          "d312d21961e6780719e163c655fa2b6f7a3240cd5af65d25067bc054025de6fd",
          "8ff06805597b7b0abfaa838e03ca052cb5925ff4592d1e85fb412edc5e793c20"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    for item in PANELS:
        source = SOURCE / item.source_name
        if not source.is_file() or sha256(source) != item.new_sha256:
            raise SystemExit(f"Source hash mismatch: {source}")
        if not item.destination.is_file() or sha256(item.destination) != item.old_sha256:
            raise SystemExit(f"Destination preflight hash mismatch: {item.destination}")

    rows: list[dict[str, str]] = []
    for item in PANELS:
        source = SOURCE / item.source_name
        with tempfile.NamedTemporaryFile(
            dir=item.destination.parent, prefix=f".{item.destination.name}.",
            suffix=".tmp", delete=False
        ) as handle:
            temporary = Path(handle.name)
        try:
            shutil.copy2(source, temporary)
            if sha256(temporary) != item.new_sha256:
                raise SystemExit(f"Temporary copy hash mismatch: {temporary}")
            os.replace(temporary, item.destination)
        finally:
            temporary.unlink(missing_ok=True)
        if sha256(item.destination) != item.new_sha256:
            raise SystemExit(f"Promoted destination hash mismatch: {item.destination}")
        rows.append({
            "callout": item.callout,
            "destination": str(item.destination.relative_to(ROOT)),
            "source": str(source.relative_to(ROOT)),
            "old_sha256": item.old_sha256,
            "new_sha256": item.new_sha256,
            "descriptive_context": "five_cohorts_or_four_stage_available_cohorts",
            "primary_inference": "GSE162694;GSE213621_only",
        })

    manifests = {
        FIG4 / "manifests/continuum_five_cohort_revision_20260824.tsv":
            [row for row in rows if row["callout"].startswith("S4") or row["callout"] == "4F"],
        FIGS3 / "manifests/continuum_five_cohort_revision_20260824.tsv":
            [row for row in rows if row["callout"].startswith("S3")],
    }
    for path, selected in manifests.items():
        if path.exists():
            raise SystemExit(f"Refusing to overwrite promotion manifest: {path}")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(selected)


if __name__ == "__main__":
    main()
