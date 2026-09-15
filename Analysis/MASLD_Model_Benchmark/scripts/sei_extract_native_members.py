#!/usr/bin/env python3
"""Extract only the five included Sei native members from its frozen archive."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import tarfile

import numpy as np


class SeiNativeMemberError(ValueError):
    """Raised when the frozen Sei archive differs from its included census."""


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def extract(
    archive: Path,
    admission_summary: Path,
    output: Path,
) -> dict[str, object]:
    if output.exists() or archive.is_symlink() or admission_summary.is_symlink():
        raise SeiNativeMemberError("Sei native-member request differs")
    summary = json.loads(admission_summary.read_text(encoding="utf-8"))
    if (
        summary.get("status") != "pass"
        or summary.get("archive_extracted")
        or summary.get("checkpoint_deserialized")
        or summary.get("pickle_loaded")
        or summary.get("next_allowed_action")
        != "isolated_no_network_tensor_only_checkpoint_inspection"
        or archive.stat().st_size != summary.get("archive_size_bytes")
        or _sha256_file(archive) != summary.get("archive_sha256")
    ):
        raise SeiNativeMemberError("Sei archive admission differs")
    expected = summary.get("expected_native_members")
    if not isinstance(expected, dict) or set(expected) != {
        "model/histone_inds.npy",
        "model/projvec_targets.npy",
        "model/sei.pth",
        "model/seqclass.names",
        "model/target.names",
    }:
        raise SeiNativeMemberError("Sei native-member census differs")
    output.mkdir(mode=0o750)
    observed: set[str] = set()
    with tarfile.open(archive, mode="r:gz") as handle:
        for member in handle:
            if member.name not in expected:
                continue
            if member.name in observed or not member.isfile():
                raise SeiNativeMemberError("Sei admitted member type differs")
            source = handle.extractfile(member)
            if source is None:
                raise SeiNativeMemberError("Sei admitted member cannot be read")
            destination = output / Path(member.name).name
            if destination.exists() or destination.is_symlink():
                raise SeiNativeMemberError("Sei output member collides")
            with source, destination.open("xb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
            contract = expected[member.name]
            if (
                destination.stat().st_size != contract["size_bytes"]
                or _sha256_file(destination) != contract["sha256"]
            ):
                raise SeiNativeMemberError("Sei extracted member bytes differ")
            observed.add(member.name)
    if observed != set(expected):
        raise SeiNativeMemberError("Sei admitted member is absent")

    projection = np.load(output / "projvec_targets.npy", allow_pickle=False)
    histone = np.load(output / "histone_inds.npy", allow_pickle=False)
    targets = (output / "target.names").read_text(encoding="utf-8").splitlines()
    classes = (output / "seqclass.names").read_text(encoding="utf-8").splitlines()
    if (
        projection.shape != (61, 21907)
        or projection.dtype.kind not in "fiu"
        or not np.isfinite(projection).all()
        or histone.shape != (10064,)
        or histone.dtype.kind not in "iu"
        or not np.isfinite(histone).all()
        or int(histone.min()) < 0
        or int(histone.max()) >= 21907
        or len(set(histone.tolist())) != 10064
        or len(targets) != 21907
        or len(set(targets)) != 21907
        or len(classes) != 40
        or len(set(classes)) != 40
    ):
        raise SeiNativeMemberError("Sei projection or target contract differs")
    receipt = {
        "schema_version": "masld-bench-sei-native-member-extraction-v1",
        "status": "pass",
        "archive_sha256": summary["archive_sha256"],
        "member_sha256": {
            path.name: _sha256_file(path) for path in sorted(output.iterdir())
        },
        "checkpoint_size_bytes": (output / "sei.pth").stat().st_size,
        "projection_shape": list(projection.shape),
        "projection_dtype": str(projection.dtype),
        "projection_basis_count": int(projection.shape[0]),
        "histone_index_shape": list(histone.shape),
        "histone_index_dtype": str(histone.dtype),
        "target_count": len(targets),
        "sequence_class_count": len(classes),
        "projection_to_sequence_class_mapping_unresolved": True,
        "tar_member_extraction_method": "exact_allowlist_stream_copy_no_extract",
        "checkpoint_deserialized": False,
        "pickle_loaded": False,
        "model_instantiated": False,
        "model_forward_executed": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
    }
    (output / "extraction_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--admission-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    extract(args.archive, args.admission_summary, args.output)


if __name__ == "__main__":
    main()
