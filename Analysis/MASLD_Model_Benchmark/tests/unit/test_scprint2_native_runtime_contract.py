from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_scprint2_runtime_is_noneditable_and_relocation_checked() -> None:
    text = (ROOT / "slurm/build_scprint2_native_runtime_v2.sbatch").read_text()
    assert "--no-editable" in text
    assert "Runtime contains a staging-path .pth reference" in text
    assert "relocated scPRINT-2 runtime import differs" in text
    assert "--partition=io" in text
    assert "--account=nslab" in text
    assert "--qos=nslab" in text


def test_scprint2_restore_is_restricted_and_outcome_blind() -> None:
    text = (ROOT / "scripts/preflight_scprint2_native_runtime.py").read_text()
    assert 'weights_only=True' in text
    assert 'mmap=True' in text
    assert '"outcomes_or_labels_read": False' in text
    assert 'model.load_state_dict(state, strict=False)' in text
    assert 'constructor.pop("_instantiator", None)' in text
    assert 'missing_keys' in text
    assert 'unexpected_keys' in text
    assert "GSE289173" not in text


def test_scprint2_exact_identity_is_bound_in_all_runtime_layers() -> None:
    checkpoint = "2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
    commit = "7ee2de5aaa3492f326f540f703d44ddeaf1f2fd4"
    preflight = (ROOT / "scripts/preflight_scprint2_native_runtime.py").read_text()
    builder = (ROOT / "slurm/build_scprint2_native_runtime_v2.sbatch").read_text()
    manifest = (ROOT / "config/artifacts/models/scprint2/checkpoints.json").read_text()
    assert checkpoint in preflight
    assert checkpoint in manifest
    assert commit in preflight
    assert commit in builder
    assert commit in manifest
