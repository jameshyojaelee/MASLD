"""v1.8.3a unit tests for LEAD_RULINGS_W4 H-13 M1: backend_torch raises the two spec-defined failures as
ModelSpecError (a RuntimeError, so fitting.py counts them as the CPU code does) and turns every other torch
RuntimeError into DevFatal (a BaseException that stops the run). No device work: device_model is replaced."""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")
os.environ.setdefault("MODEL_A_TORCH_DEVICE", "cpu")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend_torch as bt  # noqa: E402
import fitting as ft  # noqa: E402
import model_a as ma  # noqa: E402

HERE = Path(__file__).resolve().parents[1]


class Fake:
    def __init__(self, err):
        self.err = err

    def loglik3(self, *a, **k):
        raise self.err

    loglik4 = loglik3
    boundary_score = loglik3


@pytest.fixture
def torch_backend(monkeypatch):
    ma.set_backend("torch")
    yield monkeypatch
    ma.set_backend("cpu")


def test_classes():
    assert issubclass(bt.ModelSpecError, RuntimeError)
    assert issubclass(bt.DevFatal, BaseException) and not issubclass(bt.DevFatal, Exception)


def test_spec_failures_are_model_spec_errors():
    src = (HERE / "backend_torch.py").read_text()
    assert src.count('raise ModelSpecError(f"no converged local maximum for gene {g}")') == 1
    assert src.count('raise ModelSpecError(f"reference integrator: gene {g} did not agree within "') == 1
    assert "raise RuntimeError(" not in src


@pytest.mark.parametrize("name,args", [("loglik3", (1,)), ("loglik4", (1,)), ("boundary_score", ())])
def test_spec_error_passes_through_as_runtime_error(torch_backend, name, args):
    torch_backend.setattr(bt, "device_model", lambda D: Fake(bt.ModelSpecError("no converged local maximum")))
    with pytest.raises(bt.ModelSpecError):
        getattr(ma, name)(np.zeros(3), object(), *args)
    try:
        getattr(ma, name)(np.zeros(3), object(), *args)
    except RuntimeError:          # what fitting.py catches
        pass


@pytest.mark.parametrize("err", [RuntimeError("CUBLAS_STATUS_ALLOC_FAILED"), NotImplementedError("x"),
                                 RuntimeError("deterministic refusal")])
def test_other_runtime_errors_become_devfatal(torch_backend, err):
    torch_backend.setattr(bt, "device_model", lambda D: Fake(err))
    with pytest.raises(bt.DevFatal) as ei:
        ma.loglik3(np.zeros(3), object(), 1)
    assert "non-spec RuntimeError" in str(ei.value)


def test_device_errors_become_devfatal(torch_backend):
    if not bt._DEV_ERRORS:
        pytest.skip("this torch has no OutOfMemoryError/AcceleratorError")
    try:
        e = bt._DEV_ERRORS[0]("out of memory")
    except TypeError:
        pytest.skip("cannot construct the device error class here")
    torch_backend.setattr(bt, "device_model", lambda D: Fake(e))
    with pytest.raises(bt.DevFatal):
        ma.loglik3(np.zeros(3), object(), 1)


def test_non_runtime_errors_unchanged(torch_backend):
    torch_backend.setattr(bt, "device_model", lambda D: Fake(ValueError("v")))
    with pytest.raises(ValueError):
        ma.loglik3(np.zeros(3), object(), 1)


def test_devfatal_escapes_the_fit_code(torch_backend):
    """A DevFatal raised inside L-BFGS-B's objective is not turned into +inf by fitting.maximize."""
    torch_backend.setattr(bt, "device_model", lambda D: Fake(RuntimeError("boom")))
    ev = ft.Eval3(object())
    with pytest.raises(bt.DevFatal):
        ft.maximize(ev, np.zeros(3), np.full(3, -1.0), np.full(3, 1.0), np.ones(3, bool), np.ones(3))


def test_spec_error_is_a_failed_evaluation_in_the_fit_code(torch_backend):
    torch_backend.setattr(bt, "device_model", lambda D: Fake(bt.ModelSpecError("no converged local maximum")))
    ev = ft.Eval3(object())
    with pytest.raises(bt.ModelSpecError):   # maximize's final ev(x, 2) is outside the except: _run catches it
        ft.maximize(ev, np.zeros(3), np.full(3, -1.0), np.full(3, 1.0), np.ones(3, bool), np.ones(3))
