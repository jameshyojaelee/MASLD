"""Multi-covariate rank-linear partial association, extending the held-back instrument.

Stage 3 adjusts for five lineage proportions **in addition to** the other axis --
six covariates. The instrument every earlier stage shares,
``scripts/evaluate_gse267145_axis_count.py``, residualizes against exactly one
covariate: it applies a single rank-1 projection.

**That file is not edited.** Six held-back jobs hashed it into their source
manifests at ``5f4e483dc22f1143...`` -- the Stage 0b prespecification and
analysis, the derived-after addendum, Stage 0c, the Stage 1 substrate build and
the Stage 1 external check. Editing it would break the manifest of all six for a
change none of them use. So the extension lives here, imports the held-back module,
and is proved to **reduce exactly to it** when the covariate count is one. That
reduction is what preserves comparability across stages; without it, Stage 3's
numbers would sit on a quietly different instrument.

Method. Every variable is average-rank transformed by the held-back module's own
``average_ranks``. The covariate block is orthonormalized by QR, which is stable
under the strong collinearity that closure forces here -- hepatocytes sit at a
median of 0.88 to 0.91, so their rank is nearly minus the sum of the other four.
Genes and exposure are projected off that basis and rescaled to unit norm, and
the partial association is the dot product, exactly as in the one-covariate case.

The conditioning of the covariate block is measured and returned rather than
assumed, because a rank-deficient block would silently produce partial
correlations against a subspace smaller than the one named.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Sequence

import numpy as np


class MultiCovariateError(RuntimeError):
    """Raised when the covariate block cannot support the adjustment."""


#: The held-back instrument. Never edited; imported.
SEALED_INSTRUMENT = "scripts/evaluate_gse267145_axis_count.py"
SEALED_INSTRUMENT_SHA256_PREFIX = "5f4e483dc22f1143"

#: Residual norms are judged against this fraction of the original norm, never
#: against exact zero. Projecting a vector off a basis it lies in leaves a
#: rounding remnant of order 1e-16, which an exact-zero test allows as usable.
DEGENERACY_TOLERANCE = 1e-10


def load_sealed(path: Path):
    spec = importlib.util.spec_from_file_location("sealed_axis_count", path)
    if spec is None or spec.loader is None:
        raise MultiCovariateError(f"cannot import the sealed instrument at {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def covariate_basis(
    analysis, covariates: Sequence[np.ndarray], *, tolerance: float = 1e-10
) -> tuple[np.ndarray, dict[str, object]]:
    """Orthonormal basis for the centred rank space the covariates span.

    Returns the basis and a measured conditioning report. Closure makes the
    lineage ranks strongly dependent, so the rank of the block is checked
    against the number of covariates rather than presumed equal to it.
    """

    if not covariates:
        raise MultiCovariateError("at least one covariate is required")
    block = np.column_stack([analysis.centred_ranks(v) for v in covariates])
    singular = np.linalg.svd(block, compute_uv=False)
    largest = float(singular[0])
    if largest <= 0.0:
        raise MultiCovariateError("every covariate is constant")
    retained = int((singular > tolerance * largest).sum())
    basis, _ = np.linalg.qr(block)
    basis = basis[:, :retained]
    report = {
        "n_covariates": len(covariates),
        "numerical_rank": retained,
        "is_rank_deficient": retained < len(covariates),
        "singular_values": [float(s) for s in singular],
        "condition_number": (
            float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
        ),
        "why_it_is_measured": (
            "closure makes the lineage ranks strongly dependent -- hepatocytes "
            "sit near 0.9 so their rank is nearly minus the sum of the others. "
            "A rank-deficient block would produce partial correlations against "
            "a smaller subspace than the one named, silently."
        ),
    }
    return basis, report


def residualize_on_basis(target: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Project a vector or a column block off an orthonormal basis."""

    if target.ndim == 1:
        return target - basis @ (basis.T @ target)
    return target - basis @ (basis.T @ target)


def partial_associations(
    analysis,
    gene_ranks_centred: np.ndarray,
    exposure: np.ndarray,
    covariates: Sequence[np.ndarray],
) -> tuple[np.ndarray, dict[str, object]]:
    """Per-gene rank-linear partial association given several covariates.

    ``gene_ranks_centred`` is the column-centred average-rank matrix the held-back
    module produces. Genes left with no residual variance after adjustment are
    returned as zero and counted, never divided by.
    """

    basis, report = covariate_basis(analysis, covariates)

    # Degeneracy is judged on a RELATIVE tolerance, never against exact zero. A
    # vector projected off a basis it already lies in lands near 1e-16, not at
    # 0.0, so `norm > 0.0` silently allows a gene that carries no residual
    # information at all. Same lesson as the deconvolution row sums, which are
    # 1.0 to within 1.6e-15 and bitwise 1.0 on only a handful of rows.
    original = np.sqrt((gene_ranks_centred**2).sum(axis=0))
    residual_genes = residualize_on_basis(gene_ranks_centred, basis)
    norms = np.sqrt((residual_genes**2).sum(axis=0))
    usable = norms > DEGENERACY_TOLERANCE * np.maximum(original, np.finfo(float).tiny)
    scaled = np.zeros_like(residual_genes)
    scaled[:, usable] = residual_genes[:, usable] / norms[usable]

    centred_exposure = analysis.centred_ranks(exposure)
    exposure_norm = float(np.sqrt((centred_exposure**2).sum()))
    residual_exposure = residualize_on_basis(centred_exposure, basis)
    scale = float(np.sqrt((residual_exposure**2).sum()))
    if exposure_norm == 0.0 or scale <= DEGENERACY_TOLERANCE * exposure_norm:
        raise MultiCovariateError(
            "the exposure lies in the span of the covariate block; the partial "
            f"association is undefined (residual norm {scale:.3e} against "
            f"{exposure_norm:.3e})"
        )
    associations = scaled.T @ (residual_exposure / scale)
    report = {
        **report,
        "genes_with_no_residual_variance": int((~usable).sum()),
        "residual_df": gene_ranks_centred.shape[0] - 2 - report["numerical_rank"],
    }
    return associations, report


def reduces_to_sealed_single_covariate(
    analysis, gene_ranks_centred: np.ndarray, exposure: np.ndarray,
    covariate: np.ndarray, *, tolerance: float = 1e-12
) -> float:
    """Largest absolute disagreement with the held-back one-covariate path.

    This is the property that keeps Stage 3 on the same instrument as the four
    held-back stages. It is returned as a number so a caller can assert on it
    rather than trusting a docstring.
    """

    unit_covariate = analysis.unit_ranks(covariate)
    scaled, norms = analysis.unit_columns(
        analysis.residualize(gene_ranks_centred, unit_covariate)
    )
    residual = analysis.residualize(
        analysis.centred_ranks(exposure), unit_covariate
    )
    sealed = scaled.T @ (residual / np.sqrt((residual**2).sum()))
    sealed[norms == 0.0] = 0.0
    extended, _ = partial_associations(
        analysis, gene_ranks_centred, exposure, [covariate]
    )
    return float(np.max(np.abs(sealed - extended)))
