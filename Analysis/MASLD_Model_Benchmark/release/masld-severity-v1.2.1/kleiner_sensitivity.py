#!/usr/bin/env python3
"""Hypothetical classical measurement-error sensitivity calculator (masld-severity-v1.2.1).

Replaces retired kleiner.py from v1.2. The historical 0.4655--0.6433 biopsy variability
is not a qualified error range because serial biopsies mix biological change, sampling,
and reading differences. This calculator requires an explicit assumption and enforces:
  1. |r_latent| <= 1.0 feasibility guard (refuses incompatible assumptions unless
     --allow-infeasible-assumption is passed).
  2. Predictor-vs-response guard (refuses --stage-role response or --stage-response-beta,
     since 1/rho slope disattenuation applies only when stage is an error-prone predictor).
The stored severity scorer and its predictions are unchanged from v1.2.
"""
import argparse
import json
import math
import sys


def disattenuate(
    stage_var,
    *,
    sigma_e,
    assumption,
    r_obs=None,
    beta_obs=None,
    stage_role="predictor",
    stage_response_beta=None,
    allow_infeasible=False,
):
    """Assume additive independent mean-zero error in a stage predictor."""
    if stage_role != "predictor" or stage_response_beta is not None:
        raise ValueError(
            "Classical slope disattenuation (beta_obs / rho) applies ONLY when Kleiner "
            "stage is the error-prone predictor, NOT when stage is the regression response."
        )
    if not assumption or not assumption.strip():
        raise ValueError("An explicit measurement-error assumption is required")
    if not math.isfinite(stage_var) or stage_var <= 0:
        raise ValueError("stage variance must be finite and positive")
    if not math.isfinite(sigma_e) or sigma_e < 0:
        raise ValueError("assumed error SD must be finite and nonnegative")
    if sigma_e ** 2 >= stage_var:
        raise ValueError(
            "Assumed error variance must be below observed stage variance; "
            "this does not diagnose the sample"
        )
    rho = 1.0 - (sigma_e ** 2) / stage_var
    result = dict(
        evidence_state="hypothetical_sensitivity",
        assumption=assumption,
        empirically_calibrated=False,
        stage_role="predictor",
        stage_var=stage_var,
        assumed_sigma_e=sigma_e,
        conditional_reliability=rho,
        conditional_correlation_bound=math.sqrt(rho),
    )
    if r_obs is not None:
        if not math.isfinite(r_obs) or abs(r_obs) > 1.0:
            raise ValueError("observed correlation must be finite and within [-1, 1]")
        r_latent = r_obs / math.sqrt(rho)
        infeasible = abs(r_latent) > 1.0
        if infeasible and not allow_infeasible:
            raise ValueError(
                f"Infeasible measurement-error assumptions: |r_latent|={abs(r_latent):.4f} > 1.0 "
                f"(r_obs={r_obs}, sqrt(rho)={math.sqrt(rho):.4f}). Re-examine sigma_e or pass "
                f"--allow-infeasible-assumption to inspect the contradictory value."
            )
        result.update(
            r_observed=r_obs,
            conditional_latent_correlation=r_latent,
            feasible=not infeasible,
            incompatible_assumptions=infeasible,
        )
    if beta_obs is not None:
        if not math.isfinite(beta_obs):
            raise ValueError("stage-predictor slope must be finite")
        result.update(
            stage_predictor_slope=beta_obs,
            conditional_latent_slope=beta_obs / rho,
        )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-var", type=float, required=True)
    parser.add_argument("--sigma-e", type=float, required=True, help="Explicit assumed error SD in stage units")
    parser.add_argument(
        "--assumption",
        required=True,
        help="State hypothetical assumption or independent calibration source",
    )
    parser.add_argument("--r", type=float, help="Observed correlation with stage")
    parser.add_argument("--stage-predictor-beta", type=float, help="Observed slope when stage is the predictor")
    parser.add_argument(
        "--stage-response-beta",
        type=float,
        default=None,
        help="Rejected: disattenuating a slope where stage is the response is invalid under classical predictor error.",
    )
    parser.add_argument(
        "--stage-role",
        choices=("predictor", "response"),
        default="predictor",
        help="Role of stage in regression (must be predictor for 1/rho slope adjustment).",
    )
    parser.add_argument(
        "--allow-infeasible-assumption",
        action="store_true",
        help="Emit the contradictory |r_latent| > 1.0 value instead of raising on infeasibility.",
    )
    args = parser.parse_args(argv)
    try:
        result = disattenuate(
            args.stage_var,
            sigma_e=args.sigma_e,
            assumption=args.assumption,
            r_obs=args.r,
            beta_obs=args.stage_predictor_beta,
            stage_role=args.stage_role,
            stage_response_beta=args.stage_response_beta,
            allow_infeasible=args.allow_infeasible_assumption,
        )
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
