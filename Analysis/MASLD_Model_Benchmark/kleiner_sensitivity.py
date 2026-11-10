#!/usr/bin/env python3
"""Hypothetical classical measurement-error sensitivity, not empirical calibration.

The historical 0.4655--0.6433 biopsy variability is not a qualified error range.
Serial biopsies mix biological change, sampling and reading differences. This
calculator requires an explicit assumption and does not estimate its validity.
The stored severity scorer and its predictions are unchanged.
"""
import argparse
import json
import math


def disattenuate(stage_var, *, sigma_e, assumption, r_obs=None, beta_obs=None):
    """Assume additive independent mean-zero error in a stage predictor.

    The slope formula applies when stage is the error-prone predictor, not
    when stage is the response. Correlation assumes the other variable is
    measured without error. Reliability is conditional on these assumptions.
    """
    if not assumption or not assumption.strip():
        raise ValueError("An explicit measurement-error assumption is required")
    if not math.isfinite(stage_var) or stage_var <= 0:
        raise ValueError("stage variance must be finite and positive")
    if not math.isfinite(sigma_e) or sigma_e < 0:
        raise ValueError("assumed error SD must be finite and nonnegative")
    if sigma_e ** 2 >= stage_var:
        raise ValueError("Assumed error variance must be below observed variance; this does not diagnose the sample")
    rho = 1 - sigma_e ** 2 / stage_var
    result = dict(evidence_state="hypothetical_sensitivity", assumption=assumption,
                  empirically_calibrated=False, stage_var=stage_var,
                  assumed_sigma_e=sigma_e, conditional_reliability=rho,
                  conditional_correlation_bound=math.sqrt(rho))
    if r_obs is not None:
        if not math.isfinite(r_obs) or abs(r_obs) > 1:
            raise ValueError("observed correlation must be finite and within [-1, 1]")
        value = r_obs / math.sqrt(rho)
        result.update(r_observed=r_obs, conditional_latent_correlation=value,
                      incompatible_assumptions=abs(value) > 1)
    if beta_obs is not None:
        if not math.isfinite(beta_obs):
            raise ValueError("stage-predictor slope must be finite")
        result.update(stage_predictor_slope=beta_obs,
                      conditional_latent_slope=beta_obs / rho)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-var", type=float, required=True)
    parser.add_argument("--sigma-e", type=float, required=True, help="Explicit assumed error SD in stage units")
    parser.add_argument("--assumption", required=True, help="State hypothetical assumption or independent calibration source")
    parser.add_argument("--r", type=float)
    parser.add_argument("--stage-predictor-beta", type=float)
    args = parser.parse_args()
    try:
        result = disattenuate(args.stage_var, sigma_e=args.sigma_e,
                              assumption=args.assumption, r_obs=args.r,
                              beta_obs=args.stage_predictor_beta)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
