#!/usr/bin/env python3
"""Kleiner disattenuation calculator.

The Kleiner fibrosis stage is measured with error. How much error is NOT a
settled constant: it differs between the two deposits where it has been measured.

    sigma_e RANGE:  0.4655  to  0.6433  stages

| deposit / arm            |  n | sigma_e | CI95             | label ICC(1,1) | % first read F0-F1 |
|--------------------------|---:|--------:|------------------|---------------:|-------------------:|
| GSE193066                | 58 |  0.6433 | [0.517,  0.760]  |         0.1886 |               20.7 |
| PXD051911 A1, at surgery | 27 |  0.5092 | [0.304,  0.694]  |         0.6228 |               77.8 |
| PXD051911 A2, follow-up  | 30 |  0.4655 | [0.365,  0.548]  |         0.6117 |               90.0 |
| PXD051911 A3, post-surg. | 11 |  0.4264 | -                |         0.5918 |               81.8 |

A2 against GSE193066 is DISCREPANT: ratio 0.7236, CI95 [0.5366, 0.9555], which
excludes 1. A1 is consistent but underpowered (ratio CI [0.437, 1.156]).

These are NOT pooled. There are three arms, DerSimonian-Laird tau^2 is severely
biased at k=3, and the deposits are not exchangeable: their first-read stage
composition differs sharply (20.7% vs 77.8% vs 90.0% at F0-F1). For the record,
pooling was computed and not shipped: fixed-effect 0.5596, random-effects 0.5435,
I^2 58.1%, Q p = 0.092.

sigma_e itself is stage-dependent, which is the likely source of the spread:
within GSE193066 it is 0.8165 at first-read F0-F1 (n=12) and 0.5898 at F2+
(n=46). **If you have repeat reads, measure sigma_e in your own sample.**

    rho         = 1 - sigma_e^2 / Var(Kleiner stages in YOUR sample)
    r_latent    = r_observed / sqrt(rho)
    ceiling     : r_observed <= sqrt(rho)
    beta_latent = beta_observed / rho

## Which end of the range is conservative

A LARGER sigma_e gives a SMALLER rho, hence a LOWER ceiling and a LARGER
disattenuated r_latent. A lower ceiling is easier to approach and a larger
r_latent looks stronger, so **both of the usual claims are FLATTERED by the large
end**. The conservative end is therefore sigma_e = 0.4655, and it is the default
here. You must pass sigma_e explicitly to use any other value.

Usage
    python kleiner.py --demo
    python kleiner.py --stage-sd 1.2 --r 0.62
    python kleiner.py --stages stages.txt --r 0.62 --beta 0.41
    python kleiner.py --stage-var 1.50 --sigma-e 0.6433     # explicit, anti-conservative

What this does NOT license
    sigma_e is TEST-RETEST across two physically distinct biopsies. It folds
    reader variability, biopsy sampling variability and any real change over the
    interval into one term and cannot separate them. It is NOT an inter-reader
    statistic; no deposit in this project records who scored a slide.
"""
import argparse
import json
import math
import sys

SIGMA_E_LOW = 0.4654746681256314     # PXD051911 A2, n=30 -- CONSERVATIVE
SIGMA_E_HIGH = 0.6432675209026769    # GSE193066, n=58    -- anti-conservative
SIGMA_E_RANGE = (SIGMA_E_LOW, SIGMA_E_HIGH)
SIGMA_E_DEFAULT = SIGMA_E_LOW

# Measured Kleiner stage variance of the 58 GSE193066 PARTICIPANTS over their 116 reads.
# Used by --demo so the worked example sits on a real deposit rather than a round number.
# Source: executions/release-instrument-20260901T211017Z/results/job_c.json
#         diagnostics.stage_variance.all_116_reads
GSE193066_58PAIR_STAGE_VAR = 0.5091454272863568

DEPOSITS = [
    ("GSE193066", 58, 0.6432675209026769, (0.517, 0.760), 0.1886, 20.7),
    ("PXD051911 A1 at-surgery", 27, 0.5091750772173156, (0.3043, 0.6939), 0.6228, 77.8),
    ("PXD051911 A2 follow-up", 30, 0.4654746681256314, (0.3651, 0.5477), 0.6117, 90.0),
    ("PXD051911 A3 post-surgical", 11, 0.4264, None, 0.5918, 81.8),
]
POOLING_NOT_SHIPPED = dict(fixed_effect=0.5596, random_effects=0.5435, I2_pct=58.1, Q_p=0.09215,
                           k=3, why_not="k=3; DerSimonian-Laird tau^2 is severely biased at that k "
                                        "and the deposits are not exchangeable in stage composition")


class KleinerError(ValueError):
    pass


def rho_from_var(stage_var, sigma_e=SIGMA_E_DEFAULT):
    """Reliability of a SINGLE Kleiner read in a sample with this stage variance.

    sigma_e defaults to the CONSERVATIVE end of the measured range."""
    if stage_var <= 0:
        raise KleinerError(f"stage variance must be positive; got {stage_var}")
    ev = sigma_e ** 2
    if stage_var <= ev:
        raise KleinerError(
            f"stage variance {stage_var:.4f} does not exceed the error variance {ev:.4f} at "
            f"sigma_e={sigma_e:.4f}. rho would be non-positive: in a sample this homogeneous the "
            f"label carries no usable signal above its own read error and no disattenuation is "
            f"defined. That is a real answer about your sample, not a numerical failure.")
    return 1.0 - ev / stage_var


def disattenuate(stage_var, r_obs=None, beta_obs=None, sigma_e=SIGMA_E_DEFAULT):
    """Single-value entry point. sigma_e defaults to the CONSERVATIVE end; pass it
    explicitly to use any other value."""
    rho = rho_from_var(stage_var, sigma_e)
    ceiling = math.sqrt(rho)
    out = dict(sigma_e=sigma_e, is_conservative_end=bool(sigma_e == SIGMA_E_LOW),
               stage_var=stage_var, stage_sd=math.sqrt(stage_var), error_var=sigma_e ** 2,
               rho=rho, ceiling_r_observed=ceiling)
    if r_obs is not None:
        if not -1.0 <= r_obs <= 1.0:
            raise KleinerError(f"observed correlation {r_obs} is outside [-1, 1]")
        out["r_observed"] = r_obs
        out["r_latent"] = r_obs / ceiling
        out["r_latent_is_impossible"] = abs(r_obs / ceiling) > 1.0
        out["exceeds_ceiling"] = abs(r_obs) > ceiling
    if beta_obs is not None:
        out["beta_observed"] = beta_obs
        out["beta_latent"] = beta_obs / rho
        out["attenuation_factor"] = rho
    return out


def across_range(stage_var, r_obs=None, beta_obs=None):
    """Every quantity at BOTH ends of the measured sigma_e range, labelled.

    This is the entry point to prefer. A single number hides the fact that the
    measured range spans conclusions."""
    ends = {}
    for tag, s in (("conservative", SIGMA_E_LOW), ("anti_conservative", SIGMA_E_HIGH)):
        try:
            ends[tag] = disattenuate(stage_var, r_obs, beta_obs, sigma_e=s)
        except KleinerError as e:
            ends[tag] = dict(undefined=str(e), sigma_e=s)
    out = dict(sigma_e_range=list(SIGMA_E_RANGE),
               conservative_end=SIGMA_E_LOW, anti_conservative_end=SIGMA_E_HIGH,
               which_is_flattering="the LARGE end: it lowers the ceiling and inflates r_latent, so "
                                   "it flatters both 'we approach the ceiling' and 'the true "
                                   "correlation is higher than it looks'",
               ends=ends, span={})
    for key in ("rho", "ceiling_r_observed", "r_latent", "beta_latent"):
        vals = [e[key] for e in ends.values() if key in e]
        if len(vals) == 2:
            out["span"][key] = [min(vals), max(vals)]
    if r_obs is not None and all("r_latent" in e for e in ends.values()):
        imp = [e["r_latent_is_impossible"] for e in ends.values()]
        out["conclusion_depends_on_which_deposit"] = bool(imp[0] != imp[1])
    return out


def variance(values):
    n = len(values)
    if n < 2:
        raise KleinerError("need at least 2 stages to compute a variance")
    m = sum(values) / n
    return sum((v - m) ** 2 for v in values) / (n - 1)


def demo():
    var = GSE193066_58PAIR_STAGE_VAR
    a = across_range(var, r_obs=0.50)
    lo, hi = a["ends"]["conservative"], a["ends"]["anti_conservative"]
    print(f"""\
WORKED EXAMPLE — the verdict flips on which deposit you assume

Take the 58 GSE193066 PARTICIPANTS, whose Kleiner stage variance is
{var:.4f} stages^2 over their 116 reads — a measured number, not a round one — and
a reported correlation against the stage of r = 0.50. Run it at both ends of the
measured sigma_e range:

                            sigma_e      rho   ceiling   r_latent
  conservative (PXD A2)      {lo['sigma_e']:.4f}   {lo['rho']:.4f}    {lo['ceiling_r_observed']:.4f}     {lo['r_latent']:.4f}   unremarkable
  anti-cons. (GSE193066)     {hi['sigma_e']:.4f}   {hi['rho']:.4f}    {hi['ceiling_r_observed']:.4f}     {hi['r_latent']:.4f}   IMPOSSIBLE

**Same population, same reported correlation, opposite verdicts.** At sigma_e
{hi['sigma_e']:.4f} the disattenuated correlation is {hi['r_latent']:.4f}, greater than 1 and
therefore impossible. At sigma_e {lo['sigma_e']:.4f} it is {lo['r_latent']:.4f}, ordinary and
calling for no explanation at all.

That is the lesson, and it is sharper than the one this tool used to print. An
earlier version reported only the {hi['sigma_e']:.4f} end and announced "IMPOSSIBLE" in
capitals without saying that an independent deposit of the same label reverses
the verdict. The measurement error of the Kleiner stage did not replicate, so a
conclusion that depends on it is a conclusion about which deposit you picked.

If the anti-conservative row is the one you want to quote, the escapes are the
usual ones and they must be stated rather than assumed: the population may carry
more stage variance than this one, the reported r may be optimistic (in-sample,
tuned, or leaking), or the labels may have been averaged over reads, which lowers
the error variance and changes the arithmetic.

What transports:
  * report the stage variance of your sample NEXT TO any correlation against
    Kleiner, every time;
  * report the disattenuated value at BOTH ends, or measure sigma_e locally;
  * remember the LARGE end flatters you. Default to the small one.

CAUTION: the ceiling is an upper bound implied by measurement error alone. A
ceiling above your observed r means the result is DETECTABLE, not that an effect
exists.
""")
    print("Reliability of a single read, by stage SD, at both ends of the range:")
    print(f"  {'stage SD':>8} {'variance':>9} | {'rho@0.4655':>10} {'ceil':>6} "
          f"| {'rho@0.6433':>10} {'ceil':>6}")
    for s in (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5):
        v = s * s
        row = f"  {s:8.2f} {v:9.3f} |"
        for se in (SIGMA_E_LOW, SIGMA_E_HIGH):
            try:
                r = rho_from_var(v, se)
                row += f" {r:10.3f} {math.sqrt(r):6.3f} |"
            except KleinerError:
                row += f" {'undefined':>10} {'-':>6} |"
        print(row)

    print("\nMeasured deposits (not pooled):")
    print(f"  {'deposit':<28} {'n':>3} {'sigma_e':>8} {'CI95':>18} {'label ICC':>10} {'%F0-F1':>7}")
    for name, n, se, ci, icc, pct in DEPOSITS:
        cis = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "-"
        print(f"  {name:<28} {n:>3} {se:>8.4f} {cis:>18} {icc:>10.4f} {pct:>7.1f}")
    print(f"\n  pooling computed but NOT shipped: fixed {POOLING_NOT_SHIPPED['fixed_effect']}, "
          f"random {POOLING_NOT_SHIPPED['random_effects']}, I2 {POOLING_NOT_SHIPPED['I2_pct']}%")
    print(f"  {POOLING_NOT_SHIPPED['why_not']}")

    print(f"""
WHY THIS TOOL REPORTS sigma_e AND NOT A RELIABILITY COEFFICIENT

Across these deposits the Kleiner label's own ICC(1,1) moved about 3.2x
(0.1886 to 0.6228) while sigma_e moved 1.38x (0.4655 to 0.6433) — and the two
largest deposits have SIMILAR first-read stage variance (0.577 against 0.530).
So the ICC gap is a sigma_e gap, not a variance-denominator gap, and neither
number is a property of "the Kleiner stage".

Never write "Kleiner staging has ICC 0.19". Write a within-sample contrast at a
stated stage variance: "in this sample, whose stage variance was V, the label's
ICC was B."

sigma_e is TEST-RETEST across two physically distinct biopsies. It folds reader
variability, biopsy sampling variability and real change into one term. It is NOT
an inter-reader statistic.""")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Kleiner measurement-error disattenuation over a measured sigma_e RANGE.")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--stage-var", type=float, help="variance of Kleiner stages in YOUR sample")
    g.add_argument("--stage-sd", type=float, help="standard deviation of Kleiner stages")
    g.add_argument("--stages", help="text file of Kleiner stages, one per row")
    ap.add_argument("--r", type=float, help="observed correlation against the Kleiner stage")
    ap.add_argument("--beta", type=float, help="observed regression slope on the Kleiner stage")
    ap.add_argument("--sigma-e", type=float, default=None,
                    help="use ONE explicit sigma_e instead of reporting both ends. Without this the "
                         "tool reports the full measured range; the conservative default is "
                         f"{SIGMA_E_LOW:.4f}")
    ap.add_argument("--json", help="write the full result here")
    ap.add_argument("--demo", action="store_true", help="worked example and the rho tables")
    a = ap.parse_args(argv)

    if a.demo:
        return demo()
    if a.stage_var is not None:
        var = a.stage_var
    elif a.stage_sd is not None:
        var = a.stage_sd ** 2
    elif a.stages:
        vals = [float(ln) for ln in open(a.stages) if ln.strip()]
        var = variance(vals)
        print(f"{len(vals)} stages read; variance {var:.4f} (SD {math.sqrt(var):.4f})")
    else:
        ap.error("give one of --stage-var, --stage-sd, --stages, or --demo")

    print(f"stage variance   {var:.4f}  (SD {math.sqrt(var):.4f})")
    if a.sigma_e is not None:
        res = disattenuate(var, a.r, a.beta, sigma_e=a.sigma_e)
        tag = "CONSERVATIVE" if res["is_conservative_end"] else "NOT the conservative end"
        print(f"sigma_e          {a.sigma_e:.4f}  ({tag})")
        if not res["is_conservative_end"]:
            print(f"  ** you chose a sigma_e above the conservative {SIGMA_E_LOW:.4f}. A larger "
                  f"sigma_e lowers the ceiling and inflates r_latent, which FLATTERS the usual "
                  f"claims. State the value you used.")
        print(f"rho              {res['rho']:.4f}")
        print(f"ceiling on r_obs {res['ceiling_r_observed']:.4f}")
        if a.r is not None:
            print(f"r_latent         {res['r_latent']:.4f}"
                  + ("   ** IMPOSSIBLE, |r_latent| > 1" if res["r_latent_is_impossible"] else ""))
        if a.beta is not None:
            print(f"beta_latent      {res['beta_latent']:.4f}")
    else:
        res = across_range(var, a.r, a.beta)
        lo, hi = res["ends"]["conservative"], res["ends"]["anti_conservative"]
        print(f"sigma_e range    {SIGMA_E_LOW:.4f} to {SIGMA_E_HIGH:.4f} "
              f"(conservative end first; the LARGE end flatters)")
        print(f"{'':17}{'conservative':>14}{'anti-conservative':>20}")
        print(f"{'rho':17}{lo['rho']:>14.4f}{hi['rho']:>20.4f}")
        print(f"{'ceiling on r_obs':17}{lo['ceiling_r_observed']:>14.4f}"
              f"{hi['ceiling_r_observed']:>20.4f}")
        if a.r is not None:
            print(f"{'r_latent':17}{lo['r_latent']:>14.4f}{hi['r_latent']:>20.4f}")
            if res.get("conclusion_depends_on_which_deposit"):
                print("  ** the two ends DISAGREE on whether r_latent is possible. This input "
                        "cannot be interpreted without measuring sigma_e in your own sample.")
            elif hi["r_latent_is_impossible"]:
                print("  ** |r_latent| > 1 at BOTH ends: impossible however the range is read.")
        if a.beta is not None:
            print(f"{'beta_latent':17}{lo['beta_latent']:>14.4f}{hi['beta_latent']:>20.4f}")

    print("\nsigma_e is TEST-RETEST across two physically distinct biopsies. It folds reader "
          "variability,\nbiopsy sampling variability and real change into one term. It is NOT an "
          "inter-reader statistic.\nIf you have repeat reads, measure sigma_e in your own sample "
          "rather than using this range.")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(res, fh, indent=2)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KleinerError as e:
        print(f"kleiner.py: {e}", file=sys.stderr)
        sys.exit(2)
