# Codex: concurrent supplied P2/P3 execution, 2026-09-30

Root submitted P2 as job **21998518/1586** and P3 as **21998519/1587**;
both were RUNNING at the root's 23:29 UTC scheduler check. This changes
execution order only. The original job 21985999 and its saved files are unchanged.

## Why separate jobs

`check_bfix_tests_v13.py` has no role, draw-range or restart CLI. It executes all
50 P1 draws, then all 50 P2 draws and profiles, then all 50 P3 draws. Waiting
for its P1 loop would delay the remaining agreement evidence.

The new `run_bfix_role_v17.py` selects the original common initialization and
the requested P2 or P3 statements using unique AST assignment boundaries.
It uses no fixed line offsets. It checks that extracted statement ASTs are
identical to the original statements. Semicolon-separated statements are
extracted by source column, preventing duplicated initialization.

No original runner, likelihood, design, or inference helper was edited.
The preparation does not read B-MODEL's held draws, planted classes,
biological allelic results, tag counts, or restricted inputs.

## Fixed scientific behavior

Both roles independently rerun the original observed null and marginal fits
from the original converged observed fit. They require B=50 and verify the
published checksums of the actual synthetic fixture/draw inputs. They preserve:

- All prescribed warm/cold starts, original likelihood tie rule, bounds and
  projected-gradient convergence criterion through the unchanged numerical
  modules. Current source includes the separately checked constrained polish
  and memory changes; the long original job keeps its earlier loaded source.
- Complete-participant recentering, cluster-copy handling, exact supplied P2
  resample/reliability pairing, invalid-draw behavior, linear quantiles, and
  original P2 tau profiles and component formulas.
- The exact supplied P3 uniform order, conditional inverse CDF, section-3
  refit and both section-4 starts, boundary restart, exact zero LR at v=0,
  original exceedance/failure handling and Monte Carlo SE.

Before running draws, a compute-only comparison requires the rerun observed
log likelihood and kappa, plus marginal v, to differ by at most 1e-6 from
`bfix_tests_21985999/observed_{null,marginal}.json`. It saves the differences
and refuses a discrepant baseline. This comparison has **not run**.

Each role saves all unchanged per-draw outputs, observed fits, its own final
role result, complete executed local source copies, generated role source,
source/input hashes, and environment versions in a new directory:
`bfix_p2_<jobid>/` or `bfix_p3_<jobid>/`. Existing directories are refused.
All source copies and the receipt are written before numerical imports or
fitting. Sources, fixture inputs and observed reference files are rechecked
for changes during snapshot creation; changed files cause refusal. Numerical
modules import from the archived directory.

The original eventual P2/P3 outputs must remain separate. Duplicate execution
does not double B or create an independent sample. These are supplied synthetic
agreement checks, not scientific size, power, coverage or independent validation.

## Resource header and commands

Two individual jobs; no arrays or dependency on the original job. Each requests:

```text
--account=nslab --qos=nslab --partition=cpu
--cpus-per-task=1 --mem=2G --time=90:00:00
```

Reserved names are 1586 (P2) and 1587 (P3). Previous refit RSS was approximately
53 MiB; the 2 GiB request leaves room for marginal/profile work. Runtime for
the complete new roles is unmeasured; 90 hours is the authorized buffered request.

Root reviewed and submitted these two launchers, with explicit per-job log
overrides. The base launcher commands are:

```bash
sbatch /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_bfix_p2_v17.sbatch
sbatch /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/in_biopsy_allelic/rec_codex/run_bfix_p3_v17.sbatch
```

Actual logs, overriding the launchers' generic names, are
`codex-rec-20260929T142434Z/1586_21998518.{out,err}` and
`1587_21998519.{out,err}`. Actual output directories are
`bfix_p2_21998518/` and `bfix_p3_21998519/` under that execution root.
Root recorded both actual IDs in the code directory's `job_index.tsv`.
This agent performed no submission or job-index edit.

## Checks executed

Both `run_bfix_role_v17.py --role p2 --inspect` and the corresponding P3 command
passed with the existing rnaseq Python. They verify input hashes, B=50,
AST identity and complete generated Python syntax, without importing numerical
modules, fitting, reading observed numerical values, or creating result directories.
P2 selected 17 top-level statements and P3 selected 11; both kept the same
28 common statements after omitting only original output-directory setup.
`bash -n` passed for both launchers. `git diff --check` passed for the new files.

Current numerical likelihood SHA256 is
`41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6`;
original runner SHA256 is
`a4ca55a52d28dfdbfeb0fc2611c2dcb3ec412987c8f7618f59002cbf72b08c2d`.
These checks establish isolated execution equivalence of the selected source
statements; numerical fitting and draw agreement remain to be checked on compute.
Root inspected the actual archived sources after submission and confirmed both
jobs were performing observed fits at 23:29 UTC. No numerical agreement check
or supplied draw was complete at that checkpoint.
