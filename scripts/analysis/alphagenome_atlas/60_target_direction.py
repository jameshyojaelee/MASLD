#!/usr/bin/env python3
"""Step 60 (P6a): genetic direction of protection for MASLD drug targets, measured and predicted.

For every target gene with a MASLD agent (`data/external/druggability/mash_clinical_pipeline.tsv`, 45 agents) or a
MASLD development status (`data/external/drug_targets/drug_target_classification.tsv`: masld_approved / clinical /
preclinical), and for every Track 0 signal on that gene, the posterior mass is split into

  "lower expression protective"  : allele1 raises risk (gwas beta > 0) and raises expression (beta > 0), or the mirror
  "higher expression protective" : allele1 raises risk and lowers expression, or the mirror

separately for the MEASURED liver eQTL effect and for the ATLAS predicted liver RNA effect of the same allele
(`direction_variant_level.tsv.gz`, which already carries both oriented to allele1). The agent's own direction is
parsed from its mechanism string (inhibitor/antagonist/degrader/siRNA -> lowering; agonist/activator -> raising).
A target is concordant when the mass-dominant genetic direction of protection matches the agent's direction.

Scope: this is a SUPPLEMENTARY analysis and a Catalog-sidecar input. Figure 5 is molecular and physical tissue context,
not a drug-target figure, and the scope validator rejects therapeutic/calibration/top-N panels in the Figure 6 manifest
(drug content is retired to the legacy portal). No p-values: a categorical, coverage-reported table over a fixed list. Coding-mechanism targets
(PNPLA3, TM6SF2, GCKR, MTARC1, SERPINA1) are expected to be unresolved through an expression channel.

Outputs (tables/): target_direction_table.tsv, target_direction_variants.tsv.gz, target_direction_summary.json
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict

import pandas as pd

import lib_atlas as la

ROOT = la.out_root()
TABLES = ROOT / "tables"
TRACK0 = la.track0_root() / "tables"
PIPELINE = la.PROJECT / "data/external/druggability/mash_clinical_pipeline.tsv"
CLASSIFICATION = la.PROJECT / "data/external/drug_targets/drug_target_classification.tsv"
MASLD_STATUS = {"masld_approved", "masld_clinical", "masld_preclinical", "masld_discontinued"}
LOWERING = re.compile(r"inhibitor|antagonist|degrader|sirna|antisense|silenc|blocker|knockdown", re.I)
RAISING = re.compile(r"agonist|activator|analog|analogue|mimetic|replacement", re.I)
CODING_MECHANISM = {"PNPLA3", "TM6SF2", "GCKR", "MTARC1", "SERPINA1", "HSD17B13", "APOE", "ADH1B", "HFE", "CIDEB"}


def agent_direction(mechanism: str) -> str:
    low, high = bool(LOWERING.search(mechanism)), bool(RAISING.search(mechanism))
    if low and not high:
        return "lower"
    if high and not low:
        return "raise"
    return "unresolved"


MIN_SCORABLE = 10


def concordance_rate(n_concordant: int, n_discordant: int, min_scorable: int = MIN_SCORABLE) -> dict:
    """Concordance rate, or a withheld rate when too few targets can be scored.

    The funnel caps this arm at the number of targets whose agent direction is readable, and the full
    archive leaves 3 scorable. A rate on that denominator is the number that would be quoted, so it is
    withheld rather than printed beside a caveat; the deliverable is the per-target listing.
    """
    n = int(n_concordant) + int(n_discordant)
    if n < min_scorable:
        return {"rate": None, "n_scorable": n, "min_scorable": min_scorable,
                "reason": (f"withheld: {n} scorable targets is below the prespecified minimum of {min_scorable}; "
                           "the deliverable is a per-target listing and no rate is quoted")}
    return {"rate": float(n_concordant) / n, "n_scorable": n, "min_scorable": min_scorable, "reason": "reported"}


def protective_direction(gwas_beta: float, effect_beta: float) -> str:
    """Which direction of gene expression is protective, given allele1's risk effect and allele1's expression effect."""
    if not (gwas_beta == gwas_beta and effect_beta == effect_beta) or gwas_beta == 0 or effect_beta == 0:
        return "unresolved"
    return "lower" if (gwas_beta > 0) == (effect_beta > 0) else "raise"


def main() -> None:
    agents = defaultdict(list)
    for r in la.read_tsv(PIPELINE):
        agents[r["gene_symbol"].upper()].append(r)
    cls = {}
    for r in la.read_tsv(CLASSIFICATION):
        if r["drug_dev_status"] in MASLD_STATUS or r["symbol"].upper() in agents:
            cls[r["symbol"].upper()] = r
    targets = sorted(set(agents) | set(cls))
    signals = {s["signal_uid"]: s for s in la.read_tsv(TRACK0 / "eligible_signals.tsv")}
    by_gene = defaultdict(list)
    for s in signals.values():
        by_gene[s["gene"].upper()].append(s)

    var_rows = []
    with la.open_text(TRACK0 / "direction_variant_level.tsv.gz") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            s = signals.get(r["signal_uid"])
            if s is None or s["gene"].upper() not in targets:
                continue
            w = float(r["weight"])
            gb = float(r["gwas_beta_allele1"]) if r.get("gwas_beta_allele1") not in (None, "", "nan") else math.nan
            eq = float(r["eqtl_beta_allele1"]) if r.get("eqtl_beta_allele1") not in (None, "", "nan") else math.nan
            pr = float(r["pred_allele1_liver_rna"]) if r.get("pred_allele1_liver_rna") not in (None, "", "nan") else math.nan
            var_rows.append({"gene": s["gene"], "signal_uid": r["signal_uid"], "universe": s["universe"], "trait": s["trait"], "trait_class": s["trait_class"],
                             "variant_uid": r["variant_uid"], "weight": w, "gwas_beta_allele1": gb, "eqtl_beta_allele1": eq, "pred_allele1_liver_rna": pr,
                             "measured_protective_direction": protective_direction(gb, eq), "predicted_protective_direction": protective_direction(gb, pr),
                             "posterior_definition": s["posterior_definition"]})
    if var_rows:
        pd.DataFrame(var_rows).to_csv(TABLES / "target_direction_variants.tsv.gz", sep="\t", index=False)

    mass = defaultdict(lambda: defaultdict(float))
    for r in var_rows:
        mass[r["gene"].upper()][("measured", r["measured_protective_direction"])] += r["weight"]
        mass[r["gene"].upper()][("predicted", r["predicted_protective_direction"])] += r["weight"]

    rows = []
    for g in targets:
        ags = agents.get(g, [])
        dirs = {agent_direction(a["mechanism"]) for a in ags}
        agent_dir = dirs.pop() if len(dirs) == 1 else ("mixed" if len(dirs) > 1 else "no_agent")
        m = mass.get(g, {})
        row = {"gene": g, "n_agents": len(ags), "agents": ";".join(f"{a['drug']}({a['mechanism']},phase {a['max_phase']},{a['status']})" for a in ags),
               "agent_direction": agent_dir, "drug_dev_status": cls.get(g, {}).get("drug_dev_status", ""), "max_phase_masld": cls.get(g, {}).get("max_phase_masld", ""),
               "n_signals": len(by_gene.get(g, [])), "signals": ";".join(sorted({s["gwas_name"] for s in by_gene.get(g, [])})),
               "coding_mechanism_expected_unresolved": g in CODING_MECHANISM}
        for src in ("measured", "predicted"):
            lo, hi, un = m.get((src, "lower"), 0.0), m.get((src, "raise"), 0.0), m.get((src, "unresolved"), 0.0)
            tot = lo + hi + un
            row[f"{src}_mass_lower_protective"] = lo
            row[f"{src}_mass_raise_protective"] = hi
            row[f"{src}_mass_unresolved"] = un
            row[f"{src}_dominant"] = "" if tot == 0 else ("lower" if lo > max(hi, un) else "raise" if hi > max(lo, un) else "unresolved")
            row[f"{src}_verdict"] = ("no_genetic_signal" if tot == 0 else
                                     "unresolved" if row[f"{src}_dominant"] == "unresolved" or agent_dir in ("no_agent", "mixed", "unresolved") else
                                     "concordant" if row[f"{src}_dominant"] == agent_dir else "discordant")
        row["note"] = "COLOC/SuSiE posterior mass; colocalized, not causal; a coding-mechanism target cannot be resolved through an expression channel"
        rows.append(row)
    rows.sort(key=lambda r: (-r["n_agents"], -r["n_signals"], r["gene"]))
    la.write_tsv_once(TABLES / "target_direction_table.tsv", rows, list(rows[0].keys()))
    with_sig = [r for r in rows if r["measured_verdict"] != "no_genetic_signal"]
    # The denominator this table can actually speak to is targets that have BOTH an agent with a readable
    # direction AND a resolved genetic direction. Most rows carry neither: the classification file lists
    # thousands of genes while the MASH pipeline lists 45 agents. The funnel is reported so the reader sees
    # how few targets any concordance statement rests on.
    with_agent_dir = [r for r in rows if r["agent_direction"] in ("lower", "raise")]
    scorable = [r for r in with_agent_dir if r["measured_verdict"] in ("concordant", "discordant")]
    summary = {
        "funnel": {"targets_listed": len(rows),
                   "with_any_agent": sum(r["n_agents"] > 0 for r in rows),
                   "with_readable_agent_direction": len(with_agent_dir),
                   "with_a_genetic_signal": len(with_sig),
                   "with_agent_direction_AND_resolved_genetic_direction": len(scorable)},
        "n_targets": len(rows), "n_with_agent": sum(r["n_agents"] > 0 for r in rows), "n_with_genetic_signal": len(with_sig),
        "verdicts_measured": {v: sum(r["measured_verdict"] == v for r in rows) for v in ("concordant", "discordant", "unresolved", "no_genetic_signal")},
        "verdicts_predicted": {v: sum(r["predicted_verdict"] == v for r in rows) for v in ("concordant", "discordant", "unresolved", "no_genetic_signal")},
        "prediction_written": (">= 60% of targets with a readable agent direction AND a resolved genetic direction are "
                               "concordant; PNPLA3/TM6SF2/GCKR unresolved (coding). The denominator was written as "
                               "'targets with a regulatory genetic signal' in the plan, which is wrong: most such targets "
                               "have no agent and can never be scored. Corrected here before any full-run number was read."),
        "denominator_for_the_observed_concordance": "targets with a readable agent direction and a verdict of concordant or discordant",
        "cannot_conclude_if": "fewer than 10 targets reach the scorable denominator, in which case the table is a per-target listing and no rate is quoted",
        "observed_concordance_measured": concordance_rate(
            sum(r["measured_verdict"] == "concordant" for r in with_agent_dir),
            sum(r["measured_verdict"] == "discordant" for r in with_agent_dir)),
        "observed_concordance_predicted": concordance_rate(
            sum(r["predicted_verdict"] == "concordant" for r in with_agent_dir),
            sum(r["predicted_verdict"] == "discordant" for r in with_agent_dir)),
        "scorable_targets": [{"gene": r["gene"], "agent_direction": r["agent_direction"],
                              "measured_verdict": r["measured_verdict"], "predicted_verdict": r["predicted_verdict"],
                              "drug_dev_status": r["drug_dev_status"]}
                             for r in with_agent_dir if r["measured_verdict"] in ("concordant", "discordant")],
        "coding_targets_resolved": [r["gene"] for r in rows if r["coding_mechanism_expected_unresolved"] and r["measured_verdict"] in ("concordant", "discordant")],
        "family": "the fixed target list; categorical, no multiple-testing family",
    }
    json.dump(summary, (TABLES / "target_direction_summary.json").open("w"), indent=1, default=float)
    la.log(f"P6a: {summary['n_targets']} targets, {summary['n_with_genetic_signal']} with a signal; measured {summary['verdicts_measured']}")


if __name__ == "__main__":
    main()
