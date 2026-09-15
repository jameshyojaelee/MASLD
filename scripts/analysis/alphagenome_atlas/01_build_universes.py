#!/usr/bin/env python3
"""Step 01: eligible signals for Universes A (direct), C (enzyme) and B (fine-mapped)."""

from __future__ import annotations

import gzip
import csv
from collections import defaultdict

import lib_atlas as la

P = la.prespec()
OUT = la.out_root() / "tables"
PAIRS = la.PROJECT / P["universes"]["A"]["source"]
U35 = la.PROJECT / "GWAS/finemapping/runs/uniform35_v4_2026-08-04"

COLS = [
    "signal_uid", "universe", "gwas_name", "trait", "trait_class", "tier", "ancestry", "gene", "ensembl",
    "signal_pair_index", "gwas_signal", "eqtl_signal", "pp_h4", "locus_id", "cs_id", "cs_name",
    "n_cs_variants", "cs_coverage", "posterior_definition", "source_path", "source_sha256",
]


def main() -> None:
    rows = []
    # Universes A and C: one row per gene-study pair (the table repeats each pair per lineage).
    pairs = la.read_tsv(PAIRS)
    seen = set()
    sha = la.sha256_file(PAIRS)
    counts = defaultdict(int)
    for r in pairs:
        key = (r["gwas_name"], r["ensembl"], r["signal_pair_index"])
        if key in seen:
            continue
        seen.add(key)
        if r["primary_signal_pair"] != "TRUE":
            raise la.ContractError(f"non-primary pair in primary table: {key}")
        universe = {"direct_MASLD": "A_direct", "liver_enzyme": "C_enzyme"}[r["trait_class"]]
        counts[universe] += 1
        rows.append({
            "signal_uid": f"coloc:{r['gwas_name']}:{r['ensembl']}:{r['signal_pair_index']}",
            "universe": universe, "gwas_name": r["gwas_name"], "trait": r["trait"],
            "trait_class": r["trait_class"], "tier": "1" if universe == "A_direct" else "2",
            "ancestry": r["ancestry"], "gene": r["gene"], "ensembl": r["ensembl"],
            "signal_pair_index": r["signal_pair_index"], "gwas_signal": r["gwas_signal"],
            "eqtl_signal": r["eqtl_signal"], "pp_h4": r["pp_h4"], "locus_id": "", "cs_id": "",
            "cs_name": "", "n_cs_variants": "", "cs_coverage": "",
            "posterior_definition": P["universes"]["A"]["posterior_definition"],
            "source_path": str(PAIRS.relative_to(la.PROJECT)), "source_sha256": sha,
        })
    if counts["A_direct"] != P["universes"]["A"]["expected_pairs"] or counts["C_enzyme"] != P["universes"]["C"]["expected_pairs"]:
        raise la.ContractError(f"pair counts differ from prespecification: {dict(counts)}")
    studies_a = sorted({r["gwas_name"] for r in rows if r["universe"] == "A_direct"})
    if studies_a != sorted(P["direct_trait_studies_expected_in_A"]):
        raise la.ContractError(f"direct-trait study list drifted: {studies_a}")

    # Universe B: tier-1 credible sets from the uniform35 run.
    cs_path = U35 / "aggregate/credible_sets.tsv.gz"
    manifest = {r["locus_id"]: r for r in la.read_tsv(U35 / "config/locus_manifest.tsv")}
    sha_cs = la.sha256_file(cs_path)
    n_b = 0
    with gzip.open(cs_path, "rt") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            m = manifest[r["locus_id"]]
            if m["tier"] != "1":
                continue
            n_b += 1
            rows.append({
                "signal_uid": f"susie:{r['locus_id']}:{r['cs_id']}", "universe": "B_direct",
                "gwas_name": r["study_name"], "trait": m["trait"], "trait_class": "direct_MASLD", "tier": "1",
                "ancestry": m["ancestry"], "gene": "", "ensembl": "", "signal_pair_index": "",
                "gwas_signal": "", "eqtl_signal": "", "pp_h4": "", "locus_id": r["locus_id"],
                "cs_id": r["cs_id"], "cs_name": r["cs_name"], "n_cs_variants": r["n_variants"],
                "cs_coverage": r["coverage"], "posterior_definition": P["universes"]["B"]["posterior_definition"],
                "source_path": str(cs_path.relative_to(la.PROJECT)), "source_sha256": sha_cs,
            })
    la.log(f"universes: A={counts['A_direct']} C={counts['C_enzyme']} B_credible_sets={n_b}")
    n = la.write_tsv_once(OUT / "eligible_signals_provisional.tsv", rows, COLS)
    la.log(f"wrote {n} rows")


if __name__ == "__main__":
    main()
