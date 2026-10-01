#!/usr/bin/env python3
"""Write a v1-shaped input contract and expected counts for corrected inputs.

The v1 contract pins the 846-sample live bulk table and the pre-review COLOC
gene-level tables. This helper copies it, replaces the bulk and COLOC rows with
the corrected artifacts (sha256 and row counts read at run time), and writes the
frozen-contract counts that 01_freeze_and_rederive.py --expected checks, for the
same --state-rule and --untestable-state that 01 is run with.
The counts are computed here with pandas, separately from 01's csv loops, so
01 must reproduce them on the same files.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pandas as pd

from genetics_common import MISSING_STRINGS, PROJECT_ROOT, SCRIPT_DIR, inspect_table, sha256_file

UNTESTABLE_METHOD = "susie_untestable_insufficient_shared_posterior"
PRIMARY_CHROMOSOMES = {f"chr{n}" for n in range(1, 23)} | {"chrX", "chrY", "chrM"}

REPLACED = {
    "canonical_treat": "bulk",
    "tier12_genelevel": "tier12",
    "full_genelevel": "full",
    "all_gwas_coloc": "all_gwas",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bulk", type=Path, required=True, help="corrected F_five deg_results.csv")
    parser.add_argument("--tier12", type=Path, required=True, help="gene_level_coloc_tier12.csv")
    parser.add_argument("--full", type=Path, required=True, help="gene_level_coloc.csv")
    parser.add_argument("--all-gwas", type=Path, required=True, help="susie_coloc_all_gwas.csv")
    parser.add_argument(
        "--template", type=Path, default=SCRIPT_DIR / "config" / "input_contract_v1.tsv"
    )
    parser.add_argument(
        "--tier-config", type=Path,
        default=PROJECT_ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv",
    )
    parser.add_argument(
        "--gene-metadata", type=Path,
        default=PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz",
    )
    parser.add_argument("--state-rule", choices=["treat", "canonical"], default="treat")
    parser.add_argument("--untestable-state", action="store_true")
    parser.add_argument("--contract-out", type=Path, required=True)
    parser.add_argument("--expected-out", type=Path, required=True)
    return parser.parse_args()


def missing(series: pd.Series) -> pd.Series:
    return series.fillna("").str.strip().isin(MISSING_STRINGS)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.where(~missing(series)), errors="coerce")


def expected_counts(
    bulk_path: Path,
    tier12_path: Path,
    full_path: Path,
    rule: str = "treat",
    all_gwas_path: Path | None = None,
    tier_config: Path | None = None,
    gene_metadata: Path | None = None,
) -> dict[str, int]:
    bulk = pd.read_csv(bulk_path, dtype=str, keep_default_na=False)
    if rule == "treat":
        state = numeric(bulk["treat_fdr"]) < 0.05
    else:
        state = (numeric(bulk["padj"]) < 0.05) & (numeric(bulk["logFC"]).abs() > 0.5)
    base = bulk["gene"].str.strip().str.split(".").str[0]
    bulk = bulk.assign(symbol=bulk["symbol"].str.strip(), state=state, base=base)
    # Shared bulk-row rule: per symbol, a primary-assembly row first, then the
    # rule's lowest q-value, larger |t|, gene ID. The kept row gives the call.
    meta = pd.read_csv(gene_metadata, sep="\t", usecols=["gene_id", "chromosome"], dtype=str)
    primary_ids = set(meta.loc[meta["chromosome"].isin(PRIMARY_CHROMOSOMES), "gene_id"].str.split(".").str[0])
    q = numeric(bulk["treat_fdr" if rule == "treat" else "padj"]).fillna(float("inf"))
    named_bulk = bulk.assign(
        off_primary=~bulk["base"].isin(primary_ids), q=q, neg_abs_t=-numeric(bulk["t"]).fillna(0).abs()
    )[~missing(bulk["symbol"])]
    kept = named_bulk.sort_values(["symbol", "off_primary", "q", "neg_abs_t", "gene"]).drop_duplicates("symbol")
    bulk_symbols = set(kept["symbol"])
    state_symbols = set(kept.loc[kept["state"], "symbol"])
    bulk_symbol_by_base = dict(zip(named_bulk["base"], named_bulk["symbol"]))

    full = pd.read_csv(full_path, dtype=str, keep_default_na=False)
    full = full[~missing(full["ensembl"]) & ~missing(full["gene"])]
    full = full.assign(base=full["ensembl"].str.strip().str.split(".").str[0], gene=full["gene"].str.strip())
    symbol_map = full.groupby("base")["gene"].agg(lambda values: sorted(set(values)))

    # COLOC genes join bulk by Ensembl ID; their symbol is the bulk symbol when
    # the ID has a bulk row, else the COLOC symbol.
    tier = pd.read_csv(tier12_path, dtype=str, keep_default_na=False)
    tier = tier[~missing(tier["ensembl"])]
    tier = tier.assign(base=tier["ensembl"].str.strip().str.split(".").str[0])
    in_bulk = tier["base"].isin(bulk_symbol_by_base)
    tier = tier[in_bulk | tier["base"].isin(symbol_map.index)]
    in_bulk = tier["base"].isin(bulk_symbol_by_base)
    ambiguous = [b for b in tier.loc[~in_bulk, "base"] if len(symbol_map[b]) > 1]
    if ambiguous:
        raise SystemExit(f"ambiguous ensembl->symbol map: {ambiguous[:5]}")
    tier = tier.assign(symbol=[
        bulk_symbol_by_base[b] if b in bulk_symbol_by_base else symbol_map[b][0] for b in tier["base"]
    ])
    susie = numeric(tier["coloc_best_susie_pp4"])
    primary = set(tier.loc[susie > 0.5, "symbol"])
    testable = primary & bulk_symbols
    counts = {
        f"{rule}_rows": int(bulk["state"].sum()),
        f"{rule}_symbols": len(state_symbols),
        "primary_susie_named": len(primary),
        "primary_susie_bulk_testable": len(testable),
        "primary_overlap": len(testable & state_symbols),
    }
    if all_gwas_path is not None:
        tiers = pd.read_csv(tier_config, sep="\t", dtype=str)
        tier12 = set(tiers.loc[tiers["tier"].str.strip().isin(["1", "2"]), "study_name"])
        master = pd.read_csv(all_gwas_path, usecols=["gwas_name", "ensembl", "method"], dtype=str)
        hit = master[(master["method"] == UNTESTABLE_METHOD) & master["gwas_name"].isin(tier12)]
        bases = set(hit["ensembl"].dropna().str.strip().str.split(".").str[0])
        untestable_symbols = set(tier.loc[tier["base"].isin(bases), "symbol"]) - primary
        counts["genetic_untestable_named"] = len(untestable_symbols)
    return counts


def main() -> None:
    args = parse_args()
    for output in (args.contract_out, args.expected_out):
        if output.exists():
            raise SystemExit(f"refusing to overwrite {output}")
    sources = {"bulk": args.bulk, "tier12": args.tier12, "full": args.full, "all_gwas": args.all_gwas}
    with args.template.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    seen = set()
    for row in rows:
        key = REPLACED.get(row["input_id"])
        if key is None:
            continue
        seen.add(row["input_id"])
        path = sources[key].resolve()
        _, n_rows = inspect_table(path, row["format"])
        try:
            row["relative_path"] = str(path.relative_to(PROJECT_ROOT))
        except ValueError:
            row["relative_path"] = str(path)
        row["expected_sha256"] = sha256_file(path)
        row["expected_data_rows"] = str(n_rows)
    if seen != set(REPLACED):
        raise SystemExit(f"template lacks rows: {sorted(set(REPLACED) - seen)}")

    args.contract_out.parent.mkdir(parents=True, exist_ok=True)
    with args.contract_out.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = expected_counts(
        args.bulk, args.tier12, args.full, args.state_rule,
        args.all_gwas if args.untestable_state else None, args.tier_config, args.gene_metadata,
    )
    with args.expected_out.open("x", encoding="utf-8") as handle:
        json.dump(counts, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"PASS: contract={args.contract_out}; expected={json.dumps(counts, sort_keys=True)}")


if __name__ == "__main__":
    main()
