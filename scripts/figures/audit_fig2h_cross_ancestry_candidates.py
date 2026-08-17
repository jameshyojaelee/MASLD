#!/usr/bin/env python3
"""Rank provisional Fig. 2H candidates within each non-European ancestry."""

from __future__ import annotations

import csv
import math
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
COLOC = Path(os.environ.get(
    "FIG2_COLOC_INPUT",
    ROOT / "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv",
))
OUT = Path(os.environ["FIG2H_AUDIT_DIR"])
REGISTRY = ROOT / "GWAS/finemapping/config/gwas_registry.tsv"
TIER = ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv"
GTF = ROOT / "data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf"
ANCESTRIES = ("AFR", "AMR", "EAS", "SAS")
TOP_PER_ANCESTRY = 25


def number(value: str | None) -> float | None:
    try:
        result = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def protein_coding_map(symbols: set[str]) -> dict[str, set[str]]:
    result = {symbol: set() for symbol in symbols}
    with GTF.open() as handle:
        for line in handle:
            if "\tgene\t" not in line or 'gene_type "protein_coding"' not in line:
                continue
            for symbol in symbols:
                if f'gene_name "{symbol}"' in line:
                    gene_id = line.split('gene_id "', 1)[1].split('"', 1)[0].split('.', 1)[0]
                    result[symbol].add(gene_id)
    return result


def exact_pvalues(rows: list[dict[str, object]], paths: dict[str, Path]) -> dict[tuple[str, str], float]:
    wanted: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        wanted[str(row["study"])].add(str(row["key"]))
    found: dict[tuple[str, str], float] = {}
    for study, keys in wanted.items():
        path = paths[study]
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write("\n".join(sorted(keys)) + "\n")
            key_file = handle.name
        command = (
            "awk -F'\\t' 'NR==FNR{w[$1]=1;next} FNR>1 && "
            "(($1\":\"$2) in w){print $1\":\"$2\"\\t" + '"$7}' +
            f"' {key_file!r} {str(path)!r}"
        )
        try:
            output = subprocess.check_output(command, shell=True, text=True)
        finally:
            Path(key_file).unlink(missing_ok=True)
        for line in output.splitlines():
            key, pvalue = line.split("\t")
            value = number(pvalue)
            if value is not None:
                old = found.get((study, key))
                found[(study, key)] = value if old is None else min(old, value)
    return found


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tiers = {row["study_name"]: row for row in read_tsv(TIER) if row["placement"] == "main"}
    registry = {row["study_name"]: row for row in read_tsv(REGISTRY)}
    candidates: list[dict[str, object]] = []
    gene_ids: dict[str, set[str]] = defaultdict(set)
    with COLOC.open(newline="") as handle:
        for row in csv.DictReader(handle):
            pp4 = number(row.get("PP.H4.susie"))
            ancestry = row.get("ancestry")
            study = row.get("gwas_name", "")
            if study not in tiers or ancestry not in ANCESTRIES or pp4 is None or pp4 <= 0.5:
                continue
            gene = row["gene"]
            ensembl = row["ensembl"].split('.', 1)[0]
            gene_ids[gene].add(ensembl)
            top_snp = row.get("top_snp", "")
            parts = top_snp.split(":")
            if len(parts) < 2:
                continue
            strict = (
                row.get("stratum_ld_reliability") == "reliable" and
                row.get("ld_reliability") == "locus_ok"
            )
            candidates.append({
                "gene": gene, "ensembl": ensembl, "ancestry": ancestry,
                "study": study, "trait": tiers[study]["trait"],
                "trait_class": tiers[study]["tier_label"], "top_snp": top_snp,
                "key": f"{parts[0]}:{parts[1]}", "multi_signal_pp_h4": pp4,
                "top_snp_posterior": number(row.get("top_snp_PP")),
                "ld_panel": row.get("ld_panel", ""), "ld_panel_n": row.get("ld_panel_n", ""),
                "stratum_ld_reliability": row.get("stratum_ld_reliability", ""),
                "locus_ld_reliability": row.get("ld_reliability", ""),
                "strict_diagnostics": strict,
            })

    # Query exact GWAS P values for the best rows by ancestry and diagnostic tier.
    shortlist: list[dict[str, object]] = []
    for ancestry in ANCESTRIES:
        rows = [row for row in candidates if row["ancestry"] == ancestry]
        rows.sort(key=lambda row: (
            not bool(row["strict_diagnostics"]),
            -float(row["multi_signal_pp_h4"]),
            str(row["gene"]), str(row["study"]),
        ))
        shortlist.extend(rows[:TOP_PER_ANCESTRY])
    paths = {
        study: ROOT / "GWAS/finemapping" / registry[study]["sumstats_path"]
        for study in {str(row["study"]) for row in shortlist}
    }
    pvalues = exact_pvalues(shortlist, paths)
    coding = protein_coding_map({str(row["gene"]) for row in shortlist})

    for row in shortlist:
        pvalue = pvalues.get((str(row["study"]), str(row["key"])))
        ids = gene_ids[str(row["gene"])]
        row["gwas_p_at_coloc_top_snp"] = pvalue
        row["gws"] = pvalue is not None and pvalue < 5e-8
        row["unique_coloc_transcript"] = len(ids) == 1
        row["protein_coding_target"] = ids == coding[str(row["gene"])] and len(ids) == 1
        row["automated_gate"] = all((
            row["gws"], row["strict_diagnostics"],
            row["unique_coloc_transcript"], row["protein_coding_target"],
        ))
        row["release_state"] = "provisional"

    shortlist.sort(key=lambda row: (
        ANCESTRIES.index(str(row["ancestry"])),
        not bool(row["automated_gate"]), not bool(row["gws"]),
        not bool(row["strict_diagnostics"]), -float(row["multi_signal_pp_h4"]),
        str(row["gene"]),
    ))
    fields = list(shortlist[0])
    with (OUT / "Fig2H_cross_ancestry_candidate_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(shortlist)

    with (OUT / "Fig2H_cross_ancestry_shortlist.md").open("w") as handle:
        handle.write("# Provisional Figure 2H cross-ancestry shortlist\n\n")
        handle.write("Current unpromoted COLOC release. P is evaluated at the COLOC top SNP.\n\n")
        for ancestry in ANCESTRIES:
            handle.write(f"## {ancestry}\n\n")
            handle.write("| Gene–trait | GWAS | P | PP.H4 | diagnostics | automated gate |\n")
            handle.write("|---|---|---:|---:|---|---|\n")
            rows = [row for row in shortlist if row["ancestry"] == ancestry][:8]
            for row in rows:
                pvalue = row["gwas_p_at_coloc_top_snp"]
                ptext = "NA" if pvalue is None else f"{float(pvalue):.2g}"
                diagnostics = f"{row['stratum_ld_reliability']}/{row['locus_ld_reliability']}"
                handle.write(
                    f"| {row['gene']}–{row['trait']} | {row['study']} | {ptext} | "
                    f"{float(row['multi_signal_pp_h4']):.3f} | {diagnostics} | "
                    f"{'pass' if row['automated_gate'] else 'qualified/no'} |\n"
                )


if __name__ == "__main__":
    main()
