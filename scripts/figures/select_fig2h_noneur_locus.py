#!/usr/bin/env python3
"""Audit the fixed Fig. 2H non-European locus priority without promoting it.

The corrected COLOC table and non-European GWAS audit may be supplied through
FIG2_COLOC_INPUT and FIG2_NONEUR_AUDIT.  A candidate is never release-eligible
unless FIG2_RELEASE_STATE=promoted; promotion itself remains a manual,
synchronized Resource-release action.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
COLOC = Path(os.environ.get(
    "FIG2_COLOC_INPUT",
    ROOT / "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv",
))
AUDIT = Path(os.environ.get(
    "FIG2_NONEUR_AUDIT",
    ROOT / "RNA-seq/results/coloc_variant_classes/noneur_gws_audit.csv",
))
GTF = Path(os.environ.get(
    "FIG2_GENCODE_GTF",
    ROOT / "data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf",
))
OUT = Path(os.environ.get("FIG2_CANDIDATE_DIR", ROOT / "figures/candidates"))
RELEASE_STATE = os.environ.get("FIG2_RELEASE_STATE", "provisional").lower()

PRIORITY = (
    {"priority": 1, "gene": "ACADS", "trait": "GGT", "study": "BBJ_GGT", "ancestry": "EAS"},
    {"priority": 2, "gene": "EFHD1", "trait": "AST", "study": "BBJ_AST", "ancestry": "EAS"},
    {"priority": 3, "gene": "LRRC14", "trait": "ALT", "study": "MVP_ALT_AFR", "ancestry": "AFR"},
)


def finite(value: str | None) -> float | None:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def eligible_gencode_genes(path: Path, wanted: set[str]) -> dict[str, set[str]]:
    found = {gene: set() for gene in wanted}
    with path.open() as handle:
        for line in handle:
            if "\tgene\t" not in line or 'gene_type "protein_coding"' not in line:
                continue
            for gene in wanted:
                if f'gene_name "{gene}"' in line:
                    attrs = line.rstrip().split("\t", 8)[-1]
                    gene_id = attrs.split('gene_id "', 1)[1].split('"', 1)[0].split(".", 1)[0]
                    found[gene].add(gene_id)
    return found


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    coloc = read_rows(COLOC)
    audit = {row["gene"]: row for row in read_rows(AUDIT)}
    genes = {item["gene"] for item in PRIORITY}
    eligible = eligible_gencode_genes(GTF, genes)
    output: list[dict[str, object]] = []

    for item in PRIORITY:
        rows = [row for row in coloc if row.get("gene") == item["gene"]]
        exact = [row for row in rows if row.get("gwas_name") == item["study"]]
        row = max(exact, key=lambda x: finite(x.get("PP.H4.susie")) or -1) if exact else {}
        ar = audit.get(item["gene"], {})
        pp4 = finite(row.get("PP.H4.susie"))
        min_p = finite(ar.get("min_p"))
        ensembl_ids = {r.get("ensembl", "").split(".", 1)[0] for r in rows if r.get("ensembl")}
        unique_target = len(ensembl_ids) == 1 and ensembl_ids == eligible[item["gene"]]
        p_provenance = ar.get("best_gwas") == item["study"]
        locus_diag = row.get("ld_reliability") == "locus_ok"
        stratum_diag = row.get("stratum_ld_reliability") == "reliable"
        gates = {
            "non_eur_gws": bool(min_p is not None and min_p < 5e-8 and p_provenance),
            "multi_signal_pp4": bool(pp4 is not None and pp4 > 0.5),
            "reliable_finemapping": bool(locus_diag and stratum_diag),
            "unique_eligible_transcript": unique_target,
            # Readability is confirmed only after a one-page panel is rendered and reviewed.
            "readable_locus_architecture": False,
        }
        reasons = [name for name, passed in gates.items() if not passed]
        output.append({
            **item,
            "min_non_eur_gwas_p": min_p,
            "p_value_source_study": ar.get("best_gwas", ""),
            "multi_signal_pp_h4": pp4,
            "lead_variant": row.get("top_snp", ""),
            "ld_panel": row.get("ld_panel", ""),
            "ld_panel_n": row.get("ld_panel_n", ""),
            "stratum_ld_reliability": row.get("stratum_ld_reliability", ""),
            "locus_ld_reliability": row.get("ld_reliability", ""),
            "ensembl_ids": ";".join(sorted(ensembl_ids)),
            **gates,
            "automated_gates_pass": all(gates[name] for name in gates if name != "readable_locus_architecture"),
            "release_state": RELEASE_STATE,
            "selection_status": "blocked_pending_render_review" if not reasons[:-1] else "blocked_failed_gate",
            "blocking_reasons": ";".join(reasons),
        })

    fields = list(output[0])
    out_path = OUT / "Fig2H_noneur_candidate_audit.tsv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)

    leader = next((row for row in output if row["automated_gates_pass"]), None)
    status_path = OUT / "Fig2H_SELECTION_BLOCKED.txt"
    message = (
        f"Provisional leader: {leader['gene']}–{leader['trait']} ({leader['ancestry']}).\n"
        if leader else "No candidate passes the automated gates.\n"
    )
    message += (
        "2H remains blocked until corrected COLOC promotion, candidate rendering, "
        "architecture review, synchronized Resource release, and explicit approval.\n"
    )
    status_path.write_text(message)
    print(message, end="")
    print(f"Audit: {out_path}")


if __name__ == "__main__":
    main()
