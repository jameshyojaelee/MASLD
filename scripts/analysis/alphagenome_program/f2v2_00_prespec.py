#!/usr/bin/env python3
"""f2-haplotype-v2 step 0: write the prespecification and stamp its sha256 BEFORE any model API call.

Written this way, not by hand, so that (a) written_utc is the machine clock at the moment the bytes hit
disk and can never be later than the file's own mtime, which is the defect this package repairs in v1's
seven amendment files, and (b) the sha256 digest file is created in the same second as the JSON, so the
digest's mtime is itself evidence of pre-registration.

Run on the login node before submitting the scoring job. It makes no network call.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import pathlib

OUT = pathlib.Path(os.environ["F2V2_OUT_ROOT"])
HERE = pathlib.Path(__file__).resolve().parent
V1 = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/"
                  "results/alphagenome_program/f2-haplotype-20260914T230433Z")


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


SPEC = {
    "package": "f2-haplotype-v2",
    "supersedes": "f2-haplotype-20260914T230433Z",
    "executes": "section 7 (F. Local haplotypes) of scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
    "written_utc": None,
    "state_at_writing": (
        "no predict_sequence call has been made in this package and raw/ does not exist. The only thing read "
        "from v1 so far is its design metadata: the GENCODE v49 chr6 gene widths, the chr6 common-EUR panel "
        "positions, the decile edges, the two defect descriptions, and v1's own already-published numbers. "
        "No residual produced by this package exists."
    ),
    "why_this_package_exists": (
        "An independent checker found two number-level defects in f2-haplotype-20260914T230433Z. Defect 1: the "
        "primary GNMT verdict is not a matched comparison in the RNA channel, because all six GNMT-region pairs "
        "read RNA out over the 3,194 bp GENCODE v49 GNMT span while all 100 of their matched chr6 null draws read "
        "RNA out over a central 20 kb window. Defect 2: the both-variants-active stratum p-values in "
        "tables/f3_summary_n2500.json are computed against the full pooled 2,500-draw null instead of against the "
        "stratum's own matched draws, which is exactly the arm-magnitude confound v1's amendment 06 exists to avoid."
    ),
    "defect_1_repair": {
        "what_is_wrong": (
            "f2_02_score.py score_pair sets the RNA readout to the target gene's GENCODE span when that span "
            "overlaps the 1,048,576 bp scored window and to the central 20 kb otherwise. Every GNMT-region pair "
            "is centred on GNMT, so all six take the gene-span branch (3,194 bp). Every chr6 null draw carries "
            "ensembl=ENSG00000124713 but sits at a random chr6 position, so the GNMT span falls outside its window "
            "and all 100 take the central-20 kb branch. The observed and null RNA sums are taken over windows of "
            "different width, position and expression level, so the RNA arm of F2.1 is unmatched."
        ),
        "fix": (
            "draw a fresh 100-pair chr6 null in which every draw's RNA readout is a real GENCODE v49 gene span, "
            "scored through the unchanged f2_recipe / score_pair code path by assigning each draw the Ensembl id "
            "of the gene that hosts it. No scoring function is edited; only which (gene, variant pair) goes in."
        ),
        "draw_rule": {
            "chromosome": "chr6, the observed locus's chromosome",
            "host_gene": (
                "a GENCODE v49 gene row on chr6 whose span width (end - start + 1) lies in [2472, 4791] bp. "
                "4791 = 1.5 x GNMT's 3,194 bp; 2472 is the floor imposed by the separation bin below, since both "
                "variants must fit inside the span. ENSG00000124713 (GNMT) and any gene whose span overlaps "
                "chr6:42960690-42963883 are excluded so the null cannot contain the observed locus."
            ),
            "variants": (
                "two biallelic SNVs from this project's common (EUR MAF >= 0.05) 1000 Genomes EUR panel, hg19 "
                "positions lifted to hg38 (tables/null_panel/chr6.tsv.gz of v1, reused unchanged), BOTH inside "
                "the host gene's span. REF is the GRCh38 base at the lifted position and must equal one of the "
                "two panel alleles; ALT is the other. Same placement rule as v1."
            ),
            "separation": (
                "inside the GNMT primary pair's own separation decile bin, bin 6 = [2471.2, 3916.3] bp, from the "
                "199-pair observed family's decile edges. The primary pair's separation is 2,503 bp."
            ),
            "window": "the 1,048,576 bp window centred on the pair midpoint must fit inside chr6.",
            "sampling": (
                "100 distinct host genes drawn WITHOUT replacement, uniformly over the host genes that carry at "
                "least one qualifying pair, then one qualifying pair drawn uniformly inside each gene. Seed 20260915. "
                "Without replacement so no single gene's span dominates the null."
            ),
            "feasibility_measured_before_writing_this": (
                "364 chr6 genes have span width in [2472, 4791]; 189 of them host at least one qualifying common-EUR "
                "pair, 6,365 qualifying pairs in total. 189 >= 100, so sampling genes without replacement is possible."
            ),
        },
        "which_null_is_primary_in_v2": (
            "the gene-span-matched null, for ALL FOUR channels: it matches the observed readout in the RNA channel "
            "and keeps the identical +/-1 kb local readout in ATAC, DNASE and H3K27ac, so one draw set serves every "
            "channel through one loader. v1's central-20kb null is reported beside it as the superseded comparison, "
            "never instead of it."
        ),
        "what_is_NOT_rescored": (
            "the six GNMT-region observed pairs and the 199-pair P5 family with its 2,500 corrected draws are not "
            "rescored. Their rows are copied from v1's raw/scored_rows.jsonl unchanged, because defect 1 is in the "
            "null's readout, not in the observed side, and defect 2 needs no new scoring at all."
        ),
    },
    "defect_2_repair": {
        "what_is_wrong": (
            "channel_verdict in f2_03_analyse.py computes p once against the full pooled null and then masks it for "
            "the both_variants_active stratum, so the stratum's p-values are against 2,500 pooled draws rather than "
            "against the 28 (RNA), 22 (ATAC), 13 (DNASE), 20 (H3K27ac) null draws that themselves meet the stratum."
        ),
        "fix": (
            "recompute the both-active p-values with the null restricted to the null draws meeting the same stratum "
            "selection, report p_floor = 1/(n_stratum_null + 1) per channel, and run BH at q = 0.10 within the "
            "stratum family. The recipe's MIN_MATCHED_NULL = 20 gate is kept, so a channel with fewer than 20 "
            "matched draws emits no verdict; the ungated count is reported beside it as a diagnostic, labelled."
        ),
    },
    "accounting_corrections_no_new_scoring": [
        "the cancelled job 21764903's call and checkpoint accounting, recounted from its log and from the checkpoint",
        "the claim that the prespecification hash was stamped before the first API call",
        "the amendment written_utc fields, all later than their own file mtimes and out of order for 03 and 04",
        "the both-active headline, which pairs one stage's counts with another's and calls the stratum untestable "
        "while the table beside it marks three of four channels testable",
        "the 0.25 percent by-construction stratum rate, which assumes the two arms are independent",
        "the headline sentence attributing the three RNA rejections to overlapping chromatin readout windows",
    ],
    "written_predictions": {
        "P1": (
            "F2.1 re-verdict. The GNMT primary pair's residual falls inside the central 95 percent of the "
            "gene-span-matched chr6 null in all four channels, that is, the v1 verdict survives the repair."
        ),
        "P2": (
            "The gene-span-matched null's RNA q95 |residual| is more than 10 percent SMALLER than the central-20kb "
            "null's. Mechanism: a 3 kb expressed gene span carries real RNA signal, so log2 ratios of its summed "
            "signal are less noisy than over a random 20 kb window on chr6 that mostly carries none."
        ),
        "P3": (
            "The three chromatin channels' q95 |residual| agree between the two draw sets within 25 percent "
            "relative, because those channels use the identical +/-1 kb readout in both sets and only the sampled "
            "positions differ. So the ATAC, DNASE and H3K27ac verdicts do not move."
        ),
        "P4": (
            "Defect 2 recomputed inside the stratum reproduces the checker's counts: BH q<0.10 rejections RNA 5, "
            "ATAC 0, DNASE 3, H3K27ac 0, with RNA's five p-values all at the floor 1/(28+1) = 0.03448. DNASE's 13 "
            "matched draws fail the recipe's 20-draw gate, so its 3 is reported as ungated only."
        ),
        "P5": (
            "The both-active stratum's measured null rate exceeds the 0.25 percent that independent arms would give, "
            "in every channel, and the null's two arm magnitudes are positively rank-correlated."
        ),
    },
    "decision_rules": {
        "F2.1": "inside the central 95 percent (2.5th to 97.5th percentile) of the primary null's residuals, per channel",
        "empirical_p": "(count of |null| >= |observed|, plus one) / (draws plus one), two-sided, unchanged from v1",
        "bh": "q = 0.10, step-up, family stated with every count",
        "p_floor": "1/(draws + 1), reported with every p and never called significant below it",
        "unchanged_from_v1": (
            "the recipe (f2_recipe.py), the 1 Mb window, the liver ontology terms, the +/-1 kb local readout, the "
            "EPS, the log2 definition, the residual definition joint - v1 - v2, and the observed rows themselves"
        ),
    },
    "seeds": {"null_draw": 20260915, "bootstrap": 20260914},
    "api_rate_ceiling": "27.9 predict_sequence calls per minute (2.15 s between call starts); the key is shared with job 21611396",
    "claim_boundary": (
        "model prediction only. Every number is the model's response to substitutions in one GRCh38 reference "
        "sequence. Nothing here measures an interaction in cells or people and nothing here adopts a Resource "
        "number, gene membership, threshold or figure."
    ),
    "v1_numbers_this_package_may_not_quietly_change": (
        "every verdict the checker confirmed stands. RESULTS.md of v2 states explicitly, number by number, which "
        "moved and which did not."
    ),
}


def main() -> None:
    (OUT / "prespec").mkdir(parents=True, exist_ok=True)
    SPEC["written_utc"] = utc()
    spec_path = HERE / "f2v2_00_prespec.json"
    spec_path.write_text(json.dumps(SPEC, indent=1) + "\n")
    dest = OUT / "prespec" / "f2v2_00_prespec.json"
    dest.write_text(spec_path.read_text())
    digest = hashlib.sha256(dest.read_bytes()).hexdigest()
    (OUT / "prespec" / "f2v2_00_prespec.sha256").write_text(f"{digest}  f2v2_00_prespec.json\n")
    stamp = {"written_utc": SPEC["written_utc"],
             "prespec_json_mtime_utc": dt.datetime.fromtimestamp(dest.stat().st_mtime, dt.timezone.utc)
                                        .strftime("%Y-%m-%dT%H:%M:%SZ"),
             "sha256_mtime_utc": dt.datetime.fromtimestamp((OUT / "prespec" / "f2v2_00_prespec.sha256").stat().st_mtime,
                                                           dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
             "sha256": digest,
             "note": ("written_utc is the machine clock at the moment the bytes were written, so it cannot be later "
                      "than the file's own mtime. Compare with v1, where all seven amendment written_utc fields are "
                      "later than their own mtimes.")}
    (OUT / "prespec" / "WRITTEN_UTC.json").write_text(json.dumps(stamp, indent=1) + "\n")
    print(json.dumps(stamp, indent=1))
    print(f"v1 package read-only reference: {V1}")


if __name__ == "__main__":
    main()
