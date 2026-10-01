"""Collect B-Dup per-library counts into one table and a per-cohort summary.

Inputs
  --libraries   libraries.tsv from 01_regions.py (task, cohort, run, individual_id, bam, star_log)
  --restricted  directory holding per_library/<run>.tsv and markdup_stats/<run>.txt
                from 03_dup_library.sbatch
Outputs
  <restricted>/dup_per_library.tsv   one row per finished library: STAR read totals,
      gene-body subset counts, tag-exon read and duplicate counts, markdup library-size
      estimate. Per-person technical values; stays in the restricted tree.
  <exec>/dup_per_cohort.tsv          per cohort (and all_paired / all_single): libraries expected
      and finished, layout, median and quartiles of the tag-exon duplicate fraction.
      The full-20260929T134856Z output predates this and has one 'all' row pooling both
      layouts; it is superseded by 08_collect_v2.py's dup_per_cohort_role.tsv.
  <exec>/dup_missing_libraries.tsv   libraries without a finished row (empty table if none).
Rules
  Counts only. No allele, genotype or stage column is read or written. Refuses to
  overwrite an existing output.
Column meanings
  star_*                    from STAR Log.final.out; for paired-end libraries STAR counts pairs.
  subset_reads              primary MAPQ-255 reads overlapping tag-gene bodies (mates kept together).
  exon_reads                subset reads overlapping merged tag-gene exons (each read once).
  exon_dup_reads            of those, flagged duplicate by samtools markdup.
  exon_dup_fraction         exon_dup_reads / exon_reads. Pair rule for paired-end, one-end
                            rule for single-end: compare within one layout only.
  pg_programs               @PG program names; runs before 2026-09-29 read the header
                            without --no-PG, so 'samtools' there is the viewer itself.
  exon_mate_absent_reads    paired reads whose mate was outside the subset; markdup
                            treated these as single reads.
  genebody_dup_fraction     markdup DUPLICATE TOTAL / EXAMINED over the whole subset.
  estimated_library_size    markdup ESTIMATED_LIBRARY_SIZE for the subset (distinct
                            pairs the library would yield in these genes at infinite depth).
                            markdup computes it from pairs only, so it is NA for
                            single-end libraries (markdup prints 0 there).
  genebody_lw_distinct_reads  the same Lander-Waterman estimate applied to reads, so
                            every library has a value: X solving U = X * (1 - exp(-N / X)),
                            N = markdup_examined, U = N - markdup_duplicate_total. It assumes
                            uniform sampling of fragments, which RNA-seq does not satisfy
                            (expression is skewed), so read it as a within-layout ranking.
"""
import argparse
import glob
import math
import os
import sys

import pandas as pd

STAR_FIELDS = {
    "Number of input reads": "star_input_reads",
    "Uniquely mapped reads number": "star_unique_reads",
    "Number of reads mapped to multiple loci": "star_multi_reads",
}
MARKDUP_FIELDS = {
    "EXAMINED": "markdup_examined",
    "PAIRED": "markdup_paired",
    "SINGLE": "markdup_single",
    "DUPLICATE TOTAL": "markdup_duplicate_total",
    "ESTIMATED_LIBRARY_SIZE": "estimated_library_size",
}


def read_key_values(path, sep, fields):
    values = {}
    with open(path) as fh:
        for line in fh:
            if sep not in line:
                continue
            key, value = (part.strip() for part in line.split(sep, 1))
            if key in fields:
                values[fields[key]] = int(value)
    missing = set(fields.values()) - set(values)
    if missing:
        sys.exit(f"{path}: missing {sorted(missing)}")
    return values


def lander_waterman_size(n_reads, n_distinct):
    """X solving n_distinct = X * (1 - exp(-n_reads / X)) by bisection (Picard's formula).

    Returns NA when there are no duplicates (size unbounded) or no reads.
    """
    if n_reads <= 0 or n_distinct <= 0 or n_distinct >= n_reads:
        return math.nan

    def f(x):
        return x * (1 - math.exp(-n_reads / x)) - n_distinct

    lo, hi = float(n_distinct), float(n_distinct)
    while f(hi) < 0:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 0.5:
            break
    return round((lo + hi) / 2)


def summarise(group):
    frac = group["exon_dup_fraction"]
    return pd.Series({
        "n_libraries_finished": len(group),
        "n_paired": int((group["layout"] == "paired").sum()),
        "n_single": int((group["layout"] == "single").sum()),
        "n_headers_naming_a_dup_tool": int(group["pg_programs"].str.contains(
            "MarkDuplicates|markdup|dedup|UMI", case=False, regex=True).sum()),
        "preexisting_dup_flags_total": int(group["preexisting_dup_flags_in_subset"].sum()),
        "star_unique_reads_median": group["star_unique_reads"].median(),
        "exon_reads_median": group["exon_reads"].median(),
        "exon_reads_min": int(group["exon_reads"].min()),
        "exon_dup_fraction_median": frac.median(),
        "exon_dup_fraction_q25": frac.quantile(0.25),
        "exon_dup_fraction_q75": frac.quantile(0.75),
        "exon_dup_fraction_min": frac.min(),
        "exon_dup_fraction_max": frac.max(),
        "genebody_dup_fraction_median": group["genebody_dup_fraction"].median(),
        "exon_mate_absent_fraction_median": group["exon_mate_absent_fraction"].median(),
        "estimated_library_size_median": group["estimated_library_size"].median(),
        "genebody_lw_distinct_reads_median": group["genebody_lw_distinct_reads"].median(),
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--libraries", required=True)
    ap.add_argument("--restricted", required=True)
    ap.add_argument("--exec-dir", required=True)
    args = ap.parse_args()

    out_lib = os.path.join(args.restricted, "dup_per_library.tsv")
    out_cohort = os.path.join(args.exec_dir, "dup_per_cohort.tsv")
    out_missing = os.path.join(args.exec_dir, "dup_missing_libraries.tsv")
    for path in (out_lib, out_cohort, out_missing):
        if os.path.exists(path):
            sys.exit(f"refusing to overwrite {path}")

    libs = pd.read_csv(args.libraries, sep="\t", dtype=str, keep_default_na=False)
    rows_dir = os.path.join(args.restricted, "per_library")
    done = {os.path.basename(p)[:-len(".tsv")] for p in glob.glob(os.path.join(rows_dir, "*.tsv"))}
    missing = libs[~libs["run"].isin(done)]
    libs = libs[libs["run"].isin(done)]

    records = []
    for lib in libs.itertuples(index=False):
        row = pd.read_csv(os.path.join(rows_dir, f"{lib.run}.tsv"), sep="\t", dtype={"pg_programs": str},
                          keep_default_na=False).iloc[0].to_dict()
        row["individual_id"] = lib.individual_id
        row.update(read_key_values(lib.star_log, "|", STAR_FIELDS))
        row.update(read_key_values(os.path.join(args.restricted, "markdup_stats", f"{lib.run}.txt"), ":",
                                   MARKDUP_FIELDS))
        records.append(row)
    table = pd.DataFrame(records)
    table["star_mapped_reads"] = table["star_unique_reads"] + table["star_multi_reads"]
    table["layout"] = (table["markdup_paired"] > 0).map({True: "paired", False: "single"})
    table["exon_dup_fraction"] = table["exon_dup_reads"] / table["exon_reads"]
    table["exon_distinct_reads"] = table["exon_reads"] - table["exon_dup_reads"]
    table["exon_mate_absent_fraction"] = table["exon_mate_absent_reads"] / table["exon_reads"]
    table["genebody_dup_fraction"] = table["markdup_duplicate_total"] / table["markdup_examined"]
    # markdup estimates library size from pairs only and prints 0 for single-end input
    # nullable integers so counts are not written in 6-significant-digit float form
    table["estimated_library_size"] = table["estimated_library_size"].where(table["layout"] == "paired").astype("Int64")
    table["genebody_lw_distinct_reads"] = pd.array(
        [lander_waterman_size(n, n - d) for n, d in zip(table["markdup_examined"], table["markdup_duplicate_total"])],
        dtype="Int64")
    columns = ["cohort", "run", "individual_id", "layout", "pg_programs",
               "star_input_reads", "star_unique_reads", "star_multi_reads", "star_mapped_reads",
               "subset_reads", "preexisting_dup_flags_in_subset", "markdup_examined",
               "markdup_duplicate_total", "genebody_dup_fraction", "estimated_library_size",
               "genebody_lw_distinct_reads",
               "exon_reads", "exon_dup_reads", "exon_distinct_reads", "exon_dup_fraction",
               "exon_paired_reads", "exon_mate_absent_reads", "exon_mate_absent_fraction"]
    table = table[columns].sort_values(["cohort", "run"])
    table.to_csv(out_lib, sep="\t", index=False, float_format="%.6g", na_rep="NA")

    expected = pd.read_csv(args.libraries, sep="\t", dtype=str, keep_default_na=False)["cohort"].value_counts()
    per_cohort = table.groupby("cohort").apply(summarise, include_groups=False)
    # pooled rows by layout: the two layouts use different duplicate rules (pair vs one-end)
    for layout, group in table.groupby("layout"):
        per_cohort.loc[f"all_{layout}"] = summarise(group)
    expected_by_layout = table.groupby("layout").size()  # finished libraries stand in for expected here
    per_cohort.insert(0, "n_libraries_expected", [
        expected[c] if c in expected.index else expected_by_layout[c[len("all_"):]] for c in per_cohort.index])
    per_cohort.index.name = "cohort"
    per_cohort.to_csv(out_cohort, sep="\t", float_format="%.10g", na_rep="NA")
    missing[["task", "cohort", "run"]].to_csv(out_missing, sep="\t", index=False)
    print(per_cohort.to_string())
    print(f"finished {len(table)}; missing {len(missing)}")


if __name__ == "__main__":
    main()
