#!/usr/bin/env python3
"""
Plan the COLOC shared-posterior eligibility rerun (review item 2).

Only gene-study rows that are SuSiE PP.H4 > 0.5 in the adopted master need
rechecking: the GWAS SNPs in 06_susie_coloc.R are a subset of the eQTL SNPs, so
only the eQTL-side share can fail, and dropping a signal pair never changes the
PP.H4 of the pairs that remain. A negative gene therefore cannot become positive.

Selection matches Analysis/.../atac_context_v3/07_prepare_genetic_replay.py
(method == "susie", PP.H4.susie > 0.5) but over ALL 50 studies, not tier 1/2
only. Each (study, chromosome) task goes to exactly one of N batches, balanced by
gene count (longest-processing-time greedy; ties -> lowest batch index), so no
two batches write the same output file.

Writes, under a NEW --run-root (refuses to reuse one):
  plan/rerun_plan.tsv            one row per gene-study pair, with batch_id
  plan/batch_<k>.tsv             the same rows split by batch (sbatch input)
  plan/input_manifest.tsv        sha256 of the master table, runner and configs

Usage:
  python3 09_prepare_coloc_eligibility_batches.py \
      --run-root results/susie_coloc_eligibility_<UTC timestamp>
"""
import argparse, csv, hashlib, os, sys
from collections import defaultdict
from datetime import datetime, timezone

FM = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
COLUMNS = ["batch_id", "gwas_name", "chr", "gene", "ensembl", "adopted_pp_h4_susie",
           "tier", "tier_label", "placement"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", required=True,
                    help="new directory, e.g. results/susie_coloc_eligibility_<ts>")
    ap.add_argument("--master", default=os.path.join(FM, "results/susie_coloc/susie_coloc_all_gwas.csv"))
    ap.add_argument("--eqtl-dir", default=os.path.join(FM, "results/eqtl_susie_polyfun"))
    ap.add_argument("--n-batches", type=int, default=5)
    args = ap.parse_args()

    run_root = args.run_root if os.path.isabs(args.run_root) else os.path.join(FM, args.run_root)
    if os.path.lexists(run_root):
        sys.exit(f"refusing to reuse an existing run root: {run_root}")
    if os.path.realpath(run_root).startswith(os.path.realpath(os.path.join(FM, "results/susie_coloc")) + os.sep):
        sys.exit("run root must not sit inside the adopted results/susie_coloc/")

    tiers = {}
    with open(os.path.join(FM, "config/gwas_trait_tier.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tiers[r["study_name"]] = r

    rows = []
    with open(args.master) as fh:
        for r in csv.DictReader(fh):
            if r["method"] != "susie":
                continue
            try:
                pp4 = float(r["PP.H4.susie"])
            except (TypeError, ValueError):
                continue
            if not pp4 > 0.5:
                continue
            t = tiers.get(r["gwas_name"], {})
            rows.append({"gwas_name": r["gwas_name"], "chr": int(r["chr"]), "gene": r["gene"],
                         "ensembl": r["ensembl"], "adopted_pp_h4_susie": r["PP.H4.susie"],
                         "tier": t.get("tier", ""), "tier_label": t.get("tier_label", ""),
                         "placement": t.get("placement", "")})
    keys = [(r["gwas_name"], r["ensembl"]) for r in rows]
    if not rows or len(keys) != len(set(keys)):
        sys.exit("SuSiE-positive rows are empty or duplicated by (study, gene)")
    blank = [r for r in rows if not r["ensembl"].startswith("ENSG")]
    if blank:
        sys.exit(f"{len(blank)} SuSiE-positive rows without an Ensembl id")
    fits = {os.path.join(args.eqtl_dir, f"chr{r['chr']}", f"{r['ensembl']}_susie.rds") for r in rows}
    missing = sorted(p for p in fits if not os.path.isfile(p))
    if missing:
        sys.exit(f"{len(missing)} eQTL SuSiE fits missing, e.g. {missing[:3]}")

    # (study, chr) tasks, balanced by gene count.
    task_genes = defaultdict(int)
    for r in rows:
        task_genes[(r["gwas_name"], r["chr"])] += 1
    load = [0] * args.n_batches
    n_tasks = [0] * args.n_batches
    task_batch = {}
    for task in sorted(task_genes, key=lambda t: (-task_genes[t], t[0], t[1])):
        b = min(range(args.n_batches), key=lambda k: (load[k], k))
        task_batch[task] = b + 1
        load[b] += task_genes[task]
        n_tasks[b] += 1
    for r in rows:
        r["batch_id"] = task_batch[(r["gwas_name"], r["chr"])]
    rows.sort(key=lambda r: (r["batch_id"], r["gwas_name"], r["chr"], r["ensembl"]))

    t12 = [r for r in rows if r["tier"] in ("1", "2")]
    print(f"SuSiE PP.H4 > 0.5 rows: {len(rows):,} | genes: {len({r['ensembl'] for r in rows}):,} "
          f"| study x chr tasks: {len(task_genes):,}")
    print(f"  of which tier 1/2: {len(t12):,} rows | {len({r['ensembl'] for r in t12}):,} genes")
    for k in range(args.n_batches):
        print(f"  batch {k + 1}: {n_tasks[k]:>3} tasks | {load[k]:>4} gene-study rows")

    plan = os.path.join(run_root, "plan")
    os.makedirs(plan)
    write_tsv(os.path.join(plan, "rerun_plan.tsv"), COLUMNS, rows)
    for k in range(1, args.n_batches + 1):
        write_tsv(os.path.join(plan, f"batch_{k}.tsv"), COLUMNS,
                  [r for r in rows if r["batch_id"] == k])
    inputs = [("adopted_master", args.master),
              ("runner", os.path.join(FM, "src/06_susie_coloc.R")),
              ("eligibility_helper", os.path.join(FM, "src/coloc_signal_eligibility.R")),
              ("gwas_registry", os.path.join(FM, "config/gwas_registry.tsv")),
              ("gwas_trait_tier", os.path.join(FM, "config/gwas_trait_tier.tsv"))]
    write_tsv(os.path.join(plan, "input_manifest.tsv"), ["role", "path", "sha256", "created_utc"],
              [{"role": role, "path": p, "sha256": sha256(p),
                "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
               for role, p in inputs])
    print(f"wrote {plan}")


if __name__ == "__main__":
    main()
