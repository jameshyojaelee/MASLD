"""Independent cell-census, roster, and full-family BH checks for the rerun."""
import collections
import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
TASK = ROOT / "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/singlecell"
OLD = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
assert (TASK / "RERUNS_COMPLETED").is_file(), "All four rerun arms must finish before verification"


def read(path, sep="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


def write(path, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def canonical(text):
    return "".join(c if c.isascii() and c.isalnum() else "_" for c in text)


def check_bh(rows):
    ordered = sorted(range(len(rows)), key=lambda i: float(rows[i]["pvalue"]))
    running = 1.0
    for rank in range(len(ordered), 0, -1):
        row = rows[ordered[rank - 1]]
        p = float(row["pvalue"])
        assert math.isfinite(p) and 0 <= p <= 1
        running = min(running, p * len(ordered) / rank)
        assert abs(float(row["padj"]) - running) < 1e-10


libraries = read(OLD / "metadata/library_metadata_source_F.tsv")
summary = []
omissions = []
changed_calls = []
for arm, old_arm in [("binary", "pseudobulk_de_v2"),
                     ("strict_control", "pseudobulk_de_strict_control"),
                     ("coarse", "pseudobulk_de_coarse")]:
    current = TASK / "results" / arm
    old_dir = OLD / old_arm
    outputs = sorted(current.glob("*_de.csv"))
    assert outputs, arm
    for path in outputs:
        stem = path.name.removesuffix("_de.csv")
        ct = stem.split("_Steatosis_vs_Healthy")[0].split("_Steatohepatitis_vs_Steatosis")[0].split("_Cirrhosis_vs_Steatohepatitis")[0]
        census_key = [key for key in libraries[0] if canonical(key) == canonical("n_" + ct)]
        assert len(census_key) == 1
        with (ROOT / "Analysis/SingleCell/results_gpu_v2/pseudobulk" / (ct + "_pseudobulk.csv")).open() as handle:
            present = set(next(csv.reader(handle))[1:])
        expected = collections.Counter()
        for lib in libraries:
            eligible = (lib["exclude_stage_analysis"].lower() == "false" and lib["dataset"] != "GSE189600") if arm == "coarse" else lib["strict_binary_eligible" if arm == "strict_control" else "binary_eligible"].lower() == "true"
            if eligible and lib["preparation_method"] in ("nuclei", "cd45_negative", "unsorted") and lib["sample"] in present:
                expected[lib["donor"]] += int(lib[census_key[0]])
        roster = read(current / (stem + "_donor_roster.tsv"))
        assert len({row["sample"] for row in roster}) == len(roster)
        for donor in roster:
            n = int(donor["n_cells"])
            assert n == expected[donor["sample"]] and n >= 50, (arm, stem, donor, expected)
            if arm == "strict_control":
                assert donor["sample"] != "GSE202379_P30"
        old_roster_path = old_dir / (stem + "_donor_roster.tsv")
        old_roster = read(old_roster_path) if old_roster_path.exists() else []
        kept = {row["sample"] for row in roster}
        for donor in old_roster:
            if donor["sample"] not in kept:
                omissions.append(dict(arm=arm, analysis=stem, donor=donor["sample"], contributing_lineage_cells=expected[donor["sample"]]))
        rows = read(path, ",")
        check_bh(rows)
        old_path = old_dir / path.name
        previous = read(old_path, ",") if old_path.exists() else []
        calls = {row["gene"] for row in rows if float(row["padj"]) < .05}
        old_calls = {row["gene"] for row in previous if float(row["padj"]) < .05}
        assert len({row["gene"] for row in rows}) == len(rows)
        new_by_gene = {row["gene"]: row for row in rows}
        old_by_gene = {row["gene"]: row for row in previous}
        changes = [abs(float(new_by_gene[gene]["logFC"]) - float(old_by_gene[gene]["logFC"]))
                   for gene in new_by_gene.keys() & old_by_gene.keys()]
        if arm == "binary" and ct == "Hepatocytes":
            for gene in sorted(calls ^ old_calls):
                new, old = new_by_gene.get(gene, {}), old_by_gene.get(gene, {})
                changed_calls.append(dict(gene=gene, change="gained" if gene in calls else "lost",
                                          old_logFC=old.get("logFC", "not_tested"), new_logFC=new.get("logFC", "not_tested"),
                                          old_padj=old.get("padj", "not_tested"), new_padj=new.get("padj", "not_tested")))
        summary.append(dict(arm=arm, analysis=stem, old_donors=len(old_roster) if old_roster_path.exists() else "not_recorded", new_donors=len(roster),
                            old_tested=len(previous), new_tested=len(rows), old_calls=len(old_calls),
                            new_calls=len(calls), retained_calls=len(calls & old_calls),
                            lost_calls=len(old_calls - calls), gained_calls=len(calls - old_calls),
                            min_contributing_cells=min(int(row["n_cells"]) for row in roster),
                            median_abs_logFC_change=statistics.median(changes) if changes else "not_testable"))
    missing = sorted(path.name for path in old_dir.glob("*_de.csv") if not (current / path.name).exists())
    for name in missing:
        old_roster_path = old_dir / name.replace("_de.csv", "_donor_roster.tsv")
        summary.append(dict(arm=arm, analysis=name.removesuffix("_de.csv"), old_donors=len(read(old_roster_path)) if old_roster_path.exists() else "not_recorded",
                            new_donors=0, old_tested=len(read(old_dir / name, ",")), new_tested=0,
                            old_calls=sum(float(row["padj"]) < .05 for row in read(old_dir / name, ",")),
                            new_calls=0, retained_calls=0, lost_calls="not_testable", gained_calls=0,
                            min_contributing_cells="not_testable", median_abs_logFC_change="not_testable"))
write(TASK / "result_changes.tsv", summary)
write(TASK / "removed_donors.tsv", omissions)
if changed_calls:
    write(TASK / "hepatocyte_changed_calls.tsv", changed_calls)
allcell = read(TASK / "results/allcell/allcell_pseudobulk_de.csv", ",")
check_bh(allcell)
allcell_roster = read(TASK / "results/allcell/donors_in_model.tsv")
assert len({row["sample"] for row in allcell_roster}) == len(allcell_roster)
(TASK / "allcell_checks.json").write_text(json.dumps(dict(
    n_donors=len(allcell_roster), n_reference=sum(row["condition"] == "Healthy" for row in allcell_roster),
    n_cases=sum(row["condition"] == "MASLD" for row in allcell_roster), n_tested=len(allcell),
    bh_calls=sum(float(row["padj"]) < .05 for row in allcell),
    bh_and_effect_calls=sum(float(row["padj"]) < .05 and abs(float(row["lfc"])) > .5 for row in allcell),
    independent_bh_passed=True), indent=2) + "\n")
(TASK / "verification.json").write_text(json.dumps(dict(
    passed=True, source="independent CSV census against actual contributing pseudobulk columns",
    analyses_checked=len(summary), bh_method="independent monotone BH over every output's full gene family",
    min_cell_threshold=50, cell_counts_and_donor_uniqueness=True, strict_control_exclusion=True), indent=2) + "\n")
print(json.dumps(summary, indent=2))
