"""Prepare new source-label refits using existing counts and frozen programs."""
import csv
import shutil
from pathlib import Path
from qualify_target_summary import qualify

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
SRC = ROOT / "Analysis/SingleCell/scripts/diagnosis_correction"
ROSTER = OUT / "metadata/current_source_diagnosis_donor_roster.tsv"


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, rows):
    with path.open("x", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def copy_script(source, dest, replacements=()):
    text = source.read_text()
    for a, b in replacements:
        assert a in text, (source, a)
        text = text.replace(a, b)
    dest.write_text(text)


meta = {r["donor"]: r for r in read(OUT / "metadata/donor_metadata.tsv")}
for r in read(OUT / "metadata/library_metadata.tsv"):
    meta.setdefault(r["sample"], meta[r["donor"]])
scripts = OUT / "refit_scripts"
scripts.mkdir(exist_ok=True)
old_roster = "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/current_source_diagnosis_donor_roster.tsv"
follow = ROOT / "scripts/analysis/bulk_refit_followup"
copy_script(follow / "refit_qualified_program_diagnosis.R", scripts / "ols.R", [("uniqueN(p$donor)==62L", "uniqueN(p$donor)==61L")])
copy_script(follow / "refit_qualified_ambient.py", scripts / "ambient.py", [("max()==60", "max()==59")])
copy_script(follow / "refit_source_qualified_ccc.R", scripts / "ccc.R", [(old_roster, str(ROSTER.relative_to(ROOT)))])
h17 = scripts / "H17"
shutil.copytree(ROOT / "Analysis/Hypotheses/H17-source-stage-corrected-20261014_22237956/scripts", h17, dirs_exist_ok=True)
copy_script(ROOT / "Analysis/Hypotheses/H34-stage-source-corrected-20261014_22237596/h34_y_loss_qualified_stage.py", scripts / "h34.py", [
    (old_roster, str(ROSTER.relative_to(ROOT))),
    ('    md = pd.DataFrame(rows).drop_duplicates("donor")',
     '    md = pd.DataFrame(rows).drop_duplicates("donor")\n'
     '    numeric_f = pd.read_csv(ROSTER.parent / "documented_F_source.tsv", sep="\\t").set_index("donor")\n'
     '    md["F_saf"] = md.donor.map(numeric_f.F_stage_documented)'),
])
h34 = OUT / "H34"
h34.mkdir(exist_ok=True)
if not (h34 / "geo_sample_characteristics").exists():
    (h34 / "geo_sample_characteristics").symlink_to(ROOT / "Analysis/Hypotheses/H34-20261008T015507Z/geo_sample_characteristics", target_is_directory=True)
n3 = OUT / "N3"
n3.mkdir(exist_ok=True)
old = ROOT / "Analysis/Hypotheses/N3-source-diagnosis-20261014"
for kind in ("raw", "corrected"):
    dest = n3 / kind
    dest.mkdir(exist_ok=True)
    for p in (old / kind).iterdir():
        if p.name not in ("donor_meta.tsv", "donor_meta_original.tsv") and (p.is_file() or p.name == "per_library"):
            if not (dest / p.name).exists():
                (dest / p.name).symlink_to(p.resolve(), target_is_directory=p.is_dir())
    rows = read(old / kind / "donor_meta.tsv")
    for r in rows:
        m = meta[r["donor"]]
        r["disease_stage_coarse"] = m["disease_stage_coarse"]
        r["stage_ord"] = m["stage_ordinal"]
        r["exclude"] = r["exclude_stage_analysis"] = m["exclude_stage_analysis"]
        if m["exclude_stage_analysis"].lower() == "true":
            r["included"] = "False"
            r["inclusion"] = "source_diagnosis_excluded"
        if r["donor"] == "GSE202379_P98":
            r["F_documented"] = "1"
        if r["dataset"] == "GSE202379" and "F_stage_documented" in m:
            r["F_documented"] = m["F_stage_documented"]
    write(dest / "donor_meta.tsv", rows)
    assert sum(r["included"].lower() == "true" for r in rows) == 69
copy_script(old / "scripts/03_donor_model.R", scripts / "n3_model.R")
copy_script(old / "scripts/04_ambient_summary.R", scripts / "n3_summary.R", [("nrow(inc) == 70L", "nrow(inc) == 69L")])
summary_path = scripts / "n3_summary.R"
summary_path.write_text(qualify(summary_path.read_text()))
bf3 = OUT / "BF3"
(bf3 / "build").mkdir(parents=True)
old_bf = ROOT / "Analysis/Hypotheses/BF3-20261023T154841Z"
(bf3 / "build/donor_type_counts.npz").symlink_to(old_bf / "build-20261023T155101Z/donor_type_counts.npz")
rows = read(old_bf / "build-20261023T155101Z/donors.tsv")
for r in rows:
    m = meta[r["donor"]]
    r["diagnosis"] = "Unresolved" if m["exclude_stage_analysis"].lower() == "true" else m["disease_stage_coarse"]
write(bf3 / "build/donors.tsv", rows)
shutil.copytree(old_bf / "scripts", scripts / "BF3")
(OUT / "REFITS_PREPARED").touch()
print("Prepared donor-level refits; N3 has 69 eligible donors.")
