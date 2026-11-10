"""Exact-peak motif context fixtures, with no inferred RNA or causal evidence."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sqlite3

import pytest

SPEC = importlib.util.spec_from_file_location("annotate_motif_context", Path(__file__).resolve().parents[1] / "scripts/annotate_motif_context.py")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def table(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(data)


def fixture(tmp_path):
    root = tmp_path / "motifs"
    root.mkdir()
    database = root / "jaspar.sqlite"
    connection = sqlite3.connect(database)
    connection.executescript("CREATE TABLE MATRIX(ID INTEGER,BASE_ID TEXT,VERSION INTEGER,NAME TEXT);"
                             "CREATE TABLE MATRIX_ANNOTATION(ID INTEGER,TAG TEXT,VAL TEXT);"
                             "INSERT INTO MATRIX VALUES(1,'MA0001',1,'EXAMPLE_A');"
                             "INSERT INTO MATRIX VALUES(2,'MA0002',1,'EXAMPLE_B');"
                             "INSERT INTO MATRIX_ANNOTATION VALUES(1,'family','example family');")
    connection.close()
    table(root / "jaspar_database.tsv", [dict(database_path=str(database), database_sha256=module.digest(database), jaspar_package="fixture")])
    directory = root / "motif_matches/GSE244832/hepatocyte"
    table(directory / "peaks.tsv", [dict(peak_coordinate=p) for p in ["chr1:10-20", "chr1:20-30", "chr1:30-40"]])
    table(directory / "motifs.tsv", [dict(motif_id=m) for m in ["MA0001.1", "MA0002.1"]])
    (directory / "matches.mtx").write_text("%%MatrixMarket matrix coordinate pattern general\n% fixture\n3 2 3\n1 1\n1 2\n2 1\n")
    network = tmp_path / "network"
    network.mkdir()
    (network / "completion.json").write_text('{"status":"completed"}')
    link = dict(outer_fold="0", gene_id="fixture_gene", peak_id="chr1:10-20", method="distance", rank="1", score="-100", score_unit="genomic_distance_negative_bp", evidence_state="candidate_accessibility_link")
    table(network / "GSE244832_hepatocyte/links.tsv", [link, link | dict(peak_id="chr1:30-40"), link | dict(peak_id="chr1:11-21")])
    table(network / "GSE296875_hepatocyte/links.tsv", [link])
    return root, network


def test_exact_matches_and_untestable_native_peaks(tmp_path):
    root, network = fixture(tmp_path)
    out = tmp_path / "result"
    module.annotate([network], root, out)
    with gzip.open(out / "motif_peak_gene.tsv.gz", "rt") as handle:
        result = list(csv.DictReader(handle, delimiter="\t"))
    assert len(result) == 2
    assert {r["motif_id"] for r in result} == {"MA0001.1", "MA0002.1"}
    assert {r["peak_id"] for r in result} == {"chr1:10-20"}
    assert all(r["score"] == "-100" and r["outer_fold"] == "0" for r in result)
    assert all(r["evidence_state"] == "motif_compatible_accessibility_link" and r["tf_activity_inferred"] == "False" and r["causal_regulation_inferred"] == "False" for r in result)
    coverage = module.rows(out / "coverage.tsv")
    by_cohort = {r["cohort"]: r for r in coverage}
    assert by_cohort["GSE244832"]["n_exact_peak_links"] == "2"
    assert by_cohort["GSE244832"]["n_links_absent_from_motif_universe"] == "1"
    assert by_cohort["GSE296875"]["reason"] == "native_peak_motifscan_not_run"
    assert by_cohort["GSE296875"]["status"] == "untestable"
    assert json.loads((out / "completion.json").read_text())["scientific_adoption"] is False
    with pytest.raises(FileExistsError):
        module.annotate([network], root, out)


def test_matrix_dimensions_and_frozen_database_hash(tmp_path):
    root, network = fixture(tmp_path)
    directory = root / "motif_matches/GSE244832/hepatocyte"
    (directory / "matches.mtx").write_text("%%MatrixMarket matrix coordinate pattern general\n4 2 0\n")
    with pytest.raises(ValueError, match="dimensions"):
        module.selected_motif_matches(directory, {"chr1:10-20"})
    with (root / "jaspar.sqlite").open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ValueError, match="SHA256"):
        module.annotate([network], root, tmp_path / "bad_result")
