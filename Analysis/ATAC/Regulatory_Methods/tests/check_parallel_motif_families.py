"""Serial/default-gzip versus parallel/level1 exact scientific equivalence."""
import argparse
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    task = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(task / "scripts"))
    import annotate_motif_families as candidate
    import merge_motif_families as parallel
    reference_path = task / "results/20261023/full_comparison_r3/final_candidate/source_snapshot/annotate_motif_families.py"
    spec = importlib.util.spec_from_file_location("motif_serial_reference", reference_path)
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    args.output.mkdir(parents=True, exist_ok=False)
    motifs = args.output / "motifs"
    motifs.mkdir()
    database = args.output / "jaspar.sqlite"
    with sqlite3.connect(database) as connection:
        connection.executescript("CREATE TABLE MATRIX(ID INTEGER,BASE_ID TEXT,VERSION INTEGER,NAME TEXT);CREATE TABLE MATRIX_ANNOTATION(ID INTEGER,TAG TEXT,VAL TEXT);")
        connection.executemany("INSERT INTO MATRIX VALUES(?,?,?,?)", [(1,"MA1",1,"TF1"),(2,"MA2",1,"TF2"),(3,"MA3",1,"TF3")])
        connection.executemany("INSERT INTO MATRIX_ANNOTATION VALUES(?,?,?)", [(1,"family","A"),(2,"family","A;B"),(3,"family","B")])
    def table(path, fields, records):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    table(motifs / "jaspar_database.tsv", ["database_path","database_sha256","jaspar_package"], [dict(database_path=str(database.resolve()), database_sha256=candidate.context.digest(database),jaspar_package="synthetic")])
    networks = []
    for cohort, lineage in [("GSE281367","stellate"),("GSE244832","hepatocyte"),("GSE296875","hepatocyte")]:
        network = args.output / f"source_{cohort}_{lineage}"
        network.mkdir()
        (network / "completion.json").write_text(json.dumps(dict(status="completed",pilot=False,interpretation="synthetic_fixture")))
        records = [dict(outer_fold=str(fold),gene_id=gene,peak_id=peak,method=method,rank=str(rank),score=str(score),score_unit="continuous_native",evidence_state="candidate",endpoint=endpoint)
                   for fold in [0,1] for gene in ["G1","G2"] for method in ["distance","cicero"]
                   for rank,(peak,score,endpoint) in enumerate([("p1",.2,"rna"),("p1",.7,"promoter"),("p2",.1,"promoter"),("p3",.9,"rna"),("alias",.8,"rna")],1)]
        table(network / f"{cohort}_{lineage}/links.tsv", list(records[0]), records)
        directory = motifs / "motif_matches" / cohort / lineage
        table(directory / "peaks.tsv", ["peak_coordinate"], [dict(peak_coordinate=p) for p in ["p1","p2","p3"]])
        table(directory / "motifs.tsv", ["motif_id"], [dict(motif_id=p) for p in ["MA1.1","MA2.1","MA3.1"]])
        (directory / "matches.mtx").write_text("%%MatrixMarket matrix coordinate pattern general\n3 3 4\n1 1\n1 2\n2 2\n2 3\n")
        networks.append(network)
    serial = args.output / "serial"
    reference.annotate(networks, motifs, serial)
    actual = args.output / "parallel"
    parallel.parallel(networks, motifs, actual, workers=2, compression_level=1)
    merged = actual / "merged"
    with gzip.open(serial / "motif_family_gene.tsv.gz", "rb") as a, gzip.open(merged / "motif_family_gene.tsv.gz", "rb") as b:
        assert a.read() == b.read()
    assert (serial / "coverage.tsv").read_bytes() == (merged / "coverage.tsv").read_bytes()
    assert json.loads((serial / "completion.json").read_text()) == json.loads((merged / "completion.json").read_text())
    execution = json.loads((actual / "execution.json").read_text())
    assert len({entry["worker_pid"] for entry in execution["contexts"]}) == 2
    parts = sorted((actual / "parts").iterdir())
    rejected = []
    for label, inputs in [("missing_context", parts[:-1]),("duplicate_context", parts+[parts[0]])]:
        try:
            parallel.merge(networks, motifs, inputs, args.output / label, compression_level=1)
        except ValueError:
            rejected.append(label)
        else:
            raise AssertionError("Invalid context set accepted")
    (parts[0] / "INVALIDATED_fixture.json").write_text("{}")
    try:
        parallel.merge(networks, motifs, parts, args.output / "invalidated", compression_level=1)
    except ValueError:
        rejected.append("invalidated_part")
    else:
        raise AssertionError("Invalidated part accepted")
    marker = parts[1] / "completion.json"
    saved = marker.read_text()
    data = json.loads(saved)
    data["status"] = "partial"
    marker.write_text(json.dumps(data))
    # Remove only this fixture's invalidation so incomplete status is tested independently.
    (parts[0] / "INVALIDATED_fixture.json").unlink()
    try:
        parallel.merge(networks, motifs, parts, args.output / "incomplete", compression_level=1)
    except ValueError:
        rejected.append("incomplete_part")
    else:
        raise AssertionError("Incomplete part accepted")
    marker.write_text(saved)
    (args.output / "verification.json").write_text(json.dumps(dict(status="passed", workers=2,
        actual_worker_pids=sorted({entry["worker_pid"] for entry in execution["contexts"]}),
        decompressed_table_bytes_identical=True, coverage_bytes_identical=True,
        scientific_completion_identical=True, reference=str(reference_path), rejected=rejected,
        source_sha256={p.name:candidate.context.digest(p) for p in [Path(candidate.__file__),Path(parallel.__file__)]}),indent=2)+"\n")


if __name__ == "__main__":
    main()
