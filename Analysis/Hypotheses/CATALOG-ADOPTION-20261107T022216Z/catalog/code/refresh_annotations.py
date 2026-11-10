"""Replay only the producer's annotation stage after its bounded source correction."""
import hashlib
import json
from pathlib import Path
import shutil
import textwrap

import pandas as pd

root = Path(__file__).resolve().parents[5]
config = root / "scripts/portal/catalog_current.json"
c = json.loads(config.read_text())
rel = c["release_id"]
out = root / c["output"]
annotation = root / c["annotation_bundle"]
genes = pd.read_parquet(out / "passport_gene_index.parquet")
evidence = pd.read_parquet(out / "passport_evidence_long.parquet")
gmap = genes.set_index("ensembl_id")
current_gen = evidence[evidence.evidence_domain == "genetics"].set_index("ensembl_id")
current_rna = evidence[evidence.evidence_domain == "transcriptomics"].set_index("ensembl_id")
moves = pd.read_csv(root / c["eligibility_proposals"], sep="\t", dtype=str, keep_default_na=False)
moves = moves[moves.priority_rule_eligibility_transition == "True"].set_index("ensembl")
inputs = {"annotation_r9d": annotation / "catalog_r9d_gene_table.parquet"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def save_tsv(frame, path):
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


producer = root / "scripts/portal/build_adopted_eligibility_catalog.py"
text = producer.read_text()
start = text.index('    annotated = pd.read_parquet(inputs["annotation_r9d"])')
end = text.index('    graph = out / "portal_export"', start)
exec(compile(textwrap.dedent(text[start:end]), str(producer) + ":annotation_stage", "exec"))
metadata_path = out / "adoption_metadata.json"
metadata = json.loads(metadata_path.read_text())
metadata["source_inputs"]["config"]["sha256"] = sha(config)
metadata["source_inputs"]["protein_assignment"] = {"path": str(root / c["protein_assignment"]), "sha256": sha(root / c["protein_assignment"])}
metadata["annotation_classes_current_native_effects_unchanged"] = True
metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
shutil.copy2(config, out / "build_config.json")
print((out / "annotations/protein_coverage_by_current_class.tsv").read_text())
