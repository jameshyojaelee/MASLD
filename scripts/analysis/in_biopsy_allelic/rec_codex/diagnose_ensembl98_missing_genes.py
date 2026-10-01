"""Metadata-only disposition of the fixed 1010 missing modeled RNA genes."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
SOURCE = REC / "ensembl98_reference_census_21998377"
GTF = REC / "ensembl98_reference_21998355/Homo_sapiens.GRCh38.98.gtf.gz"
GTF_SHA = "7966e412a6c5864f232ebcb40d67f9ab873064af090dac28a8cab013aaea18b7"
GUARDS = {
    "missing_modeled_union_genes.tsv": "a922ce25bd465d906f27730b02a4a475d3ee56730c8a2c96fb64f72eef52cadd",
    "gtf_gene_coordinates.tsv": "594e4b78de02adc274d91d28954fbceffadba2481ca3b79e02e6da46ecc96e41",
    "gtf_transcript_gene_coordinates.tsv": "4dc4f78bd613bfb279f2e869a447f7c34625f8cccf1f9488c9ac36ca0efb7f71"}


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write(path, records, fields):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Transcript metadata census must run on compute")
    for name, expected in GUARDS.items():
        if sha(SOURCE/name) != expected:
            raise ValueError("Frozen census metadata changed: " + name)
    if sha(GTF) != GTF_SHA:
        raise ValueError("Frozen Ensembl98 GTF changed")
    missing = list(rows(SOURCE/"missing_modeled_union_genes.tsv"))
    if len(missing) != 1010 or len({r["normalized_modeled_gene_id"] for r in missing}) != 1010:
        raise ValueError("Missing modeled gene roster changed or duplicated")
    wanted = {r["normalized_modeled_gene_id"] for r in missing}
    wanted.update(r["normalized_source_gene_id"] for r in missing)
    genes, transcripts = defaultdict(list), defaultdict(list)
    all_biotypes = Counter()
    for r in rows(SOURCE/"gtf_gene_coordinates.tsv"):
        all_biotypes[r["gene_biotype"]] += 1
        if r["gene_id"] in wanted:
            genes[r["gene_id"]].append(r)
    for r in rows(SOURCE/"gtf_transcript_gene_coordinates.tsv"):
        if r["gene_id"] in wanted:
            transcripts[r["gene_id"]].append(r)
    exons = defaultdict(list)
    gtf_headers = []
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                gtf_headers.append(line.rstrip())
                continue
            cols = line.rstrip().split("\t")
            if len(cols) != 9:
                raise ValueError("Malformed GTF annotation")
            if cols[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', cols[8]))
            gene = attrs.get("gene_id", "")
            if gene in wanted:
                exons[gene].append(dict(gene_id=gene, transcript_id=attrs["transcript_id"],
                                       gene_version=attrs.get("gene_version", ""),
                                       transcript_version=attrs.get("transcript_version", ""),
                                       transcript_biotype=attrs.get("transcript_biotype", ""),
                                       seqname=cols[0], start_1based=int(cols[3]), end_inclusive=int(cols[4]),
                                       strand=cols[6], exon_number=attrs.get("exon_number", ""),
                                       exon_id=attrs.get("exon_id", ""), exon_version=attrs.get("exon_version", "")))
    disposition = []
    for r in missing:
        gene = r["normalized_modeled_gene_id"]
        gs, ts = genes[gene], transcripts[gene]
        disposition.append(dict(**r, matrix_stable_normalized_agree=gene == r["normalized_source_gene_id"],
                                gtf_gene_rows=len(gs), gtf_transcript_rows=len(ts),
                                gtf_exon_rows=len(exons[gene]),
                                transcripts_with_exons=len({e["transcript_id"] for e in exons[gene]}),
                                all_transcripts_have_exons=bool(ts) and {t["transcript_id"] for t in ts} <= {e["transcript_id"] for e in exons[gene]},
                                gtf_gene_biotypes=";".join(sorted({g["gene_biotype"] for g in gs})),
                                gtf_gene_geometry_json=json.dumps(gs, sort_keys=True),
                                disposition="gtf_gene_and_transcript_present_fasta_absent" if gs and ts else
                                "gtf_gene_present_no_transcript" if gs else "gtf_gene_absent"))
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output/"missing_1010_gene_disposition.tsv", disposition, list(disposition[0]))
    matched = [r for gene in sorted(transcripts) for r in transcripts[gene]]
    fields = list(next(iter(rows(SOURCE/"gtf_transcript_gene_coordinates.tsv"))))
    write(args.output/"missing_gene_gtf_transcripts.tsv", matched, fields)
    exon_rows = [r for gene in sorted(exons) for r in exons[gene]]
    if exon_rows:
        write(args.output/"missing_gene_gtf_exons.tsv", exon_rows, list(exon_rows[0]))
    summary = dict(source=str(SOURCE), input_sha256=GUARDS, script_sha256=sha(Path(__file__)),
                   launcher_sha256=sha(Path(__file__).with_name("run_diagnose_ensembl98_missing_genes.sbatch")),
                   slurm_job_id=os.environ["SLURM_JOB_ID"], missing_genes=len(disposition),
                   dispositions=dict(Counter(r["disposition"] for r in disposition)),
                   biotypes=dict(Counter(r["gtf_gene_biotypes"] for r in disposition)),
                   mapping_disagreements=sum(not r["matrix_stable_normalized_agree"] for r in disposition),
                   transcript_rows=len(matched), all_gtf_gene_biotype_rows=dict(all_biotypes),
                   gtf_path=str(GTF), gtf_sha256=GTF_SHA, gtf_headers=gtf_headers,
                   exon_rows=len(exon_rows), genes_all_transcripts_have_exons=sum(r["all_transcripts_have_exons"] for r in disposition),
                   missing_exon_transcripts=sorted({t["transcript_id"] for ts in transcripts.values() for t in ts} -
                                                   {e["transcript_id"] for e in exon_rows}),
                   sequencing_or_expression_values_read=False, gene_remapping_or_zero_fill=False,
                   quantification_admitted=False,
                   output_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({k:summary[k] for k in ("missing_genes", "dispositions", "biotypes", "mapping_disagreements", "transcript_rows")}, sort_keys=True))


if __name__ == "__main__":
    main()
