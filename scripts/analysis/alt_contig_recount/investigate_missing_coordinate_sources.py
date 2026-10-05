#!/usr/bin/env python3
"""Collect bounded official support for eighteen missing-coordinate annotations.

Downloads only listings, assembly tables and GRC issue pages. Never assigns a
gene disposition from missing placement, symbols, or sequence similarity.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import urllib.request
import urllib.error


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--followup", action="store_true")
    p.add_argument("--correspondence", action="store_true")
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    gca = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/"
    gcf = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/"
    requests = [
        ("gca_structure", gca+"GCA_000001405.29_GRCh38.p14_assembly_structure/"),
        ("gcf_structure", gcf+"GCF_000001405.40_GRCh38.p14_assembly_structure/"),
        ("gca_alt_alignments", gca+"GCA_000001405.29_GRCh38.p14_assembly_structure/ALT_REF_LOCI_1/alt_scaffolds/alignments/"),
        ("gcf_alt_alignments", gcf+"GCF_000001405.40_GRCh38.p14_assembly_structure/ALT_REF_LOCI_1/alt_scaffolds/alignments/"),
        ("gcf_patch_alignments", gcf+"GCF_000001405.40_GRCh38.p14_assembly_structure/PATCHES/alt_scaffolds/alignments/"),
        ("gencode49_files", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/"),
        ("gencode49_metadata", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/metadata/"),
        ("gencode49_docs", "https://www.gencodegenes.org/human/release_49.html"),
        ("ensembl115_core", "https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/"),
    ]
    requests += [(f"grc_HG-{n}", f"https://www.ncbi.nlm.nih.gov/grc/human/issues/HG-{n}") for n in range(2509,2513)]
    if args.followup:
        metadata = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/"
        core = "https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/"
        requests = [
            ("KI270825.1_CM000672.2.gff", gca+"GCA_000001405.29_GRCh38.p14_assembly_structure/ALT_REF_LOCI_1/alt_scaffolds/alignments/KI270825.1_CM000672.2.gff"),
            ("grc_exclusions.bed", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed"),
        ]
        requests += [("gencode49_"+name+".gz", metadata+"gencode.v49.metadata."+name+".gz")
                     for name in ["HGNC", "RefSeq", "EntrezGene", "Gene_source", "Annotation_remark"]]
        requests += [("ensembl115_"+name+".gz",core+name+".gz") for name in ["assembly_exception.txt", "meta.txt", "seq_region.txt"]]
    if args.correspondence:
        structure = gca+"GCA_000001405.29_GRCh38.p14_assembly_structure/"
        core = "https://ftp.ensembl.org/pub/release-115/mysql/homo_sapiens_core_115_38/"
        requests = [(f"chr{n}.agp.gz",structure+f"Primary_Assembly/assembled_chromosomes/AGP/chr{n}.agp.gz")
                    for n in [13,14,15,22]]
        requests += [("alt.scaf.agp.gz",structure+"PATCHES/alt_scaffolds/AGP/alt.scaf.agp.gz")]
        requests += [(name,core+name) for name in ["alt_allele.txt.gz", "alt_allele_attrib.txt.gz",
                     "alt_allele_group.txt.gz", "gene.txt.gz", "homo_sapiens_core_115_38.sql.gz"]]
    manifest = []
    for name, url in requests:
        row = {"name": name, "url": url}
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                data = response.read(6*1024*1024+1)
            if len(data)>6*1024*1024: raise ValueError("Source exceeds bounded six-MB limit")
            path = args.output/(name if name.endswith((".gz", ".gff", ".bed")) else name+".html")
            path.write_bytes(data)
            row.update(status="downloaded",path=str(path.resolve()),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
        except (OSError, ValueError) as exc:
            row.update(status="unavailable",reason=str(exc))
        manifest.append(row)
        print(json.dumps(row), flush=True)
    (args.output/"source_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


if __name__ == "__main__":
    main()
