"""Compute-only native protein identities and observation states, never scores."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import sys

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
CONTRACT = ROOT / "docs/technical/agent_exchange/codex_plasma_liver_source_contract_20260930"
PAIRS = CONTRACT / "codex_plasma_liver_source_contract_pairs.tsv"
RECEIPT = CONTRACT / "codex_plasma_liver_source_contract_receipt.json"
MEMBERSHIP = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
GUARDS = {PAIRS: "b1e628aa30166898c7b6eb235b932451d9e64ce973676e11dd57c496a0a73d5f",
          RECEIPT: "19f1dcc2adc859019974d0979cd216273f42acf68a01e3cfacd14ef57f35caf0",
          MEMBERSHIP: "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"}
SOURCES = {"liver": ROOT/"data/PXD051911/liver_protein_quant.txt",
           "plasma": ROOT/"data/PXD051911/plasma_protein_quant.txt"}
EXPECTED_LIVER_SHA1 = "810f33d2bbffc562fb9842df5ba75eff0bde0d15"
MISSING = {"", "NA", "NaN", "nan"}
STATES = ("missing_token", "nonfinite", "nonpositive", "finite_positive")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rows(path):
    with path.open(newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def writer(path, fields):
    handle = path.open("x", newline="")
    output = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    output.writeheader()
    return handle, output


def state(cell):
    text = str(cell).strip()
    if text in MISSING:
        return "missing_token"
    try:
        value = float(text)
    except ValueError:
        raise ValueError("Unexpected nonnumeric molecular token; value suppressed") from None
    if not math.isfinite(value):
        return "nonfinite"
    return "finite_positive" if value > 0 else "nonpositive"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Protein observation-state extraction only under SLURM")
    for path, expected in GUARDS.items():
        require(digest(path) == expected, "Frozen metadata changed: "+str(path))
    source_contract = json.loads(RECEIPT.read_text())
    pairs = list(rows(PAIRS))
    require(len(pairs)==41 and len({r["patient_name"] for r in pairs})==41, "Paired donor axis changed")
    before = {name:digest(path) for name,path in SOURCES.items()}
    require(digest(SOURCES["liver"], "sha1") == EXPECTED_LIVER_SHA1, "Published liver deposit checksum changed")
    programs = defaultdict(list)
    for r in rows(MEMBERSHIP):
        weight = float(r["original_l1_weight"])
        require(math.isfinite(weight) and weight >= 0, "Invalid frozen original weight")
        programs[r["program_uid"]].append((r["mapped_symbol"].strip().upper(),
                                            r["mapped_symbol_status"]=="gencode_v49_unique_symbol_confirmed", weight))
    require(len(programs)==117 and all(abs(sum(w for _,_,w in rs)-1)<1e-8 for rs in programs.values()),
            "Frozen117 membership/weight axis changed")
    args.output.mkdir(parents=True, exist_ok=False)
    import pandas as pd
    summary_assays = {}
    support_handle, support_writer = writer(args.output/"program_observation_support.tsv",
        ["assay","patient_name","unique_identifier","program_uid","original_members","confirmed_mapping_members",
         "source_present_members","finite_positive_members","original_l1_source_present","original_l1_finite_positive",
         "all_original_members_confirmed_source_present","all_original_members_confirmed_finite_positive"])
    try:
        for assay, path in SOURCES.items():
            filename_field = assay+"_proteomics_filename"
            selected = [r[filename_field] for r in pairs]
            require(len(set(selected))==41, "Ambiguous assay specimen ownership")
            with path.open("rb") as handle:
                header_raw = handle.readline(100001)
            require(len(header_raw)<=100000 and hashlib.sha256(header_raw).hexdigest()==
                    source_contract["headers"][path.name]["header_sha256"], "Native header changed")
            header = header_raw.decode().rstrip("\r\n").split("\t")
            require(header[:3]==["ProteinAccessions","Genes","ProteinDescriptions"] and
                    len(header)==len(set(header)) and set(selected)<=set(header), "Native feature/sample schema changed")
            fields = ["source_row_index","protein_accessions","genes","protein_description","gene_identity_state","single_gene_symbol"]
            feature_handle, feature_writer = writer(args.output/(assay+"_native_feature_identity.tsv"), fields)
            mask_handle, mask_writer = writer(args.output/(assay+"_native_group_observation_states.tsv"),
                ["source_row_index","patient_name","unique_identifier","state"])
            accession_counts, gene_groups, detected, state_counts = Counter(), Counter(), [Counter() for _ in pairs], [Counter() for _ in pairs]
            gene_counts_by_person = [Counter() for _ in pairs]
            identity_states, total = Counter(), 0
            try:
                chunks = pd.read_csv(path, sep="\t", usecols=header[:3]+selected, dtype=str,
                                     na_filter=False, chunksize=512)
                for chunk in chunks:
                    for r in chunk.to_dict("records"):
                        accession, raw_gene = r["ProteinAccessions"], r["Genes"]
                        require(accession.strip() not in MISSING, "Missing native protein accession")
                        gene = raw_gene.strip().upper()
                        if raw_gene.strip() in MISSING:
                            identity = "missing_gene_annotation"; gene = ""
                        elif any(x in gene for x in (";", ",", "|")) or not re.fullmatch(r"[A-Z0-9][A-Z0-9_.-]*", gene):
                            identity = "ambiguous_or_multigene_group"; gene = ""
                        else:
                            identity = "single_gene_annotation"
                        accession_counts[accession] += 1; identity_states[identity] += 1
                        if gene:
                            gene_groups[gene] += 1
                        feature_writer.writerow(dict(source_row_index=total, protein_accessions=accession,
                            genes=raw_gene, protein_description=r["ProteinDescriptions"], gene_identity_state=identity,
                            single_gene_symbol=gene))
                        for j,pair in enumerate(pairs):
                            s = state(r[selected[j]])
                            state_counts[j][s] += 1
                            if gene and s=="finite_positive":
                                detected[j][gene] += 1
                            if gene:
                                gene_counts_by_person[j][gene] += 1
                            mask_writer.writerow(dict(source_row_index=total, patient_name=pair["patient_name"],
                                                       unique_identifier=pair["unique_identifier"], state=s))
                        total += 1
                        require(total<=100000, "Protein metadata row cap exceeded")
            finally:
                feature_handle.close(); mask_handle.close()
            require(total>0, "Empty native protein feature axis")
            duplicate_handle, duplicate_writer = writer(args.output/(assay+"_duplicate_accession_annotations.tsv"),
                ["protein_accessions","native_rows"])
            for accession,count in sorted(accession_counts.items()):
                if count>1:
                    duplicate_writer.writerow(dict(protein_accessions=accession,native_rows=count))
            duplicate_handle.close()
            person_handle, person_writer = writer(args.output/(assay+"_person_observation_completeness.tsv"),
                ["patient_name","unique_identifier","native_groups","missing_token","nonfinite","nonpositive","finite_positive",
                 "unambiguous_source_genes","genes_any_group_finite_positive","genes_all_groups_finite_positive"])
            for j,pair in enumerate(pairs):
                person_writer.writerow(dict(patient_name=pair["patient_name"],unique_identifier=pair["unique_identifier"],
                    native_groups=total, **{s:state_counts[j][s] for s in STATES}, unambiguous_source_genes=len(gene_groups),
                    genes_any_group_finite_positive=len(detected[j]),
                    genes_all_groups_finite_positive=sum(detected[j][g]==n for g,n in gene_groups.items())))
                for uid,members in sorted(programs.items()):
                    present = [(g,ok,w) for g,ok,w in members if ok and g in gene_groups]
                    positive = [(g,ok,w) for g,ok,w in members if ok and g in detected[j]]
                    support_writer.writerow(dict(assay=assay,patient_name=pair["patient_name"],unique_identifier=pair["unique_identifier"],
                        program_uid=uid,original_members=len(members),confirmed_mapping_members=sum(ok for _,ok,_ in members),
                        source_present_members=len(present),finite_positive_members=len(positive),
                        original_l1_source_present=sum(w for _,_,w in present), original_l1_finite_positive=sum(w for _,_,w in positive),
                        all_original_members_confirmed_source_present=len(present)==len(members),
                        all_original_members_confirmed_finite_positive=len(positive)==len(members)))
            person_handle.close()
            summary_assays[assay] = dict(native_groups=total, identity_states=dict(identity_states),
                native_unambiguous_genes=len(gene_groups), duplicate_accession_ids=sum(n>1 for n in accession_counts.values()),
                states_all_paired_cells=dict(sum(state_counts,Counter())), selected_sample_columns=41)
    finally:
        support_handle.close()
    after = {name:digest(path) for name,path in SOURCES.items()}
    require(before==after, "Native source changed during extraction")
    require(all(digest(p)==s for p,s in GUARDS.items()), "Frozen metadata changed during extraction")
    report = dict(paired_people=41, frozen_programs=117, native_source_sha256_before=before, native_source_sha256_after=after,
        metadata_sha256={str(p):s for p,s in GUARDS.items()}, assays=summary_assays,
        states=dict(missing_token="source NA/NaN/nan/empty", nonfinite="other numeric nonfinite tokens",
                    nonpositive="finite<=0", finite_positive="finite>0"),
        gene_detection_rule="at least one finite-positive group with a single unambiguous native Genes annotation",
        program_support_rule="original membership and original L1 weights; uniquely confirmed symbol mapping only; no renormalization",
        intensity_values_exported=False, protein_scores_computed=False, imputation=False, target_choice_or_model_fit=False,
        source_deposited_development_only=True, prospective_spectronaut_equivalence=False, independent_validation=False,
        script_sha256=digest(Path(__file__)), launcher_sha256=digest(Path(__file__).with_name("run_codex_plasma_liver_observability.sbatch")),
        environment=dict(python=sys.version,pandas=pd.__version__,platform=platform.platform(),slurm_job_id=os.environ["SLURM_JOB_ID"]),
        output_sha256={p.name:digest(p) for p in args.output.iterdir() if p.is_file()})
    (args.output/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(paired_people=41,frozen_programs=117,assays=summary_assays),sort_keys=True))


if __name__ == "__main__":
    main()
