#!/usr/bin/env python3
"""A2. Rights, exposure and defect register for the AlphaGenome-MASLD program.

Executes IMPLEMENTATION_SPEC.md section 8 (the register deliverable of package G) and
ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md section 15 (model and Atlas usage terms).

Read-only assembly. One row per source or model artifact with training rights, evaluation
rights, redistribution class, protected-or-exposed state, prior exposure evidence, and the
exact file and line the verdict is read from. Sources with no registry entry are included and
their training rights marked unverified; a draft toml is written for each into drafts/ but is
never installed into config/datasets/.

Usage:  python3 a2_rights_exposure.py <output_dir>
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
BENCH = PROJ / "Analysis/MASLD_Model_Benchmark"
DATASETS = BENCH / "config/datasets"
MODEL_ROSTERS = BENCH / "config/models"
MODEL_ARTIFACTS = BENCH / "config/artifacts/models"
AUTHORITY_ARTIFACTS = BENCH / "config/artifacts/model-authorities"
RESOURCES = BENCH / "config/resources.toml"
DECISION = PROJ / "docs/decisions/2026-09-07-derivative-weight-release.md"
ATLAS_RUN = PROJ / "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z"
ATLAS_AUDIT = ATLAS_RUN / "provenance/exposure_audit.json"
PLAN = PROJ / "docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md"

NA = "not_recorded"

# ---------------------------------------------------------------------------
# line-addressed readers: every verdict cites file:line, never a prose summary
# ---------------------------------------------------------------------------


def toml_scalars_with_lines(path: Path) -> dict:
    """Top-level `key = value` scalars from a small toml, each with its 1-based line number.

    The registry tomls are flat enough that a real toml parser is not needed, and a line
    number is required for the register, which a parser does not give.
    """
    out = {}
    section = ""
    for i, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        line = raw.strip()
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip() + "."
            continue
        m = re.match(r'^([A-Za-z0-9_]+)\s*=\s*(.*)$', line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        v = val
        if v.startswith('"') and v.endswith('"') and len(v) >= 2:
            v = v[1:-1]
        out[section + key] = {"value": v, "line": i, "raw": val}
        if not section:
            out[key] = out[section + key]
    return out


def json_key_line(path: Path, key: str) -> int:
    """1-based line where `"key"` first appears in a JSON file, or 0."""
    needle = f'"{key}"'
    for i, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if needle in line:
            return i
    return 0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def flat(x) -> str:
    if x is None:
        return NA
    if isinstance(x, (dict, list)):
        return json.dumps(x, separators=(",", ":"))[:900]
    return str(x).replace("\t", " ").replace("\n", " ")[:900]


# ---------------------------------------------------------------------------
# verdict rules, stated once so the table is reproducible from the inputs
# ---------------------------------------------------------------------------

TRAINING_RULE = (
    "dataset: admission_blocking true and a non-empty blockers list means training is blocked "
    "outside an explicitly activated TaskSpec named in a blocker; split.outer_role train or "
    "development names the lane the activation would sit in. Verdict cites admission_blocking, "
    "split.outer_role and blockers[] by line."
)


def dataset_row(path: Path) -> dict:
    t = toml_scalars_with_lines(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    blocker_line = next((i for i, l in enumerate(lines, 1) if l.startswith("blockers")), 0)
    blockers_raw = lines[blocker_line - 1] if blocker_line else ""

    def g(key, default=NA):
        return t.get(key, {}).get("value", default)

    def gl(key):
        return t.get(key, {}).get("line", 0)

    outer = g("split.outer_role")
    label_vis = g("split.label_visibility")
    exposure = g("exposure_status")
    activated = bool(re.search(r"activat|TaskSpec frozen|is active", blockers_raw))
    lifted = "LIFTED" in blockers_raw or "AMENDMENT" in blockers_raw

    if outer in ("withheld_sealed",):
        training = "prohibited_sealed_source"
    elif outer in ("blocked", "deferred"):
        training = "blocked_source_not_admitted"
    elif activated:
        training = "allowed_only_inside_the_activated_TaskSpec_named_in_blockers"
    else:
        training = "blocked_pending_DatasetActivationContract"

    if label_vis in ("withheld_sealed", "prediction_first_outcomes_withheld"):
        evaluation = "protected_evaluation_prediction_commit_required"
    elif outer in ("blocked", "deferred"):
        evaluation = "blocked_source_not_admitted"
    elif exposure == "target_label_unexposed":
        evaluation = "evaluation_allowed_outcomes_not_yet_read_locally"
    else:
        evaluation = "development_evaluation_only_outcomes_already_exposed"

    if label_vis in ("withheld_sealed", "prediction_first_outcomes_withheld"):
        protected = "protected"
    elif label_vis.startswith("unavailable"):
        protected = "unavailable_" + label_vis.split("_", 1)[1]
    elif exposure == "target_label_unexposed":
        protected = "target_label_unexposed_locally"
    else:
        protected = "exposed"

    rel = str(path.relative_to(PROJ))
    cite = (f"{rel}:{gl('admission_blocking')} admission_blocking; "
            f"{rel}:{gl('split.outer_role')} split.outer_role; "
            f"{rel}:{blocker_line} blockers[]")
    return {
        "artifact_id": g("dataset_id", path.stem),
        "artifact_kind": "registered_dataset",
        "registry_entry": rel,
        "title": g("title"),
        "access_tier": g("access_tier"),
        "training_allowed": training,
        "evaluation_allowed": evaluation,
        "redistribution_class": g("redistribution_class"),
        "license_status": g("license_status"),
        "protected_or_exposed": protected,
        "prior_exposure_evidence": g("provenance.local_exposure"),
        "exposure_status_field": exposure,
        "split_outer_role": outer,
        "label_visibility": label_vis,
        "n_blockers_recorded": str(blockers_raw.count('", "') + 1 if blockers_raw else 0),
        "blocker_amended_or_lifted": str(lifted),
        "verdict_source": cite,
        "verdict_line_redistribution": f"{rel}:{gl('redistribution_class')}",
        "verdict_line_exposure": f"{rel}:{gl('provenance.local_exposure') or gl('exposure_status')}",
    }


def walk_json(o, pre=""):
    """Every leaf of a JSON object as (slash path, scalar value)."""
    if isinstance(o, dict):
        for k, v in o.items():
            yield from walk_json(v, pre + "/" + str(k))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from walk_json(v, pre + "/" + str(i))
    else:
        yield pre, o


LICENCE_LEAF = ("weights", "code", "algorithm", "project_code", "project_weights",
                "derivative_weights", "local_derivative_weights", "upstream_checkpoint",
                "terms_result", "terms_state", "legal_status", "release", "redistribution",
                "outputs", "training_data_and_outputs", "derivative_and_output_trained_models",
                "derivatives", "checkpoint", "candidate_data", "selected_Python_source")


def model_row(d: Path, kind: str) -> dict:
    ea = d / "exposure_audit.json"
    cw = d / "development_crosswalk.json"
    ck = d / "checkpoints.json"
    rel_ea = str(ea.relative_to(PROJ)) if ea.exists() else NA
    rel_cw = str(cw.relative_to(PROJ)) if cw.exists() else NA
    rel_ck = str(ck.relative_to(PROJ)) if ck.exists() else NA
    a = json.loads(ea.read_text()) if ea.exists() else {}
    c = json.loads(cw.read_text()) if cw.exists() else {}
    k = json.loads(ck.read_text()) if ck.exists() else {}

    # Licence evidence. checkpoints.json is the authority: 51 of 104 exposure audits carry no
    # licence field at all, while every model records one in checkpoints.json. Both are read,
    # and the cited file is whichever actually carries the verdict.
    # model-authority artifacts keep their licences in standalone files instead
    side = {}
    side_rel = NA
    for name in ("weights_license.json", "derivative_weights_license.json", "code_license.json"):
        f = d / name
        if f.exists():
            side[name] = json.loads(f.read_text())
            if side_rel == NA:
                side_rel = str(f.relative_to(PROJ))
    lic_leaves, lic_key_path, lic_file = {}, "", NA
    for src, rel in ((k, rel_ck), (side, side_rel), (a, rel_ea), (c, rel_cw)):
        for path, v in walk_json(src):
            if any(w in path.lower() for w in ("licen", "terms", "redistrib")):
                lic_leaves[f"{rel}{path}"] = v
                if not lic_key_path:
                    lic_key_path, lic_file = path, rel
        if lic_leaves:
            break
    champ = a.get("sealed_champion_eligibility", NA)

    # exposure_state lives under checkpoint_findings (62 audits) or checkpoint_finding (37)
    states, reasons = set(), []
    for path, v in walk_json(a):
        if path.endswith("/exposure_state") and v:
            states.add(str(v))
        if path.endswith("/reason") and v:
            reasons.append(str(v))
    states = sorted(states)
    blockers = a.get("remaining_blockers", []) or []

    cross = c.get("training_and_exposure_crosswalk", {}) or {}
    project_trained = cross.get("project_trained_state", NA)
    upstream = cross.get("upstream_training_corpus", NA)
    gate = (c.get("activation_gate", {}) or {}).get("current_state", NA)

    # the weight licence decides training and redistribution; the code licence alone does not
    wkeys = [p for p in lic_leaves if p.rsplit("/", 1)[-1] in
             ("weight_license", "weights", "checkpoint_license", "weight_use_status",
              "declared_license")]
    wl_specific = " ".join(str(lic_leaves[p]).lower() for p in wkeys)
    wl_all = " ".join(str(v).lower() for v in lic_leaves.values())
    wl = wl_specific or wl_all
    if not lic_leaves:
        training = "UNVERIFIED_no_licence_or_terms_field_in_checkpoints_audit_or_crosswalk"
        redis = "UNVERIFIED_no_licence_or_terms_field_exists"
    elif "gated" in wl or "noncommercial" in wl or "non-commercial" in wl:
        training = "local_inference_and_adaptation_only_under_noncommercial_model_terms"
        redis = "no_weight_redistribution_noncommercial_terms"
    elif "undeclared" in wl or "unresolved" in wl or "blocked_pending" in wl or "unknown" in wl:
        training = "UNVERIFIED_weight_terms_recorded_as_undeclared_or_blocked_pending"
        redis = "unresolved_weight_terms_undeclared"
    elif "project_owned" in wl or "not_applicable_classical" in wl:
        training = "project_owned_or_classical_no_upstream_weight_licence"
        redis = "derivative_redistribution_allowed_by_the_licence_authority_file"
    elif "not_applicable_project_trained" in wl:
        training = "project_trained_no_upstream_weight_licence_source_data_policy_governs"
        redis = "derivative_weights_releasable_owner_decision_2026_09_07_data_not_rehosted"
    elif "gpl" in wl:
        training = "permitted_under_copyleft_recorded_licence"
        redis = "copyleft_redistribution_obligations"
    elif "apache" in wl or "mit" in wl or "bsd" in wl or "cc0" in wl:
        training = "permitted_by_the_recorded_weight_licence"
        redis = "derivative_weight_redistribution_permitted_under_the_recorded_licence"
    else:
        training = "recorded_terms_present_but_not_machine_classifiable_read_the_cited_line"
        redis = "see_license_status_column"

    if gate == "deferred" or project_trained == "not_created":
        evaluation = "not_activated_no_evaluation_lane_registered"
    elif "target_label_unexposed" in states:
        evaluation = "eligible_for_protected_evaluation_target_label_unexposed"
    elif states:
        evaluation = "evaluation_restricted_" + ";".join(states)
    else:
        evaluation = NA

    protected = ";".join(states) if states else "no_exposure_state_recorded"
    lic = lic_leaves if lic_leaves else {}
    return {
        "artifact_id": d.name if kind == "model_artifact" else f"authority:{d.name}",
        "artifact_kind": kind,
        "registry_entry": rel_ea,
        "title": flat(a.get("audit_scope", {}).get("artifact_type")
                      or a.get("audit_scope", {}).get("checkpoint")
                      or (c.get("crosswalk_scope", {}) or {}).get("candidate")),
        "access_tier": NA,
        "training_allowed": training,
        "evaluation_allowed": evaluation,
        "redistribution_class": redis,
        "license_status": flat(lic) if lic else NA,
        "protected_or_exposed": protected,
        "prior_exposure_evidence": flat(reasons[0] if reasons else cross.get("exposure_rule")),
        "exposure_status_field": flat(champ),
        "split_outer_role": flat(project_trained),
        "label_visibility": flat(upstream),
        "n_blockers_recorded": str(len(blockers)),
        "blocker_amended_or_lifted": "False",
        "verdict_source": ((f"{lic_file}:{json_key_line(Path(PROJ / lic_file), lic_key_path.rsplit('/', 1)[-1])} {lic_key_path}; "
                            if lic_key_path else
                            "NO licence/terms key in checkpoints.json, exposure_audit.json or development_crosswalk.json; ")
                           + f"{rel_ea}:{json_key_line(ea, 'sealed_champion_eligibility') if ea.exists() else 0}"
                             " sealed_champion_eligibility"),
        "verdict_line_redistribution": (
            f"{lic_file}:{json_key_line(Path(PROJ / lic_file), (wkeys[0].rsplit('/', 1)[-1] if wkeys else lic_key_path.rsplit('/', 1)[-1]))}"
            if lic_key_path else "absent"),
        "verdict_line_exposure": (f"{rel_cw}:{json_key_line(cw, 'training_and_exposure_crosswalk') if cw.exists() else 0}"),
    }


# ---------------------------------------------------------------------------
# sources with NO registry entry at all: training rights unverified by construction
# ---------------------------------------------------------------------------

UNREGISTERED = [
    dict(
        artifact_id="currin_2025_caqtl",
        title="Currin et al. 2025 human liver caQTL (nominal all-pairs plus significant leads and LD proxies)",
        paths=[
            "data/external/currin_2025_caqtl/",
            "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/"
            "liver_significant_caQTL_leadVariants_1kb_analysis_with_populationAlleleFrequencies.bed.gz",
            "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/"
            "caQTL_variants_overlappingPeaks_LD-r2-0.8_withLead.bed.gz",
        ],
        program_role="C1 endpoint 1 label (beta_alt), E1 bridge endogenous side, E2 transfer target",
        note=("beta refers to the ALT allele in the nominal files; in the lead file beta is NOT "
              "re-oriented to the EA column and EA equals REF in 18,981 of 35,361 rows."),
    ),
    dict(
        artifact_id="broadaway_2024_liver_eqtl_meta",
        title="Broadaway et al. liver eQTL meta-analysis, used as COLOC eQTL side and as a benchmark truth set",
        paths=["GWAS/finemapping/results/seqfunc/broadaway_benchmark_truth.tsv"],
        program_role="E1 bridge endogenous side; the eQTL side of the adopted COLOC rerun",
        note=("named as a development source in the Atlas run exposure audit; pools GTEx v8 liver "
              "208 of 1,183, and GTEx is an AlphaGenome training source, so the eQTL side is not "
              "unseen validation for AlphaGenome."),
    ),
    dict(
        artifact_id="gtex_v8_liver_sqtl",
        title="GTEx v8 liver sQTL significant variant-phenotype pairs",
        paths=["data/external/gtex_v8_liver_sqtl/GTEx_Analysis_v8_sQTL/Liver.v8.sqtl_signifpairs.txt.gz"],
        program_role="E1 bridge endogenous side; P3a splice-channel benchmark",
        note=("GTEx v8 is a declared AlphaGenome and Borzoi training corpus, so any AlphaGenome or "
              "Borzoi score on this label is in-corpus, not held out."),
    ),
    dict(
        artifact_id="leafcutter_joint_v2_bulk_cohorts",
        title="leafcutter_joint_v2 splicing quantification over the Resource bulk cohorts",
        paths=["GWAS/finemapping/results/"],
        program_role="named as a development source in the Atlas run exposure audit",
        note=("no dataset toml exists; the five underlying bulk cohorts are registered separately "
              "under bulk_five_cohort.toml but the derived splicing product is not."),
    ),
]

DRAFT_TOML = """# DRAFT ONLY - NOT INSTALLED. Written by scripts/analysis/alphagenome_program/a2_rights_exposure.py
# This file is a proposal for config/datasets/{aid}.toml. It is deliberately left in the program
# output directory. Installing it is a user decision, and every field below marked UNVERIFIED
# requires a primary-source check first.
schema_version = "masld-bench-dataset-v1"
dataset_id = "{aid}"
title = "{title}"
role = "external_development"
status = "unregistered_in_use"
species = "human"
accession = ["UNVERIFIED"]
access_tier = "UNVERIFIED"
automatic_download = false
outcome_url = "NOT_REGISTERED"
redistribution_class = "UNVERIFIED_no_registry_entry_exists"
license_status = "UNVERIFIED_no_registry_entry_exists"
biological_unit = "UNVERIFIED"
expected_biological_units = 0
native_genome_build = "GRCh38"
native_annotation_release = "UNVERIFIED"
modalities = ["UNVERIFIED"]
pairing_levels = ["UNVERIFIED"]
modality_status_default = "observed"
exposure_status = "downstream_demo"
admission_blocking = true
blockers = ["UNVERIFIED TRAINING RIGHTS. This source is consumed by analyses on disk but has no \
registry entry, so no training, evaluation or redistribution right has ever been recorded for it. \
Resolve the source licence and any DUA before it is used as a training label."]
notes = "{note}"

[provenance]
source_url = "UNVERIFIED"
verified_date = "UNVERIFIED_never_verified_against_a_primary_source"
drafted_date = "{stamp}"
local_exposure = "Already read locally by analyses on disk: {paths}"
program_role = "{role}"

[split]
group_key = "UNVERIFIED"
outer_role = "external_development"
label_visibility = "development_visible"
"""


def main() -> None:
    out = Path(sys.argv[1])
    tables = out / "tables"
    drafts = out / "drafts"
    prov = out / "provenance"
    for d in (tables, drafts, prov):
        d.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in sorted(DATASETS.glob("*.toml")):
        rows.append(dataset_row(p))
    n_datasets = len(rows)

    for d in sorted(x for x in MODEL_ARTIFACTS.iterdir() if x.is_dir()):
        rows.append(model_row(d, "model_artifact"))
    n_models = len(rows) - n_datasets

    for d in sorted(x for x in AUTHORITY_ARTIFACTS.iterdir() if x.is_dir()):
        rows.append(model_row(d, "model_authority"))
    n_auth = len(rows) - n_datasets - n_models

    # hosted service and released weights: plan section 15 plus the Atlas run audit
    atlas = json.loads(ATLAS_AUDIT.read_text())
    rel_atlas = str(ATLAS_AUDIT.relative_to(PROJ))
    lic = atlas["license_disposition"]
    rows.append({
        "artifact_id": "alphagenome_hosted_atlas_api",
        "artifact_kind": "hosted_service_output",
        "registry_entry": rel_atlas,
        "title": atlas["audit_scope"]["artifact_type"],
        "access_tier": "hosted_api",
        "training_allowed": ("PROHIBITED: outputs are non-commercial use only and not to be used "
                             "for training other machine learning models"),
        "evaluation_allowed": "interpretation_hypothesis_generation_and_evaluation_only",
        "redistribution_class": f"unresolved: {lic['redistribution_of_api_outputs']}",
        "license_status": flat(lic),
        "protected_or_exposed": atlas["checkpoint_findings"]["alphagenome_atlas_api"]["exposure_state"],
        "prior_exposure_evidence": atlas["checkpoint_findings"]["alphagenome_atlas_api"]["reason"],
        "exposure_status_field": atlas["checkpoint_findings"]["alphagenome_atlas_api"]["exposure_confidence"],
        "split_outer_role": "evaluation_only",
        "label_visibility": "not_applicable",
        "n_blockers_recorded": str(len(atlas["remaining_blockers"])),
        "blocker_amended_or_lifted": "False",
        "verdict_source": f"{rel_atlas}:{json_key_line(ATLAS_AUDIT, 'outputs')} license_disposition.outputs",
        "verdict_line_redistribution": f"{rel_atlas}:{json_key_line(ATLAS_AUDIT, 'redistribution_of_api_outputs')}",
        "verdict_line_exposure": f"{rel_atlas}:{json_key_line(ATLAS_AUDIT, 'exposure_state')}",
    })
    rows.append({
        "artifact_id": "alphagenome_released_weights_all_folds",
        "artifact_kind": "gated_released_weights",
        "registry_entry": "Analysis/MASLD_Model_Benchmark/config/artifacts/models/alphagenome/exposure_audit.json",
        "title": "google/alphagenome-all-folds gated Hugging Face checkpoint",
        "access_tier": "gated_access_not_granted_to_this_account",
        "training_allowed": ("adaptation permitted under the accepted noncommercial model licence, "
                             "but ACCESS IS NOT GRANTED: probe 2026-09-14 returned GatedRepoError 403, "
                             "so no local weight exists to train"),
        "evaluation_allowed": "local_inference_after_terms_acceptance",
        "redistribution_class": "gated_noncommercial_no_weight_redistribution",
        "license_status": "weights gated_noncommercial; code Apache-2.0; derivative and output-trained models subject to the noncommercial model terms",
        "protected_or_exposed": "target_label_unexposed_for_gse289173_medium_confidence",
        "prior_exposure_evidence": "checkpoint_release_relation: GSE289173 was public before the checkpoint repository revision",
        "exposure_status_field": "ineligible_current_terms_and_incomplete_four_lineage_output_roster",
        "split_outer_role": "not_created_no_local_weights",
        "label_visibility": "not_applicable",
        "n_blockers_recorded": "7",
        "blocker_amended_or_lifted": "False",
        "verdict_source": "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md:20-25 (gated, 403); "
                          "docs/technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md:389 (released weights row)",
        "verdict_line_redistribution": "Analysis/MASLD_Model_Benchmark/config/artifacts/models/alphagenome/exposure_audit.json:"
                                       + str(json_key_line(MODEL_ARTIFACTS / "alphagenome/exposure_audit.json", "weights")),
        "verdict_line_exposure": "Analysis/MASLD_Model_Benchmark/config/artifacts/models/alphagenome/exposure_audit.json:"
                                 + str(json_key_line(MODEL_ARTIFACTS / "alphagenome/exposure_audit.json", "exposure_state")),
    })

    # repository-wide policy rows read off resources.toml by line
    rt = toml_scalars_with_lines(RESOURCES)
    for key, aid in (("redistribution.human_raw_data_policy", "policy_human_raw_data"),
                     ("redistribution.controlled_data_policy", "policy_controlled_data"),
                     ("redistribution.model_weight_policy", "policy_model_weights"),
                     ("redistribution.default_data_policy", "policy_default_data")):
        e = rt.get(key) or rt.get(key.split(".", 1)[1])
        if not e:
            continue
        rows.append({
            "artifact_id": aid,
            "artifact_kind": "repository_policy",
            "registry_entry": f"Analysis/MASLD_Model_Benchmark/config/resources.toml:{e['line']}",
            "title": key,
            "access_tier": "not_applicable",
            "training_allowed": "not_applicable_policy_row",
            "evaluation_allowed": "not_applicable_policy_row",
            "redistribution_class": e["value"],
            "license_status": "owner decision 2026-09-07, docs/decisions/2026-09-07-derivative-weight-release.md",
            "protected_or_exposed": "not_applicable",
            "prior_exposure_evidence": "not_applicable",
            "exposure_status_field": "in_force",
            "split_outer_role": "not_applicable",
            "label_visibility": "not_applicable",
            "n_blockers_recorded": "0",
            "blocker_amended_or_lifted": "True" if "2026_09_07" in e["value"] or "derivative" in e["value"] else "False",
            "verdict_source": f"Analysis/MASLD_Model_Benchmark/config/resources.toml:{e['line']}",
            "verdict_line_redistribution": f"Analysis/MASLD_Model_Benchmark/config/resources.toml:{e['line']}",
            "verdict_line_exposure": "not_applicable",
        })

    # unregistered sources
    stamp = out.name.split("-")[-1]
    for u in UNREGISTERED:
        rows.append({
            "artifact_id": u["artifact_id"],
            "artifact_kind": "unregistered_source_in_use",
            "registry_entry": "NONE - no config/datasets toml exists",
            "title": u["title"],
            "access_tier": "UNVERIFIED",
            "training_allowed": "UNVERIFIED - no registry entry records any training right for this source",
            "evaluation_allowed": ("in use as an evaluation label on disk; the right to do so is "
                                   "UNVERIFIED because no entry exists"),
            "redistribution_class": "UNVERIFIED",
            "license_status": "UNVERIFIED",
            "protected_or_exposed": "exposed - already read by analyses on disk",
            "prior_exposure_evidence": "; ".join(u["paths"]),
            "exposure_status_field": "unregistered",
            "split_outer_role": u["program_role"],
            "label_visibility": "development_visible",
            "n_blockers_recorded": "0",
            "blocker_amended_or_lifted": "False",
            "verdict_source": f"absence check: no file in {DATASETS.relative_to(PROJ)} carries this dataset_id",
            "verdict_line_redistribution": f"drafts/{u['artifact_id']}.toml (DRAFT, not installed)",
            "verdict_line_exposure": ATLAS_AUDIT.relative_to(PROJ).as_posix()
                                     + f":{json_key_line(ATLAS_AUDIT, 'development_sources')} development_sources",
        })
        (drafts / f"{u['artifact_id']}.toml").write_text(DRAFT_TOML.format(
            aid=u["artifact_id"], title=u["title"], note=u["note"].replace('"', "'"),
            paths="; ".join(u["paths"]), role=u["program_role"],
            stamp=f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}"))

    # which artifacts this program actually consumes, per IMPLEMENTATION_SPEC sections 3-7
    PROGRAM = {
        "gse281364": "A1 reporter reliability, E1/E2/E3 reporter side",
        "gse281367": "C1 endpoint 2 allelic imbalance, B1 measured imbalance",
        "gse244832": "C1 endpoint 2 allelic imbalance replication",
        "gse296875": "named development source in the Atlas run audit",
        "gse267145": "chromatin predictor v1.2 development source",
        "gse289173": "sealed source; touched by nothing in this program",
        "current_coloc_gwas": "Track 0 signal universe",
        "resource_atlas_current": "Track 0 signal universe",
        "alphagenome": "C1 archived zero-shot comparator only",
        "borzoi": "C1 model family 2 (liver ATAC/DNase)",
        "chrombpnet": "C1 model family 1 (project-trained adult hepatocyte ensemble)",
        "hyenadna": "C2/C3 frozen-embedding difference head",
        "caduceus": "C2/C3 frozen-embedding difference head",
        "negative_controls": "C1.3 allele-identity and position control",
        "alphagenome_hosted_atlas_api": "B1 research layer channels, C1 archived comparator, E3 model API",
        "alphagenome_released_weights_all_folds": "BLOCKED: adaptation cannot start, gated 403",
        "currin_2025_caqtl": "C1 endpoint 1 label, E1 bridge, E2 transfer target",
        "broadaway_2024_liver_eqtl_meta": "E1 bridge endogenous side",
        "gtex_v8_liver_sqtl": "E1 bridge endogenous side",
        "leafcutter_joint_v2_bulk_cohorts": "Atlas run development source, not otherwise used here",
    }
    for r in rows:
        r["program_role_in_alphagenome_program"] = PROGRAM.get(r["artifact_id"], "not_used_by_this_program")

    cols = ["artifact_id", "artifact_kind", "program_role_in_alphagenome_program",
            "registry_entry", "title", "access_tier",
            "training_allowed", "evaluation_allowed", "redistribution_class", "license_status",
            "protected_or_exposed", "prior_exposure_evidence", "exposure_status_field",
            "split_outer_role", "label_visibility", "n_blockers_recorded",
            "blocker_amended_or_lifted", "verdict_source", "verdict_line_redistribution",
            "verdict_line_exposure"]
    with (tables / "rights_exposure.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(flat(r.get(c, NA)) for c in cols) + "\n")

    counts = {
        "n_rows": len(rows),
        "n_registered_datasets": n_datasets,
        "n_model_artifacts": n_models,
        "n_model_authorities": n_auth,
        "n_hosted_or_gated_model_rows": 2,
        "n_policy_rows": sum(1 for r in rows if r["artifact_kind"] == "repository_policy"),
        "n_unregistered_sources": len(UNREGISTERED),
        "training_allowed_value_counts": {},
        "redistribution_class_value_counts": {},
        "protected_or_exposed_value_counts": {},
        "training_rule": TRAINING_RULE,
    }
    for r in rows:
        for field, key in (("training_allowed", "training_allowed_value_counts"),
                           ("redistribution_class", "redistribution_class_value_counts"),
                           ("protected_or_exposed", "protected_or_exposed_value_counts")):
            v = flat(r.get(field))[:80]
            counts[key][v] = counts[key].get(v, 0) + 1
    (tables / "rights_exposure_summary.json").write_text(json.dumps(counts, indent=1))

    # MANIFEST of inputs
    man = [("path", "sha256", "bytes", "role")]
    small = [RESOURCES, DECISION, ATLAS_AUDIT, PLAN,
             PROJ / "scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md",
             PROJ / "docs/RESULTS.md",
             PROJ / "docs/technical/alphagenome_atlas_execution_2026-09-09.md",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/50_haplotype_prespec.json",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_summary.json",
             PROJ / "scripts/analysis/alphagenome_atlas/50_haplotype_additivity.py",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p3a-benchmarks-20260909T190214Z/tables/mpra_dav_calls_S2_S5.tsv",
             PROJ / "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction/HepG2.source_reproduction_qc.json",
             PROJ / "GWAS/finemapping/results/seqfunc/mpra_benchmark/v2/source_reproduction/LX2.source_reproduction_qc.json",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p3-ase-20260909T192926Z/tables/allelic_concordance.json",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p6f-indel-rescue-20260914T135052Z/tables/indel_rescue_effects.tsv",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p6d-dossier-20260914T173952Z/tables/clinical_variant_dossier.tsv",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p0-dbsnp-orientation-20260914T121856Z/SUPERSEDED.txt",
             PROJ / "GWAS/finemapping/results/alphagenome_atlas/p0-dbsnp-orientation-20260914T174542Z/SUPERSEDED.txt",
             ]
    for p in sorted(DATASETS.glob("*.toml")):
        small.append(p)
    for p in small:
        if p.exists():
            man.append((str(p.relative_to(PROJ)), sha256(p), str(p.stat().st_size), "register input"))
    with (prov / "MANIFEST.tsv").open("w") as fh:
        for r in man:
            fh.write("\t".join(r) + "\n")

    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True)
    (prov / "pip_freeze.txt").write_text(freeze.stdout)
    (prov / "python_version.txt").write_text(sys.version + "\n" + sys.executable + "\n")
    print(json.dumps({k: v for k, v in counts.items() if not k.endswith("value_counts")}, indent=1))
    print("wrote", tables / "rights_exposure.tsv")


if __name__ == "__main__":
    main()
