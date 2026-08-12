#!/usr/bin/env python3
"""Construct isolated Plan 60 synthetic projects for fail-closed tests."""

from __future__ import annotations

import csv
import json
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Callable

from release_common import (
    ARTIFACT_FIELDS,
    CANDIDATE_ID,
    CLOSURE_FIELDS,
    PROTECTED_SCOPE_FIELDS,
    REQUIRED_SCIENTIFIC_CHECKS,
    SCIENTIFIC_REGISTRY_FIELDS,
    WORKSTREAM_LANES,
    atomic_write_tsv,
    capture_protected_baseline,
    project_relative,
    read_tsv_exact,
    sha256_file,
)
from release_products import (
    BASE_SELECTION_FIELDS,
    EXPECTED_FIGURE_TITLES,
    SOURCE_ROW_FIELDS,
)
from fibrosis_candidate_contract import (
    FIBROSIS_ANALYSIS_ID,
    FIBROSIS_BUNDLE_FILES,
    FIBROSIS_FIXTURE_N_GENES,
    FIBROSIS_PRIMARY_RELATIVE,
    FIBROSIS_PRODUCER_ARTIFACTS,
    FIBROSIS_VALIDATED_ARTIFACTS,
    TRUE_KLEINER_COHORTS,
    TRUE_KLEINER_TRANSITIONS,
    fibrosis_artifact_id,
    fibrosis_artifact_role,
    fibrosis_snapshot_relpath,
)


FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures/payloads"


class SyntheticProject:
    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.contract_dir = self.root / "fixture_contract"
        self.contract_dir.mkdir(parents=True)
        self.workstream_roots = {
            "PLAN13": self.root / "Analysis/Spatial/candidates" / CANDIDATE_ID,
            "PLAN20": self.root
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "hotspot",
            "PLAN30": self.root
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "genetics_context",
            "PLAN40": self.root
            / "Analysis/Multimodal_Program_Projection/candidates"
            / CANDIDATE_ID
            / "myojin_hlf",
            "PLAN50": self.root
            / "RNA-seq/results/evidence_passports/candidates"
            / CANDIDATE_ID,
        }
        self.source_paths: dict[str, list[Path]] = {}
        self.manifest_paths: dict[str, Path] = {}
        self.closure_path = self.contract_dir / "workstream_closure.tsv"
        self.protected_scopes_path = self.contract_dir / "protected_scopes.tsv"
        self.protected_baseline_path = (
            self.contract_dir / "protected_release_baseline.tsv"
        )
        self.base_selection_path = self.contract_dir / "base_input_selection.tsv"
        self._build_workstreams()
        self._build_protected_releases()
        self._build_base_inputs()

    @property
    def candidate_root(self) -> Path:
        return (
            self.root / "RNA-seq/results/manuscript_release/candidates" / CANDIDATE_ID
        )

    def _copy_payload(self, fixture_name: str, destination: Path) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(FIXTURE_DIR / fixture_name, destination)
        return destination

    def _artifact_row(
        self,
        workstream: str,
        artifact_id: str,
        source: Path,
        role: str,
        filename: str,
    ) -> dict[str, object]:
        relative = f"inputs/{WORKSTREAM_LANES[workstream]}/{filename}"
        return {
            "artifact_id": artifact_id,
            "source_path": project_relative(self.root, source),
            "source_sha256": sha256_file(source),
            "source_bytes": source.stat().st_size,
            "snapshot_relpath": relative,
            "artifact_role": role,
            "downstream_read_path": relative,
        }

    def _source_row(self, **updates) -> dict[str, object]:
        row: dict[str, object] = {field: "" for field in SOURCE_ROW_FIELDS[:-2]}
        row.update(
            {
                "record_id": "synthetic_record",
                "group_id": "synthetic_group",
                "artifact_id": "synthetic_artifact",
                "claim_id": "synthetic_claim",
                "section_id": "resource_interface",
                "figure_id": "Figure1",
                "panel_id": "1A",
                "panel_title": "Synthetic panel",
                "panel_role": "synthetic_fixture",
                "label": "Synthetic",
                "category": "synthetic",
                "number_role": "synthetic_value",
                "plot_role": "mark",
                "value": "0",
                "display_value": "0",
                "numerator": "0",
                "denominator": "1",
                "unit": "synthetic units",
                "biological_unit": "synthetic biological unit",
                "model_contrast": "synthetic contrast",
                "effect_unit": "synthetic effect",
                "p_value": "1",
                "q_value": "1",
                "evidence_status": "indeterminate",
                "tested_universe": "synthetic prespecified universe",
                "source_dependence": "source_dependent",
                "discovery_sources": "synthetic_source",
                "evaluation_sources": "synthetic_source",
                "reuse_detail": "Synthetic fixture reuses the same named source.",
                "independence_boundary": (
                    "Synthetic fixture only; no independent replication is claimed."
                ),
                "allowed_wording": "synthetic fixture only",
                "prohibited_wording": "biological inference",
                "claim_text": "The synthetic fixture exercises a frozen evidence role.",
                "next_experiment": "replace with a real prespecified experiment",
                "manuscript_included": "true",
                "plot_order": "1",
                "is_control": "false",
            }
        )
        row.update(updates)
        if not any(
            field in updates for field in ("discovery_sources", "evaluation_sources")
        ):
            dependence = str(row["source_dependence"])
            if dependence == "independent":
                row["discovery_sources"] = "synthetic_discovery"
                row["evaluation_sources"] = "synthetic_evaluation"
                row["reuse_detail"] = (
                    "Synthetic discovery and evaluation fixtures are distinct."
                )
                row["independence_boundary"] = "Synthetic source IDs have no overlap."
            elif dependence == "partially_dependent":
                row["discovery_sources"] = "synthetic_discovery;synthetic_shared"
                row["evaluation_sources"] = "synthetic_evaluation;synthetic_shared"
                row["reuse_detail"] = (
                    "Synthetic fixtures share one source and retain one distinct source each."
                )
                row["independence_boundary"] = (
                    "Only the distinct synthetic source components are independent."
                )
        return row

    def _write_evidence(self, path: Path, rows: list[dict[str, object]]) -> Path:
        atomic_write_tsv(path, rows, SOURCE_ROW_FIELDS[:-2])
        return path

    def _build_workstreams(self) -> None:
        payloads = {
            "PLAN13": "plan13_terminal.tsv",
            "PLAN20": "plan20_terminal.tsv",
            "PLAN30": "plan30_terminal.tsv",
            "PLAN40": "plan40_fig5_verdict.tsv",
            "PLAN50": "plan50_terminal.tsv",
        }
        closure_rows = []
        closure_contract = {
            "PLAN13": (
                "accepted_main",
                "SP_INT_FINAL",
                "ACCEPTED",
                "true",
                "true",
                "figure_4",
            ),
            "PLAN20": (
                "accepted_main",
                "HS_V2_FINAL",
                "READY",
                "true",
                "true",
                "figure_2",
            ),
            "PLAN30": (
                "accepted_main",
                "GEN_FINAL",
                "coverage_limited_terminal",
                "true",
                "true",
                "figure_3",
            ),
            "PLAN40": (
                "accepted_supplement",
                "HLF_FINAL",
                "supplement_only_nonconfirmatory",
                "false",
                "true",
                "supplementary_functional",
            ),
            "PLAN50": (
                "accepted_main",
                "PASS_FINAL",
                "ACCEPTED",
                "true",
                "true",
                "figure_5",
            ),
        }
        for workstream, fixture_name in payloads.items():
            root = self.workstream_roots[workstream]
            root.mkdir(parents=True)
            terminal = self._copy_payload(fixture_name, root / fixture_name)
            self.source_paths[workstream] = [terminal]
            artifact_rows = [
                self._artifact_row(
                    workstream,
                    "terminal_verdict",
                    terminal,
                    "terminal_verdict",
                    fixture_name,
                )
            ]
            if workstream == "PLAN20":
                ready = self._copy_payload(
                    "plan20_ready.json", root / "plan20_ready.json"
                )
                self.source_paths[workstream].append(ready)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "ready_seal",
                        ready,
                        "ready_seal",
                        "plan20_ready.json",
                    )
                )
            if workstream == "PLAN13":
                evidence = self._write_evidence(
                    root / "spatial_evidence.tsv",
                    [
                        self._source_row(
                            record_id="spatial_moran_module1",
                            group_id="spatial_module1",
                            artifact_id="spatial_evidence",
                            claim_id="claim_spatial_native",
                            section_id="physical_context",
                            figure_id="Figure4",
                            panel_id="4A",
                            panel_title="Assay-native physical context",
                            panel_role="prespecified_external_challenge",
                            label="Module 1 Moran Z",
                            category="native Moran matched-null",
                            number_role="matched_null_moran_z",
                            value="0.15",
                            display_value="0.15",
                            biological_unit="synthetic donor",
                            model_contrast="native residual Moran matched null",
                            effect_unit="matched-null Z; null SD is not sampling SE",
                            evidence_status="indeterminate",
                            source_dependence="partially_dependent",
                            allowed_wording="assay-native physical context",
                            prohibited_wording="broad validation or sampling confidence interval",
                            claim_text="The prespecified spatial challenge remains assay-native and indeterminate.",
                            plot_order="1",
                        ),
                        self._source_row(
                            record_id="vu_array_module1",
                            group_id="vu_module1",
                            artifact_id="spatial_evidence",
                            claim_id="claim_vu_coverage",
                            section_id="physical_context",
                            figure_id="Figure4",
                            panel_id="4A",
                            panel_title="Assay-native physical context",
                            panel_role="prespecified_external_challenge",
                            label="Vu array",
                            category="Vu array technical context",
                            number_role="array_effect",
                            value="-0.08",
                            display_value="-0.08",
                            unit="arrays (technical/reporting units)",
                            biological_unit="donor unknown; array technical unit",
                            model_contrast="Vu array assay-native contrast",
                            effect_unit="array-level descriptive effect",
                            evidence_status="indeterminate",
                            source_dependence="source_dependent",
                            allowed_wording="technical array-level context with donor unknown",
                            prohibited_wording="biological donor n or broad validation",
                            claim_text="The Vu array result is retained as technical context with donor n unknown.",
                            plot_order="2",
                            manuscript_included="false",
                        ),
                        self._source_row(
                            record_id="cosmx_array_module1",
                            group_id="cosmx_module1",
                            artifact_id="spatial_evidence",
                            claim_id="claim_cosmx_coverage",
                            section_id="physical_context",
                            figure_id="Figure4",
                            panel_id="4A",
                            panel_title="Assay-native physical context",
                            panel_role="prespecified_external_challenge",
                            label="CosMx array",
                            category="CosMx array technical context",
                            number_role="array_effect",
                            value="0.04",
                            display_value="0.04",
                            unit="arrays (technical/reporting units)",
                            biological_unit="donor unknown; array technical unit",
                            model_contrast="CosMx array assay-native contrast",
                            effect_unit="array-level descriptive effect",
                            evidence_status="indeterminate",
                            source_dependence="source_dependent",
                            allowed_wording="technical array-level context with donor unknown",
                            prohibited_wording="biological donor n or broad validation",
                            claim_text="The CosMx array result is retained as technical context with donor n unknown.",
                            plot_order="3",
                            manuscript_included="false",
                        ),
                        self._source_row(
                            record_id="yak_module8_median_direction",
                            group_id="yak_module8",
                            artifact_id="spatial_evidence",
                            claim_id="claim_yak_module8",
                            section_id="physical_context",
                            figure_id="Figure4",
                            panel_id="4B",
                            panel_title="Zonation-adjusted lipid challenge",
                            panel_role="prespecified_external_challenge",
                            label="Module 8 median slope",
                            category="Yakubovsky descriptive direction",
                            number_role="descriptive_median_slope_direction",
                            value="-1",
                            display_value="negative",
                            numerator="",
                            denominator="",
                            unit="signed descriptive direction",
                            biological_unit="synthetic lipid-data donor",
                            model_contrast="zonation-adjusted lipid slope",
                            effect_unit="median standardized donor-slope direction",
                            evidence_status="indeterminate",
                            source_dependence="partially_dependent",
                            allowed_wording="separate descriptive median-slope and inferential signed-Stouffer directions",
                            prohibited_wording="collapsed direction or broad validation",
                            claim_text="Module 8 has separate descriptive and inferential direction summaries.",
                            plot_order="1",
                        ),
                        self._source_row(
                            record_id="yak_module8_stouffer_direction",
                            group_id="yak_module8",
                            artifact_id="spatial_evidence",
                            claim_id="claim_yak_module8",
                            section_id="physical_context",
                            figure_id="Figure4",
                            panel_id="4B",
                            panel_title="Zonation-adjusted lipid challenge",
                            panel_role="prespecified_external_challenge",
                            label="Module 8 signed Stouffer",
                            category="Yakubovsky inferential direction",
                            number_role="inferential_signed_stouffer_direction",
                            value="1",
                            display_value="positive",
                            numerator="",
                            denominator="",
                            unit="signed inferential direction",
                            biological_unit="synthetic lipid-data donor",
                            model_contrast="equal-donor signed Stouffer",
                            effect_unit="signed-Stouffer direction",
                            evidence_status="indeterminate",
                            source_dependence="partially_dependent",
                            allowed_wording="separate descriptive median-slope and inferential signed-Stouffer directions",
                            prohibited_wording="collapsed direction or broad validation",
                            claim_text="Module 8 has separate descriptive and inferential direction summaries.",
                            plot_order="2",
                        ),
                    ],
                )
                self.source_paths[workstream].append(evidence)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "spatial_evidence",
                        evidence,
                        "candidate_evidence_rows",
                        "spatial_evidence.tsv",
                    )
                )
            if workstream == "PLAN20":
                programs = self._write_evidence(
                    root / "hotspot_programs.tsv",
                    [
                        self._source_row(
                            record_id=f"hotspot_program_{index}",
                            group_id=f"hotspot_program_{index}",
                            artifact_id="hotspot_programs",
                            claim_id=f"claim_hotspot_{index}",
                            section_id="established_state_transcriptomics",
                            figure_id="Figure2",
                            panel_id="2C",
                            panel_title="Stable established-state programs",
                            panel_role="transcriptomic_hero",
                            label=f"Program {index}",
                            category="frozen Hotspot program",
                            number_role="stage_association_beta",
                            value=value,
                            display_value=value,
                            biological_unit="synthetic donor",
                            model_contrast="Healthy=0, Steatosis=1, Steatohepatitis=2 + dataset",
                            effect_unit="donor-level beta",
                            p_value="0.4",
                            q_value="0.8",
                            evidence_status="indeterminate",
                            source_dependence="independent",
                            allowed_wording="cross-sectional established-state program",
                            prohibited_wording="longitudinal progression or discrete subtype",
                            claim_text=f"Program {index} is evaluated as a cross-sectional established-state axis.",
                            plot_order=str(index),
                        )
                        for index, value in ((1, "0.12"), (2, "-0.07"))
                    ],
                )
                self.source_paths[workstream].append(programs)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "hotspot_programs",
                        programs,
                        "candidate_evidence_rows",
                        "hotspot_programs.tsv",
                    )
                )
                nmf = self._write_evidence(
                    root / "nmf_continuous_supplement.tsv",
                    [
                        self._source_row(
                            record_id="nmf_k4_p1_continuous",
                            group_id="nmf_k4_p1",
                            artifact_id="nmf_continuous_supplement",
                            claim_id="claim_nmf_k4_p1_continuous",
                            section_id="supplementary_nmf",
                            figure_id="FigureS2",
                            panel_id="S2A",
                            panel_title="Continuous k4/k6 NMF loading distributions",
                            panel_role="supplementary_continuous_programs",
                            label="k=4 P1: synthetic continuous axis",
                            category="continuous k=4 NMF axis",
                            number_role="median_continuous_loading",
                            value="0.2",
                            display_value="median=0.200; IQR 0.100-0.300; n=1,104",
                            denominator="1104",
                            unit="nonnegative continuous NMF usage/loading",
                            biological_unit="QC-filtered human bulk RNA-seq sample in the frozen NMF fit",
                            model_contrast="descriptive sample-level continuous loading distribution",
                            effect_unit="median continuous loading; not an inferential effect or class assignment",
                            p_value="",
                            q_value="",
                            evidence_status="supported",
                            source_dependence="reused_source",
                            allowed_wording="continuous interpretive NMF axis with descriptive loading distribution",
                            prohibited_wording="subtype, hard cluster, patient class, reproducible stratification, transition, or biomarker",
                            claim_text="k=4 P1 is retained only as a continuous interpretive axis across the frozen sample set.",
                            next_experiment="prospective outcome-linked validation before any patient-stratification use",
                            plot_order="41",
                        ),
                        self._source_row(
                            record_id="nmf_k4_mean_matched_cosine",
                            group_id="nmf_k4_stability",
                            artifact_id="nmf_continuous_supplement",
                            claim_id="claim_nmf_k4_stability",
                            section_id="supplementary_nmf",
                            figure_id="FigureS2",
                            panel_id="S2B",
                            panel_title="Three-seed factor stability and hard-partition boundary",
                            panel_role="supplementary_continuous_programs",
                            label="k=4: mean Hungarian-matched factor cosine",
                            category="k=4 three-seed factor audit",
                            number_role="mean_matched_cosine",
                            value="0.976",
                            display_value="0.976",
                            denominator="3",
                            unit="mean Hungarian-matched factor cosine",
                            biological_unit="NMF factorization seed/factor comparison",
                            model_contrast="three-seed internal stability audit",
                            effect_unit="factor-stability statistic; not a sample-partition reproducibility estimate",
                            p_value="",
                            q_value="",
                            evidence_status="supported",
                            source_dependence="reused_source",
                            allowed_wording="factor-axis stability with continuous interpretation only",
                            prohibited_wording="stable subtype, reproducible hard partition, patient class, clinical stratifier, or biomarker",
                            claim_text="The k=4 factor audit is shown with an explicit boundary that factor stability does not imply stable patient partitions.",
                            next_experiment="prospective outcome-linked validation before any patient-stratification use",
                            plot_order="43",
                        ),
                    ],
                )
                self.source_paths[workstream].append(nmf)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "nmf_continuous_supplement",
                        nmf,
                        "candidate_evidence_rows",
                        "nmf_continuous_supplement.tsv",
                    )
                )
            if workstream == "PLAN30":
                genetics = self._write_evidence(
                    root / "genetics_evidence.tsv",
                    [
                        self._source_row(
                            record_id="genetics_masld_traits",
                            group_id="genetics_provenance",
                            artifact_id="genetics_evidence",
                            claim_id="claim_genetics_provenance",
                            section_id="genetics_context",
                            figure_id="Figure3",
                            panel_id="3A",
                            panel_title="Phenotype provenance",
                            panel_role="genetics_boundary",
                            label="MASLD/MASH traits",
                            category="trait provenance",
                            number_role="trait_count",
                            value="3",
                            display_value="3",
                            biological_unit="GWAS phenotype",
                            model_contrast="MASLD/MASH versus proxy provenance",
                            effect_unit="count",
                            evidence_status="context_supported",
                            source_dependence="independent",
                            allowed_wording="genetically anchored candidate",
                            prohibited_wording="causal gene or universal genetic support",
                            claim_text="Phenotype provenance separates disease traits from proxy phenotypes.",
                            plot_order="1",
                        ),
                        self._source_row(
                            record_id="genetics_proxy_traits",
                            group_id="genetics_provenance",
                            artifact_id="genetics_evidence",
                            claim_id="claim_genetics_provenance",
                            section_id="genetics_context",
                            figure_id="Figure3",
                            panel_id="3A",
                            panel_title="Phenotype provenance",
                            panel_role="genetics_boundary",
                            label="Proxy traits",
                            category="trait provenance",
                            number_role="trait_count",
                            value="4",
                            display_value="4",
                            biological_unit="GWAS phenotype",
                            model_contrast="MASLD/MASH versus proxy provenance",
                            effect_unit="count",
                            evidence_status="context_supported",
                            source_dependence="independent",
                            allowed_wording="genetically anchored candidate",
                            prohibited_wording="causal gene or universal genetic support",
                            claim_text="Phenotype provenance separates disease traits from proxy phenotypes.",
                            plot_order="2",
                        ),
                        self._source_row(
                            record_id="genetics_eqtl_coverage",
                            group_id="genetics_coverage",
                            artifact_id="genetics_evidence",
                            claim_id="claim_genetics_coverage",
                            section_id="genetics_context",
                            figure_id="Figure3",
                            panel_id="3B",
                            panel_title="eQTL observability boundary",
                            panel_role="genetics_boundary",
                            label="Adequately observed",
                            category="eQTL coverage",
                            number_role="coverage_fraction",
                            value="0.42",
                            display_value="42%",
                            biological_unit="gene",
                            model_contrast="source eGene observability audit",
                            effect_unit="fraction",
                            evidence_status="indeterminate",
                            source_dependence="source_dependent",
                            allowed_wording="coverage-limited genetic evidence",
                            prohibited_wording="negative genetic finding",
                            claim_text="Genetic interpretation is bounded by source eQTL observability.",
                            plot_order="1",
                        ),
                    ],
                )
                self.source_paths[workstream].append(genetics)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "genetics_evidence",
                        genetics,
                        "candidate_evidence_rows",
                        "genetics_evidence.tsv",
                    )
                )
            if workstream == "PLAN40":
                myojin = self._write_evidence(
                    root / "myojin_supplement.tsv",
                    [
                        self._source_row(
                            record_id="myojin_omnibus",
                            group_id="myojin_omnibus",
                            artifact_id="myojin_supplement",
                            claim_id="claim_myojin_omnibus",
                            section_id="supplementary_myojin",
                            figure_id="FigureS1",
                            panel_id="S1A",
                            panel_title="Complete Myojin evidence-class stress test",
                            panel_role="supplementary_functional",
                            label="Omnibus enrichment",
                            category="Myojin palmitate challenge",
                            number_role="omnibus_p_value",
                            value="0.537",
                            display_value="0.537",
                            unit="p value",
                            biological_unit="gene-level matched-null test",
                            model_contrast="palmitate-specific survival enrichment",
                            effect_unit="omnibus p value",
                            p_value="0.537",
                            q_value="",
                            evidence_status="indeterminate",
                            source_dependence="source_dependent",
                            allowed_wording="assay-specific non-support; nonconfirmatory result",
                            prohibited_wording="validated functional mechanism or precise null",
                            claim_text="The assay-specific non-support is nonconfirmatory and remains supplementary.",
                            plot_order="1",
                        ),
                        self._source_row(
                            record_id="myojin_program_untestable",
                            group_id="myojin_program_untestable",
                            artifact_id="myojin_supplement",
                            claim_id="claim_myojin_program_untestable",
                            section_id="supplementary_myojin",
                            figure_id="FigureS1",
                            panel_id="S1B",
                            panel_title="Complete Myojin program stress test",
                            panel_role="supplementary_functional",
                            label="Program without adequate HLF gene coverage",
                            category="Myojin palmitate challenge",
                            number_role="program_testability_indicator",
                            plot_role="mark",
                            value="1",
                            display_value="untestable",
                            unit="categorical untestable indicator; not an effect",
                            biological_unit="screened gene in one HLF cell line",
                            model_contrast="palmitate-specific program challenge",
                            effect_unit="not applicable",
                            p_value="",
                            q_value="",
                            evidence_status="untestable",
                            source_dependence="source_dependent",
                            allowed_wording="nonconfirmatory program challenge; untestable in this assay",
                            prohibited_wording="zero effect, tested negative, or broad functional validation",
                            claim_text="The prespecified program was untestable in the HLF challenge and remains supplementary.",
                            plot_order="1",
                        ),
                    ],
                )
                self.source_paths[workstream].append(myojin)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "myojin_supplement",
                        myojin,
                        "candidate_evidence_rows",
                        "myojin_supplement.tsv",
                    )
                )
            if workstream == "PLAN50":
                passports = self._write_evidence(
                    root / "passport_evidence.tsv",
                    [
                        self._source_row(
                            record_id="passport_indeterminate",
                            group_id="passport_status",
                            artifact_id="passport_evidence",
                            claim_id="claim_passport_status",
                            section_id="evidence_passports",
                            figure_id="Figure5",
                            panel_id="5A",
                            panel_title="Evidence call and testability",
                            panel_role="passport_boundary",
                            label="Indeterminate",
                            category="passport semantic state",
                            number_role="gene_count",
                            value="7",
                            display_value="7",
                            biological_unit="gene",
                            model_contrast="accepted evidence semantic audit",
                            effect_unit="count",
                            evidence_status="indeterminate",
                            source_dependence="mixed",
                            allowed_wording="provenance-preserving evidence state",
                            prohibited_wording="universal score or leaderboard",
                            claim_text="Gene Catalog entries preserve indeterminate and untestable states.",
                            plot_order="1",
                        ),
                        self._source_row(
                            record_id="passport_untestable",
                            group_id="passport_status",
                            artifact_id="passport_evidence",
                            claim_id="claim_passport_status",
                            section_id="evidence_passports",
                            figure_id="Figure5",
                            panel_id="5A",
                            panel_title="Evidence call and testability",
                            panel_role="passport_boundary",
                            label="Untestable",
                            category="passport semantic state",
                            number_role="gene_count",
                            value="5",
                            display_value="5",
                            biological_unit="gene",
                            model_contrast="accepted evidence semantic audit",
                            effect_unit="count",
                            evidence_status="untestable",
                            source_dependence="mixed",
                            allowed_wording="provenance-preserving evidence state",
                            prohibited_wording="universal score or leaderboard",
                            claim_text="Gene Catalog entries preserve indeterminate and untestable states.",
                            plot_order="2",
                        ),
                        self._source_row(
                            record_id="passport_next_experiment",
                            group_id="passport_experiment",
                            artifact_id="passport_evidence",
                            claim_id="claim_passport_experiment",
                            section_id="evidence_passports",
                            figure_id="Figure5",
                            panel_id="5B",
                            panel_title="Next discriminating experiments",
                            panel_role="passport_boundary",
                            label="Requires perturbation",
                            category="next experiment",
                            number_role="gene_count",
                            value="9",
                            display_value="9",
                            biological_unit="gene",
                            model_contrast="unresolved evidence-state audit",
                            effect_unit="count",
                            evidence_status="indeterminate",
                            source_dependence="mixed",
                            allowed_wording="next discriminating experiment",
                            prohibited_wording="treatment recommendation or leaderboard",
                            claim_text="Gene Catalog entries specify the next experiment that could change a gene call.",
                            plot_order="1",
                        ),
                    ],
                )
                self.source_paths[workstream].append(passports)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "passport_evidence",
                        passports,
                        "candidate_evidence_rows",
                        "passport_evidence.tsv",
                    )
                )
                blueprint = root / "release_blueprint.json"
                blueprint.write_text(
                    json.dumps(self._release_blueprint(), indent=2, sort_keys=True)
                    + "\n",
                    encoding="utf-8",
                )
                self.source_paths[workstream].append(blueprint)
                artifact_rows.append(
                    self._artifact_row(
                        workstream,
                        "release_blueprint",
                        blueprint,
                        "release_blueprint",
                        "release_blueprint.json",
                    )
                )
            manifest = root / "terminal_artifact_manifest.tsv"
            atomic_write_tsv(manifest, artifact_rows, ARTIFACT_FIELDS)
            self.manifest_paths[workstream] = manifest
            state, gate_id, gate_verdict, main, supplement, figure_role = (
                closure_contract[workstream]
            )
            closure_rows.append(
                {
                    "workstream_id": workstream,
                    "terminal_state": state,
                    "gate_id": gate_id,
                    "gate_verdict": gate_verdict,
                    "artifact_manifest": project_relative(self.root, manifest),
                    "manifest_sha256": sha256_file(manifest),
                    "include_main": main,
                    "include_supplement": supplement,
                    "figure_role": figure_role,
                    "allowed_wording": "synthetic fixture wording only",
                    "prohibited_wording": "no biological inference",
                    "exclusion_reason": "",
                    "adjudicator": "synthetic-release-coordinator",
                    "date": "2026-08-07",
                }
            )
        atomic_write_tsv(self.closure_path, closure_rows, CLOSURE_FIELDS)

    def _release_blueprint(self) -> dict[str, object]:
        return {
            "contract_version": 1,
            "candidate_id": CANDIDATE_ID,
            "journal_branch": "Cell Genomics Resource",
            "figure_count": 5,
            "myojin_role": "supplement_only",
            "results_sections": [
                "resource_interface",
                "established_state_transcriptomics",
                "genetics_context",
                "physical_context",
                "evidence_passports",
            ],
            "figures": [
                {
                    "figure_id": "Figure1",
                    "title": EXPECTED_FIGURE_TITLES["Figure1"],
                    "panels": [
                        {
                            "panel_id": "1A",
                            "title": "Cohort overview",
                            "role": "interface",
                        },
                        {
                            "panel_id": "1B",
                            "title": "Evidence observability",
                            "role": "interface",
                        },
                    ],
                },
                {
                    "figure_id": "Figure2",
                    "title": EXPECTED_FIGURE_TITLES["Figure2"],
                    "panels": [
                        {
                            "panel_id": "2A",
                            "title": "Cross-sectional stage-associated DEG counts",
                            "role": "transcriptomic_hero",
                        },
                        {
                            "panel_id": "2B",
                            "title": "QC-passing sample-level cell composition",
                            "role": "transcriptomic_hero",
                        },
                        {
                            "panel_id": "2C",
                            "title": "Stable established-state programs",
                            "role": "transcriptomic_hero",
                        },
                    ],
                },
                {
                    "figure_id": "Figure3",
                    "title": EXPECTED_FIGURE_TITLES["Figure3"],
                    "panels": [
                        {
                            "panel_id": "3A",
                            "title": "Phenotype provenance",
                            "role": "genetics_boundary",
                        },
                        {
                            "panel_id": "3B",
                            "title": "eQTL observability boundary",
                            "role": "genetics_boundary",
                        },
                    ],
                },
                {
                    "figure_id": "Figure4",
                    "title": EXPECTED_FIGURE_TITLES["Figure4"],
                    "panels": [
                        {
                            "panel_id": "4A",
                            "title": "Assay-native physical context",
                            "role": "prespecified_external_challenge",
                        },
                        {
                            "panel_id": "4B",
                            "title": "Zonation-adjusted lipid challenge",
                            "role": "prespecified_external_challenge",
                        },
                    ],
                },
                {
                    "figure_id": "Figure5",
                    "title": EXPECTED_FIGURE_TITLES["Figure5"],
                    "panels": [
                        {
                            "panel_id": "5A",
                            "title": "Evidence call and testability",
                            "role": "passport_boundary",
                        },
                        {
                            "panel_id": "5B",
                            "title": "Next discriminating experiments",
                            "role": "passport_boundary",
                        },
                    ],
                },
            ],
            "supplementary_figures": [
                {
                    "figure_id": "FigureS1",
                    "title": "Prespecified Myojin HLF functional stress test",
                    "panels": [
                        {
                            "panel_id": "S1A",
                            "title": "Complete Myojin evidence-class stress test",
                            "role": "supplementary_functional",
                        },
                        {
                            "panel_id": "S1B",
                            "title": "Complete Myojin program stress test",
                            "role": "supplementary_functional",
                        },
                    ],
                },
                {
                    "figure_id": "FigureS2",
                    "title": "Continuous NMF axes and factor stability",
                    "panels": [
                        {
                            "panel_id": "S2A",
                            "title": "Continuous k4/k6 NMF loading distributions",
                            "role": "supplementary_continuous_programs",
                        },
                        {
                            "panel_id": "S2B",
                            "title": "Three-seed factor stability and hard-partition boundary",
                            "role": "supplementary_continuous_programs",
                        },
                    ],
                },
            ],
            "sources": [
                {
                    "source_key": "BASE:cohort_overview",
                    "snapshot_path": "inputs/BASE/cohort_overview.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": f"BASE:{fibrosis_artifact_id(FIBROSIS_PRIMARY_RELATIVE)}",
                    "snapshot_path": fibrosis_snapshot_relpath(
                        FIBROSIS_PRIMARY_RELATIVE
                    ),
                    "kind": "fibrosis_transition_raw",
                },
                {
                    "source_key": "BASE:sample_composition",
                    "snapshot_path": "inputs/BASE/sample_composition.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": "PLAN20:hotspot_programs",
                    "snapshot_path": "inputs/HS-V2/hotspot_programs.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": "PLAN20:nmf_continuous_supplement",
                    "snapshot_path": "inputs/HS-V2/nmf_continuous_supplement.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": "PLAN30:genetics_evidence",
                    "snapshot_path": "inputs/GEN/genetics_evidence.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": "PLAN13:spatial_evidence",
                    "snapshot_path": "inputs/SP-INT/spatial_evidence.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": "PLAN50:passport_evidence",
                    "snapshot_path": "inputs/PASS/passport_evidence.tsv",
                    "kind": "evidence_rows",
                },
                {
                    "source_key": "PLAN40:myojin_supplement",
                    "snapshot_path": "inputs/HLF/myojin_supplement.tsv",
                    "kind": "evidence_rows",
                },
            ],
        }

    def _write_csv(self, path: Path, rows, fields) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t" if path.suffix == ".tsv" else ",",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _build_fibrosis_candidate_bundle(self, bundle: Path) -> Path:
        """Create a small, fully sealed, biologically valid null fixture."""

        for directory in ("results", "audits", "manifests"):
            (bundle / directory).mkdir(parents=True, exist_ok=True)

        external = self.root / "fixture_fibrosis_sources"
        external.mkdir()
        input_sources = {
            "merged_dge": external / "merged_dge.rds",
            "meta_matched": external / "meta_matched.rds",
            "sample_qc": external / "sample_qc_report.csv",
            "gse193066_source_metadata": external / "gse193066_metadata.tsv",
        }
        code_sources = {
            "producer": external / "producer.R",
            "validator": external / "validator.py",
            "slurm_wrapper": external / "wrapper.sbatch",
            "lvqw_engine": external / "de_engine_lvqw.R",
        }
        for artifact_id, source in {**input_sources, **code_sources}.items():
            source.write_text(f"synthetic {artifact_id}\n", encoding="utf-8")

        # Exact-stage donor universe: 664 unique donors. GSE193066 contributes
        # 105 eligible first biopsies and never contributes F4.
        other_cohorts = [
            cohort for cohort in TRUE_KLEINER_COHORTS if cohort != "GSE193066"
        ]
        gse_counts = {0: 20, 1: 30, 2: 30, 3: 25, 4: 0}
        other_counts = {
            0: (22, 21, 21, 21, 21),
            1: (32, 32, 31, 31, 31),
            2: (29, 29, 29, 29, 28),
            3: (22, 22, 22, 21, 21),
            4: (9, 9, 9, 9, 8),
        }
        stage_samples: dict[int, list[dict[str, object]]] = defaultdict(list)
        gse_first_rows = []
        donor_number = 0
        gse_number = 0
        for stage in range(5):
            for cohort, count in zip(other_cohorts, other_counts[stage], strict=True):
                for _ in range(count):
                    donor_number += 1
                    sample_id = f"{cohort}_S{stage}_{donor_number:04d}"
                    stage_samples[stage].append(
                        {
                            "sample_id": sample_id,
                            "donor_id": f"{cohort}::{sample_id}",
                            "donor_id_source": "single_biopsy_sample_id",
                            "dataset": cohort,
                            "source_title": sample_id,
                            "biopsy": "single deposited biopsy",
                            "fibrosis_stage": stage,
                        }
                    )
            for _ in range(gse_counts[stage]):
                gse_number += 1
                sample_id = f"SRRFIX{gse_number:06d}"
                title = f"HUnafld{gse_number:03d}" + ("_1" if gse_number <= 58 else "")
                row = {
                    "sample_id": sample_id,
                    "donor_id": f"GSE193066::HUnafld{gse_number:03d}",
                    "donor_id_source": "deposited_title_root",
                    "dataset": "GSE193066",
                    "source_title": title,
                    "biopsy": "1st biopsy",
                    "fibrosis_stage": stage,
                }
                stage_samples[stage].append(row)
                gse_first_rows.append(row)
        if (
            sum(len(rows) for rows in stage_samples.values()) != 664
            or gse_number != 105
        ):
            raise RuntimeError("synthetic true-Kleiner donor fixture census drift")

        crosswalk = []
        for index in range(1, 107):
            first_sample = (
                gse_first_rows[index - 1]["sample_id"]
                if index <= 105
                else f"SRRFIX{index:06d}"
            )
            first_title = f"HUnafld{index:03d}" + ("_1" if index <= 58 else "")
            crosswalk.append(
                {
                    "sample_id": first_sample,
                    "source_title": first_title,
                    "donor_id": f"GSE193066::HUnafld{index:03d}",
                    "biopsy": "1st biopsy",
                    "is_first_biopsy": "TRUE",
                    "included_primary": "TRUE" if index <= 105 else "FALSE",
                }
            )
            if index <= 58:
                crosswalk.append(
                    {
                        "sample_id": f"SRRFIX2{index:05d}",
                        "source_title": f"HUnafld{index:03d}_2",
                        "donor_id": f"GSE193066::HUnafld{index:03d}",
                        "biopsy": "2nd biopsy",
                        "is_first_biopsy": "FALSE",
                        "included_primary": "FALSE",
                    }
                )
        self._write_csv(
            bundle / "audits/gse193066_donor_biopsy_crosswalk.tsv",
            crosswalk,
            (
                "sample_id",
                "source_title",
                "donor_id",
                "biopsy",
                "is_first_biopsy",
                "included_primary",
            ),
        )

        sample_rows = []
        cohort_rows = []
        design_rows = []
        for contrast, spec in TRUE_KLEINER_TRANSITIONS.items():
            transition = str(spec["transition"])
            selected = [
                row
                for row in (
                    *stage_samples[int(spec["low_stage"])],
                    *stage_samples[int(spec["high_stage"])],
                )
                if row["dataset"] in spec["cohorts"]
            ]
            for row in selected:
                sample_rows.append(
                    {
                        **row,
                        "transition": transition,
                        "contrast": contrast,
                        "arm": (
                            "low"
                            if row["fibrosis_stage"] == spec["low_stage"]
                            else "high"
                        ),
                        "inferred_sex": "Female",
                        "pass_technical": "TRUE",
                    }
                )
            for cohort in sorted(spec["cohorts"]):
                cohort_subset = [row for row in selected if row["dataset"] == cohort]
                n_low = sum(
                    row["fibrosis_stage"] == spec["low_stage"] for row in cohort_subset
                )
                n_high = len(cohort_subset) - n_low
                cohort_rows.append(
                    {
                        "transition": transition,
                        "contrast": contrast,
                        "dataset": cohort,
                        "n_low": n_low,
                        "n_high": n_high,
                        "n_total": len(cohort_subset),
                        "both_arms_present": "TRUE",
                    }
                )
            design_rows.append(
                {
                    "transition": transition,
                    "contrast": contrast,
                    "formula": "~ dataset + inferred_sex + fib_group",
                    "coefficient": "fib_grouphigh",
                    "dropped_terms": "",
                    "design_rank": "7",
                    "n_design_columns": "7",
                    "design_columns": "fixture",
                    "n_genes_tested": FIBROSIS_FIXTURE_N_GENES,
                    "n_samples": spec["n_samples"],
                    "n_cohorts": len(spec["cohorts"]),
                    "cohort_set": ";".join(sorted(spec["cohorts"])),
                    "expected_cohort_set": ";".join(sorted(spec["cohorts"])),
                    "n_unique_donors": spec["n_samples"],
                }
            )
        self._write_csv(
            bundle / "audits/sample_manifest.tsv", sample_rows, tuple(sample_rows[0])
        )
        self._write_csv(
            bundle / "audits/cohort_arm_audit.tsv", cohort_rows, tuple(cohort_rows[0])
        )
        self._write_csv(
            bundle / "audits/design_audit.tsv", design_rows, tuple(design_rows[0])
        )

        result_rows = []
        transition_summary = []
        for contrast, spec in TRUE_KLEINER_TRANSITIONS.items():
            for gene_index in range(1, FIBROSIS_FIXTURE_N_GENES + 1):
                logfc = 0.1 if gene_index % 2 else -0.1
                result_rows.append(
                    {
                        "gene": f"SYNTH{gene_index:06d}",
                        "logFC": logfc,
                        "SE": "0.1",
                        "t": "1" if logfc > 0 else "-1",
                        "P.Value": "1",
                        "padj": "1",
                        "shrunk_logFC": logfc,
                        "lfsr": "1",
                        "AveExpr": "5",
                        "contrast": contrast,
                        "transition": spec["transition"],
                        "method": "limma_voom_qw_C2_true_kleiner_first_biopsy",
                        "n_cohorts": len(spec["cohorts"]),
                        "n_samples": spec["n_samples"],
                        "n_low": spec["n_low"],
                        "n_high": spec["n_high"],
                    }
                )
            transition_summary.append(
                {
                    "transition": spec["transition"],
                    "contrast": contrast,
                    "n_genes_tested": FIBROSIS_FIXTURE_N_GENES,
                    "n_deg_05": 0,
                    "n_up": 0,
                    "n_down": 0,
                    "n_low": spec["n_low"],
                    "n_high": spec["n_high"],
                    "n_samples": spec["n_samples"],
                    "n_cohorts": len(spec["cohorts"]),
                    "significance_rule": "BH padj < 0.05",
                    "lfc_gate": "none",
                }
            )
        self._write_csv(
            bundle / FIBROSIS_PRIMARY_RELATIVE,
            result_rows,
            tuple(result_rows[0]),
        )
        self._write_csv(
            bundle / "results/transition_summary.tsv",
            transition_summary,
            tuple(transition_summary[0]),
        )

        # These producer artifacts are trust-bound by both manifests; their
        # scientific detail is tested by the upstream validator, while Plan60
        # independently rederives the primary family and donor/census gates.
        simple_audits = {
            "audits/null_shuffle_results.csv": (
                "gene,contrast,P.Value,padj\nSYNTH000001,F1_vs_F0,1,1\n"
            ),
            "audits/null_shuffle_summary.tsv": ("contrast\tn_deg_05\nF1_vs_F0\t0\n"),
            "audits/model_design.tsv": (
                "sample_id\ttransition\tfib_grouphigh\nfixture\tF0_to_F1\t0\n"
            ),
            "audits/shuffle_audit.tsv": (
                "transition\tdataset\tn_changed\nF0_to_F1\tGSE130970\t2\n"
            ),
        }
        for relative, payload in simple_audits.items():
            (bundle / relative).write_text(payload, encoding="utf-8")

        manifest_fields = ("artifact_id", "path", "bytes", "sha256")
        for filename, sources in (
            ("input_manifest.tsv", input_sources),
            ("code_manifest.tsv", code_sources),
        ):
            rows = [
                {
                    "artifact_id": artifact_id,
                    "path": project_relative(self.root, source),
                    "bytes": source.stat().st_size,
                    "sha256": sha256_file(source),
                }
                for artifact_id, source in sources.items()
            ]
            atomic_write_tsv(bundle / "manifests" / filename, rows, manifest_fields)
        environment_rows = [
            {
                "component": component,
                "version": "fixture",
                "hostname": "fixture",
                "slurm_job_id": "900001",
                "slurm_partition": "cpu",
                "slurm_cpus_per_task": "1",
                "generated_at_utc": "2026-08-08 UTC",
            }
            for component in ("R", "data.table", "edgeR", "limma", "ashr")
        ]
        atomic_write_tsv(
            bundle / "manifests/environment_manifest.tsv",
            environment_rows,
            tuple(environment_rows[0]),
        )
        (bundle / "manifests/sessionInfo.txt").write_text(
            "synthetic session\n", encoding="utf-8"
        )
        run_contract = {
            "candidate_id": CANDIDATE_ID,
            "analysis_id": FIBROSIS_ANALYSIS_ID,
            "model": "~ dataset + inferred_sex + fib_group",
            "coefficient": "fib_grouphigh",
            "method": "limma_voom_qw_C2_true_kleiner_first_biopsy",
            "tested_universe": str(FIBROSIS_FIXTURE_N_GENES),
            "significance_rule": "BH padj < 0.05 per complete contrast family",
            "lfc_gate": "none",
            "allowlist": ";".join(TRUE_KLEINER_COHORTS),
            "hard_exclusions": "GSE213621;PRJNA512027",
            "gse193066_timepoint": "1st biopsy only",
            "gse193066_donor_key": "deposited !Sample_title with terminal _1/_2 removed",
            "shuffle_seed": "20260808",
            "replication_unit": "one first-biopsy biological donor/sample",
            "canonical_write": "false",
        }
        atomic_write_tsv(
            bundle / "manifests/run_contract.tsv",
            [{"key": key, "value": value} for key, value in run_contract.items()],
            ("key", "value"),
        )
        producer_status = {
            "status": "PRODUCER_COMPLETE_PENDING_VALIDATION",
            "candidate_id": CANDIDATE_ID,
            "analysis_id": FIBROSIS_ANALYSIS_ID,
            "n_contrasts": "4",
            "n_result_rows": str(FIBROSIS_FIXTURE_N_GENES * 4),
            "n_null_rows": str(FIBROSIS_FIXTURE_N_GENES * 4),
            "n_unique_primary_donors": "664",
            "canonical_write": "FALSE",
            "completed_at_utc": "2026-08-08 UTC",
        }
        atomic_write_tsv(
            bundle / "manifests/producer_status.tsv",
            [producer_status],
            tuple(producer_status),
        )

        producer_rows = []
        for relative in sorted(FIBROSIS_PRODUCER_ARTIFACTS):
            source = bundle / relative
            producer_rows.append(
                {
                    "artifact_path": relative,
                    "artifact_role": (
                        "result"
                        if relative.startswith("results/")
                        else "audit"
                        if relative.startswith("audits/")
                        else "manifest"
                    ),
                    "bytes": source.stat().st_size,
                    "sha256": sha256_file(source),
                }
            )
        atomic_write_tsv(
            bundle / "manifests/producer_artifact_manifest.tsv",
            producer_rows,
            ("artifact_path", "artifact_role", "bytes", "sha256"),
        )
        checks = [
            {
                "check_id": f"V{index:02d}_FIXTURE",
                "status": "PASS",
                "detail": "synthetic fixture",
            }
            for index in range(1, 14)
        ]
        atomic_write_tsv(
            bundle / "manifests/validation_checks.tsv",
            checks,
            ("check_id", "status", "detail"),
        )
        report_summary = {
            contrast: {
                "n_genes_tested": FIBROSIS_FIXTURE_N_GENES,
                "n_deg_05": 0,
                "n_up": 0,
                "n_down": 0,
                "n_samples": spec["n_samples"],
                "n_cohorts": len(spec["cohorts"]),
                "null_n_deg_05": 0,
            }
            for contrast, spec in TRUE_KLEINER_TRANSITIONS.items()
        }
        report = {
            "status": "READY",
            "candidate_id": CANDIDATE_ID,
            "analysis_id": FIBROSIS_ANALYSIS_ID,
            "canonical_promotion_status": "not_promoted",
            "job_id": "900001",
            "n_checks": 13,
            "summary": report_summary,
            "validated_at_utc": "2026-08-08T00:00:00+00:00",
        }
        (bundle / "manifests/validation_report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        validated_rows = []
        for relative in sorted(FIBROSIS_VALIDATED_ARTIFACTS):
            source = bundle / relative
            validated_rows.append(
                {
                    "artifact_path": relative,
                    "bytes": source.stat().st_size,
                    "sha256": sha256_file(source),
                }
            )
        atomic_write_tsv(
            bundle / "manifests/validated_artifact_manifest.tsv",
            validated_rows,
            ("artifact_path", "bytes", "sha256"),
        )
        ready = {
            "status": "READY",
            "candidate_id": CANDIDATE_ID,
            "analysis_id": FIBROSIS_ANALYSIS_ID,
            "canonical_promotion_status": "not_promoted",
            "validation_report_sha256": sha256_file(
                bundle / "manifests/validation_report.json"
            ),
            "validation_checks_sha256": sha256_file(
                bundle / "manifests/validation_checks.tsv"
            ),
            "validated_artifact_manifest_sha256": sha256_file(
                bundle / "manifests/validated_artifact_manifest.tsv"
            ),
        }
        (bundle / "READY").write_text(
            json.dumps(ready, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        observed = {
            path.relative_to(bundle).as_posix()
            for path in bundle.rglob("*")
            if path.is_file()
        }
        if observed != set(FIBROSIS_BUNDLE_FILES):
            raise RuntimeError("synthetic fibrosis sealed-bundle file universe drift")
        return bundle

    def reseal_fibrosis_bundle(self) -> None:
        """Re-sign fixture hashes after an intentional semantic mutation."""

        bundle = self.fibrosis_bundle_root
        producer_rows = []
        for relative in sorted(FIBROSIS_PRODUCER_ARTIFACTS):
            source = bundle / relative
            producer_rows.append(
                {
                    "artifact_path": relative,
                    "artifact_role": (
                        "result"
                        if relative.startswith("results/")
                        else "audit"
                        if relative.startswith("audits/")
                        else "manifest"
                    ),
                    "bytes": source.stat().st_size,
                    "sha256": sha256_file(source),
                }
            )
        atomic_write_tsv(
            bundle / "manifests/producer_artifact_manifest.tsv",
            producer_rows,
            ("artifact_path", "artifact_role", "bytes", "sha256"),
        )
        validated_rows = []
        for relative in sorted(FIBROSIS_VALIDATED_ARTIFACTS):
            source = bundle / relative
            validated_rows.append(
                {
                    "artifact_path": relative,
                    "bytes": source.stat().st_size,
                    "sha256": sha256_file(source),
                }
            )
        atomic_write_tsv(
            bundle / "manifests/validated_artifact_manifest.tsv",
            validated_rows,
            ("artifact_path", "bytes", "sha256"),
        )
        ready_path = bundle / "READY"
        ready = json.loads(ready_path.read_text(encoding="utf-8"))
        ready.update(
            {
                "validation_report_sha256": sha256_file(
                    bundle / "manifests/validation_report.json"
                ),
                "validation_checks_sha256": sha256_file(
                    bundle / "manifests/validation_checks.tsv"
                ),
                "validated_artifact_manifest_sha256": sha256_file(
                    bundle / "manifests/validated_artifact_manifest.tsv"
                ),
            }
        )
        ready_path.write_text(
            json.dumps(ready, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        selection = read_tsv_exact(self.base_selection_path, BASE_SELECTION_FIELDS)
        for row in selection:
            if row["artifact_id"] not in {
                fibrosis_artifact_id(relative) for relative in FIBROSIS_BUNDLE_FILES
            }:
                continue
            relative = row["snapshot_relpath"].removeprefix(
                "inputs/BASE/fibrosis_candidate/"
            )
            source = bundle / relative
            row["source_sha256"] = sha256_file(source)
            row["source_bytes"] = str(source.stat().st_size)
        atomic_write_tsv(self.base_selection_path, selection, BASE_SELECTION_FIELDS)

    def _build_base_inputs(self) -> None:
        source_root = self.root / "fixture_base_sources"
        source_root.mkdir()
        fibrosis_bundle = self._build_fibrosis_candidate_bundle(
            source_root / "fibrosis_candidate"
        )
        fibrosis = fibrosis_bundle / FIBROSIS_PRIMARY_RELATIVE
        composition = self._write_evidence(
            source_root / "sample_composition.tsv",
            [
                self._source_row(
                    record_id=f"composition_{order}",
                    group_id=f"composition_{celltype}",
                    artifact_id="sample_composition",
                    claim_id=f"claim_composition_{order}",
                    section_id="established_state_transcriptomics",
                    figure_id="Figure2",
                    panel_id="2B",
                    panel_title="QC-passing sample-level cell composition",
                    panel_role="transcriptomic_hero",
                    label=celltype,
                    category="cell composition",
                    number_role="composition_effect",
                    value=effect,
                    display_value=f"{float(effect):.3f}",
                    denominator="1221",
                    unit="arcsine-square-root proportion difference",
                    biological_unit="QC-passing human liver sample with MuSiC deconvolution",
                    model_contrast="Control_vs_Disease",
                    effect_unit="arcsine-square-root proportion difference",
                    p_value="0.4",
                    q_value="0.8",
                    evidence_status="indeterminate",
                    tested_universe="complete synthetic testable cell-type family",
                    source_dependence="reused_source",
                    allowed_wording="QC-filtered sample-level cross-sectional composition association",
                    prohibited_wording="donor-level, longitudinal change, lineage transition, or causal composition effect",
                    claim_text=f"{celltype} composition evidence is indeterminate in the QC-filtered synthetic fixture.",
                    next_experiment="matched tissue imaging across independently identified donors",
                    plot_order=str(order),
                )
                for order, (celltype, effect) in enumerate(
                    (
                        ("Endothelial cells", "0.025"),
                        ("Hepatocytes", "-0.053"),
                        ("Plasma cells", "0.018"),
                    ),
                    start=1,
                )
            ],
        )
        cohort = self._write_evidence(
            source_root / "cohort_overview.tsv",
            [
                self._source_row(
                    record_id="cohort_count",
                    group_id="resource_cohorts",
                    artifact_id="cohort_overview",
                    claim_id="claim_resource_cohorts",
                    section_id="resource_interface",
                    figure_id="Figure1",
                    panel_id="1A",
                    panel_title="Cohort overview",
                    panel_role="interface",
                    label="Human cohorts",
                    category="resource scale",
                    number_role="cohort_count",
                    value="9",
                    display_value="9",
                    biological_unit="human cohort",
                    model_contrast="descriptive resource inventory",
                    effect_unit="count",
                    evidence_status="descriptive",
                    source_dependence="reused_source",
                    allowed_wording="descriptive cohort inventory",
                    prohibited_wording="largest or broadest atlas",
                    claim_text="The resource integrates public human cohorts under explicit observability gates.",
                    plot_order="1",
                ),
                self._source_row(
                    record_id="observability_assays",
                    group_id="resource_observability",
                    artifact_id="cohort_overview",
                    claim_id="claim_resource_observability",
                    section_id="resource_interface",
                    figure_id="Figure1",
                    panel_id="1B",
                    panel_title="Evidence observability",
                    panel_role="interface",
                    label="Assay roles",
                    category="observability",
                    number_role="assay_role_count",
                    value="5",
                    display_value="5",
                    biological_unit="evidence role",
                    model_contrast="assay observability inventory",
                    effect_unit="count",
                    evidence_status="descriptive",
                    source_dependence="mixed",
                    allowed_wording="coverage-dependent evidence roles",
                    prohibited_wording="most modalities or universal score",
                    claim_text="Evidence roles remain distinct because assay coverage is source dependent.",
                    plot_order="1",
                ),
            ],
        )
        artifacts = [
            *[
                (
                    fibrosis_artifact_id(relative),
                    fibrosis_bundle / relative,
                    fibrosis_snapshot_relpath(relative),
                    fibrosis_artifact_role(relative),
                    (
                        "cross-sectional true-Kleiner stage-associated DEG counts"
                        if relative == FIBROSIS_PRIMARY_RELATIVE
                        else "validation/provenance artifact for the sealed true-Kleiner candidate bundle"
                    ),
                    (
                        "longitudinal progression, stale transition summary, or stand-alone unsealed coefficient table"
                        if relative == FIBROSIS_PRIMARY_RELATIVE
                        else "scientific evidence in isolation or substitution outside the sealed bundle"
                    ),
                )
                for relative in FIBROSIS_BUNDLE_FILES
            ],
            (
                "sample_composition",
                composition,
                "inputs/BASE/sample_composition.tsv",
                "figure_source_adapter",
                "QC-filtered sample-level cross-sectional composition association",
                "donor-level, cell-lineage transition, or longitudinal change",
            ),
            (
                "cohort_overview",
                cohort,
                "inputs/BASE/cohort_overview.tsv",
                "cohort_overview",
                "descriptive resource and observability inventory",
                "largest, broadest, or most-modalities claim",
            ),
        ]
        selection_rows = []
        self.base_source_paths = {}
        self.base_source_paths["fibrosis_transitions"] = fibrosis
        self.fibrosis_bundle_root = fibrosis_bundle
        for artifact_id, source, destination, role, allowed, prohibited in artifacts:
            self.base_source_paths[artifact_id] = source
            selection_rows.append(
                {
                    "selection_id": "synthetic-base-selection-v1",
                    "coordinator": "synthetic-release-coordinator",
                    "date": "2026-08-07",
                    "artifact_id": artifact_id,
                    "source_path": project_relative(self.root, source),
                    "source_sha256": sha256_file(source),
                    "source_bytes": source.stat().st_size,
                    "snapshot_relpath": destination,
                    "artifact_role": role,
                    "allowed_wording": allowed,
                    "prohibited_wording": prohibited,
                }
            )
        atomic_write_tsv(
            self.base_selection_path, selection_rows, BASE_SELECTION_FIELDS
        )

    def _build_protected_releases(self) -> None:
        release_parent = self.root / "RNA-seq/results/manuscript_release"
        release_one = release_parent / "2026-01-01-r1"
        release_two = release_parent / "2026-02-02-r2"
        release_one.mkdir(parents=True)
        release_two.mkdir(parents=True)
        (release_one / "manifest.tsv").write_text(
            "release\tvalue\nr1\tfrozen\n", encoding="utf-8"
        )
        (release_two / "manifest.tsv").write_text(
            "release\tvalue\nr2\tfrozen\n", encoding="utf-8"
        )
        v1 = self.root / "Analysis/Multimodal_Program_Projection/releases/v1"
        v1.mkdir(parents=True)
        (v1 / "frozen_programs.tsv").write_text(
            "program\tlabel\n1\tv1 frozen\n", encoding="utf-8"
        )
        scopes = [
            {
                "scope_id": "manuscript_2026-01-01-r1",
                "protection_class": "named_manuscript_release",
                "root_path": project_relative(self.root, release_one),
            },
            {
                "scope_id": "manuscript_2026-02-02-r2",
                "protection_class": "named_manuscript_release",
                "root_path": project_relative(self.root, release_two),
            },
            {
                "scope_id": "program_context_v1",
                "protection_class": "program_v1",
                "root_path": project_relative(self.root, v1),
            },
        ]
        atomic_write_tsv(self.protected_scopes_path, scopes, PROTECTED_SCOPE_FIELDS)
        capture_protected_baseline(
            self.root,
            self.protected_scopes_path,
            self.protected_baseline_path,
        )

    def closure_rows(self) -> list[dict[str, str]]:
        from release_common import read_tsv_exact

        return read_tsv_exact(self.closure_path, CLOSURE_FIELDS)

    def rewrite_closure(
        self, transform: Callable[[list[dict[str, str]]], list[dict[str, str]]]
    ) -> None:
        atomic_write_tsv(
            self.closure_path,
            transform(self.closure_rows()),
            CLOSURE_FIELDS,
        )

    def rewrite_artifact_manifest(
        self,
        workstream: str,
        transform: Callable[[list[dict[str, str]]], list[dict[str, str]]],
    ) -> None:
        from release_common import read_tsv_exact

        manifest = self.manifest_paths[workstream]
        rows = read_tsv_exact(manifest, ARTIFACT_FIELDS)
        atomic_write_tsv(manifest, transform(rows), ARTIFACT_FIELDS)

        def update_closure(closure_rows: list[dict[str, str]]) -> list[dict[str, str]]:
            for row in closure_rows:
                if row["workstream_id"] == workstream:
                    row["manifest_sha256"] = sha256_file(manifest)
            return closure_rows

        self.rewrite_closure(update_closure)

    def install_scientific_registry(
        self,
        warning_check: str | None = None,
        adjudicate_warning: bool = True,
        omit_check: str | None = None,
    ) -> Path:
        rows = []
        for check_id in REQUIRED_SCIENTIFIC_CHECKS:
            if check_id == omit_check:
                continue
            report = self.candidate_root / f"claims/validation/{check_id}.txt"
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text(
                f"synthetic validation report: {check_id}\n", encoding="utf-8"
            )
            warnings = 1 if check_id == warning_check else 0
            adjudication_path = ""
            adjudication_hash = ""
            if warnings and adjudicate_warning:
                adjudication = (
                    self.candidate_root / f"claims/adjudications/{check_id}.md"
                )
                adjudication.parent.mkdir(parents=True, exist_ok=True)
                adjudication.write_text(
                    "Synthetic warning reviewed and bounded.\n", encoding="utf-8"
                )
                adjudication_path = adjudication.relative_to(
                    self.candidate_root
                ).as_posix()
                adjudication_hash = sha256_file(adjudication)
            rows.append(
                {
                    "check_id": check_id,
                    "status": "pass",
                    "report_path": report.relative_to(self.candidate_root).as_posix(),
                    "report_sha256": sha256_file(report),
                    "warnings_count": warnings,
                    "adjudication_path": adjudication_path,
                    "adjudication_sha256": adjudication_hash,
                }
            )
        registry = self.candidate_root / "manifests/scientific_validation_registry.tsv"
        atomic_write_tsv(registry, rows, SCIENTIFIC_REGISTRY_FIELDS)
        return registry
