from __future__ import annotations

import hashlib
import json
import re
import tomllib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config"

DATASET_REQUIRED = {
    "schema_version",
    "dataset_id",
    "title",
    "role",
    "status",
    "species",
    "accession",
    "access_tier",
    "automatic_download",
    "outcome_url",
    "redistribution_class",
    "license_status",
    "biological_unit",
    "expected_biological_units",
    "native_genome_build",
    "native_annotation_release",
    "modalities",
    "pairing_levels",
    "modality_status_default",
    "exposure_status",
    "admission_blocking",
    "blockers",
    "notes",
    "provenance",
    "split",
}

MODEL_REQUIRED = {
    "model_id",
    "display_name",
    "implementation_type",
    "checkpoint_revision",
    "checkpoint_sha256",
    "license_status",
    "exposure_status",
    "status",
    "admission_blocking",
    "blockers",
    "supported_tasks",
}

TASK_REQUIRED = {
    "schema_version",
    "task_id",
    "status",
    "unit_of_inference",
    "input_modalities",
    "endpoint",
    "metrics",
    "datasets_train",
    "datasets_development",
    "datasets_sealed",
    "split_id",
    "required_pairing_levels",
    "missingness_policy",
    "admission_gates",
    "baseline_model_ids",
    "uncertainty_method",
    "resampling_units",
    "bootstrap_replicates",
    "multiplicity_family",
    "multiplicity_method",
    "promotion_gate_id",
    "promotion_gate_config_sha256",
    "primary_evaluator_id",
    "evaluator_parameters",
    "claim_gate",
}

PAIRING_LEVELS = {
    "same_molecule",
    "same_cell",
    "same_nucleus",
    "same_section",
    "adjacent_section",
    "same_sample_different_aliquot",
    "same_donor_different_tissue",
    "same_study_unpaired",
}

MODALITY_STATUSES = {
    "observed",
    "structurally_missing",
    "not_applicable",
    "below_qc",
    "unavailable_permission",
    "join_unresolved",
    "withheld_sealed",
    "derivable_not_processed",
}

TASK_IDS = {
    "cell_state_mapping",
    "variant_to_regulation",
    "rna_conditioned_atac",
    "bulk_state_transfer",
    "typed_evidence_graph",
    "perturbation_transfer",
    "gene_perturbation_response",
    "unified_model",
}

MODEL_IDS = {
    "geneformer_v1_10m",
    "geneformer_v2_104m",
    "geneformer_v2_316m",
    "scgpt_whole_human",
    "scgpt_continual",
    "uce_4l",
    "uce_33l",
    "scimilarity_v1_1",
    "cellplm_85m",
    "cellplm_85m_vae_20231027",
    "transcriptformer_tf_sapiens",
    "scprint_v1_5_medium",
    "scprint2_small_v2",
    "scprint2_medium",
    "regformer",
    "scbert",
    "sccello",
    "langcell",
    "epiagent",
    "peakvi",
    "scbasset",
    "chrombert",
    "cistopic",
    "lsi",
    "epifoundation",
    "corgi_regular",
    "borzoi_ensemble",
    "scooby_onek1k",
    "scooby_epicardioids",
    "scooby_neurips",
    "epibert",
    "epcotv2",
    "enformer",
    "sei",
    "chromdragonn",
    "crested",
    "context_borzoi",
    "epibrain",
    "corgi_plus",
    "chrombpnet",
    "bpnet",
    "mpralegnet",
    "deltasvm",
    "gkmsvm",
    "abc",
    "re2g",
    "epinformer",
    "dnabert2",
    "nucleotide_transformer",
    "hyenadna",
    "caduceus",
    "evo",
    "multivi",
    "scglue",
    "midas",
    "scmomat",
    "stabmap",
    "seurat_wnn",
    "cobolt",
    "mofaplus",
    "mira",
    "scmvp",
    "babel",
    "scbutterfly",
    "scjoint",
    "jamie",
    "monae",
    "multidgd",
    "scdiffusionx",
    "scpair",
    "gears",
    "scgen",
    "cpa",
    "cellot",
    "genepert",
    "scgpt_perturbation_head",
    "perturblib_lpm",
    "evidence_count",
    "nearest_gene",
    "node2vec",
    "bionic",
    "graphsage",
    "rgcn",
    "hgt",
    "typed_transformer",
    "late_fusion",
    "alphagenome",
    "get",
    "decima",
    "scfoundation_heads",
    "cellfm",
    "genecompass",
    "gb_cell",
    "sclong",
    "arc_state",
    "novae",
    "nichecompass",
    "spatialmeta",
    "spamosaic",
    "nicheformer",
    "graphst",
    "stagate",
    "banksy",
    "spatial_transcript_native_baselines",
    "hvg_pca_nearest_centroid",
    "hvg_pca_knn",
    "hvg_pca_logistic",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
    "harmony_logistic",
    "scvi_baseline",
    "scanvi_baseline",
    "celltypist_baseline",
    "mean_track",
    "context_only",
    "sequence_only",
    "trans_only",
    "nearest_regulatory_feature",
    "shuffled_outcome",
    "control_mean",
    "perturbed_mean",
    "ridge_perturbation",
    "bilinear_perturbation",
    "systema",
    "pca_procrustes",
    "cca",
    "mean_expression_bulk",
    "nearest_context",
    "shrunken_pseudobulk",
    "shuffled_context",
    "sequence_cnn_control",
    "sequence_transformer_control",
    "assay_native_pseudobulk",
    "graph_degree",
    "graph_source_count",
    "regularized_linear_evidence",
}


def load_toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def tables(directory: str) -> dict[str, dict]:
    return {
        path.stem: load_toml(path)
        for path in sorted((CONFIG / directory).glob("*.toml"))
    }


class RegistryFileContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.datasets_by_file = tables("datasets")
        cls.models_by_file = tables("models")
        cls.tasks_by_file = tables("tasks")
        cls.splits_by_file = tables("splits")
        cls.runtimes_by_file = tables("runtimes")
        cls.campaigns_by_file = tables("campaigns")
        cls.resources = load_toml(CONFIG / "resources.toml")
        cls.inference = load_toml(CONFIG / "evaluation" / "inference.toml")
        cls.atac_capabilities = load_toml(
            CONFIG / "evaluation" / "rna_conditioned_atac_capabilities.toml"
        )
        cls.cell_state_capabilities = load_toml(
            CONFIG / "evaluation" / "cell_state_mapping_capabilities.toml"
        )
        cls.variant_capabilities = load_toml(
            CONFIG / "evaluation" / "variant_to_regulation_capabilities.toml"
        )
        cls.variant_proxy_selection = json.loads(
            (
                CONFIG
                / "evaluation"
                / "variant_development_proxy_selection.json"
            ).read_text(encoding="utf-8")
        )
        cls.gse296875_histopathology = load_toml(
            CONFIG / "evaluation" / "gse296875_histopathology.toml"
        )
        cls.model_audit_bindings = load_toml(
            CONFIG / "evaluation" / "model_audit_bindings.toml"
        )
        cls.lead_checkpoints = load_toml(
            CONFIG / "checkpoints" / "lead_regulatory.toml"
        )

        cls.datasets = {
            value["dataset_id"]: value for value in cls.datasets_by_file.values()
        }
        cls.tasks = {value["task_id"]: value for value in cls.tasks_by_file.values()}
        cls.splits = {value["split_id"]: value for value in cls.splits_by_file.values()}

        cls.models = {}
        for family in cls.models_by_file.values():
            for model in family.get("models", []):
                model_id = model["model_id"]
                if model_id in cls.models:
                    raise AssertionError(f"duplicate model_id: {model_id}")
                cls.models[model_id] = model

    def test_all_registry_files_parse_and_have_expected_schema(self) -> None:
        self.assertGreater(len(self.datasets), 0)
        self.assertGreater(len(self.models), 0)
        self.assertEqual(set(self.tasks), TASK_IDS)

        for filename, dataset in self.datasets_by_file.items():
            self.assertEqual(DATASET_REQUIRED - dataset.keys(), set(), filename)
            self.assertEqual(dataset["schema_version"], "masld-bench-dataset-v1")
            self.assertEqual(dataset["dataset_id"], filename)
            self.assertIs(dataset["automatic_download"], False)
            self.assertNotRegex(dataset["outcome_url"], r"^https?://")
            self.assertTrue(dataset["accession"])
            self.assertTrue(dataset["native_genome_build"].strip())
            self.assertTrue(dataset["native_annotation_release"].strip())
            if "UNRESOLVED" in dataset["native_genome_build"].upper() or (
                "UNRESOLVED" in dataset["native_annotation_release"].upper()
            ):
                self.assertIs(dataset["admission_blocking"], True, filename)
            self.assertTrue(dataset["modalities"])
            self.assertTrue(dataset["pairing_levels"])
            self.assertTrue(set(dataset["pairing_levels"]) <= PAIRING_LEVELS)
            self.assertIn(dataset["modality_status_default"], MODALITY_STATUSES)
            self.assertEqual(
                {"source_url", "verified_date", "local_exposure"}
                - dataset["provenance"].keys(),
                set(),
                filename,
            )
            self.assertEqual(
                {"group_key", "outer_role", "label_visibility"}
                - dataset["split"].keys(),
                set(),
                filename,
            )

        self.assertEqual(
            self.datasets["gse249997"]["native_annotation_release"],
            "GENCODE_v29",
        )
        self.assertEqual(
            self.datasets["gse296875"]["native_genome_build"],
            "GRCh38_10x_2020-A",
        )
        self.assertIn(
            "GENCODE_v32",
            self.datasets["gse296875"]["native_annotation_release"],
        )
        self.assertEqual(
            self.datasets["gse289173"]["native_genome_build"],
            "GRCh38_10x_2020-A",
        )
        self.assertIn(
            "GENCODE_v32",
            self.datasets["gse289173"]["native_annotation_release"],
        )
        self.assertIn(
            "Cell_Ranger_3.0.2",
            self.datasets["gse244832"]["native_annotation_release"],
        )
        self.assertIn(
            "sgRNA_pseudogenes",
            self.datasets["gse281160"]["native_annotation_release"],
        )
        self.assertIn(
            "cellranger-atac_2.1.0",
            self.datasets["gse281367"]["native_genome_build"],
        )
        self.assertEqual(
            self.datasets["gse192741"]["native_genome_build"],
            "hg19_Space_Ranger_1.0",
        )
        self.assertIn(
            "Cell_Ranger_3.1.0",
            self.datasets["gse256398"]["native_genome_build"],
        )
        self.assertIn(
            "STAR_2.5.2b",
            self.datasets["gse267031"]["native_genome_build"],
        )

        for filename, family in self.models_by_file.items():
            self.assertEqual(
                {"schema_version", "family_id", "modality", "status"}
                - family.keys(),
                set(),
                filename,
            )
            self.assertEqual(family["schema_version"], "masld-bench-model-v1")
            for model in family.get("models", []):
                self.assertEqual(MODEL_REQUIRED - model.keys(), set(), model["model_id"])

        for filename, task in self.tasks_by_file.items():
            self.assertEqual(TASK_REQUIRED - task.keys(), set(), filename)
            self.assertEqual(task["schema_version"], "masld-bench-task-v1")
            self.assertEqual(task["task_id"], filename)
            self.assertTrue(set(task["required_pairing_levels"]) <= PAIRING_LEVELS)
            self.assertEqual(task["missingness_policy"], "explicit_state_and_mask")

    def test_model_roster_is_complete_unique_and_fail_closed(self) -> None:
        self.assertEqual(set(self.models), MODEL_IDS)
        sha256 = re.compile(r"^[0-9a-f]{64}$")
        for model_id, model in self.models.items():
            checkpoint = model["checkpoint_revision"]
            digest = model["checkpoint_sha256"]
            self.assertTrue(checkpoint == "UNRESOLVED" or checkpoint.strip())
            self.assertTrue(
                digest in {"UNRESOLVED", "NOT_APPLICABLE"} or sha256.fullmatch(digest)
            )
            if any(
                "UNRESOLVED" in value.upper()
                for value in (
                    checkpoint,
                    digest,
                    model["license_status"],
                    model["exposure_status"],
                )
            ):
                self.assertIs(model["admission_blocking"], True, model_id)
                self.assertTrue(model["blockers"], model_id)
            self.assertTrue(set(model["supported_tasks"]) <= TASK_IDS)
            if not model["supported_tasks"]:
                self.assertEqual(model["status"], "deferred", model_id)

        self.assertEqual(self.models["epibrain"]["license_status"], "blocked_no_license")
        self.assertEqual(self.models["epibrain"]["status"], "blocked_terms")
        self.assertEqual(self.models["corgi_plus"]["status"], "deferred")
        self.assertIn(
            "Zenodo:18630048:corgiplus_model.pt",
            self.models["corgi_plus"]["checkpoint_revision"],
        )
        self.assertIn(
            "CC-BY-4.0",
            self.models["corgi_plus"]["license_status"],
        )
        for model_id in {
            "alphagenome",
            "get",
            "decima",
            "scfoundation_heads",
            "cellfm",
            "genecompass",
            "gb_cell",
            "sclong",
            "scooby_neurips",
            "cellplm_85m_vae_20231027",
        }:
            self.assertEqual(self.models[model_id]["status"], "restricted_comparator")
            self.assertIs(self.models[model_id]["admission_blocking"], True)
        self.assertIn(
            "google/alphagenome-all-folds@a8f293a76ee73d5b57f3bf2ae146510589fcf187",
            self.models["alphagenome"]["checkpoint_revision"],
        )
        self.assertIn(
            "CC-BY-NC-4.0", self.models["decima"]["license_status"]
        )
        self.assertIn(
            "Zenodo:15092691", self.models["decima"]["checkpoint_revision"]
        )
        self.assertIn(
            "GET-Foundation/get_model@fb05735325350d83c104b5020e50b8ae8ebf13cb",
            self.models["get"]["checkpoint_revision"],
        )
        self.assertIn(
            "nontransferable", self.models["scfoundation_heads"]["license_status"]
        )
        self.assertEqual(
            self.models["scfoundation_heads"]["supported_tasks"],
            ["cell_state_mapping"],
        )

        restricted_cell_contracts = {
            "cellfm": (
                "a3ab63e580ebf779edd5d13b4b9c7cfc5795c24495f61c45e2e9f0f1eb7bfc60",
                "biomed-AI/CellFM@bfed59c0e34103231165d69b97927ecc888d623c",
                "CC-BY-NC-ND-4.0",
            ),
            "genecompass": (
                "UNRESOLVED",
                "xCompass-AI/GeneCompass@66ab82f90fa94f88613df7d89a2fb38165e2a82f",
                "NO_LICENSE_DETECTED",
            ),
            "gb_cell": (
                "4602f09951c5a2f8cb814363f68cd63d8d7bd145e02560b852f6652df78bec77",
                "genbio-ai/GB.Cell-100M@a983b388e1025860f2c413ecef60e76078253c24",
                "output_training_restricted",
            ),
            "sclong": (
                "UNRESOLVED",
                "BaiDing1234/scLong@41b72021540918e4386c4f2264351d7a6cbeed99",
                "code_Zenodo_CC-BY-4.0_GitHub_root_unlicensed_weights_outputs_derivatives_UNRESOLVED",
            ),
        }
        for model_id, (digest, revision, license_fragment) in restricted_cell_contracts.items():
            model = self.models[model_id]
            self.assertEqual(model["checkpoint_sha256"], digest, model_id)
            self.assertIn(revision, model["checkpoint_revision"], model_id)
            self.assertIn(license_fragment, model["license_status"], model_id)
            self.assertEqual(model["supported_tasks"], ["cell_state_mapping"], model_id)

        arc_state = self.models["arc_state"]
        self.assertIn(
            "arcinstitute/ST-SE-Replogle@e324967ff4cea5ec199e29bcbb5c1f00e5b9d69c",
            arc_state["checkpoint_revision"],
        )
        self.assertEqual(
            arc_state["checkpoint_sha256"],
            "9e39d5c8e8176eb80bdcab958dec9dd0bbe24aef6156fc2a5319d44d14104654",
        )
        self.assertIn("Non-Commercial", arc_state["license_status"])
        self.assertIn("wholly unseen target", arc_state["blockers"][0])
        self.assertEqual(arc_state["status"], "deferred")
        self.assertEqual(arc_state["supported_tasks"], [])

        scprint2 = self.models["scprint2_small_v2"]
        self.assertEqual(scprint2["display_name"], "scPRINT-2 small-v2")
        self.assertEqual(
            scprint2["checkpoint_revision"],
            "cantinilab/scPRINT-2@7ee2de5aaa3492f326f540f703d44ddeaf1f2fd4:1.0.3+jkobject/scPRINT@d02d4aa684a8422d8bfdc1dbb521a85c773a2772:18hebyht-final-small.ckpt",
        )
        self.assertEqual(
            scprint2["checkpoint_sha256"],
            "2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822",
        )
        self.assertEqual(scprint2["exposure_status"], "unknown")
        self.assertIs(scprint2["admission_blocking"], True)
        scprint2_medium = self.models["scprint2_medium"]
        self.assertEqual(scprint2_medium["checkpoint_revision"], "UNRESOLVED")
        self.assertEqual(scprint2_medium["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(scprint2_medium["status"], "deferred")
        self.assertIs(scprint2_medium["admission_blocking"], True)

        regformer = self.models["regformer"]
        self.assertIn(
            "BGIResearch/RegFormer@9618f1f4bb7b61c9b236d3ea8a6339c236f2da1b",
            regformer["checkpoint_revision"],
        )
        self.assertIn(
            "Figshare:28645493v2:file:57423100",
            regformer["checkpoint_revision"],
        )
        self.assertEqual(regformer["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(regformer["exposure_status"], "unknown")
        self.assertIs(regformer["admission_blocking"], True)

        exact_public_sources = {
            "scimilarity_v1_1": (
                "Zenodo:10685499:model_v1.1.tar.gz",
                "target_label_unexposed",
            ),
            "transcriptformer_tf_sapiens": (
                "s3-version:P5N5hx41EeWC7I0l8oW7yDyiJHS5InPU",
                "unknown",
            ),
        }
        for model_id, (source_id, exposure_status) in exact_public_sources.items():
            model = self.models[model_id]
            self.assertIn(source_id, model["checkpoint_revision"], model_id)
            self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED", model_id)
            self.assertEqual(
                model["exposure_status"], exposure_status, model_id
            )
            self.assertIs(model["admission_blocking"], True, model_id)

        scgpt_sources = {
            "scgpt_whole_human": "gdrive-folder:1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
            "scgpt_continual": "gdrive-folder:1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB",
        }
        for model_id, source_id in scgpt_sources.items():
            model = self.models[model_id]
            self.assertIn(source_id, model["checkpoint_revision"], model_id)
            self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED", model_id)
            self.assertEqual(
                model["exposure_status"], "target_label_unexposed", model_id
            )
            self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
            self.assertIs(model["admission_blocking"], True, model_id)

        for model in self.models_by_file["cell_foundation"].get("models", []):
            self.assertEqual(
                model["supported_tasks"],
                ["cell_state_mapping"],
                model["model_id"],
            )

        scprint_v1 = self.models["scprint_v1_5_medium"]
        self.assertEqual(
            scprint_v1["checkpoint_sha256"],
            "a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5",
        )
        for model_id in {"cellplm_85m", "scbert", "sccello", "langcell"}:
            self.assertEqual(self.models[model_id]["status"], "blocked_terms")
            self.assertIs(self.models[model_id]["admission_blocking"], True)

        cellplm = self.models["cellplm_85m"]
        cellplm_vae = self.models["cellplm_85m_vae_20231027"]
        self.assertIn("20230926_85M.best.ckpt", cellplm["checkpoint_revision"])
        self.assertIn("016213c1be94746000000029e76cf71", cellplm["checkpoint_revision"])
        self.assertEqual(cellplm["exposure_status"], "target_label_unexposed")
        self.assertIn("GMVAE", cellplm["display_name"])
        self.assertIn("20231027_85M.best.ckpt", cellplm_vae["checkpoint_revision"])
        self.assertIn("016213c1be94747000000029e76cf71", cellplm_vae["checkpoint_revision"])
        self.assertEqual(cellplm_vae["status"], "restricted_comparator")
        self.assertIn("VAE", cellplm_vae["display_name"])

        self.assertEqual(self.models["sei"]["status"], "restricted_comparator")
        self.assertIn("academic_research_only", self.models["sei"]["license_status"])
        self.assertEqual(
            self.models["sei"]["exposure_status"], "target_label_unexposed"
        )
        self.assertEqual(
            self.models["bpnet"]["checkpoint_revision"],
            "kundajelab/bpnet@f4593eedca51741f25b8c5cc0c0d647faf8c8a3c+PyPI:bpnet==2.0.0",
        )
        self.assertEqual(self.models["bpnet"]["checkpoint_sha256"], "NOT_APPLICABLE")

        # Exposure is per model and must equal that model's own audit
        # classification.  ABC ships no learned checkpoint and is consumed as a
        # reference annotation, and ENCODE-rE2G declares a clean pretraining
        # corpus, so neither is "target_label_unexposed".
        local_regulatory = {
            "mpralegnet": (
                "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39",
                "pretrained_adapter",
                "candidate",
                "target_label_unexposed",
            ),
            "deltasvm": (
                "NOT_APPLICABLE",
                "train_locally",
                "candidate",
                "target_label_unexposed",
            ),
            "gkmsvm": (
                "NOT_APPLICABLE",
                "train_locally",
                "candidate",
                "target_label_unexposed",
            ),
            "abc": (
                "NOT_APPLICABLE",
                "classical_baseline",
                "candidate",
                "reference_only",
            ),
            "re2g": (
                "UNRESOLVED",
                "pretrained_adapter",
                "candidate",
                "clean_declared",
            ),
            "epinformer": (
                "UNRESOLVED",
                "pretrained_adapter",
                "blocked_terms",
                "target_label_unexposed",
            ),
        }
        for model_id, (
            digest,
            implementation,
            status,
            exposure,
        ) in local_regulatory.items():
            model = self.models[model_id]
            self.assertEqual(model["checkpoint_sha256"], digest)
            self.assertEqual(model["implementation_type"], implementation)
            self.assertEqual(model["status"], status)
            self.assertEqual(model["exposure_status"], exposure, model_id)
        self.assertIn(
            "NO_LICENSE_DETECTED", self.models["epinformer"]["license_status"]
        )

        genomic_language_hashes = {
            "dnabert2": "7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba",
            "nucleotide_transformer": (
                "971cb721bd8cc8134d3665b5d94b195fffb18a1480ad5925d6dcc9290bf1c17f"
            ),
            "hyenadna": "deafb53209bfafb314d14bd81108546a34d473c07ed7ab7a355674376b56ad12",
            "caduceus": "a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2",
            "evo": "c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4",
        }
        for model_id, digest in genomic_language_hashes.items():
            self.assertEqual(self.models[model_id]["checkpoint_sha256"], digest)
            self.assertEqual(
                self.models[model_id]["exposure_status"], "target_label_unexposed"
            )
        self.assertEqual(
            self.models["nucleotide_transformer"]["status"],
            "restricted_comparator",
        )
        self.assertEqual(self.models["enformer"]["status"], "restricted_comparator")
        self.assertIs(self.models["enformer"]["admission_blocking"], True)
        self.assertEqual(self.models["crested"]["status"], "deferred")
        self.assertEqual(self.models["crested"]["supported_tasks"], [])
        self.assertIn("500-bp mm10", self.models["crested"]["blockers"][0])
        self.assertIs(self.models["crested"]["admission_blocking"], True)
        self.assertEqual(
            self.models["crested"]["checkpoint_sha256"],
            "8a3ddcfa29effa9e979eae769fbfc6a362de90197434a630015dccca9fecb34c",
        )
        chromdragonn = self.models["chromdragonn"]
        self.assertEqual(chromdragonn["implementation_type"], "train_locally")
        self.assertEqual(chromdragonn["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertEqual(
            chromdragonn["supported_tasks"],
            ["variant_to_regulation"],
        )
        self.assertNotIn(
            "chromdragonn",
            self.atac_capabilities["conditional_profile_candidate_models"],
        )
        epcotv2 = self.models["epcotv2"]
        self.assertEqual(
            epcotv2["checkpoint_sha256"],
            "b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d",
        )
        self.assertEqual(epcotv2["status"], "blocked_terms")
        self.assertIn("weights_UNDECLARED", epcotv2["license_status"])
        self.assertEqual(epcotv2["supported_tasks"], ["variant_to_regulation"])
        context_borzoi = self.models["context_borzoi"]
        self.assertEqual(context_borzoi["implementation_type"], "train_locally")
        self.assertEqual(context_borzoi["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertEqual(context_borzoi["status"], "deferred")
        self.assertIn(
            "context_borzoi",
            self.atac_capabilities["conditional_profile_candidate_models"],
        )
        self.assertNotIn(
            "scooby_neurips",
            self.atac_capabilities["conditional_profile_candidate_models"],
        )

        perturbation_boundaries = {
            "gears": ("deferred", ["gene_perturbation_response"], "MIT"),
            "scgen": (
                "deferred",
                ["gene_perturbation_response"],
                "repository_GPL-3.0_package_metadata_MIT_conflict",
            ),
            "cpa": (
                "deferred",
                ["gene_perturbation_response"],
                "repository_BSD-3-Clause_package_metadata_MIT_conflict",
            ),
            "cellot": ("deferred", ["gene_perturbation_response"], "BSD-3-Clause"),
            "genepert": (
                "deferred",
                ["gene_perturbation_response"],
                "code_NO_LICENSE_DETECTED_embedding_terms_UNRESOLVED",
            ),
            "scgpt_perturbation_head": (
                "deferred",
                ["gene_perturbation_response"],
                "code_MIT_weights_UNDECLARED",
            ),
            "perturblib_lpm": (
                "deferred",
                ["gene_perturbation_response"],
                "Apache-2.0",
            ),
        }
        for model_id, (status, tasks, license_status) in perturbation_boundaries.items():
            model = self.models[model_id]
            self.assertEqual(model["status"], status)
            self.assertEqual(model["supported_tasks"], tasks)
            self.assertEqual(model["license_status"], license_status)
        for model_id in {"scgen", "cpa", "cellot"}:
            self.assertIn("wholly unseen", " ".join(self.models[model_id]["blockers"]))
        self.assertIn("noncoding", " ".join(self.models["gears"]["blockers"]))
        systema = self.models["systema"]
        self.assertEqual(systema["implementation_type"], "classical_baseline")
        self.assertIn("independent_implementation", systema["license_status"])

        pyg_graph_models = {"node2vec", "graphsage", "rgcn", "hgt"}
        for model_id in pyg_graph_models:
            model = self.models[model_id]
            self.assertIn(
                "pyg-team/pytorch_geometric@76ff9c2ce18c8cebf52122b57e2aeadce9793d10",
                model["checkpoint_revision"],
            )
            self.assertEqual(model["checkpoint_sha256"], "NOT_APPLICABLE")
            self.assertEqual(model["supported_tasks"], ["typed_evidence_graph"])
        self.assertEqual(
            self.models["bionic"]["checkpoint_revision"],
            "bowang-lab/BIONIC@b74e117602a9573a72edd1f53400243d567ca373:v0.2.6",
        )
        self.assertEqual(self.models["bionic"]["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertEqual(
            self.models["typed_transformer"]["supported_tasks"],
            ["typed_evidence_graph"],
        )
        self.assertEqual(
            self.models["late_fusion"]["supported_tasks"],
            ["typed_evidence_graph", "unified_model"],
        )
        for model_id in {
            "pca_procrustes",
            "cca",
            "evidence_count",
            "regularized_linear_evidence",
        }:
            self.assertNotIn("unified_model", self.models[model_id]["supported_tasks"])
        for model_id in {
            "novae",
            "nichecompass",
            "spatialmeta",
            "spamosaic",
            "nicheformer",
            "graphst",
            "stagate",
            "banksy",
            "spatial_transcript_native_baselines",
        }:
            self.assertEqual(self.models[model_id]["supported_tasks"], [])
        harmony = self.models["harmony_logistic"]
        self.assertEqual(harmony["status"], "deferred")
        self.assertEqual(harmony["supported_tasks"], [])
        self.assertEqual(harmony["checkpoint_sha256"], "NOT_APPLICABLE")
        for model_id in {"scvi_baseline", "scanvi_baseline"}:
            self.assertEqual(
                self.models[model_id]["supported_tasks"], ["cell_state_mapping"]
            )
        celltypist = self.models["celltypist_baseline"]
        self.assertEqual(celltypist["implementation_type"], "train_locally")
        self.assertEqual(celltypist["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertEqual(celltypist["exposure_status"], "target_label_unexposed")
        self.assertNotIn(
            "harmony_logistic",
            self.tasks["cell_state_mapping"]["baseline_model_ids"],
        )
        self.assertNotIn(
            "assay_native_pseudobulk",
            self.tasks["cell_state_mapping"]["baseline_model_ids"],
        )
        self.assertEqual(
            self.models["assay_native_pseudobulk"]["supported_tasks"],
            ["rna_conditioned_atac", "bulk_state_transfer"],
        )

    def test_audited_chromatin_and_integration_task_boundaries_are_frozen(self) -> None:
        expected = {
            "scbasset": ("candidate", ["rna_conditioned_atac"]),
            "chrombert": ("blocked_terms", ["variant_to_regulation"]),
            "cistopic": ("deferred", []),
            "lsi": ("deferred", []),
            "scmvp": ("candidate", ["cell_state_mapping"]),
            "scglue": ("candidate", ["cell_state_mapping"]),
            "scmomat": ("deferred", []),
            "stabmap": ("candidate", ["cell_state_mapping"]),
            "seurat_wnn": ("deferred", []),
            "cobolt": ("candidate", ["cell_state_mapping", "rna_conditioned_atac"]),
            "mofaplus": ("deferred", []),
            "mira": ("blocked_terms", ["cell_state_mapping"]),
            "scjoint": ("blocked_terms", ["cell_state_mapping"]),
        }
        for model_id, (status, tasks) in expected.items():
            self.assertEqual(self.models[model_id]["status"], status, model_id)
            self.assertEqual(self.models[model_id]["supported_tasks"], tasks, model_id)
            self.assertNotIn("unified_model", tasks, model_id)

        self.assertIn(
            "cobolt", self.atac_capabilities["direct_profile_candidate_models"]
        )
        self.assertIn("scmvp", self.atac_capabilities["latent_only_models"])
        self.assertEqual(
            self.models["mofaplus"]["exposure_status"],
            "target_label_unexposed",
        )

    def test_every_model_audit_bundle_is_complete_and_hash_bound(self) -> None:
        audit_root = CONFIG / "artifacts" / "models"
        required = {
            "checkpoints.json",
            "development_crosswalk.json",
            "exposure_audit.json",
        }
        checkpoint_schemas = {
            "masld-bench-external-checkpoint-preflight-v1",
            "masld-bench-grouped-model-preflight-v1",
            "masld-bench-local-model-preflight-v1",
            "masld-bench-local-model-suite-preflight-v1",
            "masld-bench-upstream-checkpoint-preflight-v1",
        }
        crosswalk_schemas = {
            "masld-bench-development-crosswalk-v1",
            "masld-bench-development-exposure-crosswalk-v1",
            "masld-bench-external-model-crosswalk-v1",
            "masld-bench-grouped-model-crosswalk-v1",
            "masld-bench-local-model-crosswalk-v1",
        }
        bundle_dirs = sorted(path for path in audit_root.iterdir() if path.is_dir())
        self.assertTrue(bundle_dirs)
        for bundle in bundle_dirs:
            self.assertFalse(bundle.is_symlink(), bundle)
            members = {path.name for path in bundle.iterdir() if path.is_file()}
            self.assertTrue(required.issubset(members), bundle.name)
            checkpoint = json.loads(
                (bundle / "checkpoints.json").read_text(encoding="utf-8")
            )
            crosswalk_path = bundle / "development_crosswalk.json"
            crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
            exposure = json.loads(
                (bundle / "exposure_audit.json").read_text(encoding="utf-8")
            )
            self.assertIn(checkpoint["schema_version"], checkpoint_schemas, bundle.name)
            self.assertIn(crosswalk["schema_version"], crosswalk_schemas, bundle.name)
            self.assertEqual(
                exposure["schema_version"],
                "masld-bench-checkpoint-exposure-audit-v1",
                bundle.name,
            )
            self.assertRegex(crosswalk["verified_date"], r"^20\d{2}-\d{2}-\d{2}$")
            self.assertRegex(exposure["verified_date"], r"^20\d{2}-\d{2}-\d{2}$")
            self.assertIn("development_crosswalk_record", exposure, bundle.name)
            record = exposure["development_crosswalk_record"]
            self.assertEqual(record["path"], "development_crosswalk.json", bundle.name)
            self.assertEqual(
                record["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                bundle.name,
            )

    def test_code_policy_constants_agree_across_modules(self) -> None:
        """Policy strings are duplicated across modules; prove they agree.

        selection.py cannot import tournament.py (tournament imports selection),
        so the finalist-authority and fit-state policy strings are mirrored as
        literals. A silent divergence would let a ledger satisfy one module's
        gate and fail the other's, or worse, satisfy neither check meaningfully.
        """

        import masld_bench.conditional_model as conditional_model
        import masld_bench.firewall as firewall
        import masld_bench.selection as selection
        import masld_bench.tournament as tournament

        self.assertEqual(
            tournament.LEDGER_AUTHORITY_FINALIST,
            selection.FINALIST_LEDGER_AUTHORITY,
        )
        self.assertEqual(
            tournament._FINALIST_SEED_COUNT, selection.FINALIST_SEED_COUNT
        )
        self.assertEqual(
            tournament._ENSEMBLE_SEED_COUNT, tournament._FINALIST_SEED_COUNT
        )
        self.assertEqual(
            tournament.DEVELOPMENT_ENSEMBLE_POLICY_ID,
            selection.DEVELOPMENT_ENSEMBLE_POLICY_ID,
        )
        self.assertEqual(
            tournament.DEVELOPMENT_ENSEMBLE_POLICY_ID,
            firewall.DEVELOPMENT_ENSEMBLE_POLICY_ID,
        )
        self.assertEqual(
            tournament._DEVELOPMENT_ENSEMBLE_AGGREGATION,
            firewall.DEVELOPMENT_ENSEMBLE_AGGREGATION,
        )
        self.assertEqual(
            tournament._FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION,
            firewall.FINALIST_ENDPOINT_DISTRIBUTION_SCHEMA_VERSION,
        )
        self.assertEqual(
            tournament.FIT_STATE_ARTIFACT_POLICY,
            selection.FIT_STATE_ARTIFACT_POLICY,
        )
        self.assertEqual(
            conditional_model.STACKING_EVIDENCE_SCHEMA_VERSION,
            "masld-bench-stacking-evidence-v1",
        )

    def test_shortlist_wave_contracts_mirror_their_campaigns(self) -> None:
        """The only waves that may shortlist must match their campaign seeds."""

        import masld_bench.tournament as tournament

        expected = {
            "frozen_screen": "v1_frozen_screen",
            "full_specialist_screen": "v1_specialist_screen",
        }
        self.assertEqual(
            set(tournament.SHORTLIST_SOURCE_WAVE_CONTRACTS), set(expected)
        )
        for wave, campaign in expected.items():
            config = load_toml(CONFIG / "campaigns" / f"{campaign}.toml")
            self.assertEqual(config["wave"], wave)
            self.assertEqual(
                tuple(config["execution"]["seeds"]),
                tuple(
                    tournament.SHORTLIST_SOURCE_WAVE_CONTRACTS[wave][
                        "expected_seeds"
                    ]
                ),
                wave,
            )

    def test_every_repo_path_anchored_audit_hash_rederives(self) -> None:
        """Audit bundles pin source files by SHA-256; those pins must be live.

        Nothing rederived these bindings before, so editing a bound source (for
        example src/masld_bench/conditional_model.py, pinned by
        context_borzoi/checkpoints.json in two places) silently invalidated the
        audit while every test stayed green.
        """

        prefix = "Analysis/MASLD_Model_Benchmark/"
        checked = 0
        mismatches: list[str] = []

        def resolve(identity: str) -> Path | None:
            relative = identity[len(prefix):] if identity.startswith(prefix) else identity
            candidate = ROOT / relative
            return candidate if candidate.is_file() else None

        def walk(node: object, bundle: Path) -> None:
            nonlocal checked
            if isinstance(node, dict):
                digest_key = next(
                    (
                        key
                        for key in ("sha256", "authority_sha256")
                        if isinstance(node.get(key), str)
                        and len(node[key]) == 64
                    ),
                    None,
                )
                if digest_key is not None:
                    for path_key in ("identity", "authority_path", "path"):
                        value = node.get(path_key)
                        if not isinstance(value, str):
                            continue
                        source = resolve(value)
                        if source is None:
                            continue
                        checked += 1
                        actual = hashlib.sha256(source.read_bytes()).hexdigest()
                        if actual != node[digest_key]:
                            mismatches.append(
                                f"{bundle.parent.name}/{bundle.name}: {value} "
                                f"recorded {node[digest_key]} but is {actual}"
                            )
                for value in node.values():
                    walk(value, bundle)
            elif isinstance(node, list):
                for value in node:
                    walk(value, bundle)

        audit_root = CONFIG / "artifacts" / "models"
        for bundle in sorted(audit_root.glob("*/*.json")):
            walk(json.loads(bundle.read_text(encoding="utf-8")), bundle)

        self.assertGreater(checked, 40, "audit bundles stopped pinning sources")
        self.assertEqual(mismatches, [])

    def test_model_census_maps_exactly_to_audit_bundles(self) -> None:
        authority = self.model_audit_bindings
        self.assertEqual(
            authority["schema_version"], "masld-bench-model-audit-bindings-v1"
        )
        self.assertEqual(
            authority["default_rule"],
            "registry_model_id_equals_artifact_bundle_id",
        )
        self.assertEqual(len(self.models), authority["expected_registry_model_count"])

        grouped_by_model: dict[str, str] = {}
        grouped_bundle_ids: set[str] = set()
        for index, binding in enumerate(authority["group_bindings"]):
            self.assertEqual(set(binding), {"bundle_id", "model_ids"}, index)
            bundle_id = binding["bundle_id"]
            model_ids = binding["model_ids"]
            self.assertIsInstance(bundle_id, str, index)
            self.assertTrue(bundle_id, index)
            self.assertNotIn(bundle_id, grouped_bundle_ids, bundle_id)
            grouped_bundle_ids.add(bundle_id)
            self.assertTrue(model_ids, bundle_id)
            self.assertEqual(len(model_ids), len(set(model_ids)), bundle_id)
            for model_id in model_ids:
                self.assertIn(model_id, self.models, bundle_id)
                self.assertNotIn(model_id, grouped_by_model, model_id)
                grouped_by_model[model_id] = bundle_id

        resolved = {
            model_id: grouped_by_model.get(model_id, model_id)
            for model_id in self.models
        }
        bundle_ids = set(resolved.values())
        self.assertEqual(
            len(bundle_ids), authority["expected_artifact_bundle_count"]
        )
        audit_root = CONFIG / "artifacts" / "models"
        observed_bundle_ids = {
            path.name for path in audit_root.iterdir() if path.is_dir()
        }
        self.assertEqual(observed_bundle_ids, bundle_ids)
        required = set(authority["required_bundle_files"])
        self.assertEqual(
            required,
            {"checkpoints.json", "development_crosswalk.json", "exposure_audit.json"},
        )
        for bundle_id in sorted(bundle_ids):
            members = {
                path.name
                for path in (audit_root / bundle_id).iterdir()
                if path.is_file()
            }
            self.assertTrue(required.issubset(members), bundle_id)

    def test_restricted_cell_foundation_audits_fail_closed(self) -> None:
        cases = {
            "scfoundation": (
                "scfoundation_heads",
                "registered_scfoundation_base_and_heads",
            ),
            "cellfm": ("cellfm", "cellfm_base"),
            "genecompass": ("genecompass", "genecompass_base"),
        }
        loaded = {}
        for artifact_id, (registry_id, finding_id) in cases.items():
            audit_dir = CONFIG / "artifacts" / "models" / artifact_id
            checkpoints = json.loads(
                (audit_dir / "checkpoints.json").read_text(encoding="utf-8")
            )
            crosswalk_path = audit_dir / "development_crosswalk.json"
            crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
            exposure = json.loads(
                (audit_dir / "exposure_audit.json").read_text(encoding="utf-8")
            )
            loaded[artifact_id] = (checkpoints, crosswalk, exposure)

            self.assertEqual(
                checkpoints["registry_binding"]["artifact_bundle_id"], artifact_id
            )
            self.assertEqual(
                checkpoints["registry_binding"]["registry_model_id"], registry_id
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                {
                    finding["exposure_state"]
                    for finding in crosswalk["findings"].values()
                },
                {"unknown"},
            )
            self.assertEqual(
                crosswalk["findings"]["GSE289173"]["exposure_state"], "unknown"
            )
            self.assertEqual(
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
                "unknown",
            )
            self.assertIn("ineligible", exposure["sealed_champion_eligibility"])
            self.assertEqual(
                exposure["task_disposition"]["released_open_champion"],
                "ineligible",
            )

            model = self.models[registry_id]
            self.assertEqual(model["exposure_status"], "unknown")
            self.assertEqual(model["status"], "restricted_comparator")
            self.assertIs(model["admission_blocking"], True)
            self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])

        scfoundation = loaded["scfoundation"][0]
        self.assertEqual(
            scfoundation["registry_binding"]["registry_model_id"],
            "scfoundation_heads",
        )
        self.assertEqual(scfoundation["terms_audit"]["code_license"], "Apache-2.0")
        self.assertIs(scfoundation["terms_audit"]["open_champion_eligible"], False)
        self.assertEqual(
            scfoundation["checkpoint_contract"]["registered_sha256"], "UNRESOLVED"
        )

        cellfm = loaded["cellfm"][0]
        self.assertEqual(cellfm["checkpoint_contract"]["size_bytes"], 9659358467)
        self.assertEqual(
            cellfm["checkpoint_contract"]["sha256"],
            "a3ab63e580ebf779edd5d13b4b9c7cfc5795c24495f61c45e2e9f0f1eb7bfc60",
        )
        self.assertEqual(cellfm["terms_audit"]["weight_terms"], "CC-BY-NC-ND-4.0")
        self.assertIs(cellfm["terms_audit"]["open_champion_eligible"], False)

        genecompass, genecompass_crosswalk, _ = loaded["genecompass"]
        self.assertEqual(
            genecompass["checkpoint_contract"]["registered_sha256"], "UNRESOLVED"
        )
        self.assertIs(
            genecompass["terms_audit"]["repository_license_detected"], False
        )
        self.assertEqual(
            genecompass_crosswalk["corpus"]["declared_sum_from_component_counts"],
            101768420,
        )

    def test_graph_and_transparent_baseline_audits_freeze_exact_boundaries(
        self,
    ) -> None:
        expected_registry = {
            "node2vec": (
                "pyg-team/pytorch_geometric@76ff9c2ce18c8cebf52122b57e2aeadce9793d10:torch_geometric==2.7.0:Node2Vec",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
            "bionic": (
                "bowang-lab/BIONIC@b74e117602a9573a72edd1f53400243d567ca373:v0.2.6",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
            "graphsage": (
                "pyg-team/pytorch_geometric@76ff9c2ce18c8cebf52122b57e2aeadce9793d10:torch_geometric==2.7.0:SAGEConv",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
            "rgcn": (
                "pyg-team/pytorch_geometric@76ff9c2ce18c8cebf52122b57e2aeadce9793d10:torch_geometric==2.7.0:RGCNConv",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
            "hgt": (
                "pyg-team/pytorch_geometric@76ff9c2ce18c8cebf52122b57e2aeadce9793d10:torch_geometric==2.7.0:HGTConv",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
            "typed_transformer": ("UNRESOLVED", "UNRESOLVED", "unknown"),
            "late_fusion": (
                "masld_equal_weight_calibrated_late_fusion_v1",
                "NOT_APPLICABLE",
                "unknown",
            ),
            "evidence_count": (
                "masld_evidence_count_v1",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
            "nearest_gene": (
                "masld_nearest_gene_tss_v1:GENCODE_v49:GRCh38.p14",
                "NOT_APPLICABLE",
                "target_label_unexposed",
            ),
        }
        for model_id, (revision, digest, exposure) in expected_registry.items():
            model = self.models[model_id]
            self.assertEqual(model["checkpoint_revision"], revision, model_id)
            self.assertEqual(model["checkpoint_sha256"], digest, model_id)
            self.assertEqual(model["exposure_status"], exposure, model_id)
            self.assertEqual(model["status"], "candidate", model_id)
            self.assertIs(model["admission_blocking"], True, model_id)
            self.assertTrue(model["blockers"], model_id)

        audit_root = CONFIG / "artifacts" / "models"
        audit_ids = tuple(expected_registry)
        audits = {}
        for model_id in audit_ids:
            bundle = audit_root / model_id
            checkpoint = json.loads(
                (bundle / "checkpoints.json").read_text(encoding="utf-8")
            )
            crosswalk_path = bundle / "development_crosswalk.json"
            exposure = json.loads(
                (bundle / "exposure_audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )
            audits[model_id] = checkpoint

        self.assertEqual(
            audits["node2vec"]["implementation"]["upstream_checkpoint"],
            "NOT_APPLICABLE_NO_UPSTREAM_WEIGHTS",
        )
        self.assertEqual(
            audits["node2vec"]["native_assumptions_and_limits"]["closed_world"],
            "The Embedding parameter has one row per training node. A node absent from the retained graph has no learned representation; adding it after fitting is impossible without a new model.",
        )
        self.assertEqual(
            audits["bionic"]["weight_license"],
            "NOT_APPLICABLE_NO_UPSTREAM_WEIGHTS",
        )
        self.assertEqual(
            audits["graphsage"]["implementation"]["upstream_checkpoint"],
            "NOT_APPLICABLE_NO_UPSTREAM_WEIGHTS",
        )
        self.assertIn(
            "root-feature embedding",
            audits["graphsage"]["feature_contract"]["inductive_query"],
        )
        self.assertIn(
            "cannot be used to claim recovery",
            audits["rgcn"]["held_source_output_contract"][
                "held_relation_parameter_rule"
            ],
        )
        self.assertIn(
            "cannot be used to score",
            audits["hgt"]["held_source_output_contract"]["held_edge_type_rule"],
        )
        self.assertEqual(
            audits["typed_transformer"]["architecture_contract"][
                "architecture_status"
            ],
            "UNRESOLVED_BLOCKING",
        )
        self.assertEqual(
            audits["late_fusion"]["algorithm_contract"]["fitting"],
            "No fusion weight is fitted. Equal component weights are fixed a priori; task-specific component calibration occurs upstream inside development training folds.",
        )
        self.assertIn(
            "distinct remaining-source evidence units",
            audits["evidence_count"]["algorithm_contract"]["definition"],
        )
        self.assertEqual(
            audits["nearest_gene"]["output_contract"]["allele_effect"],
            "none; REF and ALT at one anchor receive the same proximity score",
        )
        self.assertIn(
            "signed_eQTL_effect",
            audits["nearest_gene"]["output_contract"]["not_supported"],
        )

        source_bound_registries = {
            "typed_transformer": "evidence_graph.toml",
            "late_fusion": "evidence_graph.toml",
            "evidence_count": "mandatory_baselines.toml",
            "nearest_gene": "mandatory_baselines.toml",
        }
        for model_id, registry_name in source_bound_registries.items():
            registry_path = CONFIG / "models" / registry_name
            identity = (
                "Analysis/MASLD_Model_Benchmark/config/models/" + registry_name
            )
            record = next(
                item
                for item in audits[model_id]["source_records"]
                if item["identity"] == identity
            )
            self.assertEqual(
                record["sha256"],
                hashlib.sha256(registry_path.read_bytes()).hexdigest(),
                model_id,
            )
            self.assertEqual(record["size_bytes"], registry_path.stat().st_size)

    def test_new_direct_translation_audits_fail_closed(self) -> None:
        scpair = self.models["scpair"]
        self.assertEqual(
            scpair["supported_tasks"],
            ["cell_state_mapping", "rna_conditioned_atac"],
        )
        self.assertIn(
            "scpair", self.atac_capabilities["direct_profile_candidate_models"]
        )
        scpair_preflight = json.loads(
            (CONFIG / "artifacts/models/scpair/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            scpair_preflight["output_contract"]["primary_atac_output"].split(",")[0],
            "Use decoder output_px",
        )

        jamie = self.models["jamie"]
        self.assertEqual(jamie["license_status"], "GPL-3.0-only")
        self.assertEqual(jamie["supported_tasks"], ["cell_state_mapping"])
        self.assertNotIn(
            "jamie", self.atac_capabilities["direct_profile_candidate_models"]
        )
        jamie_preflight = json.loads(
            (CONFIG / "artifacts/models/jamie/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn(
            "unbounded continuous vector",
            jamie_preflight["output_contract"]["profile_limitation"],
        )

        monae = self.models["monae"]
        self.assertEqual(monae["status"], "deferred")
        self.assertEqual(monae["supported_tasks"], [])
        self.assertNotIn(
            "monae", self.atac_capabilities["direct_profile_candidate_models"]
        )
        monae_preflight = json.loads(
            (CONFIG / "artifacts/models/monae/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            monae_preflight["architecture_contract"]["release_fidelity_status"],
            "BLOCKED_PAPER_CODE_MISMATCH",
        )

        multidgd = self.models["multidgd"]
        self.assertEqual(multidgd["status"], "deferred")
        self.assertEqual(multidgd["supported_tasks"], [])
        self.assertNotIn(
            "multidgd", self.atac_capabilities["direct_profile_candidate_models"]
        )
        multidgd_preflight = json.loads(
            (CONFIG / "artifacts/models/multidgd/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIs(
            multidgd_preflight["inductive_inference_contract"][
                "frozen_raw_rna_to_atac_path"
            ],
            False,
        )

        scdiffusionx = self.models["scdiffusionx"]
        self.assertEqual(
            scdiffusionx["license_status"],
            "root_MIT_nested_package_BSD-3-Clause",
        )
        self.assertIn(
            "scdiffusionx",
            self.atac_capabilities["direct_profile_candidate_models"],
        )
        scdiffusionx_preflight = json.loads(
            (CONFIG / "artifacts/models/scdiffusionx/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIs(
            scdiffusionx_preflight["inductive_inference_contract"][
                "frozen_raw_rna_to_atac_path"
            ],
            True,
        )
        self.assertEqual(
            scdiffusionx_preflight["inductive_inference_contract"]["sampler_rule"].split(".")[0],
            "Use conditional DDPM only",
        )

    def test_regulatory_component_audits_preserve_endpoint_boundaries(self) -> None:
        mpralegnet = json.loads(
            (CONFIG / "artifacts/models/mpralegnet/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            mpralegnet["checkpoint"]["lfs_oid_sha256"],
            self.models["mpralegnet"]["checkpoint_sha256"],
        )
        self.assertIn(
            "ineligible",
            mpralegnet["representation_and_task_contract"][
                "standalone_primary_endpoint_status"
            ],
        )

        abc = json.loads(
            (CONFIG / "artifacts/models/abc/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            abc["implementation"]["repository_revision"],
            "d5976ee642edf94de54b1d694c29477983f88df2",
        )
        self.assertIn(
            "signed_expression_effect",
            abc["output_and_task_contract"]["not_supported"],
        )
        self.assertEqual(
            abc["output_and_task_contract"]["variant_overlap_helper"].split(";")[0],
            "prohibited",
        )

        re2g = json.loads(
            (CONFIG / "artifacts/models/re2g/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            re2g["checkpoint"]["feature_table_sha256"],
            "fc4d831be3692162bc45dc20558cd3cc0e3546d1acddadbc6003f31efbfb7d5c",
        )
        self.assertIn("atac_megamap", self.models["re2g"]["checkpoint_revision"])
        self.assertEqual(re2g["checkpoint"]["sha256"], None)

        deltasvm = json.loads(
            (CONFIG / "artifacts/models/deltasvm/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        gkmsvm = json.loads(
            (CONFIG / "artifacts/models/gkmsvm/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn(
            "same artifacts",
            gkmsvm["algorithm_contract"]["model_family_relationship"],
        )
        self.assertIn(
            "not an independently trained model family",
            deltasvm["algorithm_contract"]["model_family_relationship"],
        )
        self.assertIn(
            "ALT_minus_REF",
            deltasvm["score_contract"]["canonical_alt_minus_ref"],
        )
        self.assertIn(
            "ALT_minus_REF",
            gkmsvm["output_and_task_contract"]["project_variant_delta"],
        )

        epinformer = json.loads(
            (CONFIG / "artifacts/models/epinformer/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(epinformer["checkpoint_set"]["member_count"], 24)
        self.assertEqual(len(epinformer["enhancer_encoder_checkpoints"]), 12)
        self.assertEqual(len(epinformer["expression_checkpoints"]), 12)
        self.assertEqual(
            epinformer["output_and_task_contract"]["task_role"].split()[0],
            "link_only",
        )
        self.assertEqual(self.models["epinformer"]["status"], "blocked_terms")

        gates = " ".join(
            self.tasks["variant_to_regulation"]["admission_gates"]
        )
        self.assertIn(
            "deltaSVM and direct gkm-SVM consume the same", gates
        )
        self.assertIn("one fitted-model family", gates)

        for model_id in ("deltasvm", "gkmsvm", "epinformer"):
            crosswalk_path = (
                CONFIG
                / f"artifacts/models/{model_id}/development_crosswalk.json"
            )
            exposure = json.loads(
                (
                    CONFIG
                    / f"artifacts/models/{model_id}/exposure_audit.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )

    def test_restricted_sequence_context_audits_fail_closed(self) -> None:
        alphagenome = json.loads(
            (CONFIG / "artifacts/models/alphagenome/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn(
            "cholangiocyte",
            alphagenome["output_contract"]["sealed_lineage_coverage"],
        )
        self.assertIn(
            "no_matching_RNA_SEQ_track",
            alphagenome["output_contract"]["sealed_lineage_coverage"][
                "cholangiocyte"
            ],
        )
        self.assertIn(
            "no_hepatic_stellate_RNA_SEQ_track",
            alphagenome["output_contract"]["sealed_lineage_coverage"][
                "stellate_cell"
            ],
        )

        get_model = json.loads(
            (CONFIG / "artifacts/models/get/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIs(get_model["input_contract"]["observed_atac_required"], True)
        self.assertIn(
            "still requires an observed cell-context accessible-region universe",
            get_model["input_contract"]["binary_atac_mode"],
        )

        decima = json.loads(
            (CONFIG / "artifacts/models/decima/development_crosswalk.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            decima["sealed_output_roster"]["roster_status"],
            "UNRESOLVED_BLOCKING",
        )
        self.assertEqual(
            decima["findings"]["current_coloc_gwas"]["resolved_alias"],
            "QTS000038=OneK1K=GSE196830",
        )
        self.assertIs(
            decima["findings"]["current_coloc_gwas"]["sealed_overlap"], False
        )

        for model_id in ("alphagenome", "get", "decima"):
            crosswalk_path = (
                CONFIG
                / f"artifacts/models/{model_id}/development_crosswalk.json"
            )
            exposure = json.loads(
                (
                    CONFIG / f"artifacts/models/{model_id}/exposure_audit.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )

    def test_sequence_component_audits_preserve_native_outputs(self) -> None:
        enformer = json.loads(
            (CONFIG / "artifacts/models/enformer/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            len(enformer["native_checkpoint"]["checkpoint_objects"]), 3
        )
        self.assertEqual(
            enformer["output_contract"]["human_track_counts"]["CAGE"], 638
        )
        self.assertEqual(self.models["enformer"]["status"], "restricted_comparator")

        sei = json.loads(
            (CONFIG / "artifacts/models/sei/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            sei["native_checkpoint"]["archive"]["md5"],
            "4297aafb711aec4ecccb645b8928ea26",
        )
        self.assertEqual(
            sei["output_contract"]["sequence_class_manifest"]["rows"], 40
        )

        chromdragonn = json.loads(
            (CONFIG / "artifacts/models/chromdragonn/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIs(
            chromdragonn["output_contract"]["multinomial_profile"], False
        )
        self.assertEqual(
            chromdragonn["task_fit"]["rna_conditioned_atac"].split(";")[0],
            "supported_only_as_secondary_binary_region_accessibility_classification",
        )

        chrombert = json.loads(
            (CONFIG / "artifacts/models/chrombert/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            chrombert["registered_checkpoint"]["sha256"],
            self.models["chrombert"]["checkpoint_sha256"],
        )
        self.assertIn(
            "same native base representation",
            chrombert["input_contract"]["variant_limit"],
        )
        self.assertEqual(
            chrombert["post_census_alternate"]["status"],
            "separate_post_seal_identity_not_admitted",
        )
        variant_records = {
            record["model_id"]: record
            for record in self.variant_capabilities["models"]
        }
        self.assertEqual(
            variant_records["chrombert"]["allowed_endpoints"],
            ["eqtl_retrieval"],
        )

        epibrain = json.loads(
            (CONFIG / "artifacts/models/epibrain/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            epibrain["architecture_contract"][
                "paper_cell_type_embedding_dimension"
            ],
            8,
        )
        self.assertEqual(
            epibrain["architecture_contract"][
                "checkpoint_receipt_cell_type_embedding_dimension"
            ],
            6,
        )
        self.assertEqual(self.models["epibrain"]["exposure_status"], "unknown")
        self.assertEqual(
            epibrain["task_fit"]["variant_to_regulation"],
            "brain_domain_secondary_diagnostic_only",
        )

        crested = json.loads(
            (CONFIG / "artifacts/models/crested/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(crested["architecture_contract"]["native_reference"], "mm10")
        self.assertEqual(crested["architecture_contract"]["output_dim"], 82)
        self.assertEqual(
            crested["registered_converted_checkpoint"][
                "declared_sha256_from_pinned_pooch_registry"
            ],
            self.models["crested"]["checkpoint_sha256"],
        )
        self.assertEqual(self.models["crested"]["supported_tasks"], [])

        for model_id in ("chrombert", "epibrain", "crested"):
            crosswalk_path = (
                CONFIG
                / f"artifacts/models/{model_id}/development_crosswalk.json"
            )
            exposure = json.loads(
                (
                    CONFIG
                    / f"artifacts/models/{model_id}/exposure_audit.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )

    def test_scmvp_mira_and_mofaplus_inductive_boundaries_are_frozen(self) -> None:
        scmvp = json.loads(
            (CONFIG / "artifacts/models/scmvp/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(scmvp["output_contract"]["profile_status"], "latent_only")
        self.assertIn(
            "observed ATAC",
            scmvp["inductive_inference_contract"]["upstream_profile_audit"],
        )

        mira = json.loads(
            (CONFIG / "artifacts/models/mira/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(mira["output_contract"]["profile_status"], "latent_only")
        self.assertEqual(
            mira["terms_audit"]["result"],
            "blocked_pending_complete_license_text_or_written_confirmation",
        )

        mofaplus = json.loads(
            (CONFIG / "artifacts/models/mofaplus/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(mofaplus["output_contract"]["profile_status"], "latent_only")
        self.assertIn(
            "no amortized encoder",
            mofaplus["architecture_contract"]["inductive_limitation"],
        )
        mofaplus_exposure = json.loads(
            (CONFIG / "artifacts/models/mofaplus/exposure_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            mofaplus_exposure["checkpoint_findings"]["mofaplus_local"][
                "exposure_state"
            ],
            "continual_seen",
        )

        for model_id in {"scmvp", "mira", "mofaplus"}:
            crosswalk_path = (
                CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            )
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )

    def test_scooby_context_and_output_boundaries_are_frozen(self) -> None:
        expected = {
            "scooby_onek1k": (10, 2, "onek1k"),
            "scooby_epicardioids": (50, 3, "epicardioids"),
            "scooby_neurips": (14, 3, "neurips"),
        }
        for model_id, (context_dim, tracks, finding_id) in expected.items():
            checkpoint = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/checkpoints.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                checkpoint["architecture_contract"]["cell_decoder"].split()[1],
                f"{context_dim}-dimensional",
                model_id,
            )
            self.assertEqual(
                checkpoint["architecture_contract"]["decoder_tracks"],
                tracks,
                model_id,
            )
            self.assertIn(
                "NOT_BUNDLED",
                checkpoint["input_contract"]["cell_context_encoder_status"],
                model_id,
            )
            crosswalk_path = (
                CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            )
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )
            self.assertEqual(
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
                "target_label_unexposed",
                model_id,
            )

        self.assertEqual(
            self.models["scooby_onek1k"]["supported_tasks"],
            ["variant_to_regulation"],
        )
        self.assertEqual(
            self.models["scooby_epicardioids"]["supported_tasks"],
            ["variant_to_regulation"],
        )
        self.assertEqual(
            self.models["scooby_neurips"]["supported_tasks"],
            ["variant_to_regulation"],
        )

    def test_primary_evaluator_rosters_are_frozen_and_authoritative(self) -> None:
        expected = {
            "cell_state_mapping": (
                "cell_donor_balanced_macro_f1_v1",
                "class_roster",
                [
                    "cholangiocyte",
                    "endothelial",
                    "hepatocyte",
                    "immune",
                    "mesenchymal_stromal",
                ],
            ),
            "variant_to_regulation": (
                "variant_ld_block_fisher_z_spearman_gain_v1",
                "strata",
                [
                    "cholangiocyte",
                    "endothelial_cell",
                    "hepatocyte",
                    "stellate_cell",
                ],
            ),
            "rna_conditioned_atac": (
                "rna_atac_two_way_deviance_reduction_v1",
                "strata",
                [
                    "cholangiocyte",
                    "fibroblast",
                    "hepatocyte",
                    "macrophage",
                    "t_cell",
                ],
            ),
        }
        for task_id, (evaluator_id, roster_field, roster) in expected.items():
            task = self.tasks[task_id]
            self.assertEqual(task["primary_evaluator_id"], evaluator_id)
            parameters = task["evaluator_parameters"]
            self.assertEqual(parameters[roster_field], roster)
            authority = parameters["roster_authority"]
            self.assertEqual(
                authority["role"], f"task_evaluator_roster:{task_id}"
            )
            authority_path = CONFIG / authority["path"]
            self.assertEqual(authority_path.stat().st_size, authority["size_bytes"])
            self.assertEqual(
                hashlib.sha256(authority_path.read_bytes()).hexdigest(),
                authority["sha256"],
            )
            self.assertEqual(
                json.loads(authority_path.read_text(encoding="utf-8")),
                {
                    "schema_version": "masld-bench-evaluator-roster-v1",
                    "task_id": task_id,
                    "primary_evaluator_id": evaluator_id,
                    "roster_field": roster_field,
                    "roster": roster,
                },
            )

    def test_geneformer_and_corgi_upstream_preflights_are_frozen(self) -> None:
        geneformer = json.loads(
            (CONFIG / "artifacts/models/geneformer/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            geneformer["repository_revision"],
            "04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5",
        )
        observed_geneformer = {
            model_id: (
                self.models[model_id]["checkpoint_sha256"],
                self.models[model_id]["checkpoint_revision"],
            )
            for model_id in (
                "geneformer_v1_10m",
                "geneformer_v2_104m",
                "geneformer_v2_316m",
            )
        }
        self.assertEqual(
            [value[0] for value in observed_geneformer.values()],
            [artifact["sha256"] for artifact in geneformer["artifacts"]],
        )
        for _, revision in observed_geneformer.values():
            self.assertIn(geneformer["repository_revision"], revision)
        self.assertEqual(
            geneformer["auxiliary_artifacts"]["V1"][-1]["sha256"],
            "ab9dc40973fa5224d77b793e2fd114cacf3d08423ed9c4c49caf0ba9c7f218f1",
        )
        self.assertEqual(
            geneformer["auxiliary_artifacts"]["V2_shared"][-1]["sha256"],
            "67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f",
        )
        self.assertEqual(
            geneformer["input_contract"],
            {
                "V1": {"model_input_size": 2048, "special_token": False},
                "V2": {"model_input_size": 4096, "special_token": True},
                "required_h5ad_fields": {
                    "obs": ["n_counts"],
                    "var": ["ensembl_id"],
                },
                "required_matrix": "raw_unselected_counts",
            },
        )
        self.assertEqual(
            {
                model_id: contract["max_position_embeddings"]
                for model_id, contract in geneformer[
                    "architecture_contract"
                ].items()
            },
            {
                "geneformer_v1_10m": 2048,
                "geneformer_v2_104m": 4096,
                "geneformer_v2_316m": 4096,
            },
        )
        self.assertEqual(
            geneformer["architecture_contract"]["geneformer_v2_316m"][
                "num_hidden_layers"
            ],
            18,
        )
        exposure = json.loads(
            (CONFIG / "artifacts/models/geneformer/exposure_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            set(exposure["checkpoint_findings"]), set(observed_geneformer)
        )
        for model_id, finding in exposure["checkpoint_findings"].items():
            self.assertEqual(finding["exposure_state"], "target_label_unexposed")
            self.assertEqual(
                self.models[model_id]["supported_tasks"], ["cell_state_mapping"]
            )
            self.assertEqual(
                self.models[model_id]["exposure_status"],
                finding["exposure_state"],
            )
            self.assertEqual(
                finding["corpus_annotation_state"],
                "unlabeled_rank_value_encodings",
            )
        self.assertEqual(
            exposure["audit_scope"]["development_encoder_overlap"],
            "classified_for_all_registered_cell_state_task_datasets_with_unknown_retained",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "eligible_on_exposure_only_other_admission_gates_still_apply",
        )
        crosswalk = json.loads(
            (crosswalk_path := (
                CONFIG / "artifacts/models/geneformer/development_crosswalk.json"
            )).read_text(encoding="utf-8")
        )
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            set(crosswalk["findings"]),
            {
                "GSE136103",
                "GSE174748",
                "GSE185477",
                "GSE189600",
                "GSE202379",
                "GSE244832",
                "GSE256398",
                "GSE281367",
                "GSE289173",
                "GSE296875",
                "Liver_Atlas",
            },
        )

        v1_seen = {
            dataset_id
            for dataset_id, finding in crosswalk["findings"].items()
            if finding["corpus_findings"]["Genecorpus-30M"][
                "exposure_state"
            ]
            == "encoder_seen"
        }
        v2_seen = {
            dataset_id
            for dataset_id, finding in crosswalk["findings"].items()
            if finding["corpus_findings"]["Genecorpus-104M"][
                "exposure_state"
            ]
            == "encoder_seen"
        }
        self.assertEqual(v1_seen, {"GSE136103", "Liver_Atlas"})
        self.assertEqual(
            v2_seen, {"GSE136103", "GSE185477", "Liver_Atlas"}
        )
        liver_atlas = crosswalk["findings"]["Liver_Atlas"]["corpus_findings"]
        self.assertEqual(
            sum(
                record["cells_passed_filtering"]
                for record in liver_atlas["Genecorpus-30M"]["matched_records"]
            ),
            768462,
        )
        self.assertEqual(
            sum(
                record["cells_passed_filtering"]
                for record in liver_atlas["Genecorpus-104M"]["matched_records"]
            ),
            782916,
        )
        v1_clean = {
            dataset_id
            for dataset_id, finding in crosswalk["findings"].items()
            if finding["corpus_findings"]["Genecorpus-30M"][
                "exposure_state"
            ]
            == "clean_declared"
        }
        v2_clean = {
            dataset_id
            for dataset_id, finding in crosswalk["findings"].items()
            if finding["corpus_findings"]["Genecorpus-104M"][
                "exposure_state"
            ]
            == "clean_declared"
        }
        v2_unknown = {
            dataset_id
            for dataset_id, finding in crosswalk["findings"].items()
            if finding["corpus_findings"]["Genecorpus-104M"][
                "exposure_state"
            ]
            == "unknown"
        }
        self.assertEqual(
            v1_clean,
            {
                "GSE174748",
                "GSE185477",
                "GSE189600",
                "GSE202379",
                "GSE244832",
                "GSE256398",
                "GSE281367",
                "GSE296875",
            },
        )
        self.assertEqual(v2_clean, {"GSE256398", "GSE281367", "GSE296875"})
        self.assertEqual(
            v2_unknown,
            {"GSE174748", "GSE189600", "GSE202379", "GSE244832"},
        )
        self.assertTrue(
            all(
                crosswalk["findings"]["GSE289173"]["corpus_findings"][corpus][
                    "exposure_state"
                ]
                == "target_label_unexposed"
                for corpus in ("Genecorpus-30M", "Genecorpus-104M")
            )
        )

        corgi = json.loads(
            (CONFIG / "artifacts/models/corgi/release.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(corgi["weight"]["zenodo_doi"], "10.5281/zenodo.18630048")
        self.assertIn("Zenodo:18630048", self.models["corgi_regular"]["checkpoint_revision"])
        adaptation = load_toml(CONFIG / "adaptation" / "corgi.toml")
        self.assertEqual(
            adaptation["trans_input"]["released_inference_contract"],
            "order_2891_trans_regulators_then_rank_quantile_normalize_to_packaged_tf_reference",
        )
        self.assertEqual(
            adaptation["trans_input"]["candidate_mappers"],
            [
                "released_rank_quantile_from_counts",
                "length_adjusted_tpm_then_released_rank_quantile",
                "learned_count_to_rank_context_encoder",
            ],
        )
        self.assertEqual(
            adaptation["trans_input"]["redundant_mappers_excluded"],
            [
                "cpm_before_released_rank_quantile",
                "log_cpm_before_released_rank_quantile",
                "custom_monotone_quantile_before_released_rank_quantile",
            ],
        )

    def test_scgpt_upstream_preflight_and_exposure_are_frozen(self) -> None:
        checkpoints = json.loads(
            (CONFIG / "artifacts/models/scgpt/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            checkpoints["repository_revision"],
            "cebd6fae655b9c585a4807daa3ac31bb764f06b4",
        )
        self.assertEqual(checkpoints["weight_license"], "UNDECLARED")
        self.assertIs(checkpoints["weight_content_downloaded"], False)
        self.assertEqual(
            checkpoints["architecture_contract"],
            {
                "d_hid": 512,
                "dropout": 0.2,
                "embedding_size": 512,
                "fast_transformer": True,
                "max_sequence_length": 1200,
                "mvc": True,
                "n_expression_bins": 51,
                "n_heads": 8,
                "n_layers": 12,
                "n_layers_cls": 3,
                "pad_token": "<pad>",
                "pad_value": -2,
            },
        )
        expected = {
            "scgpt_whole_human": {
                "weight_bytes": 205385258,
                "args_sha256": "c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60",
                "vocab_sha256": "acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f",
            },
            "scgpt_continual": {
                "weight_bytes": 207861754,
                "args_sha256": "77fec83c32306225f37cc023321b716a49ce9e778039f2c1c4b80e1e4e7008bd",
                "vocab_sha256": "ee2b2c90158eedb97c2318e49abaaed0a02c6fdf7e3f7ca6a821906413c4d2a4",
            },
        }
        for model_id, values in expected.items():
            artifacts = {
                item["path"]: item
                for item in checkpoints["bundles"][model_id]["artifacts"]
            }
            self.assertEqual(
                artifacts["best_model.pt"]["size_bytes"], values["weight_bytes"]
            )
            self.assertEqual(artifacts["best_model.pt"]["sha256"], "UNRESOLVED")
            self.assertEqual(
                artifacts["args.json"]["sha256"], values["args_sha256"]
            )
            self.assertEqual(
                artifacts["vocab.json"]["sha256"], values["vocab_sha256"]
            )
            self.assertEqual(
                self.models[model_id]["exposure_status"],
                "target_label_unexposed",
            )
            self.assertEqual(
                self.models[model_id]["supported_tasks"], ["cell_state_mapping"]
            )
        self.assertEqual(
            checkpoints["bundles"]["scgpt_whole_human"]["training_contract"][
                "cellxgene_census_version"
            ],
            "2023-05-08",
        )
        self.assertIs(
            checkpoints["bundles"]["scgpt_continual"]["training_contract"][
                "supervised_cell_type_objective"
            ],
            True,
        )
        self.assertIn(
            "per-cell sampling seeds",
            checkpoints["input_contract"]["sampling_warning"],
        )

        exposure = json.loads(
            (CONFIG / "artifacts/models/scgpt/exposure_audit.json").read_text(
                encoding="utf-8"
            )
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/scgpt/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(set(exposure["checkpoint_findings"]), set(expected))
        self.assertEqual(
            set(crosswalk["findings"]),
            {
                "GSE136103",
                "GSE174748",
                "GSE185477",
                "GSE189600",
                "GSE202379",
                "GSE244832",
                "GSE256398",
                "GSE281367",
                "GSE289173",
                "GSE296875",
                "Liver_Atlas",
            },
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "unknown"
            },
            {
                "GSE136103",
                "GSE174748",
                "GSE185477",
                "GSE189600",
                "Liver_Atlas",
            },
        )
        acquisition = (
            ROOT / "slurm" / "acquire_scgpt_checkpoints.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=8G", acquisition)
        self.assertIn("#SBATCH --time=04:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn("No checkpoint content was deserialized", acquisition)

    def test_uce_upstream_preflight_and_exposure_are_frozen(self) -> None:
        checkpoints = json.loads(
            (CONFIG / "artifacts/models/uce/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            checkpoints["repository_revision"],
            "9c416007be15ad6753dc84af4468c1dc10421ab9",
        )
        self.assertEqual(checkpoints["weight_license"], "CC-BY-4.0")
        self.assertIs(checkpoints["weight_content_downloaded"], False)
        expected = {
            "uce_4l": {
                "layers": 4,
                "source": "Figshare:24320806v4:file:42706576",
                "bytes": 3403514339,
                "md5": "3e2f59d6da6eaa5396aa297edfcec08f",
            },
            "uce_33l": {
                "layers": 33,
                "source": "Figshare:24320806v5:file:43423236",
                "bytes": 5686283829,
                "md5": "eed62d4fe51873c87bcf9660faa07707",
            },
        }
        for model_id, values in expected.items():
            model = self.models[model_id]
            weight = checkpoints["weight_artifacts"][model_id]
            self.assertIn(values["source"], model["checkpoint_revision"])
            self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED")
            self.assertEqual(model["exposure_status"], "target_label_unexposed")
            self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
            self.assertEqual(
                checkpoints["architecture_contract"][model_id]["n_layers"],
                values["layers"],
            )
            self.assertEqual(weight["size_bytes"], values["bytes"])
            self.assertEqual(weight["source_md5"], values["md5"])
        self.assertEqual(
            {item["figshare_file_id"] for item in checkpoints["auxiliary_artifacts"]},
            {42706555, 42706558, 42706585, 42715213},
        )
        self.assertIn(
            "per-cell seed",
            checkpoints["input_contract"]["sampling_warning"],
        )

        exposure = json.loads(
            (CONFIG / "artifacts/models/uce/exposure_audit.json").read_text(
                encoding="utf-8"
            )
        )
        crosswalk_path = CONFIG / "artifacts/models/uce/development_crosswalk.json"
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(set(exposure["checkpoint_findings"]), set(expected))
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "encoder_seen"
            },
            {"GSE185477", "Liver_Atlas"},
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )
        acquisition = (
            ROOT / "slurm" / "acquire_uce_checkpoints.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=16G", acquisition)
        self.assertIn("#SBATCH --time=12:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn("No checkpoint content was deserialized", acquisition)

    def test_scimilarity_preflight_and_exposure_are_frozen(self) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/scimilarity/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["scimilarity_v1_1"]
        self.assertEqual(
            checkpoints["repository_revision"],
            "3ce3ec81e122a115eb3e667706725120f7ee252d",
        )
        self.assertEqual(checkpoints["code_license"], "Apache-2.0")
        self.assertEqual(checkpoints["weight_license"], "CC-BY-SA-4.0")
        self.assertIs(checkpoints["weight_content_downloaded"], False)
        archive = checkpoints["archive"]
        self.assertEqual(archive["size_bytes"], 30310810843)
        self.assertEqual(
            archive["source_md5"], "546251b7c435f3b1dbe38e2e420ad57f"
        )
        self.assertEqual(archive["sha256"], "UNRESOLVED")
        self.assertIn("Zenodo:10685499", model["checkpoint_revision"])
        self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(
            model["license_status"],
            "code_Apache-2.0_weights_CC-BY-SA-4.0",
        )
        self.assertEqual(model["exposure_status"], "target_label_unexposed")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(
            checkpoints["architecture_contract"],
            {
                "hidden_dimensions": [1024, 1024, 1024],
                "hidden_dropout": 0.5,
                "input_dimension": 28231,
                "input_dropout": 0.4,
                "latent_dimension": 128,
                "output_normalization": "unit_L2_hypersphere",
            },
        )
        self.assertEqual(
            checkpoints["input_contract"]["normalization"],
            "per_cell_total_10000_then_natural_log1p",
        )
        self.assertEqual(
            checkpoints["input_contract"]["minimum_gene_overlap"], 5000
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "released_loader_uses_weights_only_false"
            ],
            True,
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/scimilarity/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/scimilarity/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["scimilarity_v1_1"][
                "exposure_state"
            ],
            "target_label_unexposed",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "eligible_on_exposure_and_terms_other_admission_gates_still_apply",
        )
        self.assertEqual(crosswalk["corpus"]["annotated_training_studies"], 56)
        self.assertEqual(
            crosswalk["corpus"]["annotated_training_table_rows"], 52
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "encoder_seen"
            },
            {"GSE185477"},
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "reference_only"
            },
            {"GSE136103"},
        )
        self.assertEqual(
            crosswalk["findings"]["GSE136103"][
                "matched_reference_samples"
            ],
            26,
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )

        acquisition = (
            ROOT / "slurm" / "acquire_scimilarity_checkpoint.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=8G", acquisition)
        self.assertIn("#SBATCH --time=12:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn("md5sum -c -", acquisition)
        self.assertIn("The archive was not extracted", acquisition)

    def test_scprint_v1_preflight_rejects_silent_replacement(self) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/scprint/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["scprint_v1_5_medium"]
        self.assertEqual(
            checkpoints["repository_revision"],
            "6bbe1a7cfebb557dc745e23d1ddf5be6e38daf3f",
        )
        self.assertEqual(checkpoints["release"], "1.6.4")
        self.assertEqual(checkpoints["code_license"], "MIT")
        self.assertEqual(checkpoints["weight_license"], "Apache-2.0")

        original = checkpoints["checkpoint"]
        self.assertEqual(
            original["hf_commit"],
            "f60c87152e83fedb932a84f8ba850fbf489a5a86",
        )
        self.assertEqual(original["name"], "v2-medium.ckpt")
        self.assertEqual(original["size_bytes"], 221775592)
        self.assertEqual(
            original["sha256"],
            "a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5",
        )
        self.assertEqual(model["checkpoint_sha256"], original["sha256"])
        self.assertIn(original["hf_commit"], model["checkpoint_revision"])
        self.assertEqual(
            model["license_status"], "code_MIT_weights_Apache-2.0"
        )
        self.assertEqual(model["exposure_status"], "target_label_unexposed")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])

        alias = checkpoints["historical_alias"]
        self.assertEqual(alias["operation"], "rename_only")
        self.assertEqual(alias["sha256"], original["sha256"])
        self.assertEqual(alias["size_bytes"], original["size_bytes"])
        replacement = checkpoints["replacement_warning"]
        self.assertEqual(replacement["eligibility"], "ineligible")
        self.assertEqual(replacement["name"], "medium-v1.5.ckpt")
        self.assertNotEqual(replacement["sha256"], original["sha256"])
        self.assertNotEqual(replacement["size_bytes"], original["size_bytes"])
        self.assertNotIn(
            replacement["hf_commit"], model["checkpoint_revision"]
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "checkpoint_content_downloaded"
            ],
            False,
        )
        self.assertEqual(
            checkpoints["runtime_contract"]["declared_torch"], "2.2.0"
        )
        self.assertIn(
            "row ID",
            checkpoints["embedding_contract"]["row_identity_rule"],
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/scprint/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/scprint/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["scprint_v1_5_medium"][
                "exposure_state"
            ],
            "target_label_unexposed",
        )
        self.assertEqual(
            exposure["rejected_checkpoint"]["exposure_state"], "unknown"
        )
        self.assertEqual(
            exposure["rejected_checkpoint"]["sealed_champion_eligibility"],
            "ineligible",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "eligible_on_exposure_and_terms_other_admission_gates_still_apply",
        )
        self.assertEqual(crosswalk["corpus"]["post_qc_cells"], 54084961)
        self.assertEqual(crosswalk["corpus"]["post_qc_datasets"], 548)
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "unknown"
            },
            {
                "GSE136103",
                "GSE174748",
                "GSE185477",
                "GSE189600",
                "Liver_Atlas",
            },
        )
        self.assertEqual(
            crosswalk["findings"]["GSE185477"][
                "historical_archive_records"
            ][0]["dataset_id"],
            "ddb22b3d-a75c-4dd1-9730-dff7fc8ca530",
        )

        acquisition = (
            ROOT / "slurm" / "acquire_scprint_checkpoint.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=4G", acquisition)
        self.assertIn("#SBATCH --time=02:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn(original["hf_commit"], acquisition)
        self.assertIn(original["sha256"], acquisition)
        self.assertNotIn(replacement["sha256"], acquisition)
        self.assertIn("No checkpoint content was deserialized", acquisition)

    def test_scprint2_preflight_fails_closed_on_exposure_and_neighbors(
        self,
    ) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/scprint2/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["scprint2_small_v2"]
        self.assertEqual(
            checkpoints["repository_revision"],
            "7ee2de5aaa3492f326f540f703d44ddeaf1f2fd4",
        )
        self.assertEqual(checkpoints["release"], "1.0.3")
        self.assertEqual(checkpoints["code_license"], "GPL-3.0-or-later")
        self.assertEqual(checkpoints["weight_license"], "Apache-2.0")
        original = checkpoints["checkpoint"]
        self.assertEqual(original["size_bytes"], 889706878)
        self.assertEqual(original["sha256"], model["checkpoint_sha256"])
        self.assertIn(original["hf_commit"], model["checkpoint_revision"])
        self.assertEqual(
            checkpoints["historical_alias"]["sha256"], original["sha256"]
        )
        self.assertEqual(model["exposure_status"], "unknown")
        self.assertIn(
            "outer fold",
            checkpoints["embedding_contract"][
                "common_lane_neighbor_policy"
            ],
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/scprint2/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/scprint2/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["scprint2_small_v2"][
                "exposure_state"
            ],
            "unknown",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "ineligible_while_exposure_unknown",
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "unknown",
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "clean_declared"
            },
            {"GSE281367"},
        )

        acquisition = (
            ROOT / "slurm" / "acquire_scprint2_checkpoint.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=4G", acquisition)
        self.assertIn("#SBATCH --time=04:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn(original["hf_commit"], acquisition)
        self.assertIn(original["sha256"], acquisition)
        self.assertIn("No checkpoint content was deserialized", acquisition)

    def test_regformer_preflight_rejects_ambiguous_weight_and_label_leakage(
        self,
    ) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/regformer/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["regformer"]
        canonical = checkpoints["canonical_checkpoint"]
        self.assertEqual(
            checkpoints["publication_archive"]["repository_mirror_revision"],
            "9618f1f4bb7b61c9b236d3ea8a6339c236f2da1b",
        )
        self.assertEqual(checkpoints["code_license"], "MIT")
        self.assertEqual(checkpoints["weight_license"], "CC-BY-4.0")
        self.assertEqual(canonical["figshare_file_id"], 57423100)
        self.assertEqual(canonical["size_bytes"], 197677201)
        self.assertEqual(
            canonical["md5"], "6f0dea8ae1508077d1851cdcda6fce42"
        )
        self.assertEqual(canonical["sha256"], "UNRESOLVED")
        self.assertIn("file:57423100", model["checkpoint_revision"])
        self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(model["license_status"], "code_MIT_weights_CC-BY-4.0")
        self.assertEqual(model["exposure_status"], "unknown")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(
            checkpoints["architecture_contract"]["maximum_genes_per_cell"],
            1200,
        )
        self.assertEqual(
            checkpoints["architecture_contract"]["state_space_layers"], 10
        )
        self.assertIn(
            "true cell-type labels",
            checkpoints["embedding_contract"]["forbidden_postprocessing"],
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "checkpoint_content_downloaded"
            ],
            False,
        )

        noncanonical = {
            item["figshare_file_id"]: item
            for item in checkpoints["noncanonical_artifacts"]
        }
        self.assertEqual(
            noncanonical[57423085]["eligibility"],
            "not_admitted_ambiguous_provenance",
        )
        self.assertEqual(
            noncanonical[53164916]["eligibility"],
            "legacy_preprint_variant_not_primary_census_checkpoint",
        )
        self.assertNotIn("57423085", model["checkpoint_revision"])
        self.assertNotIn("53164916", model["checkpoint_revision"])

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/regformer/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/regformer/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["regformer"]["exposure_state"],
            "unknown",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "ineligible_while_exposure_unknown",
        )
        self.assertEqual(
            exposure["label_leakage_incident"]["benchmark_action"],
            "prohibited",
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "unknown",
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "clean_declared"
            },
            {"GSE281367", "GSE296875"},
        )

        acquisition = (
            ROOT / "slurm" / "acquire_regformer_checkpoint.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=4G", acquisition)
        self.assertIn("#SBATCH --time=02:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn("fetch 57423100 best_model.pt", acquisition)
        self.assertNotIn("57423085", acquisition)
        self.assertNotIn("53164916", acquisition)
        self.assertIn(canonical["md5"], acquisition)
        self.assertIn("The checkpoint was not deserialized", acquisition)

    def test_scbert_preflight_fails_closed_on_weight_terms_and_identity(
        self,
    ) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/scbert/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["scbert"]
        checkpoint = checkpoints["checkpoint"]
        self.assertEqual(
            checkpoints["repository_revision"],
            "262fd4b91f3f1c21a6e595d03d4ef423e16ffc99",
        )
        self.assertEqual(checkpoints["code_license"], "GPL-3.0-or-later")
        self.assertEqual(checkpoints["weight_license"], "UNDECLARED")
        self.assertEqual(checkpoint["share_id"], "AJEAIQdfAAoUxhXE7r")
        self.assertEqual(
            checkpoint["file_id"],
            "s.1970325010981265.642406688Tub_f.642406714080m",
        )
        self.assertEqual(checkpoint["name"], "panglao_pretrain.pth")
        self.assertEqual(checkpoint["size_bytes"], 87631538)
        self.assertEqual(checkpoint["sha256"], "UNRESOLVED")
        self.assertEqual(checkpoint["terms"], "UNDECLARED")
        self.assertIn(checkpoint["file_id"], model["checkpoint_revision"])
        self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(model["license_status"], "code_GPL-3.0_weights_UNDECLARED")
        self.assertEqual(model["status"], "blocked_terms")
        self.assertEqual(model["exposure_status"], "unknown")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(
            checkpoints["architecture_contract"]["gene_context_length"],
            16906,
        )
        self.assertEqual(
            checkpoints["embedding_contract"]["common_lane_dimension"],
            200,
        )
        self.assertIn(
            "not treated as a pretrained CLS token",
            checkpoints["embedding_contract"]["selection_rule"],
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "require_explicit_weight_terms_before_download"
            ],
            True,
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "checkpoint_content_downloaded"
            ],
            False,
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/scbert/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/scbert/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        finding = exposure["checkpoint_findings"]["scbert"]
        self.assertEqual(
            finding["historical_paper_checkpoint_exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(finding["registered_object_exposure_state"], "unknown")
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "ineligible_while_weight_terms_and_object_identity_are_unresolved",
        )
        self.assertEqual(crosswalk["corpus"]["cells"], 1126580)
        self.assertEqual(crosswalk["corpus"]["datasets"], 209)
        self.assertEqual(crosswalk["corpus"]["tissues"], 74)
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            crosswalk["findings"]["Liver_Atlas"]["exposure_state"],
            "downstream_demo",
        )
        self.assertFalse(
            (ROOT / "slurm" / "acquire_scbert_checkpoint.sbatch").exists()
        )

    def test_sccello_preflight_is_clean_on_exposure_but_blocked_on_terms(
        self,
    ) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/sccello/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["sccello"]
        checkpoint = checkpoints["checkpoint"]
        self.assertEqual(
            checkpoints["repository_revision"],
            "767585b6cc2ef5a3e03234164c665e751a432194",
        )
        self.assertEqual(checkpoints["code_license"], "NO_LICENSE_DETECTED")
        self.assertEqual(checkpoints["weight_license"], "NO_LICENSE_DETECTED")
        self.assertEqual(
            checkpoint["hf_commit"],
            "e59304e5f49d661b7bd142e9b926737e34ae222d",
        )
        self.assertEqual(checkpoint["size_bytes"], 44049097)
        self.assertEqual(
            checkpoint["sha256"],
            "b30ddf0aadfbcca32f1a92c3b77743091e433e4f4a38dcd77c87a89236c3e676",
        )
        self.assertEqual(model["checkpoint_sha256"], checkpoint["sha256"])
        self.assertIn(checkpoint["hf_commit"], model["checkpoint_revision"])
        self.assertEqual(
            model["license_status"],
            "blocked_no_detected_code_or_weight_license",
        )
        self.assertEqual(model["status"], "blocked_terms")
        self.assertEqual(model["exposure_status"], "target_label_unexposed")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(
            checkpoints["architecture_contract"]["hidden_dimension"], 256
        )
        self.assertEqual(
            checkpoints["architecture_contract"]["maximum_tokens"], 2048
        )
        self.assertEqual(
            checkpoints["embedding_contract"]["output_dimension"], 256
        )
        self.assertIn(
            "raw_last_layer_cls_embedding",
            checkpoints["embedding_contract"]["common_lane_candidate"],
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "require_explicit_code_and_weight_terms_before_download"
            ],
            True,
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/sccello/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/sccello/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["sccello"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "ineligible_while_code_and_weight_terms_are_unresolved",
        )
        self.assertEqual(crosswalk["corpus"]["cells"], 22316072)
        self.assertEqual(crosswalk["corpus"]["hf_shards"], 1050)
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertFalse(
            (ROOT / "slurm" / "acquire_sccello_checkpoint.sbatch").exists()
        )

    def test_langcell_preflight_resolves_bundle_but_fails_closed_on_terms(
        self,
    ) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/langcell/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["langcell"]
        modules = {
            item["module"]: item
            for item in checkpoints["checkpoint_bundle"]["modules"]
        }
        self.assertEqual(
            checkpoints["repository_revision"],
            "69e41ef2fae485b67703294b50f3150bb1d5bb9b",
        )
        self.assertEqual(checkpoints["code_license"], "MIT")
        self.assertEqual(checkpoints["weight_license"], "UNDECLARED")
        self.assertEqual(len(modules), 8)
        self.assertEqual(
            modules["cell_bert"]["file_id"],
            "119ic1oiZr6Lxa5Bguv240L-OL5lPGEAq",
        )
        self.assertEqual(modules["cell_bert"]["size_bytes"], 157274017)
        self.assertEqual(modules["cell_projection"]["size_bytes"], 526383)
        self.assertEqual(modules["text_bert"]["size_bytes"], 532629833)
        self.assertEqual(modules["text_projection"]["size_bytes"], 788527)
        self.assertEqual(
            modules["cell_bert_config"]["sha256"],
            "aea9ca5f130ef3370d4ad40e173578929232fa4cb7058e1900ac23f84d719ef6",
        )
        self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED")
        self.assertIn(modules["cell_bert"]["file_id"], model["checkpoint_revision"])
        self.assertIn(
            modules["cell_projection"]["file_id"], model["checkpoint_revision"]
        )
        self.assertEqual(model["license_status"], "code_MIT_weights_UNDECLARED")
        self.assertEqual(model["status"], "blocked_terms")
        self.assertEqual(model["exposure_status"], "target_label_unexposed")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(
            checkpoints["architecture_contract"]["cell_encoder"][
                "hidden_dimension"
            ],
            512,
        )
        self.assertEqual(
            checkpoints["embedding_contract"]["common_lane_dimension"], 512
        )
        self.assertIn(
            "never loads the text encoder",
            checkpoints["embedding_contract"]["selection_rule"],
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "require_explicit_weight_terms_before_download"
            ],
            True,
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "checkpoint_content_downloaded"
            ],
            False,
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/langcell/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/langcell/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["langcell"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "ineligible_while_weight_terms_are_unresolved",
        )
        self.assertEqual(crosswalk["corpus"]["cells_approximate"], 27500000)
        self.assertEqual(
            crosswalk["findings"]["GSE136103"]["exposure_state"],
            "encoder_seen",
        )
        self.assertEqual(
            crosswalk["findings"]["Liver_Atlas"]["exposure_state"],
            "encoder_seen",
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertFalse(
            (ROOT / "slurm" / "acquire_langcell_checkpoint.sbatch").exists()
        )

    def test_chrombpnet_preflight_is_local_sequence_only_and_donor_safe(
        self,
    ) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/chrombpnet/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["chrombpnet"]
        implementation = checkpoints["implementation"]
        self.assertEqual(
            implementation["repository_revision"],
            "eaa0fe58b6a43da62ea23b75cfc2bef4ecd3550c",
        )
        self.assertIn(
            implementation["repository_revision"], model["checkpoint_revision"]
        )
        self.assertEqual(implementation["pypi_release"]["version"], "1.0.1")
        self.assertEqual(
            implementation["pypi_release"]["sha256"],
            "fc2986a04d9cb0c0e4737c823b056ddd97edc77fac9b495fe0db9963531aad59",
        )
        self.assertEqual(model["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertEqual(model["license_status"], "MIT")
        self.assertEqual(model["exposure_status"], "target_label_unexposed")
        self.assertEqual(
            set(model["supported_tasks"]),
            {"rna_conditioned_atac", "variant_to_regulation"},
        )
        architecture = checkpoints["architecture_contract"]["assay_model"]
        self.assertEqual(architecture["input_length_bp"], 2114)
        self.assertEqual(architecture["output_length_bp"], 1000)
        self.assertEqual(architecture["filters"], 512)
        self.assertEqual(architecture["dilated_convolution_layers"], 8)
        inference = checkpoints["inference_contract"]
        self.assertIs(inference["atac_available_at_inference"], False)
        self.assertIn("sequence_only", inference["donor_context"])
        self.assertIn("ALT_minus_REF", inference["variant_primary_score"])
        self.assertIn(
            "not enhancer-gene links",
            inference["variant_to_gene_limit"],
        )
        self.assertIn(
            "+4/-4", checkpoints["input_contract"]["atac_shift"]
        )
        self.assertIs(
            checkpoints["safe_loading_contract"][
                "external_assay_checkpoint_downloaded"
            ],
            False,
        )

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/chrombpnet/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG / "artifacts/models/chrombpnet/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["chrombpnet"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            crosswalk["findings"]["gse296875"]["exposure_state"],
            "continual_seen",
        )
        self.assertEqual(
            crosswalk["findings"]["gse289173"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            crosswalk["local_model_roster"]["seeds"], [11, 29, 47, 71, 101]
        )

    def test_bpnet_preflight_is_matched_local_sequence_control(self) -> None:
        checkpoints = json.loads(
            (
                CONFIG / "artifacts/models/bpnet/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["bpnet"]
        implementation = checkpoints["implementation"]
        self.assertEqual(
            implementation["repository_revision"],
            "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
        )
        self.assertEqual(
            implementation["pypi_release"]["sha256"],
            "dccff1bfe7c1d726967e95e270b712a9eacd580756d80f6c19ae55ea1be43633",
        )
        self.assertEqual(model["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertEqual(model["license_status"], "MIT")
        self.assertEqual(model["exposure_status"], "target_label_unexposed")
        architecture = checkpoints["architecture_contract"]
        self.assertEqual(architecture["input_length_bp"], 2114)
        self.assertEqual(architecture["output_length_bp"], 1000)
        self.assertEqual(architecture["motif_module"]["filters"], [64])
        self.assertEqual(architecture["syntax_module"]["layers"], 8)
        self.assertEqual(
            checkpoints["runtime_contract"]["declared_tensorflow"], "2.4.1"
        )
        inference = checkpoints["inference_contract"]
        self.assertIs(inference["atac_available_at_inference"], False)
        self.assertEqual(
            inference["bias_inputs"],
            "none_for_the_prespecified_sequence_only_baseline",
        )
        self.assertIn("ALT_minus_REF", inference["variant_primary_score"])

        exposure = json.loads(
            (
                CONFIG / "artifacts/models/bpnet/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = CONFIG / "artifacts/models/bpnet/development_crosswalk.json"
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["bpnet"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            crosswalk["local_model_roster"]["seeds"], [11, 29, 47, 71, 101]
        )

    def test_cellplm_release_split_and_nonleaking_context_are_frozen(self) -> None:
        checkpoints = json.loads(
            (CONFIG / "artifacts/models/cellplm/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        documented = checkpoints["documented_checkpoint"]
        variant = checkpoints["registered_checkpoint"]
        self.assertEqual(documented["name"], "20230926_85M.best.ckpt")
        self.assertEqual(documented["configuration"]["latent_mod"], "gmvae")
        self.assertEqual(documented["size_bytes"], 907020069)
        self.assertEqual(variant["name"], "20231027_85M.best.ckpt")
        self.assertEqual(variant["configuration"]["latent_mod"], "vae")
        self.assertEqual(variant["size_bytes"], 906500243)
        self.assertEqual(
            checkpoints["architecture_contract"]["declared_parameter_count"],
            82402543,
        )
        self.assertEqual(checkpoints["architecture_contract"]["latent_dimension"], 512)
        self.assertEqual(checkpoints["embedding_contract"]["common_lane_dimension"], 1024)
        self.assertIn(
            "within one biological donor",
            checkpoints["embedding_contract"]["common_lane_context"],
        )
        self.assertIs(
            checkpoints["safe_loading_contract"]["checkpoint_content_downloaded"],
            False,
        )

        crosswalk_path = CONFIG / "artifacts/models/cellplm/development_crosswalk.json"
        exposure = json.loads(
            (CONFIG / "artifacts/models/cellplm/exposure_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["cellplm_20230926_documented"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertIn(
            "terms_identity_and_census_release",
            exposure["sealed_champion_eligibility"],
        )

    def test_chromatin_preflights_match_active_task_inputs(self) -> None:
        epiagent = json.loads(
            (CONFIG / "artifacts/models/epiagent/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(epiagent["checkpoint"]["size_bytes"], 5801822387)
        self.assertEqual(epiagent["weight_license"], "UNRESOLVED")
        self.assertEqual(epiagent["task_fit"]["active_tasks"], [])
        self.assertIn("RNA_only", epiagent["input_contract"]["unsupported_inputs"])
        self.assertEqual(self.models["epiagent"]["exposure_status"], "unknown")

        peakvi = json.loads(
            (CONFIG / "artifacts/models/peakvi/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            peakvi["implementation"]["repository_revision"],
            "56520c713eb1d2b245c72a0f60bc393b74198c91",
        )
        self.assertEqual(peakvi["task_fit"]["active_tasks"], [])
        self.assertEqual(self.models["peakvi"]["status"], "deferred")
        self.assertIn("1.5.0.post1", self.models["peakvi"]["checkpoint_revision"])

        scbasset = json.loads(
            (CONFIG / "artifacts/models/scbasset/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(scbasset["architecture_contract"]["input_length_bp"], 1344)
        self.assertIs(scbasset["inference_contract"]["atac_available_at_inference"], False)
        self.assertEqual(
            scbasset["inference_contract"]["held_cell_embedding"],
            "structurally_unavailable",
        )
        self.assertEqual(self.models["scbasset"]["supported_tasks"], ["rna_conditioned_atac"])

        for model_id, finding_id in {
            "epiagent": "epiagent",
            "peakvi": "peakvi",
            "scbasset": "scbasset_local",
        }.items():
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertIn(
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
                {"unknown", "target_label_unexposed"},
            )

    def test_borzoi_preflight_freezes_canonical_same_split_ensemble(self) -> None:
        checkpoints = json.loads(
            (CONFIG / "artifacts/models/borzoi/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        ensemble = checkpoints["canonical_ensemble"]
        self.assertEqual(ensemble["ensemble_size"], 4)
        self.assertIn("same genomic split", ensemble["member_relationship"])
        self.assertEqual(ensemble["released_total_size_bytes"], 2976449872)
        self.assertEqual(checkpoints["architecture_contract"]["input_length_bp"], 524288)
        self.assertEqual(checkpoints["architecture_contract"]["human_output_tracks"], 7611)
        self.assertEqual(checkpoints["weight_license"], "UNDECLARED")
        self.assertIs(checkpoints["safe_loading_contract"]["checkpoint_content_downloaded"], False)

        crosswalk_path = CONFIG / "artifacts/models/borzoi/development_crosswalk.json"
        exposure = json.loads(
            (CONFIG / "artifacts/models/borzoi/exposure_audit.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"]["borzoi_ensemble"]["exposure_state"],
            "target_label_unexposed",
        )
        self.assertEqual(
            exposure["license_disposition"]["open_champion_eligibility"],
            "ineligible_until_weight_terms_are_explicit",
        )

    def test_multivi_scglue_and_babel_pairing_and_inference_are_frozen(self) -> None:
        multivi = json.loads(
            (CONFIG / "artifacts/models/multivi/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            multivi["implementation"]["repository_revision"],
            "56520c713eb1d2b245c72a0f60bc393b74198c91",
        )
        self.assertIn("1.5.0.post1", self.models["multivi"]["checkpoint_revision"])
        self.assertIn(
            "BLOCKED_PENDING_EXPLICIT_MASK_ADAPTER",
            multivi["missingness_contract"]["admission_status"],
        )
        self.assertIn(
            "same_sample_different_aliquot_as_cell_pairs",
            multivi["input_contract"]["unsupported_pairing_topologies"],
        )
        self.assertEqual(
            set(self.models["multivi"]["supported_tasks"]),
            {"cell_state_mapping", "rna_conditioned_atac"},
        )

        scglue = json.loads(
            (CONFIG / "artifacts/models/scglue/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            scglue["implementation"]["repository_revision"],
            "091bdbbdb3dd6d652ed4677200b034ee3c85dae1",
        )
        self.assertIn("EXPERIMENTAL", scglue["inductive_inference_contract"]["experimental_cross_decode"])
        self.assertIn("dill", scglue["safe_loading_contract"]["upstream_risk"])
        self.assertEqual(self.models["scglue"]["supported_tasks"], ["cell_state_mapping"])

        babel = json.loads(
            (CONFIG / "artifacts/models/babel/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            babel["implementation"]["checkpoint_identity"]["google_drive_file_id"],
            "1uJDbiDrBb5M0d9I5hjj2Ext-N08CXESS",
        )
        self.assertEqual(babel["code_license"], "NO_LICENSE_FILE_OR_REPOSITORY_LICENSE_DETECTED_BLOCKED")
        self.assertIs(babel["safe_loading_contract"]["archive_downloaded"], False)
        self.assertIn(
            "same_sample_different_aliquot_for_training",
            babel["input_contract"]["unsupported_pairing_topologies"],
        )
        self.assertEqual(self.models["babel"]["implementation_type"], "pretrained_adapter")
        self.assertEqual(self.models["babel"]["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(self.models["babel"]["supported_tasks"], ["rna_conditioned_atac"])

        for model_id, finding_id in {
            "multivi": "multivi_local",
            "scglue": "scglue_local",
            "babel": "babel_human_v1_1",
        }.items():
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertIn(
                "target_label_unexposed",
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
            )

    def test_epibert_and_epcotv2_require_observed_atac_context(self) -> None:
        epibert = json.loads(
            (CONFIG / "artifacts/models/epibert/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(epibert["release"]["revision"], "dae61b434f885991c29e2fde3fbf7f70d0ff9b80")
        self.assertEqual(epibert["checkpoint_release"]["archive"]["size_bytes"], 3738627947)
        self.assertEqual(epibert["weight_license"], "CC-BY-4.0")
        self.assertIs(epibert["input_contract"]["accessibility"]["required"], True)
        self.assertIn("not_supported", epibert["task_fit"]["rna_conditioned_atac"])
        self.assertIn("dae61b434", self.models["epibert"]["checkpoint_revision"])
        self.assertEqual(self.models["epibert"]["exposure_status"], "target_label_unexposed")

        epcot = json.loads(
            (CONFIG / "artifacts/models/epcotv2/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        checkpoint = epcot["canonical_checkpoint"]
        self.assertEqual(checkpoint["canonical_origin"], "Paper-linked Hugging Face Space luosanj/EPCOTv2")
        self.assertEqual(checkpoint["file"]["size_bytes"], 468200299)
        self.assertEqual(
            checkpoint["file"]["git_lfs_oid_sha256"],
            "b824a80ed238e64e15c5a6eeae91329d24f04825c62d2a9dd328658822aba45d",
        )
        self.assertIs(epcot["input_contract"]["accessibility"]["required"], True)
        self.assertIn("not_supported", epcot["task_fit"]["rna_conditioned_atac"])
        self.assertIn(
            "no frozen paper eQTL classifier checkpoint was found",
            epcot["native_scoring_contract"]["published_eqtl_path"],
        )
        self.assertIn("luosanj/EPCOTv2", self.models["epcotv2"]["checkpoint_revision"])
        self.assertEqual(self.models["epcotv2"]["status"], "blocked_terms")

        for model_id in ("epibert", "epcotv2"):
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertIn("input", exposure["sealed_champion_eligibility"])

    def test_genomic_language_encoders_freeze_sequence_and_label_exposure(self) -> None:
        expected = {
            "dnabert2": ("7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba", 468354983),
            "nucleotide_transformer": ("971cb721bd8cc8134d3665b5d94b195fffb18a1480ad5925d6dcc9290bf1c17f", 1942934337),
            "hyenadna": ("deafb53209bfafb314d14bd81108546a34d473c07ed7ab7a355674376b56ad12", 218252304),
            "caduceus": ("a3e6976fe90460ff5d90d371b457d0263dc65e156ea5651b4a452e98228daef2", 30937760),
            "evo": ("c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4", 13766621200),
        }
        for model_id, (digest, size_bytes) in expected.items():
            checkpoints = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/checkpoints.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(checkpoints["checkpoint"]["sha256"], digest)
            self.assertEqual(checkpoints["checkpoint"]["size_bytes"], size_bytes)
            self.assertEqual(self.models[model_id]["checkpoint_sha256"], digest)
            self.assertEqual(
                checkpoints["representation_and_task_contract"]["admitted_task"],
                "variant_to_regulation",
            )
            # The guarantee, not one wording: a supervised signed-effect head
            # needs a development label source, and it may not be the seal.
            # dnabert2 names GSE289173 explicitly; the other four say
            # "nonsealed".  Both express the same limit.
            signed_effect_limit = checkpoints[
                "representation_and_task_contract"
            ]["signed_effect_limit"]
            self.assertTrue(
                "nonsealed" in signed_effect_limit
                or "GSE289173 cannot be used for fitting or calibration"
                in signed_effect_limit,
                f"{model_id} signed_effect_limit must exclude the sealed "
                f"dataset as a fitting source: {signed_effect_limit!r}",
            )
            self.assertIs(
                checkpoints["safe_loading_contract"]["checkpoint_content_downloaded"],
                False,
            )

            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            finding = exposure["checkpoint_finding"]
            self.assertEqual(finding["exposure_state"], "target_label_unexposed")
            self.assertIn("encoder_seen", finding["reference_sequence_exposure"])
            self.assertIn("unexposed", finding["target_effect_label_exposure"])

        self.assertEqual(
            self.models["nucleotide_transformer"]["status"],
            "restricted_comparator",
        )
        self.assertIn(
            "UNRESOLVED",
            json.loads(
                (CONFIG / "artifacts/models/evo/checkpoints.json").read_text(
                    encoding="utf-8"
                )
            )["code"]["exact_inference_engine_revision"],
        )
        gates = " ".join(self.tasks["variant_to_regulation"]["admission_gates"])
        self.assertIn("No nonsealed signed eQTL/ieQTL development source", gates)
        self.assertIn("eQTL-specific head fitting and calibration are disabled", gates)
        self.assertIn("cross-assay MPRA transfer remain eligible", gates)

    def test_midas_scmomat_and_stabmap_query_contracts_are_frozen(self) -> None:
        midas = json.loads(
            (CONFIG / "artifacts/models/midas/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(midas["implementation"]["repository_tag"], "v0.3.0")
        self.assertEqual(
            midas["implementation"]["pypi_release"]["wheel"]["sha256"],
            "24df81372decc6d3edf6c2481d7b0c58f4685386d2ac938564946916f1bd9522",
        )
        self.assertIn("weights_only=False", midas["safe_loading_contract"]["upstream_risk"])
        self.assertIn(
            "without validating ordered cell identifiers",
            midas["input_contract"]["warning"],
        )
        self.assertIn("RNA-only to ATAC", midas["inductive_inference_contract"]["designed_output"])
        self.assertEqual(
            set(self.models["midas"]["supported_tasks"]),
            {"cell_state_mapping", "rna_conditioned_atac"},
        )

        scmomat = json.loads(
            (CONFIG / "artifacts/models/scmomat/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIs(scmomat["inductive_inference_contract"]["frozen_query_transform"], False)
        self.assertIn(
            "optimizable row for every cell",
            scmomat["inductive_inference_contract"]["reason"],
        )
        self.assertEqual(self.models["scmomat"]["status"], "deferred")
        self.assertEqual(self.models["scmomat"]["supported_tasks"], [])

        stabmap = json.loads(
            (CONFIG / "artifacts/models/stabmap/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(stabmap["implementation"]["bioconductor_release"]["version"], "1.6.0")
        self.assertIn("terminal query leaf", stabmap["inductive_inference_contract"]["allowed_wrapper"])
        self.assertIn("training-reference-only", stabmap["inductive_inference_contract"]["combined_centering_fix"])
        self.assertEqual(stabmap["output_contract"]["profile_status"], "latent-only under the frozen RNA-conditioned ATAC capability registry")
        self.assertEqual(self.models["stabmap"]["supported_tasks"], ["cell_state_mapping"])

        for model_id, finding_id in {
            "midas": "midas_local",
            "scmomat": "scmomat_local",
            "stabmap": "stabmap_local",
        }.items():
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertIn(
                "target_label_unexposed",
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
            )

    def test_corgi_family_loaders_and_input_modes_fail_closed(self) -> None:
        release_path = CONFIG / "artifacts/models/corgi/release.json"
        self.assertEqual(
            hashlib.sha256(release_path.read_bytes()).hexdigest(),
            "604c1d368770a4ab02472bd718c75b05c16e32b7f273c0f5839342561bdd35be",
        )
        corgi = json.loads(
            (CONFIG / "artifacts/models/corgi/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(corgi["architecture_contract"]["declared_parameter_count"], 196000000)
        self.assertEqual(corgi["architecture_contract"]["checkpoint_state_numel"], 195870252)
        self.assertEqual(corgi["architecture_contract"]["trans_context_features"], 2891)
        self.assertEqual(corgi["output_contract"]["count"], 22)
        self.assertIn("strict=False", corgi["safe_loading_contract"]["current_loader"])
        self.assertEqual(corgi["checkpoint"]["weight_downloaded"], True)
        self.assertEqual(
            corgi["checkpoint"]["sha256"],
            self.models["corgi_regular"]["checkpoint_sha256"],
        )
        self.assertEqual(corgi["checkpoint"]["state_key_count"], 177)
        self.assertIs(
            corgi["verified_artifacts"]["strict_architecture_restore_performed"],
            True,
        )
        self.assertIs(
            corgi["verified_artifacts"]["full_524288_bp_l40s_forward_performed"],
            True,
        )
        self.assertIs(
            corgi["verified_artifacts"]["native_numeric_parity_performed"],
            False,
        )
        self.assertEqual(
            corgi["verified_artifacts"]["schema_artifacts_sha256"],
            "312bb3c90ae6da3eff20f43d920b03fbf6a5b05b2b480bbc98d24aaad71d868e",
        )
        self.assertEqual(self.models["corgi_regular"]["exposure_status"], "target_label_unexposed")

        corgi_plus = json.loads(
            (CONFIG / "artifacts/models/corgi_plus/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            corgi_plus["identity_reconciliation"]["checkpoint_mode"],
            "FINAL_26_CHANNEL_AUXILIARY_WIDTH_CONFIRMED_BY_STATE_SHAPE",
        )
        self.assertEqual(corgi_plus["architecture_contract"]["published_final_auxiliary_input"].split()[-2:], ["26", "channels"])
        self.assertIn("[128,26,1]", corgi_plus["architecture_contract"]["checkpoint_auxiliary_input"])
        self.assertEqual(corgi_plus["architecture_contract"]["checkpoint_state_numel"], 198664365)
        self.assertIn("NotImplementedError", corgi_plus["loader_audit"]["high_level_api"])
        self.assertIn("structurally missing", corgi_plus["input_contract"]["warning"])
        self.assertEqual(corgi_plus["checkpoint"]["weight_downloaded"], True)
        self.assertIs(
            corgi_plus["verified_artifacts"][
                "strict_architecture_restore_performed"
            ],
            True,
        )
        self.assertIs(
            corgi_plus["verified_artifacts"]["released_input_bundle_complete"],
            False,
        )
        self.assertEqual(
            corgi_plus["checkpoint"]["sha256"],
            self.models["corgi_plus"]["checkpoint_sha256"],
        )
        self.assertEqual(self.models["corgi_plus"]["supported_tasks"], ["rna_conditioned_atac"])
        self.assertEqual(self.models["corgi_plus"]["exposure_status"], "target_label_unexposed")

        for model_id, finding_id in {
            "corgi": "regular_corgi",
            "corgi_plus": "corgi_plus",
        }.items():
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
                "target_label_unexposed",
            )

    def test_perturbation_models_reject_incompatible_held_targets(self) -> None:
        gears = json.loads(
            (CONFIG / "artifacts/models/gears/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("noncoding variant loci", gears["query_contract"]["reason"])
        self.assertIn("pickle.loads", gears["safe_loading_contract"]["finding"])
        self.assertEqual(self.models["gears"]["supported_tasks"], ["gene_perturbation_response"])

        scgen = json.loads(
            (CONFIG / "artifacts/models/scgen/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("observed cells", scgen["query_contract"]["finding"])
        self.assertIn("held target", scgen["query_contract"]["test_outcome_rule"])
        self.assertEqual(self.models["scgen"]["supported_tasks"], ["gene_perturbation_response"])

        cpa = json.loads(
            (CONFIG / "artifacts/models/cpa/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("integer embeddings", cpa["query_contract"]["finding"])
        self.assertIn("reject_unknown_test_atoms", cpa["query_contract"]["identity_leakage_rule"])
        self.assertEqual(self.models["cpa"]["supported_tasks"], ["gene_perturbation_response"])

        cellot = json.loads(
            (CONFIG / "artifacts/models/cellot/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertIn("observed target distribution", cellot["query_contract"]["finding"])
        self.assertIn("not_a_same_cell_pair", cellot["architecture_contract"]["pairing_semantics"])
        self.assertEqual(self.models["cellot"]["supported_tasks"], ["gene_perturbation_response"])

        for model_id in ("gears", "scgen", "cpa", "cellot"):
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                exposure["checkpoint_finding"]["declared_algorithm_exposure_state"],
                "target_label_unexposed",
            )
            self.assertIn("ineligible", exposure["sealed_champion_eligibility"])

        gates = " ".join(self.tasks["perturbation_transfer"]["admission_gates"])
        self.assertIn("holds out noncoding CRISPRi target loci", gates)
        self.assertIn("deferred and unschedulable", gates)
        self.assertIn("cannot be used by pretending libraries are cells", gates)

    def test_remaining_perturbation_models_reject_unseen_noncoding_targets(self) -> None:
        expected = {
            "genepert": (
                "blocked_target_representation_mismatch",
                "nearest_gene_substitution",
            ),
            "scgpt_perturbation_head": (
                "blocked_target_representation_mismatch",
                "map_locus_or_guide_to_pad",
            ),
            "perturblib_lpm": (
                "blocked_unseen_symbol_mismatch",
                "unknown_perturbation_to_control",
            ),
        }
        for model_id, (result, prohibited) in expected.items():
            checkpoint = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/checkpoints.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                checkpoint["query_contract"]["current_task_result"],
                result,
                model_id,
            )
            self.assertIn(
                prohibited,
                checkpoint["query_contract"]["prohibited_behaviors"],
                model_id,
            )
            crosswalk_path = (
                CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            )
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )
            self.assertIn("ineligible", exposure["sealed_champion_eligibility"])
            self.assertEqual(
                self.models[model_id]["supported_tasks"],
                ["gene_perturbation_response"],
            )

        gates = " ".join(self.tasks["perturbation_transfer"]["admission_gates"])
        self.assertIn("All seven are deferred and unschedulable", gates)
        self.assertIn("random vector", gates)

    def test_scbutterfly_cobolt_and_scjoint_profile_boundaries_are_frozen(self) -> None:
        butterfly = json.loads(
            (CONFIG / "artifacts/models/scbutterfly/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            butterfly["architecture_contract"]["released_variant_for_common_screen"],
            "scButterfly-B with no cell-type, cluster, or MultiVI augmentation",
        )
        self.assertIn("requires both RNA and ATAC", butterfly["inductive_inference_contract"]["upstream_api_audit"])
        self.assertIn("must never be derived from held ATAC", butterfly["inductive_inference_contract"]["inverse_tfidf_warning"])
        self.assertIn(
            "continuous frozen-peak decoder scores",
            butterfly["task_fit"]["rna_conditioned_atac"],
        )

        cobolt = json.loads(
            (CONFIG / "artifacts/models/cobolt/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(cobolt["implementation"]["release_tag"], "v1.0.1")
        self.assertIn("unrelated", cobolt["implementation"]["warning"])
        self.assertIn("disables dataset-specific", cobolt["inductive_inference_contract"]["new_dataset_rule"])
        self.assertIn("not a fragment count", cobolt["inductive_inference_contract"]["profile_validity"])
        self.assertIn("metadata_warning", cobolt["terms_audit"]["result"])

        scjoint = json.loads(
            (CONFIG / "artifacts/models/scjoint/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(scjoint["terms_audit"]["result"], "blocked_no_detected_license")
        self.assertIn("co-trains", scjoint["inductive_inference_contract"]["native_behavior"])
        self.assertEqual(scjoint["output_contract"]["profile_status"], "latent_only")
        self.assertIn("no peak decoder", scjoint["output_contract"]["not_admitted"])

        for model_id, finding_id in {
            "scbutterfly": "scbutterfly_local",
            "cobolt": "cobolt_local",
            "scjoint": "scjoint_repository_and_demo_checkpoint",
        }.items():
            crosswalk_path = CONFIG / f"artifacts/models/{model_id}/development_crosswalk.json"
            exposure = json.loads(
                (CONFIG / f"artifacts/models/{model_id}/exposure_audit.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
            )
            self.assertIn(
                exposure["checkpoint_findings"][finding_id]["exposure_state"],
                {"target_label_unexposed", "unknown"},
            )

    def test_transcriptformer_preflight_fails_closed_on_exposure(self) -> None:
        checkpoints = json.loads(
            (
                CONFIG
                / "artifacts/models/transcriptformer/checkpoints.json"
            ).read_text(encoding="utf-8")
        )
        model = self.models["transcriptformer_tf_sapiens"]
        self.assertEqual(
            checkpoints["repository_revision"],
            "c943a89a4de9511a8ab1010715bf1cfa8827cc99",
        )
        self.assertEqual(checkpoints["release"], "v0.6.0")
        self.assertEqual(checkpoints["weight_license"], "MIT")
        self.assertIs(checkpoints["weight_content_downloaded"], False)
        archive = checkpoints["archive"]
        self.assertEqual(archive["size_bytes"], 1819975447)
        self.assertEqual(
            archive["etag"], "a841d3efcfc21eb14c41c6151bc39b73-217"
        )
        self.assertEqual(
            archive["s3_version_id"],
            "P5N5hx41EeWC7I0l8oW7yDyiJHS5InPU",
        )
        self.assertEqual(archive["sha256"], "UNRESOLVED")
        self.assertIn(archive["s3_version_id"], model["checkpoint_revision"])
        self.assertEqual(model["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(model["license_status"], "MIT")
        self.assertEqual(model["exposure_status"], "unknown")
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(
            checkpoints["architecture_contract"],
            {
                "assay_token_count": 1,
                "attention_heads": 16,
                "context_genes": 2047,
                "hidden_dimension": 2048,
                "transformer_layers": 12,
                "trainable_parameters": 368000000,
                "frozen_embedding_parameters": 61000000,
                "vocabulary_size": 23829,
            },
        )
        self.assertEqual(checkpoints["input_contract"]["count_clip"], 30)
        self.assertEqual(
            checkpoints["input_contract"]["feature_identifier"],
            "Ensembl_gene_stable_id",
        )
        self.assertIs(
            checkpoints["input_contract"]["inference_default_randomize_genes"],
            False,
        )
        self.assertIs(
            checkpoints["extraction_contract"][
                "upstream_cli_is_admission_safe"
            ],
            False,
        )

        exposure = json.loads(
            (
                CONFIG
                / "artifacts/models/transcriptformer/exposure_audit.json"
            ).read_text(encoding="utf-8")
        )
        crosswalk_path = (
            CONFIG
            / "artifacts/models/transcriptformer/development_crosswalk.json"
        )
        crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
        self.assertEqual(
            exposure["development_crosswalk_record"]["sha256"],
            hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            exposure["checkpoint_findings"][
                "transcriptformer_tf_sapiens"
            ]["exposure_state"],
            "unknown",
        )
        self.assertEqual(
            exposure["sealed_champion_eligibility"],
            "ineligible_while_exposure_unknown",
        )
        self.assertEqual(
            set(crosswalk["findings"]),
            {
                "GSE136103",
                "GSE174748",
                "GSE185477",
                "GSE189600",
                "GSE202379",
                "GSE244832",
                "GSE256398",
                "GSE281367",
                "GSE289173",
                "GSE296875",
                "Liver_Atlas",
            },
        )
        self.assertEqual(
            crosswalk["findings"]["GSE289173"]["exposure_state"],
            "unknown",
        )
        self.assertEqual(
            {
                dataset_id
                for dataset_id, finding in crosswalk["findings"].items()
                if finding["exposure_state"] == "clean_declared"
            },
            {"GSE281367", "GSE296875"},
        )
        self.assertEqual(
            crosswalk["corpus"]["source_record_sha256"], "UNRESOLVED"
        )
        self.assertEqual(
            crosswalk["corroborating_live_catalog"][
                "gse289173_accession_or_title_matches"
            ],
            0,
        )

        acquisition = (
            ROOT / "slurm" / "acquire_transcriptformer_checkpoint.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --mem=8G", acquisition)
        self.assertIn("#SBATCH --time=04:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertIn("versionId=${version_id}", acquisition)
        self.assertIn("The archive was not extracted", acquisition)
        self.assertNotIn("tar -", acquisition)

    def test_first_geneformer_smoke_lane_is_exact_and_fail_closed(self) -> None:
        campaign = self.campaigns_by_file["v1_cell_geneformer_smoke"]
        self.assertEqual(campaign["campaign_id"], "v1-cell-geneformer-smoke")
        self.assertEqual(campaign["task_ids"], ["cell_state_mapping"])
        self.assertEqual(
            campaign["selection"]["model_ids"],
            [
                "geneformer_v1_10m",
                "geneformer_v2_104m",
                "geneformer_v2_316m",
            ],
        )
        self.assertIs(campaign["selection"]["include_mandatory_baselines"], True)
        self.assertEqual(campaign["execution"]["cell_budget"], 1000)
        self.assertEqual(campaign["execution"]["seeds"], [1103])
        self.assertEqual(
            campaign["execution"]["common_head_ids"],
            ["linear", "two_layer_mlp"],
        )
        self.assertEqual(
            campaign["execution"]["hyperparameters"]["head_policy"],
            "separate_training_only_linear_and_two_layer_mlp_runs",
        )
        self.assertEqual(
            campaign["execution"]["hyperparameters"]["common_head_contract_id"],
            "cell_state_common_heads_v1",
        )
        self.assertEqual(
            campaign["dataset_view_ids"],
            ["resource_atlas_geneformer_smoke_1000_v1"],
        )
        self.assertEqual(
            campaign["status"], "blocked_environment_and_checkpoints"
        )
        self.assertIs(campaign["submit_enabled"], False)
        self.assertTrue(campaign["status"].startswith("blocked_"))
        runtime = self.runtimes_by_file["gpu_geneformer"]
        self.assertIs(runtime["admission_blocking"], True)
        self.assertEqual(runtime["environment_lock"], "UNRESOLVED")

        acquisition = (
            ROOT / "slurm" / "acquire_geneformer_checkpoints.sbatch"
        ).read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=io", acquisition)
        self.assertIn("#SBATCH --cpus-per-task=1", acquisition)
        self.assertIn("#SBATCH --mem=16G", acquisition)
        self.assertIn("#SBATCH --time=08:00:00", acquisition)
        self.assertNotIn("#SBATCH --array", acquisition)
        self.assertNotRegex(acquisition, r"\bsbatch\b")
        geneformer = json.loads(
            (CONFIG / "artifacts/models/geneformer/checkpoints.json").read_text(
                encoding="utf-8"
            )
        )
        expected_artifacts = {
            item["path"]: item["sha256"]
            for item in geneformer["artifacts"]
        }
        for roster in geneformer["auxiliary_artifacts"].values():
            for item in roster:
                expected_artifacts[item["path"]] = item["sha256"]
        observed_artifacts = {
            path: digest
            for digest, path in re.findall(
                r"^([0-9a-f]{64})  (\S+)$", acquisition, flags=re.MULTILINE
            )
        }
        self.assertEqual(observed_artifacts, expected_artifacts)
        for model_id in (
            "geneformer_v1_10m",
            "geneformer_v2_104m",
            "geneformer_v2_316m",
        ):
            self.assertIn(model_id, acquisition)

        frozen_screen = self.campaigns_by_file["v1_frozen_screen"]
        self.assertEqual(
            frozen_screen["execution"]["common_head_ids"],
            ["linear", "two_layer_mlp"],
        )

    def test_first_real_cell_baseline_is_exact_and_submission_guarded(self) -> None:
        campaign = self.campaigns_by_file["v1_cell_baseline_smoke"]
        self.assertEqual(campaign["campaign_id"], "v1-cell-baseline-smoke")
        self.assertEqual(campaign["status"], "blocked_completed")
        self.assertIs(campaign["submit_enabled"], False)
        self.assertIs(campaign["allow_arrays"], False)
        self.assertEqual(campaign["selection"]["model_ids"], ["hvg_pca_logistic"])
        self.assertEqual(campaign["execution"]["seeds"], [1103])
        self.assertEqual(campaign["execution"]["folds"], [0])
        self.assertEqual(campaign["execution"]["cell_budget"], 1000)
        self.assertEqual(campaign["execution"]["hyperparameters"]["outer_folds"], 5)
        self.assertEqual(
            campaign["dataset_view_ids"],
            ["resource_atlas_geneformer_smoke_1000_v1"],
        )

        model = self.models["hvg_pca_logistic"]
        self.assertIs(model["admission_blocking"], False)
        self.assertEqual(model["supported_tasks"], ["cell_state_mapping"])
        self.assertEqual(model["checkpoint_sha256"], "NOT_APPLICABLE")
        self.assertIs(model["execution"]["ready"], True)
        self.assertEqual(
            model["execution"]["supported_actions"],
            ["prepare", "fit", "predict"],
        )
        self.assertEqual(
            model["execution"]["supported_adaptation_regimes"],
            ["native_lane"],
        )
        self.assertEqual(
            model["execution"]["environment_artifact"]["sha256"],
            "8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298",
        )
        self.assertNotIn(
            "hvg_pca_logistic",
            self.tasks["bulk_state_transfer"]["baseline_model_ids"],
        )

    def test_classical_cell_baseline_screen_is_exact_and_submission_guarded(self) -> None:
        campaign = self.campaigns_by_file["v1_cell_classical_baselines_smoke"]
        self.assertEqual(
            campaign["campaign_id"], "v1-cell-classical-baselines-smoke"
        )
        self.assertEqual(campaign["status"], "ready")
        self.assertIs(campaign["submit_enabled"], True)
        self.assertIs(campaign["allow_arrays"], False)
        model_ids = [
            "hvg_pca_nearest_centroid",
            "hvg_pca_knn",
            "hvg_pca_logistic",
            "hvg_pca_elastic_net",
            "hvg_pca_linear_svm",
        ]
        self.assertEqual(campaign["selection"]["model_ids"], model_ids)
        self.assertEqual(campaign["execution"]["seeds"], [1103])
        self.assertEqual(campaign["execution"]["folds"], [0])
        self.assertEqual(campaign["execution"]["cell_budget"], 1000)
        self.assertEqual(
            campaign["execution"]["hyperparameters"]["svm_calibration_folds"],
            3,
        )
        for model_id in model_ids:
            model = self.models[model_id]
            self.assertIs(model["admission_blocking"], False)
            expected_supported = ["cell_state_mapping"]
            if model_id in {
                "hvg_pca_nearest_centroid",
                "hvg_pca_elastic_net",
                "hvg_pca_linear_svm",
            }:
                expected_supported = ["bulk_state_transfer", "cell_state_mapping"]
            self.assertEqual(model["supported_tasks"], expected_supported)
            self.assertEqual(model["checkpoint_sha256"], "NOT_APPLICABLE")
            self.assertIs(model["execution"]["ready"], True)
            self.assertEqual(
                model["execution"]["executable_tasks"], ["cell_state_mapping"]
            )
            self.assertEqual(
                model["execution"]["supported_actions"],
                ["prepare", "fit", "predict"],
            )

    def test_gse289173_preflight_and_seal_incident_are_recorded(self) -> None:
        preflight = json.loads(
            (
                CONFIG
                / "artifacts/datasets/gse289173/zenodo_central_directory_preflight.json"
            ).read_text(encoding="utf-8")
        )
        self.assertIs(preflight["content_opened"], False)
        self.assertEqual(preflight["inventory"]["donor_matrix_directories"], 48)
        self.assertEqual(
            preflight["central_directory"]["sha256"],
            "60691c678ae499509506ed37b0bd92a205baf54fa99689a77941db710abc9b35",
        )
        incident = json.loads(
            (
                CONFIG
                / "artifacts/incidents/gse289173_geo_metadata_20260821.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(incident["affected_primary_task_outcomes"], [])
        self.assertIs(incident["primary_cell_state_labels_opened"], False)
        self.assertIs(incident["primary_eqtl_or_ieqtl_outcomes_opened"], False)
        self.assertIs(
            incident["containment"]["disease_control_secondary_sealed_eligible"],
            False,
        )
        protocol = load_toml(CONFIG / "evaluation" / "sealed_protocol.toml")
        self.assertIs(
            protocol["gse289173"]["disease_control_secondary_sealed_eligible"],
            False,
        )

    def test_no_floating_latest_tags(self) -> None:
        for path in sorted(CONFIG.rglob("*.toml")):
            text = path.read_text(encoding="utf-8")
            self.assertIsNone(
                re.search(r"(?i)(?:^|[/:@._-])latest(?:$|[/:@._-])", text),
                path,
            )

    def test_every_campaign_binds_the_frozen_model_census(self) -> None:
        expected = hashlib.sha256(
            json.dumps(
                sorted(MODEL_IDS),
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
            ).encode("utf-8")
        ).hexdigest()
        for path in sorted((CONFIG / "campaigns").glob("*.toml")):
            campaign = load_toml(path)
            self.assertEqual(
                campaign["selection"]["census_model_ids_sha256"],
                expected,
                path.name,
            )

    def test_sealed_dataset_is_never_planning_visible(self) -> None:
        sealed = self.datasets["gse289173"]
        self.assertEqual(sealed["role"], "withheld_sealed")
        self.assertEqual(sealed["status"], "withheld_sealed")
        self.assertEqual(sealed["modality_status_default"], "withheld_sealed")
        self.assertEqual(sealed["outcome_url"], "WITHHELD_SEALED")
        self.assertIs(sealed["automatic_download"], False)
        self.assertIs(sealed["admission_blocking"], True)
        self.assertEqual(sealed["split"]["outer_role"], "withheld_sealed")
        self.assertEqual(sealed["split"]["label_visibility"], "withheld_sealed")
        self.assertIn("one co-unblinding event", sealed["notes"])
        self.assertIn("10.5281/zenodo.14586466", " ".join(sealed["blockers"]))

        for task in self.tasks.values():
            self.assertNotIn("gse289173", task["datasets_train"])
            self.assertNotIn("gse289173", task["datasets_development"])

        conditional_bulk = self.datasets["gse267031"]
        self.assertEqual(conditional_bulk["role"], "withheld_sealed")
        self.assertEqual(conditional_bulk["outcome_url"], "WITHHELD_SEALED")
        self.assertIs(conditional_bulk["automatic_download"], False)
        self.assertIs(conditional_bulk["admission_blocking"], True)
        for task in self.tasks.values():
            self.assertNotIn("gse267031", task["datasets_train"])
            self.assertNotIn("gse267031", task["datasets_development"])

        fnih = self.datasets["fnih_86_liver"]
        self.assertEqual(fnih["role"], "deferred")
        self.assertEqual(fnih["expected_biological_units"], 86)
        self.assertIs(fnih["admission_blocking"], True)

    def test_fnih_activation_path_remains_fail_closed(self) -> None:
        activation = self.campaigns_by_file["v1_fnih_activation"]
        self.assertEqual(activation["campaign_id"], "v1-fnih-activation")
        self.assertEqual(activation["status"], "blocked_data_access")
        self.assertEqual(activation["activation_manifest_sha256"], "UNRESOLVED")
        self.assertEqual(activation["selection_lock_id"], "UNRESOLVED")
        self.assertIs(activation["allow_sealed_features"], False)
        self.assertIs(activation["allow_sealed_labels"], False)
        self.assertIs(activation["submit_enabled"], False)
        self.assertEqual(
            set(activation["task_ids"]), {"rna_conditioned_atac", "unified_model"}
        )

        sealed = self.campaigns_by_file["v1_sealed_inference"]
        self.assertEqual(sealed["selection_lock_id"], "UNRESOLVED")
        self.assertEqual(sealed["fnih_activation_manifest_sha256"], "UNRESOLVED")
        self.assertIn("rna_conditioned_atac", sealed["task_ids"])
        self.assertIn("unified_model", sealed["task_ids"])
        self.assertIs(sealed["submit_enabled"], False)

    def test_intervention_stress_uses_three_deposited_batches(self) -> None:
        stress = self.datasets["gse313774"]
        self.assertEqual(stress["biological_unit"], "independent_experimental_batch")
        self.assertEqual(stress["expected_biological_units"], 3)
        self.assertEqual(stress["split"]["group_key"], "cohort_family_id+batch_id")
        contract_text = " ".join(stress["blockers"] + [stress["notes"]])
        self.assertIn("33 deposited libraries", contract_text)
        self.assertIn("three independent batches", contract_text)
        self.assertIn("no resmetirom", contract_text)
        self.assertIn("CXCL1", contract_text)
        self.assertIn("IL8", contract_text)
        policy = stress["prediction_first_outcome"]
        self.assertEqual(policy["scope"], "non_champion_stress_only")
        self.assertEqual(
            policy["outcome_visibility"],
            "withheld_until_prediction_commit",
        )
        self.assertEqual(policy["prediction_scope"], "genome_wide")
        self.assertEqual(policy["prediction_universe_sha256"], "UNRESOLVED")
        self.assertEqual(
            policy["prediction_universe_role"],
            "prediction_universe:gse313774",
        )
        self.assertEqual(policy["excluded_endpoint_ids"], ["CXCL1", "IL8"])
        self.assertIs(policy["forbid_fit"], True)
        self.assertIs(policy["forbid_model_selection"], True)
        self.assertIs(policy["require_prediction_commit"], True)
        self.assertIn("gse313774", self.tasks["perturbation_transfer"]["datasets_development"])
        task_text = " ".join(
            self.tasks["perturbation_transfer"]["admission_gates"]
            + [self.tasks["perturbation_transfer"]["claim_gate"]]
        )
        self.assertIn("three-batch", task_text)
        self.assertIn("no deposited resmetirom", task_text)

    def test_gse296875_histopathology_is_development_only_and_fail_closed(self) -> None:
        substrate = self.datasets["gse296875"]
        self.assertEqual(substrate["role"], "train_development")
        self.assertEqual(substrate["status"], "available")
        self.assertIs(substrate["admission_blocking"], True)
        self.assertEqual(substrate["expected_biological_units"], 39)
        self.assertIn("pathology_scores", substrate["modalities"])
        self.assertIn("clinical_covariates", substrate["modalities"])
        self.assertIn("person_key", substrate["split"]["group_key"])
        self.assertNotIn("well_id", substrate["split"]["group_key"])
        contract_text = " ".join(substrate["blockers"] + [substrate["notes"]])
        self.assertIn("barcode-to-donor", contract_text)
        self.assertIn("leave-one-well-out", contract_text)
        self.assertIn("steatosis is observed for 38", contract_text)
        self.assertIn("fibrosis for 37", contract_text)
        self.assertIn("fibrosis_categorical", contract_text)
        self.assertIn("steatosis_numeric", contract_text)
        self.assertIn("BMI", contract_text)
        self.assertIn("no explicit MASLD or MASH adjudication", contract_text)
        self.assertIn("do not make a MASLD diagnosis", contract_text)
        self.assertIn("10.1016/j.ajhg.2025.11.009", contract_text)
        self.assertIn("mmc2.xlsx", contract_text)
        for task_id in {"cell_state_mapping", "rna_conditioned_atac", "unified_model"}:
            task_text = " ".join(self.tasks[task_id]["admission_gates"])
            self.assertIn("barcode-to-donor", task_text, task_id)
            self.assertIn("leave-one-well-out", task_text, task_id)

        endpoint = self.gse296875_histopathology
        self.assertEqual(endpoint["role"], "secondary_development_only")
        self.assertEqual(endpoint["unit_of_inference"], "donor")
        self.assertEqual(endpoint["analyzed_donors"], 39)
        self.assertEqual(endpoint["technical_batches"], 8)
        self.assertEqual(endpoint["donors_under_18"], 5)
        self.assertEqual(
            (endpoint["minimum_age_years"], endpoint["maximum_age_years"]),
            (13, 75),
        )
        self.assertEqual(
            endpoint["source_supplement_sha256"],
            "8356c2002059073faca89831181c8dc60d87857ebd408fe1021a7fabc15f73eb",
        )
        self.assertEqual(endpoint["source_supplement_size_bytes"], 5884547)
        self.assertEqual(endpoint["steatosis"]["observed_donors"], 38)
        self.assertEqual(endpoint["fibrosis"]["observed_donors"], 37)
        self.assertEqual(endpoint["fibrosis"]["negative_source_values"], ["None"])
        self.assertEqual(
            endpoint["fibrosis"]["positive_source_values"], ["mild", "yes"]
        )
        self.assertIs(endpoint["fibrosis"]["ordinal_stage_inference"], False)
        self.assertIs(endpoint["external_or_sealed"], False)
        self.assertIs(endpoint["champion_eligible"], False)
        self.assertIs(endpoint["admission_blocking"], True)
        self.assertIn("adult-only", endpoint["covariates"]["age_scope_policy"])
        self.assertEqual(endpoint["inference"]["bootstrap_replicates"], 10000)
        self.assertIn("MASLD diagnosis", endpoint["claims"]["forbidden"])

    def test_atac_profile_and_latent_capabilities_are_disjoint(self) -> None:
        capability = self.atac_capabilities
        self.assertEqual(
            capability["schema_version"], "masld-bench-atac-capability-v2"
        )
        self.assertEqual(capability["task_id"], "rna_conditioned_atac")
        direct_models = set(capability["direct_profile_candidate_models"])
        conditional_models = set(
            capability["conditional_profile_candidate_models"]
        )
        baseline_models = set(capability["profile_baseline_models"])
        profile_models = direct_models | conditional_models | baseline_models
        latent_models = set(capability["latent_only_models"])
        task_models = {
            model["model_id"]
            for family in self.models_by_file.values()
            for model in family["models"]
            if "rna_conditioned_atac" in model["supported_tasks"]
        }
        self.assertFalse(direct_models & conditional_models)
        self.assertFalse(direct_models & baseline_models)
        self.assertFalse(conditional_models & baseline_models)
        self.assertFalse(profile_models & latent_models)
        self.assertEqual(profile_models, task_models)
        self.assertEqual(
            self.tasks["rna_conditioned_atac"]["baseline_model_ids"],
            capability["profile_baseline_models"],
        )
        self.assertEqual(
            capability["profile_fixture_passed_models"],
            [
                "assay_native_pseudobulk",
                "cobolt",
                "mean_track",
                "nearest_context",
                "scpair",
                "shrunken_pseudobulk",
                "shuffled_context",
                "trans_only",
            ],
        )
        for model_id in latent_models:
            self.assertNotIn(
                "rna_conditioned_atac", self.models[model_id]["supported_tasks"]
            )
        passed_models = set(capability["profile_fixture_passed_models"])
        for model_id in profile_models - passed_models:
            self.assertIs(self.models[model_id]["admission_blocking"], True)
        for model_id in passed_models:
            self.assertIs(self.models[model_id]["admission_blocking"], False)
        self.assertEqual(
            set(capability["profile_metrics"]),
            {
                "conditional_multinomial_profile_deviance_skill",
                "donor_peak_auprc",
                "donor_profile_spearman",
            },
        )
        self.assertEqual(capability["retrieval_metrics"], ["donor_cell_state_retrieval"])
        self.assertNotIn(
            "donor_cell_state_retrieval",
            self.tasks["rna_conditioned_atac"]["metrics"],
        )
        self.assertIs(capability["admission_blocking"], True)
        task_text = " ".join(
            self.tasks["rna_conditioned_atac"]["admission_gates"]
            + [self.tasks["rna_conditioned_atac"]["claim_gate"]]
        )
        self.assertIn("Latent-only integration models are not scheduled", task_text)
        self.assertIn("cannot support an ATAC-profile winner claim", task_text)

    def test_variant_capabilities_are_endpoint_safe_and_complete(self) -> None:
        capability = self.variant_capabilities
        self.assertEqual(
            capability["schema_version"], "masld-bench-variant-capability-v1"
        )
        self.assertEqual(capability["task_id"], "variant_to_regulation")
        self.assertEqual(
            capability["sealed_available_context"],
            ["dna_sequence", "single_nucleus_rna"],
        )
        self.assertEqual(
            capability["sealed_unavailable_context"],
            [
                "cell_type_eqtl_summary",
                "interaction_eqtl_summary",
                "single_nucleus_atac",
            ],
        )
        records = capability["models"]
        record_ids = [record["model_id"] for record in records]
        self.assertEqual(record_ids, sorted(record_ids))
        task_models = sorted(
            model_id
            for model_id, model in self.models.items()
            if "variant_to_regulation" in model["supported_tasks"]
        )
        self.assertEqual(record_ids, task_models)
        by_model = {record["model_id"]: record for record in records}
        mandatory = {
            model_id
            for model_id, record in by_model.items()
            if record["is_mandatory_baseline"]
        }
        self.assertEqual(
            mandatory,
            set(self.tasks["variant_to_regulation"]["baseline_model_ids"]),
        )
        primary = {
            model_id
            for model_id, record in by_model.items()
            if record["primary_eligible"]
        }
        self.assertEqual(
            primary,
            {
                "borzoi_ensemble",
                "context_borzoi",
                "corgi_regular",
                "sequence_only",
                "shuffled_context",
            },
        )
        observed_context_only = {
            model_id
            for model_id, record in by_model.items()
            if record["role"] == "observed_context_only"
        }
        self.assertEqual(
            observed_context_only,
            {
                "epibert",
                "epcotv2",
                "get",
                "scooby_epicardioids",
                "scooby_neurips",
            },
        )
        for model_id in observed_context_only:
            record = by_model[model_id]
            self.assertIs(record["requires_observed_target_context"], True)
            self.assertIs(record["primary_eligible"], False)
        for model_id in {"abc", "epinformer", "re2g"}:
            record = by_model[model_id]
            self.assertEqual(record["role"], "link_only")

        for model_id in {"alphagenome", "decima"}:
            record = by_model[model_id]
            self.assertEqual(record["role"], "secondary_only")
            self.assertIs(record["primary_eligible"], False)
            self.assertNotIn("signed_cell_type_eqtl_effect", record["allowed_endpoints"])
            # A secondary-only model may use any registered endpoint EXCEPT
            # the primary one.  Freezing a literal whitelist here went stale
            # when the endpoint roster grew; bind it to the registry instead so
            # the real guarantee (no primary endpoint, nothing off-roster) is
            # what is enforced.
            self.assertTrue(
                set(record["allowed_endpoints"])
                <= set(capability["endpoint_ids"])
                - {capability["primary_endpoint_id"]},
                f"{model_id} allows an endpoint outside the non-primary roster",
            )
        for model_id in {
            "caduceus",
            "chrombert",
            "dnabert2",
            "evo",
            "hyenadna",
            "nucleotide_transformer",
        }:
            record = by_model[model_id]
            self.assertEqual(record["role"], "representation_only")
            self.assertIs(record["requires_fitted_head"], True)
            self.assertNotIn(
                "signed_cell_type_eqtl_effect", record["allowed_endpoints"]
            )
        task_text = " ".join(
            self.tasks["variant_to_regulation"]["admission_gates"]
            + [self.tasks["variant_to_regulation"]["claim_gate"]]
        )
        self.assertIn("frozen variant capability registry is binding", task_text)
        self.assertIn("cannot be treated as signed-effect baselines", task_text)
        self.assertIn("Link-only", task_text)
        self.assertIs(capability["admission_blocking"], True)

    def test_variant_development_proxies_cannot_be_relabelled_as_eqtl(self) -> None:
        policy = self.variant_proxy_selection
        self.assertEqual(
            policy["schema_version"],
            "masld-bench-variant-proxy-selection-v1",
        )
        self.assertEqual(
            policy["confirmation_endpoint_id"],
            "signed_cell_type_eqtl_effect",
        )
        self.assertEqual(policy["confirmation_dataset_id"], "gse289173")
        self.assertIs(policy["eqtl_specific_head_fitting_allowed"], False)
        self.assertIs(policy["weighted_composite_allowed"], False)
        self.assertIs(policy["sealed_results_used"], False)
        by_endpoint = {
            item["endpoint_id"]: item for item in policy["proxy_endpoints"]
        }
        self.assertEqual(
            set(by_endpoint),
            {"mpra_allelic_direction", "rna_atac_regulatory_profile"},
        )
        self.assertEqual(
            by_endpoint["mpra_allelic_direction"]["dataset_ids"],
            ["gse281364"],
        )
        self.assertEqual(
            by_endpoint["rna_atac_regulatory_profile"]["dataset_ids"],
            ["gse244832", "gse296875"],
        )
        self.assertEqual(policy["diagnostic_only_dataset_ids"], ["gse281160"])
        self.assertIn("not eQTL performance", policy["claim_boundary"])

    def test_cell_state_sealed_modality_boundary_is_frozen(self) -> None:
        capability = self.cell_state_capabilities
        self.assertEqual(
            capability["schema_version"],
            "masld-bench-cell-state-capability-v1",
        )
        self.assertEqual(capability["sealed_dataset_id"], "gse289173")
        self.assertEqual(
            capability["sealed_input_modalities"], ["single_nucleus_rna"]
        )
        self.assertEqual(
            self.tasks["cell_state_mapping"]["input_modalities"],
            ["single_cell_rna", "single_nucleus_rna"],
        )
        atac_only = set(capability["atac_only_development_models"])
        self.assertEqual(
            atac_only,
            {"cistopic", "epiagent", "epifoundation", "lsi", "peakvi"},
        )
        for model_id in atac_only:
            self.assertEqual(self.models[model_id]["supported_tasks"], [])
            self.assertEqual(self.models[model_id]["status"], "deferred")
            self.assertIs(self.models[model_id]["admission_blocking"], True)
        self.assertIs(capability["admission_blocking"], True)

    def test_variant_training_routes_are_assay_native(self) -> None:
        task = self.tasks["variant_to_regulation"]
        self.assertEqual(
            task["datasets_train"],
            ["current_coloc_gwas", "gse244832", "gse281364", "gse296875"],
        )
        self.assertEqual(task["datasets_development"], ["gse281160"])
        gates = " ".join(task["admission_gates"])
        self.assertIn("locus-cross-fitted training source", gates)
        self.assertIn("guides and cells are not biological replicates", gates)
        self.assertIn("donor-held GSE296875/GSE244832", gates)

    def test_bulk_baselines_are_assay_native_and_inductive(self) -> None:
        baselines = set(self.tasks["bulk_state_transfer"]["baseline_model_ids"])
        self.assertEqual(
            baselines,
            {
                "assay_native_pseudobulk",
                "hvg_pca_elastic_net",
                "hvg_pca_linear_svm",
                "hvg_pca_nearest_centroid",
                "mean_expression_bulk",
            },
        )
        self.assertFalse(
            {"harmony_logistic", "scanvi_baseline", "scvi_baseline"} & baselines
        )

    def test_adaptation_campaign_uses_the_frozen_regime_vocabulary(self) -> None:
        campaign_regimes = self.campaigns_by_file["v1_adaptation"]["execution"][
            "adaptation_regimes"
        ]
        native_regimes = load_toml(
            CONFIG / "adaptation" / "common_and_native.toml"
        )["native_lane"]["regimes"]
        self.assertEqual(campaign_regimes, native_regimes)
        self.assertNotIn("peft", campaign_regimes)
        common = load_toml(
            CONFIG / "adaptation" / "common_and_native.toml"
        )["common_lane"]["cell_state_head_contract"]
        self.assertEqual(common["contract_id"], "cell_state_common_heads_v1")
        self.assertEqual(common["head_ids"], ["linear", "two_layer_mlp"])
        self.assertEqual(common["biological_group"], "donor")
        self.assertEqual(common["inner_folds"], 5)
        self.assertTrue(
            common["probability_calibration"].startswith("one_positive_temperature")
        )

    def test_gse249997_is_early_masld_domain_transport_only(self) -> None:
        dataset = self.datasets["gse249997"]
        self.assertEqual(dataset["expected_biological_units"], 77)
        self.assertEqual(dataset["biological_unit"], "participant")
        self.assertEqual(dataset["role"], "external_development")
        self.assertIs(dataset["admission_blocking"], True)
        text = " ".join(dataset["blockers"] + [dataset["notes"]])
        self.assertIn("77", text)
        self.assertIn("70", text)
        self.assertIn("no stage endpoint", text)
        self.assertIn(
            "gse249997", self.tasks["bulk_state_transfer"]["datasets_development"]
        )

    def test_conditional_bulk_power_freezes_user_specified_effect(self) -> None:
        claim_gate = self.tasks["bulk_state_transfer"]["claim_gate"]
        self.assertIn("0.05", claim_gate)
        self.assertIn("paired development variance", claim_gate)
        self.assertIn("at least 80% power", claim_gate)
        self.assertIn("downgrade", claim_gate)
        self.assertIn("Do not assume", claim_gate)
        power = self.inference["power"]
        self.assertEqual(power["prospective_target"], 0.80)
        self.assertEqual(
            power["effect_estimator"], "cross_fitted_or_shrunk_development_gain"
        )
        self.assertEqual(power["bulk_state_transfer_minimum_effect"], 0.05)
        self.assertIn("exactly_0.05", power["effect_size_rule"])
        self.assertEqual(power["underpowered_action"], "downgrade_without_unblinding")
        self.assertEqual(power["reachability_assumption"], "none")

    def test_borzoi_canonical_ensemble_is_official_and_same_split(self) -> None:
        artifact = next(
            artifact
            for artifact in self.lead_checkpoints["artifacts"]
            if artifact["artifact_id"] == "borzoi_calico_official_human_ensemble"
        )
        self.assertEqual(
            [member["path"] for member in artifact["members"]],
            [f"borzoi/f{index}/model0_best.h5" for index in range(4)],
        )
        self.assertEqual(
            [member["gcs_generation"] for member in artifact["members"]],
            [
                "1691539358988138",
                "1691539419241923",
                "1691539441571542",
                "1691539462862235",
            ],
        )
        self.assertIn("same_fold3_test_fold4_validation", artifact["member_relationship"])
        self.assertIn("not_cross_validation_folds", artifact["member_relationship"])
        self.assertEqual(artifact["weight_license"], "UNDECLARED_BLOCKED")
        self.assertEqual(self.models["borzoi_ensemble"]["checkpoint_sha256"], "UNRESOLVED")
        self.assertEqual(
            self.models["borzoi_ensemble"]["exposure_status"],
            "target_label_unexposed",
        )

        converted = next(
            candidate
            for candidate in self.lead_checkpoints["artifacts"]
            if candidate["artifact_id"] == "borzoi_genentech_grelu_conversion"
        )
        self.assertEqual(
            converted["ordered_composite_sha256"],
            "75bf9ee19c41575747850b6c5aae6836a206b04985fe3da0f39f4ed813db1921",
        )

    def test_spatial_and_pathology_assets_fail_closed(self) -> None:
        for dataset_id in {"gse292268", "gse312698", "vu2025"}:
            dataset = self.datasets[dataset_id]
            self.assertEqual(dataset["role"], "blocked")
            self.assertEqual(dataset["status"], "blocked")
            self.assertIs(dataset["admission_blocking"], True)

        for dataset_id in {"hra007511", "steatosite"}:
            dataset = self.datasets[dataset_id]
            self.assertEqual(dataset["role"], "deferred")
            self.assertEqual(dataset["status"], "deferred")
            self.assertIs(dataset["admission_blocking"], True)

        spatial_family = self.models_by_file["spatial_deferred"]
        self.assertEqual(spatial_family["status"], "deferred")
        self.assertEqual(
            {model["model_id"] for model in spatial_family["models"]},
            {
                "novae",
                "nichecompass",
                "spatialmeta",
                "spamosaic",
                "nicheformer",
                "graphst",
                "stagate",
                "banksy",
                "spatial_transcript_native_baselines",
            },
        )
        self.assertTrue(
            all(model["admission_blocking"] for model in spatial_family["models"])
        )
        self.assertNotIn("spatial", self.tasks)

    def test_spatial_audit_bundles_match_deferred_registry(self) -> None:
        expected_registry = {
            "novae": (
                "prism-oncology/novae@b71c2234c2d1c8e88981297dfdcb8b61b42c03c0:v1.1.1+prism-oncology/novae-human-0@b8c0a5d7612bac6bc719ab57ed3cd16ad814728c:model.safetensors",
                "7d1e4a051469d14e796994faa34a41556183e78757f5dc70a99dbd9f6a19ce87",
                "code_BSD-3-Clause_checkpoint_terms_UNRESOLVED",
            ),
            "nichecompass": (
                "Lotfollahi-lab/nichecompass@5bf7d2a09e70bd490c59e978e37f4df30f812994:0.3.3+PyPI:nichecompass==0.3.3",
                "NOT_APPLICABLE",
                "code_BSD-3-Clause_prior_and_training_asset_terms_source_specific",
            ),
            "spatialmeta": (
                "WanluLiuLab/SpatialMETA@909844e03c2ddffa331e8c45212830f8ffc85553+PyPI:spatialMETA==0.0.3.0+Zenodo:16750012:SpatialMETA-0.0.3.zip",
                "NOT_APPLICABLE",
                "top_level_BSD-3-Clause_vendored_STalign_GPL-3.0_patent_and_asset_terms_UNRESOLVED",
            ),
            "spamosaic": (
                "JinmiaoChenLab/SpaMosaic@cc1336755de8d1ddd1337b3cd67c9d41a61c2cb4:v1.0.3+PyPI:spamosaic==1.0.3",
                "NOT_APPLICABLE",
                "MIT",
            ),
            "nicheformer": (
                "theislab/nicheformer@485cadbc5caa15119adfd54228f8a8af835fcabc:0.0.1+theislab/Nicheformer@0ba4ba162f9ad3d6c8b713761954c31be2fa69ed:model.safetensors",
                "65fff02ac9bdd978cbb06b158a60d94ec59a96f30e706aed2c56d213d582e603",
                "source_BSD-3-Clause_HuggingFace_code_and_weights_MIT_training_assets_source_specific",
            ),
            "graphst": (
                "JinmiaoChenLab/GraphST@d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74+PyPI:GraphST==1.1.1",
                "NOT_APPLICABLE",
                "repository_AGPL-3.0_setup_and_PyPI_MIT_conflict_UNRESOLVED",
            ),
            "stagate": (
                "QIFEIDKN/STAGATE_pyG@ae1158ca8cf1eb6bb8ee198298552d44c9ac21db:1.0.0",
                "NOT_APPLICABLE",
                "MIT",
            ),
            "banksy": (
                "prabhakarlab/Banksy_py@9278996c39e376277d57ef95278000447ba6c57c:1.3.5",
                "NOT_APPLICABLE",
                "GPL-3.0_selected_Python_source",
            ),
            "spatial_transcript_native_baselines": (
                "PyPI:scanpy==1.12.3+squidpy==1.8.3+scikit-learn==1.9.0+leidenalg==0.12.0+igraph==1.0.0+Bioconductor:BayesSpace==1.22.0",
                "NOT_APPLICABLE",
                "BSD-3-Clause_components_plus_GPL_Leiden_and_igraph_plus_BayesSpace_MIT_plus_file",
            ),
        }

        loaded = {}
        for model_id, (revision, digest, license_status) in expected_registry.items():
            model = self.models[model_id]
            self.assertEqual(model["checkpoint_revision"], revision, model_id)
            self.assertEqual(model["checkpoint_sha256"], digest, model_id)
            self.assertEqual(model["license_status"], license_status, model_id)
            self.assertEqual(model["exposure_status"], "unknown", model_id)
            self.assertEqual(model["status"], "deferred", model_id)
            self.assertIs(model["admission_blocking"], True, model_id)
            self.assertEqual(model["supported_tasks"], [], model_id)

            audit_dir = CONFIG / "artifacts" / "models" / model_id
            checkpoints = json.loads(
                (audit_dir / "checkpoints.json").read_text(encoding="utf-8")
            )
            crosswalk_path = audit_dir / "development_crosswalk.json"
            crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
            exposure = json.loads(
                (audit_dir / "exposure_audit.json").read_text(encoding="utf-8")
            )
            loaded[model_id] = checkpoints
            self.assertEqual(
                checkpoints["activation_disposition"]["admission_status"],
                "deferred",
                model_id,
            )
            self.assertTrue(
                any(
                    key.startswith("v1_") and value == "prohibited"
                    for key, value in checkpoints["activation_disposition"].items()
                ),
                model_id,
            )
            self.assertEqual(
                exposure["development_crosswalk_record"]["sha256"],
                hashlib.sha256(crosswalk_path.read_bytes()).hexdigest(),
                model_id,
            )
            self.assertIn(
                "ineligible", exposure["sealed_champion_eligibility"], model_id
            )

        novae = loaded["novae"]
        self.assertEqual(
            novae["registered_checkpoint"]["sha256"],
            expected_registry["novae"][1],
        )
        self.assertIn(
            "Transductive", novae["held_donor_capability"]["native_zero_shot_lane"]
        )
        self.assertIn(
            "UNRESOLVED", novae["terms_and_release_contract"]["checkpoint"]
        )

        nichecompass = loaded["nichecompass"]
        self.assertEqual(
            nichecompass["checkpoint_contract"][
                "admitted_upstream_biological_checkpoint"
            ],
            "none",
        )
        self.assertEqual(
            nichecompass["spatial_graph_contract"]["cross_donor_edges"],
            "forbidden",
        )

        spatialmeta = loaded["spatialmeta"]
        self.assertIn("2025110969999", spatialmeta["terms_and_release_contract"]["patent"])
        self.assertIn(
            "adjacent_section", spatialmeta["input_and_missingness_contract"]["pairing"]
        )

        spamosaic = loaded["spamosaic"]
        self.assertIn(
            "transductive",
            spamosaic["held_donor_capability"]["query_mapping_lane"],
        )
        self.assertIn(
            "predeclared_transductive",
            spamosaic["graph_and_leakage_contract"]["cross_batch_mnn_edges"],
        )

        nicheformer = loaded["nicheformer"]
        self.assertEqual(
            nicheformer["registered_checkpoint"]["model_file"]["sha256"],
            expected_registry["nicheformer"][1],
        )
        self.assertEqual(nicheformer["checkpoint_code_alignment"]["status"], "UNRESOLVED")
        self.assertIn(
            "does not prove spatial structure transfer",
            nicheformer["held_donor_and_spatial_contract"]["base_limit"],
        )

        graphst = loaded["graphst"]
        self.assertIn("UNRESOLVED", graphst["source_terms_contract"]["disposition"])
        self.assertIn(
            "dense all-pairs", graphst["graph_and_section_contract"]["memory_risk"]
        )

        stagate = loaded["stagate"]
        self.assertIn(
            "prohibited", stagate["graph_and_topology_contract"]["three_dimensional_path"]
        )

        banksy = loaded["banksy"]
        self.assertEqual(
            banksy["terms_and_release_contract"]["selected_Python_source"],
            "GPL-3.0",
        )
        self.assertEqual(
            banksy["package_artifacts"]["admission_status"], "not_selected"
        )

        baselines = loaded["spatial_transcript_native_baselines"]
        self.assertIn("must beat", baselines["suite_contract"]["mandatory_comparison"])
        self.assertIn(
            "not an inductive predictor", baselines["suite_contract"]["unsupported_claim"]
        )

    def test_pairing_topology_prevents_false_pairs(self) -> None:
        self.assertEqual(self.datasets["gse296875"]["pairing_levels"], ["same_nucleus"])
        self.assertEqual(
            self.datasets["gse244832"]["pairing_levels"],
            ["same_sample_different_aliquot"],
        )
        self.assertNotIn("same_cell", self.datasets["gse244832"]["pairing_levels"])
        self.assertEqual(self.datasets["gse192741"]["pairing_levels"], ["same_section"])
        self.assertEqual(
            self.datasets["pxd051911"]["pairing_levels"],
            ["same_donor_different_tissue"],
        )
        bulk_task = self.tasks["bulk_state_transfer"]
        self.assertIn("pxd051911", bulk_task["datasets_development"])
        self.assertIn("protein_program_rank_transportability", bulk_task["metrics"])
        self.assertIn(
            "never paired RNA-protein learning", " ".join(bulk_task["admission_gates"])
        )
        self.assertIn("adjacent_section", self.datasets["hra007511"]["pairing_levels"])
        self.assertNotIn("same_cell", self.datasets["vu2025"]["pairing_levels"])

        semantics = self.resources["semantics"]
        self.assertEqual(set(semantics["pairing_levels"]), PAIRING_LEVELS)
        self.assertEqual(set(semantics["modality_statuses"]), MODALITY_STATUSES)
        self.assertEqual(
            set(semantics["exposure_statuses"]),
            {
                "clean_declared",
                "target_label_unexposed",
                "encoder_seen",
                "continual_seen",
                "reference_only",
                "downstream_demo",
                "unknown",
            },
        )
        self.assertIn("never a biological zero", semantics["missingness_policy"])

    def test_reference_paths_hashes_and_coordinate_contract_are_frozen(self) -> None:
        reference = self.resources["reference"]
        self.assertEqual(reference["bundle_id"], "grch38p14_gencode_v49")
        self.assertEqual(reference["assembly"], "GRCh38")
        self.assertEqual(reference["assembly_patch"], "p14")
        self.assertEqual(reference["annotation_release"], "GENCODE_v49")
        self.assertEqual(
            reference["genome"]["path"],
            "/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz",
        )
        self.assertEqual(
            reference["genome"]["sha256"],
            "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215",
        )
        self.assertEqual(
            reference["annotation"]["sha256"],
            "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
        )
        prohibited = set(self.resources["prohibited_references"]["paths"])
        self.assertNotIn(reference["genome"]["path"], prohibited)
        self.assertNotIn(reference["annotation"]["path"], prohibited)
        self.assertEqual(
            self.resources["coordinate_contract"]["bed_system"],
            "zero_based_half_open",
        )
        self.assertEqual(
            set(self.resources["coordinate_contract"]["analysis_contigs"]),
            {f"chr{index}" for index in range(1, 23)} | {"chrX", "chrY", "chrM"},
        )
        self.assertEqual(
            self.resources["coordinate_contract"]["excluded_contig_policy"],
            "reject",
        )
        self.assertIn(
            "native_build_before_liftover",
            self.resources["coordinate_contract"]["native_ld_policy"],
        )

    def test_task_dataset_and_split_references_resolve(self) -> None:
        for task_id, task in self.tasks.items():
            referenced_datasets = (
                task["datasets_train"]
                + task["datasets_development"]
                + task["datasets_sealed"]
            )
            self.assertTrue(set(referenced_datasets) <= set(self.datasets), task_id)
            self.assertIn(task["split_id"], self.splits, task_id)
            self.assertTrue(task["metrics"], task_id)
            self.assertTrue(task["admission_gates"], task_id)

        self.assertEqual(self.tasks["typed_evidence_graph"]["status"], "blocked")
        self.assertEqual(self.tasks["unified_model"]["status"], "blocked")

    def test_split_contract_uses_biological_units(self) -> None:
        required = {
            "schema_version",
            "split_id",
            "entity",
            "group_keys",
            "roles",
            "seed",
            "locus_grouping",
            "label_policy",
            "exposure_policy",
            "admission_gates",
        }
        for filename, split in self.splits_by_file.items():
            self.assertEqual(required - split.keys(), set(), filename)
            self.assertEqual(split["schema_version"], "masld-bench-split-v1")
            self.assertEqual(split["split_id"], filename)
            self.assertNotIn(split["entity"], {"cell", "spot", "section", "array", "run"})
            self.assertIn("withheld_sealed", split["roles"])
            self.assertEqual(split["seed"], 20260821)

        locus_controls = self.splits["locus_outer"]["locus_controls"]
        self.assertEqual(
            locus_controls,
            {
                "schema_version": "masld-bench-locus-split-controls-v1",
                "block_merge_rule": "within_distance_or_ld_threshold",
                "merge_distance_bp": 1_000_000,
                "ld_r2_threshold": 0.8,
                "ld_ancestry_rule": "any_represented_ancestry",
                "sequence_identity_grouping": [
                    "exact_forward_identity",
                    "exact_reverse_complement_identity",
                ],
                "boundary_buffer_policy": (
                    "exclude_each_outer_boundary_by_largest_admitted_model_receptive_field"
                ),
            },
        )

    def test_runtime_profiles_are_bounded_and_unresolved_runs_block(self) -> None:
        profiles = self.resources["profiles"]
        profile_required = {
            "partition",
            "cpus",
            "memory_gb",
            "wall_time",
            "gpus",
            "accelerator",
            "qos",
            "admission_blocking",
            "blockers",
        }
        self.assertIn("cpu_contract", profiles)
        self.assertEqual(profiles["gpu_single"]["accelerator"], "l40s")
        self.assertIs(profiles["gpu_single"]["admission_blocking"], False)
        for profile_id, profile in profiles.items():
            self.assertEqual(profile_required - profile.keys(), set(), profile_id)
            self.assertLessEqual(profile["cpus"], 16, profile_id)
            if profile["partition"] in {"cpu", "io", "gpu"}:
                self.assertLessEqual(profile["memory_gb"], 210, profile_id)
            if profile["partition"] == "bigmem":
                self.assertGreater(profile["memory_gb"], 210, profile_id)
            if profile["accelerator"] == "UNRESOLVED":
                self.assertIs(profile["admission_blocking"], True, profile_id)
                self.assertTrue(profile["blockers"], profile_id)

        runtime_required = {
            "schema_version",
            "runtime_id",
            "status",
            "purpose",
            "scheduler",
            "resource_profile",
            "environment_lock",
            "container_digest",
            "modules",
            "admission_blocking",
            "blockers",
        }
        for filename, runtime in self.runtimes_by_file.items():
            self.assertEqual(runtime_required - runtime.keys(), set(), filename)
            self.assertEqual(runtime["schema_version"], "masld-bench-runtime-v1")
            self.assertEqual(runtime["runtime_id"], filename)
            self.assertIn(runtime["resource_profile"], profiles)
            if "UNRESOLVED" in {
                runtime["environment_lock"], runtime["container_digest"]
            } or "UNRESOLVED" in runtime["modules"]:
                self.assertIs(runtime["admission_blocking"], True, filename)
                self.assertTrue(runtime["blockers"], filename)

        self.assertEqual(
            self.runtimes_by_file["cpu_baseline_smoke"]["modules"],
            ["python/3.11.5-GCCcore-13.2.0"],
        )


if __name__ == "__main__":
    unittest.main()
