from __future__ import annotations

import ast
import csv
import json
from pathlib import Path
import pickle
import tempfile
import unittest

from scripts import build_uce_frozen_screen_50000 as builder
from scripts import extract_uce_frozen_screen_50000 as extractor
from scripts import uce_frozen_screen_50000_common as common
from scripts import fit_predict_common_cell_heads_study_50000 as common_heads


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/artifacts/models/uce/frozen_screen_50000.json"
CHECKPOINTS = ROOT / "config/artifacts/models/uce/checkpoints.json"
CPU_WRAPPER = ROOT / "slurm/preflight_uce_frozen_screen_50000_cpu.sbatch"
GPU_WRAPPER = ROOT / "slurm/run_uce_frozen_screen_50000.sbatch"


class UCEFrozenScreen50000Test(unittest.TestCase):
    def _synthetic_inputs(self, root: Path, *, crossing_study: bool = False):
        import anndata
        import numpy as np
        import pandas as pd
        from scipy import sparse
        import torch

        row_ids = [f"cell-{index}" for index in range(5)]
        donors = [f"donor-{index}" for index in range(5)]
        studies = [f"study-{index}" for index in range(5)]
        obs = pd.DataFrame(
            {
                "row_id": row_ids,
                "donor_id": donors,
                "dataset_id": studies,
                "outer_fold": np.arange(5, dtype=np.int8),
                "assay": ["unknown"] * 5,
            },
            index=pd.Index(row_ids, name="row_id_index"),
        )
        ensembl = [f"ENSG{index:06d}" for index in range(1, 5)]
        symbols = ["GENEA", "GENEB", "GENEC", "GENED"]
        matrix = sparse.csr_matrix(
            np.asarray(
                [
                    [5, 2, 0, 0],
                    [0, 3, 4, 0],
                    [0, 0, 5, 1],
                    [2, 0, 0, 6],
                    [3, 1, 0, 0],
                ],
                dtype=np.float32,
            )
        )
        count_source = root / "counts.h5ad"
        anndata.AnnData(
            X=matrix,
            obs=obs,
            var=pd.DataFrame(index=pd.Index(ensembl, name="ensembl_id")),
        ).write_h5ad(count_source)
        symbol_source = root / "symbol_source.h5ad"
        symbol_var = pd.DataFrame(
            {"ensembl_id": ensembl, "source_feature_id": symbols},
            index=pd.Index(ensembl, name="feature_index"),
        )
        symbol_obs = pd.DataFrame(
            {"broad_label": ["must_not_be_read"] * 5}, index=pd.Index(row_ids)
        )
        anndata.AnnData(X=matrix, obs=symbol_obs, var=symbol_var).write_h5ad(
            symbol_source
        )
        split = root / "split.tsv"
        with split.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["row_id", "donor_id", "dataset", "outer_fold"])
            for index, (row_id, donor, study) in enumerate(
                zip(row_ids, donors, studies, strict=True)
            ):
                if crossing_study and index == 4:
                    study = studies[0]
                writer.writerow([row_id, donor, study, index])
        protein = root / "human.pt"
        torch.save(
            {symbol: torch.ones(4, dtype=torch.float32) * offset for offset, symbol in enumerate(symbols, 1)},
            protein,
        )
        offsets = root / "offsets.pkl"
        offsets.write_bytes(pickle.dumps({"human": 10}, protocol=4))
        chrom = root / "chrom.csv"
        chrom.write_text(
            "gene_symbol,chromosome,start,species\n"
            "GENEA,1,10,human\n"
            "GENEB,1,20,human\n"
            "GENEC,2,30,human\n"
            "GENED,2,40,human\n",
            encoding="utf-8",
        )
        return count_source, symbol_source, split, protein, offsets, chrom

    def test_exact_checkpoints_licenses_runtime_and_exposure_are_frozen(self) -> None:
        contract = common.load_activation_contract(CONTRACT, ROOT)
        checkpoints = json.loads(CHECKPOINTS.read_text(encoding="utf-8"))
        self.assertEqual(contract["licenses"]["code"], "MIT")
        self.assertEqual(contract["licenses"]["weights"], "CC-BY-4.0")
        self.assertEqual(contract["runtime"]["torch"], "2.5.1")
        self.assertEqual(
            contract["runtime"]["artifacts_sha256"],
            "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb",
        )
        self.assertEqual(
            contract["checkpoints"]["uce_4l"]["sha256"],
            "acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560",
        )
        self.assertEqual(
            contract["checkpoints"]["uce_33l"]["sha256"],
            "3f458726196308e171611ed28394b55865708749f82af68cfd7771d5dfee661e",
        )
        self.assertEqual(
            {record["sha256"] for record in checkpoints["weight_artifacts"].values()},
            {
                contract["checkpoints"]["uce_4l"]["sha256"],
                contract["checkpoints"]["uce_33l"]["sha256"],
            },
        )
        for record in contract["checkpoints"].values():
            self.assertEqual(record["exposure_status"], "target_label_unexposed")
            self.assertTrue(record["sealed_champion_eligible"])
        self.assertEqual(
            contract["exposure"]["known_encoder_seen_development_studies"],
            ["GSE185477", "Liver_Atlas"],
        )

    def test_builder_preserves_rows_and_never_reads_symbol_source_observations(self) -> None:
        import numpy as np

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = self._synthetic_inputs(root)
            output = root / "fixture"
            receipt = builder.build(
                *inputs,
                output,
                expected_rows=5,
                expected_donors=5,
                expected_studies=5,
                gene_min_cells=1,
                cell_min_genes=1,
            )
            self.assertTrue(receipt["row_roster_preserved"])
            self.assertEqual(receipt["cells_removed"], 0)
            self.assertEqual(receipt["cells_truncated"], 0)
            self.assertEqual(receipt["evaluation_label_columns_read"], [])
            self.assertEqual(receipt["symbol_source_observation_fields_read"], [])
            self.assertFalse(receipt["downstream_head_fit"])
            self.assertEqual(np.load(output / "tokens.npy").shape, (5, 1536))
            fixture_root = root / "campaign"
            fixture_root.mkdir()
            output.rename(fixture_root / "fixture")
            observed = extractor.load_fixture(
                fixture_root,
                expected_rows=5,
                expected_donors=5,
                expected_studies=5,
            )
            self.assertEqual(observed["row_ids"], [f"cell-{index}" for index in range(5)])

    def test_cell_seed_is_order_invariant_and_overlength_fails(self) -> None:
        import numpy as np

        counts = np.asarray([5.0, 2.0, 1.0])
        token_rows = np.asarray([10, 11, 12])
        chromosome_codes = np.asarray([0, 0, 1])
        starts = np.asarray([10, 20, 30])
        first = {
            row_id: common.tokenize_cell(
                counts,
                token_rows,
                chromosome_codes,
                starts,
                row_id=row_id,
                run_seed=20260824,
            )[0]
            for row_id in ("a", "b")
        }
        second = {
            row_id: common.tokenize_cell(
                counts,
                token_rows,
                chromosome_codes,
                starts,
                row_id=row_id,
                run_seed=20260824,
            )[0]
            for row_id in ("b", "a")
        }
        self.assertTrue(np.array_equal(first["a"], second["a"]))
        self.assertTrue(np.array_equal(first["b"], second["b"]))
        with self.assertRaisesRegex(common.UCEActivationError, "truncation is forbidden"):
            common.tokenize_cell(
                np.ones(600),
                np.arange(600) + 10,
                np.arange(600),
                np.arange(600),
                row_id="overlength",
                run_seed=20260824,
            )

    def test_study_firewall_rejects_cross_fold_study(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = self._synthetic_inputs(root, crossing_study=True)
            with self.assertRaisesRegex(common.UCEActivationError, "metadata join"):
                builder.build(
                    *inputs,
                    root / "fixture",
                    expected_rows=5,
                    expected_donors=5,
                    expected_studies=5,
                    gene_min_cells=1,
                    cell_min_genes=1,
                )

    def test_each_model_root_satisfies_hardened_common_head_identity(self) -> None:
        from masld_bench.artifacts import freeze_tree

        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            for model_id in ("uce_4l", "uce_33l"):
                root = temporary_root / model_id
                embedding_dir = root / "embeddings"
                embedding_dir.mkdir(parents=True)
                # Identity derivation is deliberately independent of loading the
                # large matrix; shard loading separately enforces exact NPZ fields.
                (embedding_dir / "common_embeddings.npz").write_bytes(b"identity fixture")
                record = contract["checkpoints"][model_id]
                (embedding_dir / "receipt.json").write_text(
                    json.dumps(
                        {
                            "model_id": model_id,
                            "checkpoint_sha256": record["sha256"],
                            "exposure_status": "target_label_unexposed",
                            "policies": {"common": "test"},
                            "rows": 50_000,
                            "embedding_width": 1280,
                            "evaluation_label_columns_read": [],
                            "sealed_outcomes_read": False,
                        }
                    ),
                    encoding="utf-8",
                )
                freeze_tree(
                    root,
                    {
                        "artifact_class": "uce_frozen_screen_raw_embeddings",
                        "model_id": model_id,
                        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
                        "development_rows": 50_000,
                        "exposure_status": "target_label_unexposed",
                        "sealed_champion_eligible": True,
                        "evaluation_labels_read": False,
                        "sealed_outcomes_read": False,
                    },
                )
                identity, path = common_heads.derive_embedding_identity(
                    root, Path("embeddings/common_embeddings.npz")
                )
                self.assertEqual(identity["model_id"], model_id)
                self.assertEqual(identity["embedding_width"], 1280)
                self.assertTrue(identity["sealed_champion_eligible"])
                self.assertEqual(path.name, "common_embeddings.npz")

    def test_activation_code_has_no_head_fitting_surface(self) -> None:
        forbidden_functions = {"fit", "predict", "build_head", "fit_head", "train_head"}
        for path in (
            ROOT / "scripts/uce_frozen_screen_50000_common.py",
            ROOT / "scripts/build_uce_frozen_screen_50000.py",
            ROOT / "scripts/extract_uce_frozen_screen_50000.py",
        ):
            tree = ast.parse(path.read_text(encoding="utf-8"))
            names = {
                node.name
                for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            }
            self.assertTrue(names.isdisjoint(forbidden_functions), path.name)
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("head_epochs", text)
            self.assertNotIn("optimizer.step", text)
        extraction = (ROOT / "scripts/extract_uce_frozen_screen_50000.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('"common_embeddings.npz"', extraction)
        self.assertIn('set(bundle.files) != {"embeddings", "outer_folds", "row_ids"}', extraction)

    def test_wrappers_are_generic_nslab_bundled_and_do_not_submit(self) -> None:
        cpu = CPU_WRAPPER.read_text(encoding="utf-8")
        gpu = GPU_WRAPPER.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --partition=cpu", cpu)
        self.assertIn("#SBATCH --cpus-per-task=8", cpu)
        self.assertIn("#SBATCH --mem=96G", cpu)
        self.assertIn("#SBATCH --time=12:00:00", cpu)
        self.assertIn("#SBATCH --partition=gpu", gpu)
        self.assertIn("#SBATCH --account=nslab", gpu)
        self.assertIn("#SBATCH --qos=nslab", gpu)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", gpu)
        self.assertIn("#SBATCH --cpus-per-task=12", gpu)
        self.assertIn("#SBATCH --mem=128G", gpu)
        self.assertIn("#SBATCH --time=120:00:00", gpu)
        self.assertIn('MODEL_IDS=${MODEL_IDS:-"uce_4l uce_33l"}', gpu)
        self.assertNotIn("innovation", gpu.lower())
        self.assertNotIn("#SBATCH --array", gpu)
        self.assertNotIn("#SBATCH --array", cpu)
        self.assertFalse(any(line.strip().startswith("sbatch ") for line in gpu.splitlines()))
        self.assertNotIn("common_cell_head", gpu)
        for text in (cpu, gpu):
            job_name = next(
                line.split("=", 1)[1] for line in text.splitlines() if "--job-name=" in line
            )
            self.assertNotIn("uce", job_name.lower())
            self.assertNotIn("masld", job_name.lower())


if __name__ == "__main__":
    unittest.main()
