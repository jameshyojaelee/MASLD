from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from scripts.plan_gse296875_observed_multiome_masks import (
    ObservedMultiomeMaskPlanError,
    build_plan,
    read_identifier_axes,
    validate_config,
    write_plan,
)
from masld_bench.adapters.rna_atac_classical import fold_index


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_mask_plan_20260825.json"


def synthetic_axes() -> dict[str, list[object]]:
    donors_by_fold: dict[int, str] = {}
    index = 0
    while len(donors_by_fold) < 5:
        donor = f"synthetic_donor_{index}"
        donors_by_fold.setdefault(
            fold_index(donor, seed=20260821, outer_folds=5), donor
        )
        index += 1
    donors: list[str] = []
    lineages: list[str] = []
    for donor in donors_by_fold.values():
        for lineage in ("hepatocyte", "macrophage"):
            donors.append(donor)
            lineages.append(lineage)
    peak_ids: list[str] = []
    chromosomes: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    for chromosome_index, chromosome in enumerate(
        ("chr1", "chr2", "chr3", "chr4", "chr5")
    ):
        for peak_index in range(8):
            peak_ids.append(f"peak_{chromosome}_{peak_index}")
            chromosomes.append(chromosome)
            start = chromosome_index * 100000 + peak_index * 1000
            starts.append(start)
            ends.append(start + 500)
    return {
        "donors": donors,
        "lineages": lineages,
        "peak_ids": peak_ids,
        "chromosomes": chromosomes,
        "starts": starts,
        "ends": ends,
    }


class GSE296875ObservedMultiomeMaskPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_config_passes_with_closed_matrix_firewall(self) -> None:
        resolved = validate_config(ROOT, deepcopy(self.config))
        self.assertTrue(resolved["input_h5"].is_file())
        self.assertFalse(self.config["firewall"]["atac_count_values_read"])
        self.assertFalse(
            self.config["disposition"]["biological_matrix_execution_authorized"]
        )

    def test_identifier_reader_does_not_require_count_groups(self) -> None:
        import h5py

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "identifiers_only.h5"
            with h5py.File(path, "x") as handle:
                handle.attrs["schema_version"] = "masld-bench-multimodal-h5-v1"
                handle.attrs["view_id"] = "gse296875_rna_atac_smoke_1000_v1"
                handle.attrs["pairing_level"] = "same_nucleus"
                handle.attrs["bed_coordinate_system"] = "0_based_half_open"
                obs = handle.create_group("obs")
                obs.create_dataset("donor_id", data=[b"donor1", b"donor2"])
                obs.create_dataset("broad_label", data=[b"hepatocyte", b"macrophage"])
                atac = handle.create_group("atac")
                atac.create_dataset("peak_id", data=[b"p1", b"p2"])
                atac.create_dataset("chromosome", data=[b"chr1", b"chr2"])
                atac.create_dataset("bed_start_0based", data=[0, 100])
                atac.create_dataset("bed_end_half_open", data=[50, 150])
            axes = read_identifier_axes(path)
            self.assertEqual(axes["donors"], ["donor1", "donor2"])
            self.assertNotIn("counts_csr", axes)

    def test_whole_chromosome_masks_are_disjoint_and_raw_ids_stay_private(self) -> None:
        axes = synthetic_axes()
        plan = build_plan(
            axes,
            donor_folds=5,
            donor_seed=20260821,
            genomic_folds=5,
            genomic_seed=20260825,
            canonical_contigs=("chr1", "chr2", "chr3", "chr4", "chr5"),
            targets_per_fold=2,
            input_peaks_per_fold=3,
            donor_namespace="synthetic:donor",
            peak_namespace="synthetic:peak",
        )
        targets: dict[int, set[str]] = defaultdict(set)
        inputs: dict[int, set[str]] = defaultdict(set)
        for row in plan["target_rows"]:
            targets[int(row["genomic_fold"])].add(str(row["chromosome"]))
        for row in plan["input_rows"]:
            inputs[int(row["held_genomic_fold"])].add(str(row["chromosome"]))
        for fold in range(5):
            self.assertTrue(targets[fold].isdisjoint(inputs[fold]))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            write_plan(output, plan)
            tokens = {
                token
                for path in output.glob("*.tsv")
                for line in path.read_text().splitlines()
                for token in line.split("\t")
            }
            self.assertTrue(set(plan["raw_donor_ids"]).isdisjoint(tokens))

    def test_open_count_firewall_is_rejected(self) -> None:
        config = deepcopy(self.config)
        config["firewall"]["atac_count_values_read"] = True
        with self.assertRaises(ObservedMultiomeMaskPlanError):
            validate_config(ROOT, config)

    def test_zero_peak_canonical_contig_is_retained(self) -> None:
        axes = synthetic_axes()
        plan = build_plan(
            axes,
            donor_folds=5,
            donor_seed=20260821,
            genomic_folds=5,
            genomic_seed=20260825,
            canonical_contigs=("chr1", "chr2", "chr3", "chr4", "chr5", "chrY"),
            targets_per_fold=2,
            input_peaks_per_fold=3,
            donor_namespace="synthetic:donor",
            peak_namespace="synthetic:peak",
        )
        chr_y = [row for row in plan["genomic_rows"] if row["chromosome"] == "chrY"]
        self.assertEqual(chr_y[0]["eligible_peaks"], 0)


if __name__ == "__main__":
    unittest.main()
