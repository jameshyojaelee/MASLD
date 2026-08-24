from __future__ import annotations

import unittest

from masld_bench.reference import CoordinateContract, ReferenceError
from masld_bench.sequence import (
    centered_window,
    require_ref_allele,
    reverse_complement,
    validate_bed_interval,
    vcf_to_bed,
)


def coordinate_contract() -> CoordinateContract:
    return CoordinateContract(
        bed_system="zero_based_half_open",
        contig_policy="fixture_primary_only",
        analysis_contigs=("chr1", "chrX"),
        excluded_contig_policy="reject",
        variant_normalization="reference_match_then_left_normalize",
        liftover_policy="unique_only",
        native_ld_policy="native_first",
    )


class SequenceContractTests(unittest.TestCase):
    def test_vcf_to_bed_is_one_based_to_zero_based_half_open(self) -> None:
        interval = vcf_to_bed(
            "chr1",
            1,
            "AC",
            "A",
            coordinate_contract=coordinate_contract(),
            contig_lengths={"chr1": 10},
        )
        self.assertEqual((interval.contig, interval.start, interval.end), ("chr1", 0, 2))

    def test_bad_coordinate_and_excluded_contig_fail_closed(self) -> None:
        contract = coordinate_contract()
        with self.assertRaisesRegex(ReferenceError, "one-based"):
            vcf_to_bed("chr1", 0, "A", "G", coordinate_contract=contract)
        with self.assertRaisesRegex(ReferenceError, "excluded"):
            validate_bed_interval(
                "chr1_KI270706v1_random",
                0,
                1,
                coordinate_contract=contract,
            )

    def test_ref_alt_and_reverse_complement_checks(self) -> None:
        self.assertEqual(require_ref_allele("ac", "AC"), "AC")
        with self.assertRaisesRegex(ReferenceError, "REF allele mismatch"):
            require_ref_allele("AC", "AG")
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")
        self.assertEqual(reverse_complement(reverse_complement("AACGTN")), "AACGTN")
        with self.assertRaisesRegex(ReferenceError, "unsupported"):
            reverse_complement("ACGU")

    def test_window_boundary_and_exact_width(self) -> None:
        contract = coordinate_contract()
        interval = centered_window(
            "chr1",
            2,
            4,
            coordinate_contract=contract,
            contig_lengths={"chr1": 10},
        )
        self.assertEqual((interval.start, interval.end, interval.length), (0, 4, 4))
        with self.assertRaisesRegex(ReferenceError, "0 <= start"):
            centered_window(
                "chr1",
                1,
                4,
                coordinate_contract=contract,
                contig_lengths={"chr1": 10},
            )
        with self.assertRaisesRegex(ReferenceError, "boundary"):
            centered_window(
                "chr1",
                9,
                4,
                coordinate_contract=contract,
                contig_lengths={"chr1": 10},
            )


if __name__ == "__main__":
    unittest.main()
