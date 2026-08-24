from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

import pyfaidx

from scripts import filter_training_peaks as peaks


class TrainingPeakFilterTests(unittest.TestCase):
    def test_jittered_sequence_and_signal_windows_are_centered_on_summit(self) -> None:
        sequence, signal = peaks.peak_windows(1000, 250)
        self.assertEqual(sequence, (-307, 2807))
        self.assertEqual(signal, (250, 2250))

    def test_half_open_blacklist_overlap(self) -> None:
        intervals = peaks.merge_intervals(((100, 200), (180, 250), (300, 400)))
        self.assertEqual(intervals, ((100, 250), (300, 400)))
        self.assertTrue(peaks.overlaps(intervals, 249, 251))
        self.assertFalse(peaks.overlaps(intervals, 250, 300))

    def test_peak_selection_uses_only_three_genomic_training_folds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "fixture.fa"
            with fasta.open("x", encoding="utf-8") as handle:
                for fold in range(5):
                    handle.write(f">chr{fold + 1}\n{'ACGT' * 1250}\n")
            reference = pyfaidx.Fasta(str(fasta), rebuild=True)
            reference.close()
            narrowpeak = root / "fixture.narrowPeak"
            with narrowpeak.open("x", encoding="utf-8") as handle:
                for fold in range(5):
                    handle.write(
                        f"chr{fold + 1}\t2490\t2510\tpeak_{fold}\t100\t.\t10\t5\t5\t10\n"
                    )
            blacklist = root / "blacklist.bed.gz"
            with gzip.open(blacklist, "wt", encoding="utf-8"):
                pass
            output = root / "selected.bed"
            summary_path = root / "summary.json"
            result = peaks.build(
                narrowpeak=narrowpeak,
                blacklist=blacklist,
                fasta=fasta,
                fold_by_contig={f"chr{fold + 1}": fold for fold in range(5)},
                genomic_test_fold=4,
                genomic_valid_fold=3,
                output=output,
                summary_path=summary_path,
                peaks_per_allowed_fold=1,
                minimum_total_peaks=1,
            )
            self.assertEqual(result["held_genomic_test"], 1)
            self.assertEqual(result["held_genomic_valid"], 1)
            self.assertEqual(result["allowed_genomic_folds"], [0, 1, 2])
            self.assertEqual(result["retained"], 3)
            self.assertEqual(sum(1 for _line in output.open()), 3)


if __name__ == "__main__":
    unittest.main()
