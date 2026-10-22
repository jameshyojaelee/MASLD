"""Exercise the workflow's actual SAF command and fragment boundary counts."""
import csv
import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


WORKFLOW = Path(__file__).resolve().parents[1] / 'Snakefile'


def convert_bed(bed):
    source = Path(os.environ.get('MASLD_ATAC_WORKFLOW', WORKFLOW)).read_text()
    command = re.search(r"awk '([^']+)' \{input.peaks\}", source).group(1)
    # Snakemake formats braces after Python decodes the shell string.
    command = command.replace('{{', '{').replace('}}', '}').replace('\\\\', '\\')
    return subprocess.check_output(['awk', command], input=bed, text=True)


class BedToSafTests(unittest.TestCase):
    def test_coordinates_and_interval_lengths(self):
        saf = list(csv.reader(convert_bed('chr1\t100\t200\nchr2\t0\t1\n').splitlines(), delimiter='\t'))
        self.assertEqual(saf, [['peak_1', 'chr1', '101', '200', '.'],
                               ['peak_2', 'chr2', '1', '1', '.']])
        self.assertEqual([int(row[3]) - int(row[2]) + 1 for row in saf], [100, 1])

    def test_fragment_boundaries_with_featurecounts(self):
        binary = os.environ.get('FEATURECOUNTS') or shutil.which('featureCounts')
        self.assertTrue(binary, 'Set FEATURECOUNTS to the installed featureCounts executable')
        with tempfile.TemporaryDirectory(prefix='masld_atac_boundary_') as folder:
            folder = Path(folder)
            annotation = folder / 'peaks.saf'
            annotation.write_text(convert_bed('chr1\t100\t200\n'))
            sam = folder / 'pairs.sam'
            records = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:chr1\tLN:1000']
            # One pair ends at BED's excluded left boundary, two overlap,
            # and one begins just after the included right boundary.
            positions = [(71, 91), (111, 131), (181, 191), (201, 221)]
            for index, (left, right) in enumerate(positions):
                length = right + 10 - left
                records += [f'r{index}\t99\tchr1\t{left}\t60\t10M\t=\t{right}\t{length}\tAAAAAAAAAA\tIIIIIIIIII',
                            f'r{index}\t147\tchr1\t{right}\t60\t10M\t=\t{left}\t{-length}\tAAAAAAAAAA\tIIIIIIIIII']
            sam.write_text('\n'.join(records) + '\n')
            counts = folder / 'counts.txt'
            subprocess.run([binary, '-T', '1', '-p', '--countReadPairs', '-F', 'SAF',
                            '-a', str(annotation), '-o', str(counts), str(sam)],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            with counts.open() as handle:
                rows = list(csv.reader((line for line in handle if not line.startswith('#')), delimiter='\t'))
            # Independent BED overlap oracle on each mate's aligned interval.
            expected = sum(any(start - 1 < 200 and start - 1 + 10 > 100 for start in pair)
                           for pair in positions)
            self.assertEqual(expected, 2)
            self.assertEqual(int(rows[1][6]), expected)
            self.assertEqual(int(rows[1][5]), 100)


if __name__ == '__main__':
    unittest.main()
