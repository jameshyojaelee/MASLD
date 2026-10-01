#!/usr/bin/env python3
"""Check the user's compute amendment without submitting a scheduler job."""
import tempfile
from pathlib import Path
import unittest
from campaign_budget import apply_authorized_limits, check_compute_limits, parse_header


class ResourceAmendmentTests(unittest.TestCase):
    def test_amendment_retains_storage_history_and_accounting(self):
        old = dict(gpu_hours=300, cpu_core_hours=1500, simultaneous_gpus=4, new_bytes=500_000_000_000)
        journal = dict(limits=old, jobs=[dict(state="RUNNING", gpus=8, gpu_hours=1000, cpu_core_hours=5000)])
        apply_authorized_limits(journal)
        self.assertEqual(journal["resource_authorizations"][0]["previous_limits"], old)
        self.assertEqual(journal["limits"]["new_bytes"], old["new_bytes"])
        check_compute_limits(journal, dict(gpus=8, gpu_hours=1000, cpu_core_hours=5000))
        self.assertEqual(len(journal["jobs"]), 1)
        apply_authorized_limits(journal)
        self.assertEqual(len(journal["resource_authorizations"]), 1)

    def test_explicit_hardware_and_io_wall_limit(self):
        base = "#SBATCH --account=nslab\n#SBATCH --qos=nslab\n#SBATCH --cpus-per-task=4\n#SBATCH --mem=32G\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.sbatch"
            path.write_text(base + "#SBATCH --partition=gpu\n#SBATCH --gres=gpu:l40s:8\n#SBATCH --time=06:00:00\n")
            _, resource = parse_header(path)
            self.assertEqual(resource["gpus"], 8)
            self.assertEqual(resource["gpu_hours"], 48)
            path.write_text(base + "#SBATCH --partition=io\n#SBATCH --time=72:00:00\n")
            self.assertEqual(parse_header(path)[1]["gpus"], 0)
            path.write_text(base + "#SBATCH --partition=io\n#SBATCH --time=90:00:00\n")
            with self.assertRaises(ValueError):
                parse_header(path)


if __name__ == "__main__":
    unittest.main()
