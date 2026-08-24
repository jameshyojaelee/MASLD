from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from masld_bench.slurm import (
    ResourceProfile,
    SlurmError,
    parse_duration,
    render_sbatch,
    resource_totals,
)


class SlurmTests(unittest.TestCase):
    def test_partition_limits_and_memory_routing(self) -> None:
        self.assertEqual(parse_duration("3-00:00:00"), 259200)
        with self.assertRaises(SlurmError):
            ResourceProfile("bad", "io", 8, 64, "90:00:00").validate()
        with self.assertRaises(SlurmError):
            ResourceProfile("bad", "cpu", 16, 211, "90:00:00").validate()
        with self.assertRaises(SlurmError):
            ResourceProfile("bad", "cpu", 8, 32, "02:00:00", gpus=1).validate()
        with self.assertRaises(SlurmError):
            ResourceProfile("bad", "cpu", True, 32, "02:00:00").validate()
        with self.assertRaises(SlurmError):
            ResourceProfile("bad", "cpu", 1, 4, "00:00:00").validate()

    def test_exact_render_and_totals(self) -> None:
        profile = ResourceProfile("l40s_smoke", "gpu", 8, 96, "02:00:00", 1, "L40S")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = render_sbatch(
                job_name="fixture",
                profile=profile,
                command=("python", "-m", "fixture"),
                stdout_path=root / "logs" / "out",
                stderr_path=root / "logs" / "err",
            )
        self.assertIn("#SBATCH --partition=gpu", script)
        self.assertIn("#SBATCH --gres=gpu:l40s:1", script)
        self.assertNotIn("--array", script)
        totals = resource_totals({"l40s_smoke": profile}, ["l40s_smoke", "l40s_smoke"])
        self.assertEqual(totals["jobs"], 2)
        self.assertEqual(totals["gpu_hours_requested"], 4.0)
        self.assertEqual(totals["cpu_hours_requested"], 32.0)

    def test_directive_injection_and_lossy_profile_coercion_are_rejected(self) -> None:
        profile = ResourceProfile("cpu", "cpu", 1, 4, "01:00:00")
        with self.assertRaisesRegex(SlurmError, "stdout_path"):
            render_sbatch(
                job_name="fixture",
                profile=profile,
                command=("true",),
                stdout_path="relative.out\n#SBATCH --array=1-9",
                stderr_path="/tmp/err",
            )
        with self.assertRaisesRegex(SlurmError, "must be an integer"):
            ResourceProfile.from_mapping(
                "bad", {"partition": "cpu", "cpus": 1.5, "memory_gb": 4, "wall_time": "01:00:00"}
            )

    def test_serialized_profile_round_trip_is_exact(self) -> None:
        profile = ResourceProfile("cpu", "cpu", 2, 8, "01:00:00")
        self.assertEqual(
            ResourceProfile.from_serialized("cpu", profile.as_dict()),
            profile,
        )
        forged = profile.as_dict()
        forged["profile_id"] = "different"
        with self.assertRaisesRegex(SlurmError, "differs"):
            ResourceProfile.from_serialized("cpu", forged)


if __name__ == "__main__":
    unittest.main()
