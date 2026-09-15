"""A wrapper must not pre-create a path it later hands to a script as --output.

Scripts in this campaign create their own output directory and then refuse to
run if it already exists, which is how they guarantee they never overwrite an
immutable receipt. A wrapper that creates that same directory claims ownership
a second time and makes the guard fire on the wrapper's own scaffolding. The
job then cannot succeed on any run, on any node.
"""

from __future__ import annotations

from pathlib import Path
import re
import unittest


SLURM = Path(__file__).resolve().parents[2] / "slurm"
CREATES = re.compile(r'^\s*(?:install -d[^"]*|mkdir(?: -p)?)\s+"([^"]+)"', re.M)
OUTPUT = re.compile(r'--output(?:=|\s+)"([^"]+)"')


MKDIR = re.compile(r"output\.mkdir\(([^)]*)\)")


def script_owns_its_output(script: Path) -> bool:
    """True if the script insists on creating its own output directory.

    Pre-creating is only a defect when the script claims ownership: it either
    makes the directory without ``exist_ok`` or refuses to run when it already
    exists. A script that uses ``exist_ok=True`` and has no such guard is
    indifferent, and flagging it would make this test noise -- and a test
    people learn to ignore has become vacuous by a slower route.
    """

    try:
        source = script.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    made = MKDIR.search(source)
    tolerant = bool(made) and "exist_ok" in made.group(1)
    guarded = "output.exists()" in source
    return (bool(made) and not tolerant) or guarded


def conflicts(text: str, root: Path | None = None) -> list[str]:
    """Paths this wrapper both creates and passes as --output.

    Only counted when the invoked script owns its output directory.
    """

    if root is not None:
        scripts = [root / m for m in re.findall(r"(scripts/[\w/\-]+\.py)", text)]
        if not any(script_owns_its_output(s) for s in scripts):
            return []
    created = {m.group(1) for m in CREATES.finditer(text)}
    # A wrapper that creates a path and then removes it again before handing it
    # over has not claimed ownership -- that pair is a working idiom, not a
    # defect. Only an un-undone creation conflicts.
    removed = {m.group(1) for m in re.finditer(r'^\s*rm(?:dir| -rf?)\s+"([^"]+)"', text, re.M)}
    created -= removed
    # An sbatch --output=... directive is a log file, not a script argument.
    passed = {
        m.group(1)
        for m in OUTPUT.finditer(text)
        if not text[: m.start()].rstrip().endswith("#SBATCH")
    }
    return sorted(created & passed)


class WrapperOutputOwnershipTests(unittest.TestCase):
    def test_no_wrapper_precreates_a_path_it_passes_as_output(self) -> None:
        offenders = {}
        for wrapper in sorted(SLURM.glob("*.sbatch")):
            found = conflicts(
                wrapper.read_text(encoding="utf-8", errors="replace"),
                root=SLURM.parent,
            )
            if found:
                offenders[wrapper.name] = found
        self.assertEqual(
            offenders,
            {},
            "these wrappers create a directory they then pass as --output; the "
            "script owns that directory and will refuse to run",
        )

    def test_the_detector_fires_on_the_shape_it_guards(self) -> None:
        """The bug this was written for, as a fixture."""

        bad = 'install -d -m 0750 "${STAGE}/audit"\npython x.py --output "${STAGE}/audit" \\\n'
        self.assertEqual(conflicts(bad), ["${STAGE}/audit"])

    def test_creating_only_the_parent_is_allowed(self) -> None:
        good = 'install -d -m 0750 "${STAGE}"\npython x.py --output "${STAGE}/audit" \\\n'
        self.assertEqual(conflicts(good), [])

    def test_an_sbatch_output_directive_is_not_a_script_argument(self) -> None:
        directive = '#SBATCH --output=/logs/%x-%j.out\ninstall -d "${STAGE}"\n'
        self.assertEqual(conflicts(directive), [])


if __name__ == "__main__":
    unittest.main()
