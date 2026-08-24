from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "scglue_fit_predict_retrieval_smoke.py"
SPEC = importlib.util.spec_from_file_location("scglue_retrieval_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
retrieval = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(retrieval)


class SCGLUERetrievalTests(unittest.TestCase):
    def test_model_roster_contains_family_native_and_linear_baseline(self) -> None:
        self.assertEqual(retrieval.MODELS, ("scglue", "paired_scglue", "linear_cca"))

    def test_hidden_pair_map_is_not_an_adapter_argument(self) -> None:
        parameter_names = retrieval.run.__code__.co_varnames[
            : retrieval.run.__code__.co_argcount
        ]
        self.assertEqual(parameter_names, ("inputs", "output", "seed"))


if __name__ == "__main__":
    unittest.main()
