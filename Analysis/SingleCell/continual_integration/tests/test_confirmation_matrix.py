from __future__ import annotations

import unittest

from masld_cl.confirmation_matrix import _expected_keys, confirmation_settings


class TestConfirmationMatrix(unittest.TestCase):
    def test_duplicate_roles_share_one_fit_matrix(self):
        selection = {
            "selected": {"ewc_lambda": 1.0, "replay_fraction": 0.2},
            "pareto_runner_up": {"ewc_lambda": 100.0, "replay_fraction": 0.2},
        }
        settings = confirmation_settings(selection)
        self.assertEqual(len(settings), 2)
        self.assertEqual(settings[(100.0, 0.2)], ["paper", "pareto_runner_up"])
        config = {"screen": {"confirmation_seeds": [17, 41]}, "lineages": ["H", "M"]}
        self.assertEqual(len(_expected_keys(config, settings)), 12)


if __name__ == "__main__":
    unittest.main()
