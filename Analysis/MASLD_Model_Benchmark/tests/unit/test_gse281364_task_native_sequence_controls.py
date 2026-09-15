from __future__ import annotations

import torch
import unittest

from scripts.gse281364_task_native_sequence_controls import (
    MODEL_IDS,
    allele_pair_activity,
    build_model,
    joint_orientation_augmentation,
    reverse_complement_one_hot,
)


class GSE281364TaskNativeSequenceControlTests(unittest.TestCase):
    def test_reverse_complement_is_an_involution(self) -> None:
        values = torch.zeros((2, 4, 230), dtype=torch.float32)
        values[:, 0, :] = 1.0
        self.assertTrue(
            torch.equal(
                reverse_complement_one_hot(reverse_complement_one_hot(values)),
                values,
            )
        )

    def test_models_have_two_context_activity_outputs(self) -> None:
        values = torch.zeros((2, 4, 230), dtype=torch.float32)
        values[:, 0, :] = 1.0
        for model_id in MODEL_IDS:
            model = build_model(model_id, seed=1103)
            self.assertEqual(tuple(model(values).shape), (2, 2))

    def test_ref_alt_swap_negates_shared_encoder_delta(self) -> None:
        reference = torch.zeros((2, 4, 230), dtype=torch.float32)
        alternative = reference.clone()
        reference[:, 0, :] = 1.0
        alternative[:, 3, :] = 1.0
        for model_id in MODEL_IDS:
            model = build_model(model_id, seed=1103).eval()
            with torch.inference_mode():
                forward = allele_pair_activity(
                    model, reference, alternative, reverse_complement_average=True
                )["alt_minus_ref"]
                swapped = allele_pair_activity(
                    model, alternative, reference, reverse_complement_average=True
                )["alt_minus_ref"]
            self.assertTrue(torch.equal(forward, -swapped))

    def test_orientation_augmentation_is_joint(self) -> None:
        reference = torch.zeros((2, 4, 230), dtype=torch.float32)
        alternative = torch.ones((2, 4, 230), dtype=torch.float32)
        mask = torch.tensor([True, False], dtype=torch.bool)
        observed_ref, observed_alt = joint_orientation_augmentation(
            reference, alternative, mask
        )
        self.assertTrue(torch.equal(observed_ref[1], reference[1]))
        self.assertTrue(torch.equal(observed_alt[1], alternative[1]))
        self.assertTrue(
            torch.equal(observed_ref[0], reverse_complement_one_hot(reference)[0])
        )
        self.assertTrue(
            torch.equal(observed_alt[0], reverse_complement_one_hot(alternative)[0])
        )

    def test_five_seed_initializations_are_distinct(self) -> None:
        seeds = (1103, 2909, 4721, 6673, 8111)
        for model_id in MODEL_IDS:
            first_values = []
            for seed in seeds:
                model = build_model(model_id, seed=seed)
                first_values.append(next(model.parameters()).detach().flatten()[0].item())
            self.assertEqual(len(set(first_values)), 5)


if __name__ == "__main__":
    unittest.main()
