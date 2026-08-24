from __future__ import annotations

import unittest

from masld_bench.context_preconditions import (
    ContextMapperPrecondition,
    ContextPreconditionError,
    ContextRow,
    ContextVector,
    apply_film_lora_precondition,
    require_precondition_only_campaign,
    stable_rank_quantile,
)


def vector(*values: float | None) -> ContextVector:
    return ContextVector(
        values=values,
        observed_mask=tuple(value is not None for value in values),
        states=tuple("observed" if value is not None else "structurally_missing" for value in values),
    )


def row(
    donor: str,
    block: str,
    *,
    role: str = "train",
    values: tuple[float | None, ...] = (1.0, 2.0, 3.0),
    dataset: str = "GSE296875",
    pairing: str = "same_nucleus",
    atac_input: bool = False,
) -> ContextRow:
    return ContextRow(
        row_id=f"{donor}:{block}",
        dataset_id=dataset,
        donor_id=donor,
        study_id=dataset,
        split_role=role,
        genomic_block=block,
        pairing_level=pairing,
        rna=vector(*values),
        observed_atac_input=atac_input,
    )


class ContextPreconditionTests(unittest.TestCase):
    def mapper(self, *, arm: str = "released_rank_quantile_from_counts") -> ContextMapperPrecondition:
        return ContextMapperPrecondition.fit(
            [row("train-1", "train-block")],
            ordered_gene_ids=("g1", "g2", "g3"),
            reference_values=(10.0, 20.0, 30.0),
            mapper_arm=arm,
            held_donor_ids=("held-1",),
            held_genomic_blocks=("held-block",),
            receptive_field_bp=524_288,
            boundary_buffer_bp=524_288,
        )

    def test_released_rank_mapper_matches_stable_upstream_contract(self) -> None:
        expected = stable_rank_quantile((3.0, 1.0, 1.0), (10.0, 20.0, 30.0))
        self.assertEqual(expected, (30.0, 10.0, 20.0))
        self.assertEqual(
            expected,
            stable_rank_quantile((30.0, 10.0, 10.0), (10.0, 20.0, 30.0)),
        )

    def test_fit_and_query_cross_donor_and_genomic_block(self) -> None:
        mapped = self.mapper().transform_query(
            row("held-1", "held-block", role="test", values=(3.0, 1.0, 2.0))
        )
        self.assertEqual(mapped, (30.0, 10.0, 20.0))
        with self.assertRaisesRegex(ContextPreconditionError, "cross both"):
            self.mapper().transform_query(
                row("train-1", "held-block", role="test")
            )

    def test_held_rows_and_inadequate_boundary_cannot_enter_fit(self) -> None:
        kwargs = dict(
            ordered_gene_ids=("g1", "g2", "g3"),
            reference_values=(10.0, 20.0, 30.0),
            mapper_arm="released_rank_quantile_from_counts",
            held_donor_ids=("held-1",),
            held_genomic_blocks=("held-block",),
            receptive_field_bp=524_288,
            boundary_buffer_bp=524_288,
        )
        with self.assertRaisesRegex(ContextPreconditionError, "nontraining"):
            ContextMapperPrecondition.fit(
                [row("held-1", "held-block", role="test")], **kwargs
            )
        kwargs["boundary_buffer_bp"] = 1000
        with self.assertRaisesRegex(ContextPreconditionError, "boundary buffer"):
            ContextMapperPrecondition.fit(
                [row("train-1", "train-block")], **kwargs
            )

    def test_missing_context_stays_null_and_routes_sequence_only(self) -> None:
        missing = row(
            "held-1", "held-block", role="test", values=(1.0, None, 3.0)
        )
        self.assertIsNone(self.mapper().transform_query(missing))
        self.assertIsNone(missing.rna.values[1])
        self.assertFalse(missing.rna.observed_mask[1])
        with self.assertRaisesRegex(ContextPreconditionError, "never zero"):
            ContextVector(
                values=(1.0, 0.0, 3.0),
                observed_mask=(True, False, True),
                states=("observed", "structurally_missing", "observed"),
            ).validate(3)

    def test_target_atac_sealed_outcomes_and_false_pairing_are_rejected(self) -> None:
        with self.assertRaisesRegex(ContextPreconditionError, "ATAC"):
            self.mapper().transform_query(
                row("held-1", "held-block", role="test", atac_input=True)
            )
        false_pair = row(
            "held-1",
            "held-block",
            role="test",
            dataset="GSE244832",
            pairing="same_nucleus",
        )
        with self.assertRaisesRegex(ContextPreconditionError, "same-donor"):
            self.mapper().transform_query(false_pair)

    def test_length_adjusted_arm_uses_reference_gene_lengths(self) -> None:
        mapped = self.mapper(
            arm="length_adjusted_tpm_then_released_rank_quantile"
        ).transform_query(
            row("held-1", "held-block", role="test", values=(10.0, 10.0, 10.0)),
            gene_lengths=(10.0, 5.0, 2.0),
        )
        self.assertEqual(mapped, (10.0, 20.0, 30.0))

    def test_neutral_film_and_zero_initialized_lora_are_identity(self) -> None:
        features = ((1.0, 2.0, 3.0), (4.0, 5.0, 6.0))
        observed = apply_film_lora_precondition(
            features,
            gamma_delta=(0.0, 0.0, 0.0),
            beta=(0.0, 0.0, 0.0),
            lora_a=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            lora_b=((0.0, 0.0), (0.0, 0.0), (0.0, 0.0)),
            alpha=2.0,
            context_observed=True,
            gate=1.0,
        )
        self.assertEqual(observed, features)
        missing = apply_film_lora_precondition(
            features,
            gamma_delta=(9.0, 9.0, 9.0),
            beta=(9.0, 9.0, 9.0),
            lora_a=((1.0, 0.0, 0.0),),
            lora_b=((1.0,), (1.0,), (1.0,)),
            alpha=1.0,
            context_observed=False,
            gate=0.0,
        )
        self.assertEqual(missing, features)

    def test_precondition_fixture_cannot_authorize_campaign(self) -> None:
        require_precondition_only_campaign({
            "wave": "conditional_model",
            "status": "blocked_trigger",
            "submit_enabled": False,
            "conditional_model_spec_path": "UNRESOLVED",
        })
        with self.assertRaises(ContextPreconditionError):
            require_precondition_only_campaign({
                "wave": "conditional_model",
                "status": "ready",
                "submit_enabled": True,
                "conditional_model_spec_path": "/tmp/spec",
            })


if __name__ == "__main__":
    unittest.main()
