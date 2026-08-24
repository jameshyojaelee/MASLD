from __future__ import annotations

import json
import unittest

from scripts.audit_gse267145_authoritative_join import (
    GSE267145JoinError,
    build_join_rows,
    participant_fold,
)


def sample(series: str, accession: str, participant: str, *, h3: bool) -> dict[str, str]:
    characteristics = ["tissue: Liver"]
    title = participant
    if h3:
        title = f"{participant}_nor_H3K27ac_cutrun"
        characteristics = [
            "tissue: Liver",
            "steatosis: 0",
            "balloning: 0",
            "lobular inflammation: 0",
            "lobular necrosis: 0",
            "fibrosis nas_tidy: 0",
            "Stage: NOR",
            "Sex: F",
            "chip antibody: H3K27ac",
        ]
    return {
        "series": series,
        "accession": accession,
        "title": title,
        "source_name": "Liver",
        "candidate_participant_id": participant,
        "characteristics_json": json.dumps(characteristics),
        "data_processing_json": "[]",
        "supplementary_files_json": "[]",
    }


class GSE267145JoinTests(unittest.TestCase):
    def test_fold_is_stage_balanced_rank_modulo_five(self) -> None:
        self.assertEqual(participant_fold("ABC", "NOR", 7), 2)

    def test_join_rejects_any_axis_other_than_99(self) -> None:
        with self.assertRaises(GSE267145JoinError):
            build_join_rows(
                [sample("GSE267119", "GSM1", "ABC", h3=True)],
                [sample("GSE269412", "GSM2", "ABC", h3=False)],
            )


if __name__ == "__main__":
    unittest.main()
