from collections import defaultdict

from sequence_overlap_patch import SequenceReference, Transcript, sequence_far_enough


def row(*intervals):
    return {"_actual_intervals": tuple(intervals)}


def test_plus_projection():
    tx = Transcript("g", "chr1", "+", [(101, 110), (201, 220)])
    assert tx.genomic_intervals(7, 6) == (("chr1", 108, 110), ("chr1", 201, 203))


def test_minus_projection():
    tx = Transcript("g", "chr1", "-", [(101, 110), (201, 220)])
    assert tx.genomic_intervals(18, 6) == (("chr1", 201, 202), ("chr1", 107, 110))


def test_rejects_one_nt_overlap():
    assert not sequence_far_enough(row(("chr1", 122, 144)), [row(("chr1", 100, 122))], 23)


def test_accepts_adjacent_sites():
    assert sequence_far_enough(row(("chr1", 123, 145)), [row(("chr1", 100, 122))], 23)


def test_rejects_same_site_across_isoforms():
    assert not sequence_far_enough(row(("chr2", 500, 522)), [row(("chr2", 500, 522))], 23)


def test_different_chromosomes():
    assert sequence_far_enough(row(("chr2", 100, 122)), [row(("chr1", 100, 122))], 23)


def test_unmapped_fails_closed():
    assert not sequence_far_enough({}, [], 23)


def test_isoform_count_comes_from_exact_sequence_matches():
    ref = SequenceReference.__new__(SequenceReference)
    ref.by_gene = defaultdict(list, {"g": ["tx1", "tx2", "tx3"]})
    ref.transcripts = {
        "tx1": Transcript("g", "chr1", "+", [(101, 150)], "A" * 10 + "C" * 23 + "A" * 17),
        "tx2": Transcript("g", "chr1", "+", [(201, 250)], "A" * 5 + "C" * 23 + "A" * 22),
        "tx3": Transcript("g", "chr1", "+", [(301, 350)], "A" * 50),
    }
    ref.candidate_rows = ref.candidate_rows_located = ref.candidate_rows_unlocated = 0
    rows = [{"gene_id_base": "g", "target_seq": "C" * 23,
             "n_isoforms_targeted": 99}]
    found = ref.annotate_rows(rows)
    assert len(found) == 1
    assert found[0]["n_isoforms_targeted"] == 2
