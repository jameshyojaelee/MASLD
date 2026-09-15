"""
Unit tests for the GSE202379 SAF recovery artifact.

Every expected value below is written by hand from the deposited GEO record.
Nothing here is generated from the artifact it validates.

The gate is deliberately NOT "agree with 343 everywhere". Agreement with a
harvest that is known to be defective cannot be a correctness criterion at the
donors where that harvest is the thing that is wrong. The gate is:
SAF-derived F equals 343's F on every SAF-graded donor EXCEPT where 343's
label-map pass overrides a present SAF string, with the exception set
enumerated below.
"""
from __future__ import annotations
import importlib.util, json, os
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
H343 = ROOT/"Analysis/SingleCell/scripts/343_harvest_documented_fstage.py"
ART = Path(os.environ["SAF_ARTIFACT_DIR"])

# ---- hand-pinned expectations -------------------------------------------
N_GSM                 = 59
N_DONORS              = 47
N_SAF_GRADED          = 40
SAMPLES_NO_LIVER_LOBE = ["GSM6112243", "GSM6112244"]      # P30, P98
F_EXCEPTION_SET       = ["P98"]                            # overwritten real grade
AGREES_BY_LUCK        = ["P30"]                            # same short-circuit, F0 either way
INVENTED_NO_SAF       = 7                                  # 5 'end stage' + 2 'healthy control'
EXPR_JOIN_N           = 39                                 # 40 - P70
CCC_JOIN_N            = 37                                 # 39 - P62 - P67
CCC_DROPS             = ["P62", "P67", "P70"]

MARGINALS_N40 = {
    "steatosis": {"0": 1, "1": 17, "2": 18, "3": 4},
    "activity":  {"0": 1, "1": 8,  "2": 4,  "3": 18, "4": 9},
    "fibrosis":  {"0": 3, "1": 9,  "2": 12, "3": 12, "4": 4},
}
MARGINALS_N38 = {
    "steatosis": {"1": 16, "2": 18, "3": 4},
    "activity":  {"0": 1,  "1": 6,  "2": 4,  "3": 18, "4": 9},
    "fibrosis":  {"0": 2,  "1": 8,  "2": 12, "3": 12, "4": 4},
}


@pytest.fixture(scope="module")
def donor():
    return pd.read_csv(ART/"donor_saf_grades.tsv", sep="\t")

@pytest.fixture(scope="module")
def marg():
    return json.loads((ART/"marginals.json").read_text())

@pytest.fixture(scope="module")
def defect():
    return json.loads((ART/"label_map_defect.json").read_text())

@pytest.fixture(scope="module")
def h343():
    spec = importlib.util.spec_from_file_location("h343", H343)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---- cardinality ---------------------------------------------------------
def test_donor_cardinality(donor):
    assert len(donor) == N_DONORS
    assert donor["donor_id"].nunique() == N_DONORS
    assert int(donor["n_gsm"].sum()) == N_GSM          # rows vs distinct donors
    assert int(donor["saf_recovered"].sum()) == N_SAF_GRADED


# ---- structural: the thing that is actually wrong with the GEO record ----
def test_two_samples_lack_liver_lobe_field(donor):
    """Replaces the earlier 38->40 'repair' assertion. There is no repair:
    P30/P98 simply lack 'liver lobe', shifting later fields up one row.
    Key-based parsing recovers all 40 directly."""
    flagged = sorted(donor.loc[donor["lacks_liver_lobe_field"], "donor_id"])
    assert flagged == ["P30", "P98"]
    gsms = sorted(
        g for row in donor.loc[donor["lacks_liver_lobe_field"], "gsm_ids"]
        for g in str(row).split(";"))
    assert gsms == SAMPLES_NO_LIVER_LOBE

def test_all_forty_saf_strings_recovered_by_key(donor):
    g = donor[donor["saf_recovered"]]
    assert len(g) == N_SAF_GRADED
    assert g["steatosis_S"].notna().all()
    assert g["activity_A"].notna().all()
    assert g["fibrosis_F"].notna().all()


# ---- marginals: BOTH sets pinned ----------------------------------------
def test_marginals_n40(marg):
    m = marg["n40_all_saf_graded"]
    for k, exp in MARGINALS_N40.items():
        assert m[k]["n"] == 40, k
        assert m[k]["levels"] == exp, k

def test_marginals_n38(marg):
    m = marg["n38_excluding_P30_P98"]
    for k, exp in MARGINALS_N38.items():
        assert m[k]["n"] == 38, k
        assert m[k]["levels"] == exp, k

def test_steatosis_has_four_levels_at_n40(marg):
    """Consequence that breaks downstream code assuming a 1-3 scale."""
    assert set(marg["n40_all_saf_graded"]["steatosis"]["levels"]) == {"0","1","2","3"}
    assert set(marg["n38_excluding_P30_P98"]["steatosis"]["levels"]) == {"1","2","3"}


# ---- THE GATE ------------------------------------------------------------
# NOTE 2026-08-30: these three tests validate the SAF recovery ARTIFACT, whose
# F_343 column was computed with the PRE-correction extract_fstage(). They are
# pinned to that historical artifact and still describe it correctly. If the
# recovery is ever re-run against the corrected source, F_343 for P98 becomes
# 1, the exception set becomes empty, and these expectations must be restated
# as "no disagreement remains, because the defect was fixed at source".
def test_f_agreement_gate_exception_set_is_exactly_P98(donor):
    g = donor[donor["saf_recovered"]]
    disagree = sorted(g.loc[g["f_agrees_with_343"] == False, "donor_id"])
    assert disagree == F_EXCEPTION_SET

def test_every_exception_is_a_label_map_override(donor):
    """The exception set is allowed ONLY because the label-map overrode a
    present SAF string. Any other cause is a real parse bug."""
    g = donor[donor["saf_recovered"]]
    exc = g[g["f_agrees_with_343"] == False]
    assert bool(exc["label_map_fired"].all())
    assert bool(exc["saf_recovered"].all())


# ---- THE MECHANISM (the cause, not the symptom) -------------------------
def test_regex_pass_now_precedes_label_map_pass(h343):
    """CORRECTED 2026-08-30. This test used to assert the defect; it now
    asserts the fix. A row carrying BOTH a mapped disease-status label AND a
    real SAF string must return the MEASURED value. Under the pre-fix source
    the first assertion returns 0 and this test fails."""
    both = pd.Series({"char_0": "disease status: Healthy control",
                      "char_1": "saf score: S1A1F1"})
    assert h343.extract_fstage(both) == 1        # measurement wins

    saf_only = pd.Series({"char_0": "disease status: NASH w/o cirrhosis",
                          "char_1": "saf score: S1A1F1"})
    assert h343.extract_fstage(saf_only) == 1    # unchanged where no label fires

def test_end_stage_label_no_longer_short_circuits(h343):
    both = pd.Series({"char_0": "disease status: end stage",
                      "char_1": "saf score: S2A3F2"})
    assert h343.extract_fstage(both) == 2        # pre-fix: 4

def test_label_map_survives_as_the_last_resort(h343):
    """Demotion, not deletion: a label with no measurement still supplies the
    only available value for the 7 donors that have nothing else."""
    assert h343.extract_fstage(pd.Series({"char_0": "saf score: end stage"})) == 4
    assert h343.extract_fstage(pd.Series({"char_0": "saf score: healthy control"})) == 0

def test_P30_agrees_by_luck_and_is_still_flagged(donor):
    """A disagreement-only audit would miss P30. It hits the identical
    short-circuit; its SAF is F0 so the wrong path yields the right answer."""
    r = donor.set_index("donor_id").loc["P30"]
    assert bool(r["label_map_fired"])
    assert bool(r["saf_recovered"])
    assert r["fibrosis_F"] == 0 and r["F_343"] == 0
    assert bool(r["f_agrees_with_343"])
    luck = donor[(donor["saf_recovered"]) & (donor["label_map_fired"])
                 & (donor["f_agrees_with_343"] == True)]
    assert sorted(luck["donor_id"]) == AGREES_BY_LUCK

def test_seven_stages_invented_with_no_saf_string(donor, defect):
    inv = sorted(donor.loc[(~donor["saf_recovered"]) & donor["label_map_fired"], "donor_id"])
    assert len(inv) == INVENTED_NO_SAF
    bd = defect["invented_breakdown"]
    assert len(bd["end stage"]) == 5
    assert len(bd["healthy control"]) == 2

def test_one_mechanism_three_symptoms(defect):
    s = defect["one_mechanism_three_symptoms"]
    assert s["overwritten_real_grade"]["donors"] == F_EXCEPTION_SET
    assert s["agrees_by_luck"]["donors"] == AGREES_BY_LUCK
    assert len(s["invented_stage_no_saf_string"]["donors"]) == INVENTED_NO_SAF


# ---- joins: two columns, and the reason must not read as missing data ----
def test_expression_and_ccc_joins_are_separate_and_differ(donor):
    g = donor[donor["saf_recovered"]]
    assert int(g["usable_for_expression_join"].sum()) == EXPR_JOIN_N
    assert int(g["usable_for_ccc_join"].sum()) == CCC_JOIN_N
    assert EXPR_JOIN_N != CCC_JOIN_N

def test_ccc_drops_are_the_expected_three(donor):
    g = donor[donor["saf_recovered"]]
    assert sorted(g.loc[~g["usable_for_ccc_join"], "donor_id"]) == CCC_DROPS

def test_P62_P67_drop_for_depth_not_phenotype(donor):
    d = donor.set_index("donor_id")
    for p in ["P62", "P67"]:
        assert bool(d.loc[p, "usable_for_expression_join"])   # phenotype+expression fine
        assert not bool(d.loc[p, "usable_for_ccc_join"])
        assert "DEPTH_not_missing_phenotype" in d.loc[p, "ccc_exclusion_reason"]

def test_P70_absent_from_expression_substrate(donor):
    r = donor.set_index("donor_id").loc["P70"]
    assert bool(r["saf_recovered"])
    assert not bool(r["usable_for_expression_join"])
    assert not bool(r["usable_for_ccc_join"])


# ---- the artifact must not have touched 343's OUTPUT ---------------------
# Digests taken from the deposited/committed state before this artifact was
# built. If the output changes, the recovery stopped being additive.
#
# 343's SCRIPT digest has been repinned TWICE, both times deliberately:
#   2026-08-30 (a) SAF regex fixed at source so S and A stop being discarded
#                  (SAF_FULL_PATTERN + extract_saf).
#   2026-08-30 (b) LABEL-MAP PRECEDENCE fixed: extract_fstage() now runs
#                  SAF_FULL_PATTERN, then the remaining FSTAGE patterns, then
#                  SAF_LABEL_MAP last. This CHANGES a frozen value -- GSE202379
#                  P98 (GSM6112244) moves F0 -> F1 -- and was authorised
#                  explicitly by the user.
# A bare digest swap would turn this test into a rubber stamp, so
# test_343_script_byte_identical also asserts the properties both fixes were
# for: a repin that reverted either one still fails.
SHA_343_SCRIPT_PRE_SAF_FIX  = "b9aa938afefad21f8323d40cbdd14baa16c84bc355d3319422ae2475612854a9"
SHA_343_SCRIPT_PRE_PREC_FIX = "601956711202dee9bb44d70fbea08e25d1f98fad5c5960465e6a730de652060a"
SHA_343_SCRIPT = "f2d421419fa55a3aaa0c9256913381d10e5ccfec3ea43b856644a554aa2a05a2"
# The PRE-correction release. Left in place and never overwritten; the
# corrected harvest is written beside it as a new timestamped artifact.
SHA_343_OUTPUT = "9548296ff43888a583e6d39552e87f77b37520ca96bfca798daa8f07485735c6"
SHA_SERIES     = "2994ba6b7cddad4b5f44626fbc3e12f495d6c617db7b70a72770a75e49f060cc"


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def test_343_script_byte_identical():
    """Pinned to the POST-fix state, plus the properties both fixes deliver."""
    got = _sha256(H343)
    assert got != SHA_343_SCRIPT_PRE_SAF_FIX, \
        "343 reverted to the pre-SAF-fix source; S and A are discarded again"
    assert got != SHA_343_SCRIPT_PRE_PREC_FIX, \
        "343 reverted to the pre-precedence-fix source; the label map is " \
        "overwriting deposited SAF grades again (P98 back to F0)"
    assert got == SHA_343_SCRIPT

def test_343_source_fix_recovers_S_and_A(h343):
    """The teeth behind the repin: the source must capture all three groups
    and hand back S and A, not just F."""
    h = h343
    assert h.SAF_FULL_PATTERN.groups == 3
    row = pd.Series({"char_0": "saf score: S0A1F0", "sample": "GSMX"})
    assert h.extract_saf(row) == (0, 1, 0)
    assert h.extract_fstage(row) == 0        # F, not S

def test_343_output_byte_identical():
    """The NAMED release must never be overwritten, including by the
    precedence correction. The corrected harvest lives beside it under
    donor_fstage_documented_labelmap_corrected_<UTC stamp>.tsv."""
    out = ROOT/"Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
    assert _sha256(out) == SHA_343_OUTPUT

def test_source_series_matrix_byte_identical():
    src = ROOT/"Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/geo_cache/GSE202379_series_matrix.txt.gz"
    assert _sha256(src) == SHA_SERIES

def test_recovery_wrote_nothing_into_the_singlecell_tree():
    """The overlay must live outside the consumed tree entirely."""
    sc = (ROOT/"Analysis/SingleCell").resolve()
    assert not str(ART.resolve()).startswith(str(sc))
    assert not (ART/"donor_fstage_documented.tsv").exists()
