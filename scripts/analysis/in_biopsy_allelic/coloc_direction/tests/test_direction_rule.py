"""Synthetic self-tests of the direction rule: every state reachable, every sign flip seen once."""
import sys
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coloc_direction as cd  # noqa: E402


def call(beta_a1, a1, a2, ea, nea, beta_ea, same_lead=True, r_study=np.nan, r_eur=np.nan, pop="EUR",
         f1=(np.nan, np.nan), f2=(np.nan, np.nan)):
    r, block = cd.linkage(same_lead, a1, a2, ea, nea, r_study, r_eur, pop)
    blocks = [cd.palindrome_block(a1, a2, *f1), cd.palindrome_block(ea, nea, *f2), block]
    return cd.direction_sign(beta_a1, r, beta_ea, blocks)


def test_shared_lead_same_letters():
    assert call(0.2, "A", "G", "A", "G", 0.3)[0] == 1
    assert call(-0.2, "A", "G", "A", "G", 0.3)[0] == -1      # strong reverse GWAS effect prints decreases
    assert call(0.2, "A", "G", "A", "G", -0.3)[0] == -1
    assert call(-0.2, "A", "G", "A", "G", -0.3)[0] == 1


def test_shared_lead_swapped_letters_flips_once():
    assert call(0.2, "A", "G", "G", "A", 0.3)[0] == -1
    assert call(-0.2, "A", "G", "G", "A", 0.3)[0] == 1


def test_shared_lead_allele_mismatch_blocks():
    s, why = call(0.2, "A", "G", "A", "C", 0.3)
    assert s == 0 and "alleles differ" in why


def test_different_leads_use_ld_sign():
    assert call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=0.9)[0] == 1
    assert call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=-0.9)[0] == -1
    assert call(-0.2, "A", "G", "C", "T", -0.3, same_lead=False, r_study=-0.9)[0] == -1


def test_weak_or_missing_ld_blocks():
    s, why = call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=0.6)
    assert s == 0 and "weak LD" in why
    s, why = call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=np.nan)
    assert s == 0 and "no LD estimate" in why


def test_non_eur_needs_eur_agreement():
    assert call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=0.9, r_eur=0.8, pop="AFR")[0] == 1
    s, why = call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=0.9, r_eur=-0.8, pop="AFR")
    assert s == 0 and "differs" in why
    s, why = call(0.2, "A", "G", "C", "T", 0.3, same_lead=False, r_study=0.9, r_eur=0.3, pop="EAS")
    assert s == 0 and "EUR" in why


def test_palindromes():
    s, why = call(0.2, "A", "T", "A", "T", 0.3, f1=(0.45, 0.44))
    assert s == 0 and "MAF > 0.4" in why
    assert call(0.2, "A", "T", "A", "T", 0.3, f1=(0.2, 0.22), f2=(0.21, 0.2))[0] == 1
    s, why = call(0.2, "C", "G", "C", "G", 0.3, f1=(0.2, 0.8))
    assert s == 0 and "opposite strand" in why
    s, why = call(0.2, "C", "G", "C", "G", 0.3, f1=(0.2, 0.5))
    assert s == 0 and "MAF > 0.4" in why               # either frequency above 0.4 blocks
    s, why = call(0.2, "C", "G", "C", "G", 0.3, f1=(0.05, 0.3))
    assert s == 0 and "disagree" in why
    s, why = call(0.2, "C", "G", "C", "G", 0.3)
    assert s == 0 and "without a 1kGP frequency" in why
    assert cd.palindrome_block("A", "G", np.nan, np.nan) is None


def test_palindrome_without_file_frequency_needs_the_convention():
    assert cd.palindrome_block("A", "T", np.nan, 0.2, "GWAS lead", True) is None
    assert "breaks" in cd.palindrome_block("A", "T", np.nan, 0.2, "GWAS lead", False)
    assert "no allele-order convention" in cd.palindrome_block("A", "T", np.nan, 0.2, "GWAS lead", None)
    assert "MAF > 0.4" in cd.palindrome_block("A", "T", np.nan, 0.45, "GWAS lead", True)
    assert "without a 1kGP frequency" in cd.palindrome_block("A", "T", 0.2, np.nan, "GWAS lead", True)
    # a file frequency is used before the convention
    assert "opposite strand" in cd.palindrome_block("A", "T", 0.2, 0.8, "GWAS lead", True)


def test_convention_followed():
    assert cd.convention_followed("allele1_is_alt", "allele1", "a=ALT") is True
    assert cd.convention_followed("allele1_is_alt", "allele1", "a=REF") is False
    assert cd.convention_followed("allele1_is_ref", "allele1", "a=REF") is True
    assert cd.convention_followed("allele1_is_alt", "allele2", "a=REF") is True     # effect allele2 = REF, so allele1 = ALT
    assert cd.convention_followed("mixed", "allele1", "a=ALT") is None
    assert cd.convention_followed("allele1_is_alt", "allele1", "allele pair not in 1kGP panel") is None


def test_frequency_block():
    assert "1kGP record may not be" in cd.frequency_block("A", "G", 0.1, 0.5, "GWAS lead")
    assert cd.frequency_block("A", "G", 0.4, 0.5) is None
    assert cd.frequency_block("A", "T", 0.1, 0.9) is None           # palindromes have their own rule
    assert cd.frequency_block("AT", "A", np.nan, 0.9) is None


def test_marginal_floor():
    assert cd.marginal_block(1e-4, "hit1", "GWAS") is None
    assert "p >= 0.001 at hit1" in cd.marginal_block(1e-3, "hit1", "GWAS")
    assert cd.marginal_block(np.nan, "hit2", "Broadaway") is not None
    assert cd.marginal_floor_label(1e-9, 1e-9) == "pass"
    assert cd.marginal_floor_label(0.5, 1e-9) == "fail: GWAS"
    assert cd.marginal_floor_label(0.5, 0.2) == "fail: GWAS and eQTL"


def test_excluded_span_is_bed_half_open():
    spans = [("chr6", 28_500_000, 33_500_000, "MHC_A1_definition")]
    assert cd.excluded_span(spans, "chr6", 30_000_000) == "MHC_A1_definition"
    assert cd.excluded_span(spans, "chr6", 28_500_000) == ""
    assert cd.excluded_span(spans, "chr6", 28_500_001) == "MHC_A1_definition"
    assert cd.excluded_span(spans, "chr6", 33_500_000) == "MHC_A1_definition"
    assert cd.excluded_span(spans, "chr7", 30_000_000) == ""


def test_zero_beta_blocks():
    assert call(0.0, "A", "G", "A", "G", 0.3)[0] == 0
    assert call(0.2, "A", "G", "A", "G", np.nan)[0] == 0


def test_labels():
    assert cd.direction_label("ALT", 1) == "liver-enzyme-raising haplotype increases expression"
    assert cd.direction_label("GGT", -1) == "liver-enzyme-raising haplotype decreases expression"
    assert cd.direction_label("NAFLD", -1) == "NAFLD-risk-raising haplotype decreases expression"
    assert cd.direction_label("GGT", 0) == "not_directional"


def test_combine_rows_agreement_and_disagreement():
    d = pd.DataFrame(dict(k=["a", "a", "b", "b", "c"], rank=[0.9, 0.6, 0.8, 0.7, 0.9], t=[1, 2, 1, 2, 1],
                          direction_sign=[0, 1, 1, -1, 0], not_directional_reason=["x", "", "", "", "y"],
                          eqtl_beta_of_hit2_allele_on_trait_raising_haplotype=[np.nan, 0.2, 0.3, -0.1, np.nan]))
    top = cd.combine_rows(d, ["k"], "rank", ["t"], "rows").set_index("k")
    assert top.loc["a", "direction_sign"] == 1 and top.loc["a", "rank"] == 0.6     # directional row reported
    assert top.loc["b", "direction_sign"] == 0 and top.loc["b", "not_directional_reason"] == "rows disagree"
    assert np.isnan(top.loc["b", "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype"])
    assert top.loc["c", "direction_sign"] == 0 and top.loc["c", "not_directional_reason"] == "y"
    assert [top.loc[k, "n_up"] for k in "abc"] == [1, 1, 0]
    assert [top.loc[k, "n_down"] for k in "abc"] == [0, 1, 0]


# ------------------------------------------------------------------ GWAS effect-allele audit

def test_gwas_effect_alleles():
    assert cd.gwas_effect_alleles("A", "G", 0.3, "allele1") == ("A", "G", 0.3)
    ea, oa, f = cd.gwas_effect_alleles("A", "G", 0.3, "allele2")
    assert (ea, oa) == ("G", "A") and f == pytest.approx(0.7)


def test_effect_allele_column():
    assert cd.effect_allele_column(20, 0.95, [1]) == "allele1"
    assert cd.effect_allele_column(20, 0.05, [-1, -1]) == "allele2"
    assert cd.effect_allele_column(20, 0.95, [-1]) == "unresolved"     # control contradicts the audit
    assert cd.effect_allele_column(20, 0.05, [1]) == "unresolved"
    assert cd.effect_allele_column(9, 1.0, []) == "unresolved"         # too few loci
    assert cd.effect_allele_column(50, 0.5, []) == "unresolved"        # neither convention


def _snv(chrom, pos, a1, a2, beta, p):
    return pd.DataFrame(dict(chromosome=chrom, position=pos, allele1=a1, allele2=a2, beta_num=beta, pval_num=p))


def test_load_snv_extract_keeps_nonpalindromic_snvs(tmp_path):
    f = tmp_path / "x.tsv"
    f.write_text("\t".join(cd.GWAS_COLS) + "\n"
                 "1\t100\tA\tG\t0.1\t0.01\t1e-9\t0.2\n"
                 "1\t200\tA\tT\t0.1\t0.01\t1e-9\t0.2\n"      # palindrome
                 "1\t300\tAT\tA\t0.1\t0.01\t1e-9\t0.2\n"     # indel
                 "1\t400\tC\tT\t-0.2\t0.01\t1e-3\tNA\n")
    d = cd.load_snv_extract(f)
    assert d.position.tolist() == [100, 400] and d.beta_num.tolist() == [0.1, -0.2]
    assert cd.load_snv_extract(tmp_path.joinpath("x.tsv")).pval_num.tolist() == [1e-9, 1e-3]


def test_orientation_concordance_matches_letters_both_orders_and_thins():
    pos = [i * 2_000_000 for i in range(1, 13)]
    ref = _snv(["1"] * 12, pos, ["A"] * 12, ["G"] * 12, [0.1] * 12, [1e-10] * 12)
    same = _snv(["1"] * 12, pos, ["A"] * 12, ["G"] * 12, [0.2] * 12, [1e-3] * 12)
    swapped = _snv(["1"] * 12, pos, ["G"] * 12, ["A"] * 12, [-0.2] * 12, [1e-3] * 12)   # same effect, letters swapped
    allele2 = _snv(["1"] * 12, pos, ["A"] * 12, ["G"] * 12, [-0.2] * 12, [1e-3] * 12)   # beta of G under A/G
    assert cd.orientation_concordance(ref, same) == (12, 1.0)
    assert cd.orientation_concordance(ref, swapped) == (12, 1.0)
    assert cd.orientation_concordance(ref, allele2) == (12, 0.0)
    # two reference hits in one 1 Mb window count once; a study row with other letters is ignored
    ref2 = pd.concat([ref, _snv(["1"], [pos[0] + 10], ["C"], ["T"], [0.1], [1e-9])])
    other = pd.concat([same, _snv(["1"], [pos[0] + 10], ["C"], ["A"], [-0.3], [1e-4])])
    assert cd.orientation_concordance(ref2, other) == (12, 1.0)
    weak = same.assign(pval_num=0.5)
    assert cd.orientation_concordance(ref, weak)[0] == 0
    # hits of the audited study count too: a reference with no genome-wide hit still audits it
    small_ref = ref.assign(pval_num=1e-3)
    big_study = allele2.assign(pval_num=1e-12)
    assert cd.orientation_concordance(small_ref, big_study) == (12, 0.0)
    assert cd.orientation_concordance(small_ref, small_ref)[0] == 0      # neither file genome-wide significant


# ------------------------------------------------------------------ haplotype panel and pair evaluation

def make_panel():
    # 60 EUR + 60 AFR samples; site 100 REF=A ALT=G, site 200 REF=C ALT=T
    rng = np.random.default_rng(1)
    n = 120
    h100 = rng.integers(0, 2, 2 * n).astype(np.int8)
    h200 = (1 - h100).astype(np.int8)            # T at 200 sits on the A haplotype at 100
    samples = [f"s{i}" for i in range(n)]
    pops = {s: ("EUR" if i < 60 else "AFR") for i, s in enumerate(samples)}
    haps = {"chr1": {100: [("A", "G", h100)], 200: [("C", "T", h200)]}}
    return cd.Panel(haps, samples, pops), h100


def test_panel_indicator_both_orders_and_strand():
    panel, h100 = make_panel()
    x_alt, o = panel.indicator("chr1", 100, "+", "G", "A")
    assert o == "a=ALT" and (x_alt == h100).all()
    x_ref, o = panel.indicator("chr1", 100, "+", "A", "G")
    assert o == "a=REF" and (x_ref == 1 - h100).all()
    x_minus, o = panel.indicator("chr1", 100, "-", "C", "T")     # complement of G/A
    assert o == "a=ALT" and (x_minus == h100).all()
    assert panel.indicator("chr1", 100, "+", "G", "C")[0] is None
    y, _ = panel.indicator("chr1", 200, "+", "T", "C")
    assert panel.r(x_ref, y, "EUR") == pytest.approx(1.0)
    assert panel.r(x_alt, y, "EUR") == pytest.approx(-1.0)


def test_resolve_hit1_follows_coloc_merge():
    g = pd.DataFrame(dict(allele1=["A", "C"], allele2=["AT", "T"], beta_num=[0.1, 0.2], pval_num=[1e-9, 1e-8],
                          af_num=[np.nan, np.nan]))
    e = pd.DataFrame(dict(EA=["T"], NEA=["C"], beta_num=[0.5], p_num=[1e-5], eaf_num=[0.3]))
    row, erow, why, rank = cd.resolve_hit1(g, e)
    assert why == "" and rank == 2 and row.allele1 == "C" and erow.EA == "T"


Row = namedtuple("Row", "study ensembl hit1 hit2")


def pair_inputs(gwas_beta=0.3):
    gw = {("S", "1", 90): pd.DataFrame(dict(allele1=["A"], allele2=["G"], beta_num=[gwas_beta], pval_num=[1e-20],
                                            af_num=[np.nan]))}
    eq = {("E", "1", 90): pd.DataFrame(dict(EA=["G"], NEA=["A"], beta_num=[-0.1], p_num=[1e-3], eaf_num=[0.5])),
          ("E", "1", 190): pd.DataFrame(dict(EA=["T"], NEA=["C"], beta_num=[0.4], p_num=[1e-30], eaf_num=[0.5]))}
    lifted = {("1", 90): ("chr1", 100, "+"), ("1", 190): ("chr1", 200, "+")}
    return gw, eq, lifted


def test_evaluate_pair_end_to_end():
    panel, _ = make_panel()
    gw, eq, lifted = pair_inputs()
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {})
    # allele1 A raises the trait; A at hit1 sits with T at hit2 (r(A,T) = +1); T raises expression.
    assert out["ld_r_effect_ea_study"] == pytest.approx(1.0)
    assert out["trait_raising_allele"] == "A" and out["ld_r_trait_raising_allele_vs_ea"] == pytest.approx(1.0)
    assert out["direction_sign"] == 1 and out["not_directional_reason"] == ""
    # the Broadaway marginal effect at hit1 of A is +0.1 (EA=G has -0.1): the check agrees
    assert out["check_eqtl_at_hit1_sign"] == 1
    # reverse the GWAS effect: the same haplotype now lowers the trait
    gw, _, _ = pair_inputs(gwas_beta=-0.3)
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {})
    assert out["direction_sign"] == -1 and out["trait_raising_allele"] == "G"
    assert out["ld_r_trait_raising_allele_vs_ea"] == pytest.approx(-1.0)


def test_evaluate_pair_allele2_study_flips_once():
    """A file whose beta belongs to allele2: G (not A) raises the trait, and the sign flips exactly once."""
    panel, _ = make_panel()
    gw, eq, lifted = pair_inputs()
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {}, "allele2")
    assert out["gwas_effect_allele"] == "G" and out["trait_raising_allele"] == "G"
    assert out["ld_r_effect_ea_study"] == pytest.approx(-1.0)
    assert out["ld_r_trait_raising_allele_vs_ea"] == pytest.approx(-1.0)
    assert out["direction_sign"] == -1
    assert out["check_eqtl_at_hit1_sign"] == -1                  # G carries the -0.1 expression effect


def test_evaluate_pair_blocks():
    panel, _ = make_panel()
    gw, eq, lifted = pair_inputs()
    row = Row("S", "E", "1:90", "1:190")
    out, _, _ = cd.evaluate_pair(row, gw, eq, lifted, panel, "EUR", {}, "unresolved")
    assert out["direction_sign"] == 0 and "unresolved" in out["not_directional_reason"]
    assert out["trait_raising_allele"] == ""
    out, _, _ = cd.evaluate_pair(row, gw, eq, lifted, panel, "EUR", {}, "allele1", conditional=False)
    assert out["direction_sign"] == 0 and "SuSiE component" in out["not_directional_reason"]
    out, _, _ = cd.evaluate_pair(row, gw, eq, lifted, panel, "EUR", {}, "allele1", conditional=True)
    assert out["direction_sign"] == 1
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:191"), gw, eq, lifted, panel, "EUR", {})
    assert out["direction_sign"] == 0 and "not in Broadaway" in out["not_directional_reason"]


def test_evaluate_pair_gtex_check_reads_alt_slope():
    panel, _ = make_panel()
    gw, eq, lifted = pair_inputs()
    gtex = {("E", "chr1", 200): [("C", "T", 0.5)]}                # GTEx slope is the ALT (T) effect
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", gtex)
    assert out["check_gtex_liver_same_sign_hit2"] == 1
    gtex = {("E", "chr1", 200): [("T", "C", 0.5)]}                # ALT is C: T effect is -0.5
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", gtex)
    assert out["check_gtex_liver_same_sign_hit2"] == 0


def test_evaluate_pair_marginal_floor_keeps_technical_sign():
    panel, _ = make_panel()
    gw, eq, lifted = pair_inputs()
    gw[("S", "1", 90)] = gw[("S", "1", 90)].assign(pval_num=0.5)
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {})
    assert out["direction_sign"] == 0 and "GWAS marginal p" in out["not_directional_reason"]
    assert out["sign_before_evidence_rules"] == 1 and out["marginal_floor"] == "fail: GWAS"
    gw, eq, lifted = pair_inputs()
    eq[("E", "1", 190)] = eq[("E", "1", 190)].assign(p_num=0.01)
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {}, conditional=False)
    assert out["direction_sign"] == 0 and "Broadaway marginal p" in out["not_directional_reason"]
    assert "SuSiE component" in out["all_block_reasons"] and out["sign_before_evidence_rules"] == 1


def test_evaluate_pair_excluded_span_and_frequency_blocks():
    panel, _ = make_panel()
    gw, eq, lifted = pair_inputs()
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {},
                                 spans=[("chr1", 150, 250, "IG_locus")])
    assert out["direction_sign"] == 0 and out["not_directional_reason"] == "hit2 in excluded span (IG_locus)"
    assert out["sign_before_evidence_rules"] == 0
    eq[("E", "1", 190)] = eq[("E", "1", 190)].assign(eaf_num=0.97)          # panel frequency of T is about 0.5
    out, _, _ = cd.evaluate_pair(Row("S", "E", "1:90", "1:190"), gw, eq, lifted, panel, "EUR", {})
    assert out["direction_sign"] == 0 and "eQTL lead: source and 1kGP" in out["not_directional_reason"]


def palindrome_panel():
    n = 120
    h100 = np.array([1 if i % 10 < 3 else 0 for i in range(2 * n)], dtype=np.int8)   # ALT T at 0.3
    h200 = (1 - h100).astype(np.int8)
    samples = [f"s{i}" for i in range(n)]
    haps = {"chr1": {100: [("A", "T", h100)], 200: [("C", "T", h200)]}}
    return cd.Panel(haps, samples, {s: "EUR" for s in samples})


def test_evaluate_pair_palindromic_gwas_lead_uses_the_file_convention():
    panel = palindrome_panel()
    _, eq, lifted = pair_inputs()
    eq[("E", "1", 90)] = pd.DataFrame(dict(EA=["T"], NEA=["A"], beta_num=[-0.1], p_num=[1e-3], eaf_num=[0.3]))
    eq[("E", "1", 190)] = eq[("E", "1", 190)].assign(eaf_num=0.7)
    row = Row("S", "E", "1:90", "1:190")
    gw = {("S", "1", 90): pd.DataFrame(dict(allele1=["T"], allele2=["A"], beta_num=[0.3], pval_num=[1e-20],
                                            af_num=[np.nan]))}
    out, _, _ = cd.evaluate_pair(row, gw, eq, lifted, panel, "EUR", {}, convention="allele1_is_alt")
    # T at 100 sits with C at 200 (r(T, T@200) = -1); T at 200 raises expression: -1
    assert out["gwas_allele1_follows_file_convention"] is True
    assert out["direction_sign"] == -1 and out["not_directional_reason"] == ""
    out, _, _ = cd.evaluate_pair(row, gw, eq, lifted, panel, "EUR", {}, convention="allele1_is_ref")
    assert out["direction_sign"] == 0 and "breaks" in out["not_directional_reason"]
    out, _, _ = cd.evaluate_pair(row, gw, eq, lifted, panel, "EUR", {}, convention="mixed")
    assert out["direction_sign"] == 0 and "no allele-order convention" in out["not_directional_reason"]


def test_blocked_pair_counts():
    floor = "GWAS marginal p >= 0.001 at hit1"
    p = pd.DataFrame(dict(k=["a", "a", "a", "b", "c", "c"], direction_sign=[1, 0, 0, 0, -1, 0],
                          sign_before_evidence_rules=[1, -1, 0, 1, -1, 1],
                          evidence_blocks=["", floor, "", cd.COMPONENT_REASON, "", cd.COMPONENT_REASON]))
    calls = pd.DataFrame(dict(k=["a", "b", "c"], direction_sign=[1, 0, -1]))
    c = cd.blocked_pair_counts(p, calls, ["k"], "direction_sign").set_index("k")
    assert c.loc["a", "n_pairs_blocked"] == 2 and c.loc["a", "n_pairs_evidence_blocked_opposing_call"] == 1
    assert c.loc["a", "n_pairs_evidence_blocked_agreeing_call"] == 0
    assert c.loc["a", "n_pairs_component_blocked_opposing_call"] == 0     # floor-blocked: its sign is not evidence
    assert c.loc["b", "n_pairs_evidence_blocked_opposing_call"] == 0     # no call to oppose
    assert c.loc["c", "n_pairs_component_blocked_opposing_call"] == 1


def test_opposing_component_pair_block():
    t = pd.DataFrame(dict(k=["a", "c", "d"], direction_sign=[1, -1, 0], not_directional_reason=["", "", "x"],
                          eqtl_beta_of_hit2_allele_on_trait_raising_haplotype=[0.2, -0.3, np.nan],
                          n_pairs_component_blocked_opposing_call=[0, 1, 1])).set_index("k")
    t = cd.opposing_component_pair_block(t)
    assert t.loc["a", "direction_sign"] == 1
    assert t.loc["c", "direction_sign"] == 0 and t.loc["c", "not_directional_reason"] == cd.OPPOSING_COMPONENT_REASON
    assert np.isnan(t.loc["c", "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype"])
    assert t.loc["d", "not_directional_reason"] == "x"                  # already not directional: reason kept


def test_tag_flags_labels():
    tag_ld = pd.DataFrame(dict(release="r", ensembl=["g1", "g1", "g2", "g2", "g3"], hit2=["1:1", "1:2", "1:3", "1:4", "1:5"],
                               tag_pos=[10, 10, 20, 20, 30], r2_eur=[0.9, np.nan, 0.5, np.nan, np.nan], in_L=True))
    units = pd.DataFrame(dict(release="r", ensembl=["g1", "g2", "g3", "g4", "g5"]))
    f = cd.tag_flags(tag_ld, units, ["release", "ensembl"], {"g1", "g2", "g3", "g5"}, "L", "r2_eur", "gene_") \
        .set_index("ensembl")["gene_tag_r2_ge_0.8__L"]
    assert f["g1"] == "yes"
    assert f["g2"] == "no among estimable (1 of 2 leads not estimable)"
    assert f["g3"] == "not estimable"
    assert f["g4"] == "no Model A tag"
    assert f["g5"] == "not estimable"          # tags but no colocalized lead record
