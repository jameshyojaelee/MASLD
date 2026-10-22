#!/usr/bin/env python3
"""Coordinate and biological-object invariants for annotation correspondence."""
from pathlib import Path
import tempfile
import unittest
from resolve_annotation_correspondence import (Gene, alignment_blocks, project, union,
    verified_exon_segments, allele_identity_candidates, added_patch_component,
    reference_exclusion, write_grouped_annotation, exact_transcript_correspondence,
    resolve_exact_RNA_transcripts, hold_unresolved, DISPOSITIONS)
import gzip


class CoordinateTests(unittest.TestCase):
    def alignment(self, strand, target, gap, start=101, end=112):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"alignment.gff"
            p.write_text(f'chr1\tNCBI\tmatch\t{start}\t{end}\t.\t{strand}\t.\tTarget=ALT {target};Gap={gap}\n')
            return alignment_blocks(p)[1]

    def test_deletion_does_not_project_an_alt_base_into_reference_gap(self):
        blocks = self.alignment("+", "1 10 +", "M5 D2 M5")
        mapped, covered, directions = project([(0,10)], blocks)
        self.assertEqual(mapped, [(100,105),(107,112)])
        self.assertEqual(covered,10)
        self.assertEqual(directions,{1})

    def test_insertion_is_unresolved_for_exonic_bases_in_insert(self):
        blocks = self.alignment("+", "1 14 +", "M5 I2 M7")
        self.assertEqual(project([(0,14)],blocks)[1],12)

    def test_partial_equivalent_rescue_excludes_unaligned_inserted_bases(self):
        blocks = self.alignment("+", "1 14 +", "M5 I2 M7")
        self.assertEqual(verified_exon_segments([(0,14)], blocks), [(0,5),(7,14)])

    def test_reverse_orientation_preserves_half_open_coordinates(self):
        blocks = self.alignment("-", "1 12 +", "M12")
        self.assertEqual(project([(0,3)],blocks), ([(109,112)],3,{-1}))

    def test_overlapping_isoform_exons_are_not_counted_twice(self):
        self.assertEqual(union([(0,10),(5,15),(0,10)]),[(0,15)])

    def test_malformed_gap_cannot_license_a_correspondence(self):
        with self.assertRaises(ValueError):
            self.alignment("+", "1 12 +", "M10")

    def test_exact_isoform_can_match_when_gene_union_has_other_isoforms(self):
        blocks = self.alignment("+", "1 12 +", "M12")
        self.assertTrue(exact_transcript_correspondence([(0,3),(9,12)], [(100,103),(109,112)], "+", "+", blocks))

    def test_shared_transcript_identity_does_not_rescue_incomplete_exon_coordinates(self):
        blocks = self.alignment("+", "1 14 +", "M5 I2 M7")
        self.assertFalse(exact_transcript_correspondence([(0,14)], [(100,112)], "+", "+", blocks))

    def test_exact_transcript_requires_matching_strand(self):
        blocks = self.alignment("+", "1 12 +", "M12")
        self.assertFalse(exact_transcript_correspondence([(0,12)], [(100,112)], "+", "-", blocks))


class ExactRNAIdentityTests(unittest.TestCase):
    def resolve(self, duplicate_accession=False, flags=None, conflict=False, complete=True):
        with tempfile.TemporaryDirectory() as tmp:
            gtf, metadata = Path(tmp)/"genes.gtf.gz", Path(tmp)/"RefSeq.gz"
            lines = []
            for gene, transcript, chrom, start in [("A", "AT", "ALT", 1), ("P", "PT", "chr1", 101)]:
                for feature in ("transcript", "exon"):
                    lines.append(f'{chrom}\tfixture\t{feature}\t{start}\t{start+9}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{transcript}";\n')
            if duplicate_accession:
                lines.append('chr2\tfixture\ttranscript\t201\t210\t.\t+\t.\tgene_id "P2"; transcript_id "P2T";\n')
            with gzip.open(gtf, "wt") as handle: handle.writelines(lines)
            with gzip.open(metadata, "wt") as handle:
                handle.write("AT\tNR_123456.1\nPT\tNR_123456.1\n")
                if duplicate_accession: handle.write("P2T\tNR_123456.1\n")
            genes = {"A": Gene("A", "", "lncRNA", "ALT", "+", 0, 10),
                     "P": Gene("P", "", "lncRNA", "chr1", "+", 100, 110),
                     "P2": Gene("P2", "", "lncRNA", "chr2", "+", 200, 210)}
            rows = [dict(gene_id="A", disposition="unresolved", counting_gene_id="A", primary_gene_candidates="P"),
                    dict(gene_id="P", disposition="primary_gene", counting_gene_id="P")]
            groups, memberships = {}, {}
            if flags or conflict:
                groups = {"G": [("A", set(flags or [])), ("P2" if conflict else "P", set())]}
                memberships = {"A": ["G"]}
            blocks = [(0,10 if complete else 9,100,110 if complete else 109,1)]
            resolve_exact_RNA_transcripts(gtf, metadata, genes, {"ALT": ("chr1", blocks)}, rows, groups, memberships)
            return rows[0]

    def test_unique_curated_RNA_and_complete_isoform_coordinates_support_identity(self):
        self.assertEqual(self.resolve()["counting_gene_id"], "P")

    def test_xref_alone_never_resolves_incomplete_coordinates(self):
        self.assertEqual(self.resolve(complete=False)["disposition"], "unresolved")

    def test_RNA_accession_shared_by_primary_genes_remains_ambiguous(self):
        self.assertEqual(self.resolve(duplicate_accession=True)["disposition"], "unresolved")

    def test_declared_gene_group_conflict_cannot_be_overruled(self):
        self.assertEqual(self.resolve(conflict=True)["disposition"], "unresolved")

    def test_invalid_reference_gene_flags_prevent_xref_rescue(self):
        self.assertEqual(self.resolve(flags=["IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY"])["disposition"], "unresolved")


class ReferenceExclusionTests(unittest.TestCase):
    def setUp(self):
        self.gene = Gene("G", "", "lncRNA", "ALTACC", "+", 10, 100, exons=[(10,20),(90,100)])
        seq = {"genbank": "ALTACC", "ucsc": "chr_alt"}
        self.report = {"ALTACC": seq, "chr_alt": seq}

    def test_whole_gene_exclusion_accepts_verified_assembly_alias(self):
        self.assertEqual(reference_exclusion(self.gene, self.report,
            [("chr_alt", 0, 100, "false_duplication")]), "false_duplication")

    def test_intron_only_overlap_retains_expression_features(self):
        self.assertIsNone(reference_exclusion(self.gene, self.report,
            [("chr_alt", 30, 80, "false_duplication")]))

    def test_partial_gene_overlap_does_not_license_complete_exclusion(self):
        self.assertIsNone(reference_exclusion(self.gene, self.report,
            [("chr_alt", 0, 99, "contamination")]))

    def test_unknown_alias_cannot_license_exclusion(self):
        self.assertIsNone(reference_exclusion(self.gene, self.report,
            [("unknown", 0, 100, "contamination")]))

    def test_excluded_record_survives_correspondence_but_has_no_counting_exon(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dest = Path(tmp)/"src.gtf.gz", Path(tmp)/"out.gtf.gz"
            with gzip.open(src, "wt") as handle:
                handle.write('chr1\tfixture\texon\t1\t10\t.\t+\t.\tgene_id "valid";\n')
                handle.write('ALTACC\tfixture\texon\t11\t20\t.\t+\t.\tgene_id "excluded";\n')
            rows = [dict(gene_id="valid", counting_gene_id="valid", disposition="primary_gene"),
                    dict(gene_id="excluded", counting_gene_id="", disposition="excluded_reference")]
            write_grouped_annotation(src, dest, rows, {})
            with gzip.open(dest, "rt") as handle:
                result = handle.read()
            self.assertIn('gene_id "valid"', result)
            self.assertNotIn('excluded', result)
            self.assertEqual(len(rows), 2)


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.alt = Gene("ALT", "", "lncRNA", "ALT1", "+", 0, 10)
        self.primary = Gene("P", "", "lncRNA", "chr1", "+", 100, 110)
        self.genes = {g.gene_id: g for g in [self.alt, self.primary]}
        self.groups = {"G": [("ALT", set()), ("P", set())]}

    def candidates(self, bodies=None):
        return allele_identity_candidates(self.alt, [self.primary] if bodies is None else bodies,
                                           self.genes, self.groups, {"ALT": ["G"]})

    def test_source_group_and_independent_locus_both_required(self):
        self.assertEqual(self.candidates(), [self.primary])
        self.assertEqual(self.candidates([]), [])

    def test_conflicting_HGNC_cannot_be_overruled_by_allele_group(self):
        self.alt.hgnc = "HGNC:1"; self.primary.hgnc = "HGNC:2"
        self.assertEqual(self.candidates(), [])

    def test_same_haplotype_copies_remain_distinct(self):
        self.genes["ALT2"] = Gene("ALT2", "", "lncRNA", "ALT1", "+", 20, 30)
        self.groups["G"].append(("ALT2", set()))
        self.assertEqual(self.candidates(), [])

    def test_missing_exact_GTF_member_cannot_look_like_singleton_group(self):
        self.groups["G"].append(("P.other_version", set()))
        self.assertEqual(self.candidates(), [])

    def test_added_patch_requires_positive_gap_and_component_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with gzip.open(root/"alt.scaf.agp.gz", "wt") as out:
                out.write("ALTACC\t1\t100\t1\tW\tNEWACC\t1\t100\t+\n")
            with gzip.open(root/"chr13.agp.gz", "wt") as out:
                out.write("PRIMARY\t1\t1000\t1\tN\t1000\tshort_arm\tno\tna\n")
            seq = {"unit": "PATCHES", "genbank": "ALTACC"}
            placement = {"parent_acc": "PRIMARY", "parent_name": "13", "parent_start": "200", "parent_stop": "299"}
            self.assertTrue(added_patch_component(self.alt, seq, placement, root))
            self.alt.end = 101
            self.assertFalse(added_patch_component(self.alt, seq, placement, root))

    def test_multiple_group_memberships_keep_conflicting_primary_targets(self):
        second = Gene("P2", "", "lncRNA", "chr1", "+", 100, 110)
        self.genes["P2"] = second
        self.groups["G2"] = [("ALT", set()), ("P2", set())]
        selected = allele_identity_candidates(self.alt, [self.primary, second], self.genes,
                                             self.groups, {"ALT": ["G", "G2"]})
        self.assertEqual({p.gene_id for p in selected}, {"P", "P2"})

    def test_held_annotation_keeps_its_own_counting_target(self):
        rows = [dict(gene_id="P", counting_gene_id="P", disposition="primary_gene", reason=""),
                dict(gene_id="U", counting_gene_id="", disposition="unresolved", reason="no_coordinate_alignment")]
        held = hold_unresolved(rows)
        self.assertEqual([r["gene_id"] for r in held], ["U"])
        self.assertEqual((rows[1]["disposition"], rows[1]["counting_gene_id"], rows[1]["reason"]),
                         ("held_untestable", "U", "no_coordinate_alignment"))
        self.assertEqual(rows[0]["disposition"], "primary_gene")
        self.assertIn("held_untestable", DISPOSITIONS)


if __name__ == "__main__": unittest.main()
