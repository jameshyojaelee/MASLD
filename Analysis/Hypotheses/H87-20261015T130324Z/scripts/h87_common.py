"""Shared constants and helpers for H87 (prespec PRESPEC_H87.md, sha256 58085a41...)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

H87 = Path(__file__).resolve().parent.parent
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
LIB = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
REVIEW = ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z"
ROUND4 = REVIEW / "round4-20261015T112148Z"

PRESPEC_COPY = H87 / "PRESPEC_H87.md"
PRESPEC_ORIG = REVIEW / "plans/PRESPEC_H87.md"
PRESPEC_SHA = "58085a41ce3f71ba38e8991b5da2863cf0b898e9806cb9b79039fba08f1665dc"

# Inputs table of the prespec (role, path, bytes, sha256), in prespec order.
INPUTS = {
    "members": (ROOT / "GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/credible_set_members_annotated.tsv",
                13625617, "408e0cfe94d9a22549deff04a7d1d7621ed149ac5678ec388aeb76965d0a053f"),
    "repo_vep": (LIB / "GWAS/finemapping/results/candidates/noncoding-dna-vep-2026-08-12/run_19782908/credible_set_variant_consequences_complete.tsv",
                 3088477, "0859df37d2de0516870232e10c22c2fed74834d326b43f6a77ef2ab360f9c01d"),
    "liftover": (ROOT / "GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv",
                 719936, "33697ea30d8a320475ec2f486b0d88c2203f030735c69b4189eab527cd955bb2"),
    "gtf": (Path("/gpfs/commons/home/jameslee/reference_genome/grch38_primary_chr_gencode_v49/gencode.v49.primary_chr.gtf"),
            3323462848, "ff32fd55c6799b3b94fe10aa17b2b5d4da952fa1de12fe44afadf32e949ec914"),
    "catalog_reviewed": (ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet",
                         468456, "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b"),
    "catalog_completion": (ROOT / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/passport_gene_index.parquet",
                           659749, "ca9326d259da3c427f7bd21b225fe58d1844061a6afbb3fd24a12dfe53f22223"),
    "tier_table": (ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv",
                   2296, "865f108de0782a708469b7a4e42943075de79b5ae062dc452f4050034bc15fc7"),
    "fig6c": (ROOT / "figures/main/fig6_gene_catalog/panels/fig6c_worked_entries.pdf",
              16426, "4132b889695b33cf44a8f38d43be0173b7a16b790f9830bf7093931e145e0fcd"),
    "gate_count_script": (ROUND4 / "gate/R15-01/count_coding_genes.py",
                          3657, "19e9ac80d4e0a3305545f8fefe9619ea8dcada37f4dbf44ea938efe996cf54cd"),
    "gate_merge_script": (ROUND4 / "gate/R15-01/merge_vep.py",
                          2632, "07ed9a154f67109ed1aef271a3e671788040a14808076cc035c8c212dccf9f60"),
    "gate_gene_summary": (ROUND4 / "gate/R15-01/out/coding_gene_summary_full.tsv",
                          5608, "a3a5af2aa1694e180f4707e13b1c2aa911e94e0d27b735eec06ff19d1a9cd68a"),
    "gate_member_map": (ROUND4 / "gate/R15-01/out/coding_member_gene_map.tsv",
                        74346, "a0494344bb6696803b63976db8ba3fef61d94ffc4640372ad75cdcd79bc69e2a"),
    "count_script": (ROUND4 / "write/H87/count_n.sh",
                     1700, "65e63bf140c1dce728e7a1ac54208a1a48d6bdb4edb1a3350763c681a2c0461d"),
    "count_output": (ROUND4 / "write/H87/count_n.out",
                     949, "f8b60bd02fd2721ad8cdce30262de1577f2ffcf924227f31b508bf836ef8a37f"),
    "repo_splice_rule": (LIB / "scripts/manuscript/noncoding_resource/noncoding_genetics.py",
                         37675, "d4300fc06ab0a4641684fb8254eea0dbca54ab81faa58aa3b4810435722daff3"),
    "repo_vep_producer": (LIB / "scripts/manuscript/noncoding_resource/annotate_credible_set_vep.py",
                          16411, "cdb3425dcdecdf66bd7a9672c473c86334eb851102ce8fad94938a35d136e41b"),
    "catalog_evidence": (ROOT / "figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_evidence_long.parquet",
                         3194057, "a2d4131bd38464501e595a6554ae1b9956d9ae1fd503f57f240d87120d0c392b"),
    "gate_eppk1_vep": (ROUND4 / "gate/R15-01/tmp/vep/8_144942477.json",
                       1633, "f27c08164cb66a85eb2e432f15ab199d602636ad5ce2f69bccfaca44fb0279fd"),
}

# A2 sets (annotate_credible_set_vep.py:219-232 and noncoding_genetics.py:377-389).
PA = {
    "transcript_ablation", "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant",
}
SPL = {"splice_donor_variant", "splice_acceptor_variant"}
T = PA | SPL
CODING_CATS = {"protein_altering", "canonical_splice"}

SEVERITY_ORDER = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification",
    "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant",
]

SEED = 20260923
ENSEMBL37 = "https://grch37.rest.ensembl.org"
ENSEMBL_MAIN = "https://rest.ensembl.org"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify_prespec() -> str:
    got = sha256_file(PRESPEC_COPY)
    if got != PRESPEC_SHA:
        raise SystemExit(f"S1 STOP: prespec copy sha256 {got} != {PRESPEC_SHA}")
    got_orig = sha256_file(PRESPEC_ORIG)
    if got_orig != PRESPEC_SHA:
        raise SystemExit(f"S1 STOP: prespec original sha256 {got_orig} != {PRESPEC_SHA}")
    return got


def check_inputs(roles=None) -> list[dict]:
    """Return one record per input with observed size and sha256; ok flag per row."""
    rows = []
    for role, (path, size, sha) in INPUTS.items():
        if roles is not None and role not in roles:
            continue
        if not path.exists():
            rows.append(dict(role=role, path=str(path), expected_bytes=size, observed_bytes=-1,
                             expected_sha256=sha, observed_sha256="MISSING", ok=False))
            continue
        obs_size = path.stat().st_size
        obs_sha = sha256_file(path)
        rows.append(dict(role=role, path=str(path), expected_bytes=size, observed_bytes=obs_size,
                         expected_sha256=sha, observed_sha256=obs_sha,
                         ok=(obs_size == size and obs_sha == sha)))
    return rows


def sort_terms(terms) -> list[str]:
    rank = {t: i for i, t in enumerate(SEVERITY_ORDER)}
    return sorted(set(terms), key=lambda t: (rank.get(t, 10_000), t))


def write_json(path: Path, obj) -> None:
    with open(path, "x") as handle:
        json.dump(obj, handle, indent=1, sort_keys=True, default=str)
