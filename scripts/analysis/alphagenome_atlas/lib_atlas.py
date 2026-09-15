"""Shared helpers for the AlphaGenome Atlas extension (protocol §3–§5, §10).

Reuses, by copy, the 1-Mb single-linkage block rule from
Analysis/Multimodal_Program_Projection/scripts/source_independent_risk_state_relay/
02_build_corrected_registry_v3.py and the key/liver-track conventions of
GWAS/finemapping/src/65_alphagenome_score.py. Nothing here touches the network.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import pathlib
import re
import sys
from collections import defaultdict
from typing import Callable, Iterable

PROJECT = pathlib.Path(
    os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
).resolve()
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
PRESPEC_PATH = SCRIPT_DIR / "00_prespecification.json"

LOCUS_WINDOW_BP = 1_000_000
LIVER_ONTOLOGY = {"UBERON:0002107", "CL:0000182"}
HEPATOCYTE_ONTOLOGY = {"CL:0000182"}
LIVER_NAME_RE = re.compile(r"\bliver\b|hepatocyte", re.IGNORECASE)
HEPG2_RE = re.compile(r"hepg2", re.IGNORECASE)
STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


class ContractError(RuntimeError):
    """Raised when a table violates a prespecified invariant."""


# ------------------------------------------------------------------ io helpers
def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def prespec() -> dict:
    with PRESPEC_PATH.open() as handle:
        return json.load(handle)


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: pathlib.Path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv_once(path: pathlib.Path, rows: Iterable[dict], columns: list[str]) -> int:
    path = pathlib.Path(path)
    if path.exists():
        raise ContractError(f"refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row.get(c)) for c in columns})
            n += 1
    return n


TRACK0_DEFAULT = "GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z"


def track0_root() -> pathlib.Path:
    """Root of the Track 0 deposits a step script reads.

    Defaults to the named run root; AGA_TRACK0_ROOT overrides it. Several arms were scored against the
    preview (direct-only) archive while the full run was still building, so which archive a number came from
    has to be a stated switch that lands in the provenance, not an edit to a constant.
    """
    v = os.environ.get("AGA_TRACK0_ROOT", "")
    return pathlib.Path(v).resolve() if v.strip() else (PROJECT / TRACK0_DEFAULT)


def leafcutter_intron_key(junction_start: int, junction_end: int) -> tuple:
    """Atlas (half-open) junction coordinates -> LeafCutter (closed) intron coordinates.

    Measured on both LeafCutter substrates, not assumed. GTEx v8 liver sQTL: start matches at offset 0
    (795 hits, next best 5) and end at +1 (814, next best 4). joint_v2 disease splicing: start at offset 0
    (17,159 hits, next best 11) and end at +1 (19,893, next best 133). A +/-1 tolerance admits the -1 and -2
    near-misses and can match a neighbouring junction inside a dense cluster.
    """
    return (int(junction_start), int(junction_end) + 1)


def marginal_expected_concordance(pred, meas) -> float:
    """Sign concordance expected from the two marginal sign skews alone.

    If predictor and measurement lean the same way for a reason unrelated to the variant, sign agreement
    follows with no real signal behind it. An observed concordance only means something above this number.
    """
    import numpy as _np

    pp = float(_np.mean(_np.asarray(pred) > 0))
    pm = float(_np.mean(_np.asarray(meas) > 0))
    return pp * pm + (1.0 - pp) * (1.0 - pm)


def out_root() -> pathlib.Path:
    value = os.environ.get("AGA_OUT_ROOT", "")
    if not value:
        raise ContractError("AGA_OUT_ROOT is unset")
    return pathlib.Path(value).resolve()


def load_api_key() -> tuple[str, str]:
    k = os.environ.get("ALPHA_GENOME_API_KEY", "")
    if k.strip():
        return k.strip(), "env:ALPHA_GENOME_API_KEY"
    kf = os.path.expanduser("~/.alphagenome_key")
    if os.path.exists(kf):
        with open(kf) as handle:
            v = handle.read().strip()
        if v:
            return v, "keyfile:~/.alphagenome_key"
    raise ContractError("no API key in ALPHA_GENOME_API_KEY or ~/.alphagenome_key")


# ------------------------------------------------------------------ variants
def variant_uid(chrom: str, position_1based: int, ref: str, alt: str) -> str:
    if not str(chrom).startswith("chr"):
        raise ValueError(f"chromosome must carry the chr prefix: {chrom}")
    return f"{chrom}:{int(position_1based)}:{ref.upper()}:{alt.upper()}"


def reverse_complement(value: str) -> str:
    return "".join(COMP[b] for b in reversed(value.upper()))


SOURCE_ID_RE = re.compile(r"^(?P<asm>[A-Za-z0-9]+):(?:chr)?(?P<chrom>[0-9]{1,2}|[XYM]|MT):(?P<pos>\d+):"
                          r"(?P<ref>[ACGTNacgtn]+):(?P<alt>[ACGTNacgtn]+)$")


def normalise_source_id(source_id: str) -> str:
    """One canonical form for a source variant id, so `hg19:15:...` and `hg19:chr15:...` are the same key."""
    m = SOURCE_ID_RE.match(str(source_id).strip())
    if m is None:
        raise ContractError(f"source id is not assembly:chrom:pos:ref:alt: {source_id!r}")
    return f"{m['asm']}:chr{m['chrom']}:{int(m['pos'])}:{m['ref'].upper()}:{m['alt'].upper()}"


def posterior_key(variant_uid: str, source_variant_id: str) -> str:
    """The identity a posterior entry is keyed on. NEVER key on variant_uid alone.

    A variant the Atlas defers or fails to lift has an EMPTY variant_uid: step 04 records the exclusion
    reason instead. Keying on that column collapses every such variant of a signal into one entry whose
    weight is their sum, removing exactly the class a coverage analysis exists to count.
    """
    uid = str(variant_uid or "").strip()
    if uid:
        return uid
    src = str(source_variant_id or "").strip()
    if not src:
        raise ContractError("a posterior entry with neither an hg38 uid nor a source id cannot be identified")
    return normalise_source_id(src)


def assert_unique_identities(pairs: Iterable[tuple]) -> None:
    """Refuse a posterior whose entries share an identity, or carry none."""
    seen = set()
    for group, key in pairs:
        if not str(key or "").strip():
            raise ContractError(f"{group}: a posterior entry has a blank identity")
        if (group, key) in seen:
            raise ContractError(f"{group}: two posterior entries share the identity {key}")
        seen.add((group, key))


def resolve_indel_alleles(chrom: str, position_1based: int, allele1: str, allele2: str,
                          fetch: Callable[[str, int, int], str]) -> dict:
    """hg38 reference/alternate alleles for an indel whose liftover was 1:1, validated against the FASTA.

    The Atlas point-query API refuses indels, so these stay `excluded` from the query set. What they gain is
    an identity: the coordinate was always known, only the alleles and the uid were never written, and
    without a uid every deferred indel of a signal shares the empty key.

    Forward strand only. A multi-base allele's strand flip is not a base-wise complement of the VCF-anchored
    representation, so matching a reverse complement would invent a variant that is not at that site.

    `orientation_ambiguous` is True when both allele orders are representable against the reference, which
    for prefix-anchored indels means the long allele is present. 87.7% of this Resource's deferred indels
    are in that state, so it is a property of the class and not an edge case: for those, the reported
    reference allele rests on the source's allele order, not on the reference.
    """
    row = {"hg38_ref": None, "hg38_alt": None, "orientation": None, "allele_swap": None,
           "variant_uid": None, "mapping_status": "excluded", "exclusion_reason": "indel_deferred",
           "orientation_ambiguous": False}
    a1, a2 = str(allele1).upper(), str(allele2).upper()
    at_ref = [ref for ref in (a1, a2)
              if fetch(chrom, position_1based - 1, position_1based - 1 + len(ref)).upper() == ref]
    # Both readings match whenever the long allele is present, because a VCF indel's short allele is a
    # prefix of its long one. The reference therefore cannot orient those; the source's own allele order is
    # the only other information, so it is tried first and the ambiguity is flagged rather than decided.
    row["orientation_ambiguous"] = len(at_ref) > 1
    for ref, alt, swap in ((a1, a2, "none"), (a2, a1, "swapped")):
        if ref in at_ref:
            row.update(hg38_ref=ref, hg38_alt=alt, orientation="forward", allele_swap=swap,
                       variant_uid=variant_uid(chrom, position_1based, ref, alt))
            return row
    row["exclusion_reason"] = "indel_deferred_ref_mismatch"
    return row


def crosswalk_variant(
    chrom: str,
    position_1based: int,
    allele1: str,
    allele2: str,
    fetch: Callable[[str, int, int], str],
    n_liftover_mappings: int,
) -> dict:
    """Validate one hg38 site against the reference. Never repairs silently.

    ``fetch`` uses zero-based half-open coordinates. Returns a crosswalk row with
    mapping_status ∈ {mapped, excluded}; alleles are oriented to reference→alternate
    and every reorientation is recorded (allele_swap, orientation, palindromic).
    """
    row = {
        "hg38_chrom": chrom,
        "hg38_position_1based": int(position_1based),
        "source_alleles": f"{allele1}/{allele2}",
        "n_liftover_mappings": int(n_liftover_mappings),
        "hg38_ref": None,
        "hg38_alt": None,
        "orientation": None,
        "allele_swap": None,
        "palindromic": False,
        "is_snv": len(allele1) == 1 and len(allele2) == 1,
        "variant_uid": None,
        "mapping_status": "excluded",
        "exclusion_reason": None,
    }
    if n_liftover_mappings == 0:
        row["exclusion_reason"] = "liftover_failed"
        return row
    if n_liftover_mappings > 1:
        row["exclusion_reason"] = "multimapped"
        return row
    if chrom not in STANDARD_CHROMS:
        row["exclusion_reason"] = "nonstandard_chrom"
        return row
    if not row["is_snv"]:
        # The Atlas point-query API refuses indels, so this stays out of the query set. It still gets an
        # identity: without a uid every deferred indel of a signal shares the empty key and a posterior
        # dict merges them.
        row.update(resolve_indel_alleles(chrom, position_1based, allele1, allele2, fetch))
        return row
    a1, a2 = allele1.upper(), allele2.upper()
    base = fetch(chrom, position_1based - 1, position_1based).upper()
    row["palindromic"] = a2 == COMP.get(a1)
    if base == a1:
        row.update(hg38_ref=a1, hg38_alt=a2, orientation="forward", allele_swap="none")
    elif base == a2:
        row.update(hg38_ref=a2, hg38_alt=a1, orientation="forward", allele_swap="swapped")
    elif base == COMP.get(a1):
        row.update(hg38_ref=COMP[a1], hg38_alt=COMP[a2], orientation="reverse_complement", allele_swap="none")
    elif base == COMP.get(a2):
        row.update(hg38_ref=COMP[a2], hg38_alt=COMP[a1], orientation="reverse_complement", allele_swap="swapped")
    else:
        row["exclusion_reason"] = "hg38_allele_validation_failed"
        return row
    row["mapping_status"] = "mapped"
    row["variant_uid"] = variant_uid(chrom, position_1based, row["hg38_ref"], row["hg38_alt"])
    return row


# ------------------------------------------------------------------ blocks
def coarse_blocks(anchors: Iterable[tuple[str, str, int]]) -> dict[str, str]:
    """1-Mb single-linkage components per chromosome (copied rule, registry v3)."""
    by_chrom: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for key, chrom, position in anchors:
        by_chrom[str(chrom).replace("chr", "")].append((int(position), str(key)))
    out: dict[str, str] = {}
    for chrom, items in by_chrom.items():
        items.sort()
        component = 1
        previous: int | None = None
        for position, key in items:
            if previous is not None and position - previous > LOCUS_WINDOW_BP:
                component += 1
            out[key] = f"chr{chrom}:component{component:04d}"
            previous = position
    return out


# ------------------------------------------------------------------ posterior rules
def select_query_variants(weights: dict[str, float], floor: float, cumulative: float, cap: int) -> tuple[set[str], float]:
    """Prespecified query floor: keep w>=floor ∪ smallest set reaching `cumulative`, capped."""
    ordered = sorted(weights.items(), key=lambda kv: (-kv[1], kv[0]))
    selected: set[str] = set()
    running = 0.0
    for key, w in ordered:
        if len(selected) >= cap:
            break
        if w >= floor or running < cumulative:
            selected.add(key)
            running += w
        else:
            break
    excluded = sum(w for k, w in weights.items() if k not in selected)
    return selected, excluded


def assert_signal_mass(rows: Iterable[tuple[str, str, float]], tolerance: float) -> dict[str, float]:
    """rows = (signal_uid, variant_uid, weight). Duplicated joins inflate mass → raise."""
    seen: set[tuple[str, str]] = set()
    mass: dict[str, float] = defaultdict(float)
    for signal, variant, w in rows:
        if (signal, variant) in seen:
            raise ContractError(f"duplicated variant within signal: {signal} {variant}")
        seen.add((signal, variant))
        mass[signal] += float(w)
    for signal, total in mass.items():
        if total > 1 + tolerance:
            raise ContractError(f"posterior mass exceeds one for {signal}: {total}")
    return dict(mass)


def posterior_summary(weights: dict[str, float], scores: dict[str, float], is_signed: bool, complete_cutoff: float = 0.95) -> dict:
    """§5: coverage C, signed S (NaN when unsigned), magnitude A, missing mass."""
    coverage = signed = magnitude = 0.0
    for key, w in weights.items():
        if key in scores and scores[key] is not None and not math.isnan(scores[key]):
            coverage += w
            signed += w * scores[key]
            magnitude += w * abs(scores[key])
    total = sum(weights.values())
    return {
        "coverage": coverage,
        "signed": signed if is_signed else math.nan,
        "magnitude": magnitude,
        "signed_over_coverage": (signed / coverage) if (is_signed and coverage > 0) else math.nan,
        "magnitude_over_coverage": (magnitude / coverage) if coverage > 0 else math.nan,
        "missing_mass": total - coverage,
        "complete_0.95": coverage >= complete_cutoff,
        "conditional": coverage < total - 1e-12,
    }


def require_signed(scorer: str, is_signed: bool) -> None:
    if not is_signed:
        raise ContractError(f"unsigned scorer passed to a directional analysis: {scorer}")


# ------------------------------------------------------------------ tracks
def track_class(ontology_curie: str, biosample_name: str, biosample_type: str) -> str:
    name = str(biosample_name)
    if HEPG2_RE.search(name):
        return "HepG2"
    if ontology_curie in HEPATOCYTE_ONTOLOGY or re.search(r"hepatocyte", name, re.IGNORECASE):
        return "hepatocyte"
    if ontology_curie in LIVER_ONTOLOGY or LIVER_NAME_RE.search(name):
        return "primary_liver"
    return "other"


# ------------------------------------------------------------------ the archived 211-variant AlphaGenome gate set
def gate_set(truth_rows, ag_rows, fetch) -> list[dict]:
    """Rebuild the 65b gate set with hg38 alleles from the 65-series score file, FASTA-validated, oriented to the eQTL effect allele.

    truth_rows: broadaway_benchmark_truth.tsv rows (hg19-oriented ref/alt; effect_allele; eqtl_sign).
    ag_rows: alphagenome_eqtl_scores.tsv rows (hg38_ref/hg38_alt resolved by script 65 against the FASTA).
    fetch(chrom, start0, end0) returns the reference base. Indels are excluded (the Atlas gate is SNV-only);
    a lead without an hg38 row, a FASTA mismatch, or an effect allele outside {hg38_ref, hg38_alt} raises.
    """
    def strip(v):
        return str(v).split(".")[0]
    ag = {(r["variant_id"], strip(r["ensembl"])): r for r in ag_rows}
    out = []
    for r in truth_rows:
        if r.get("is_signal_lead") != "TRUE" or r.get("strand_ambiguous") == "TRUE":
            continue
        key = (r["variant_id"], strip(r["ensembl"]))
        a = ag.get(key)
        if a is None:
            raise ContractError(f"gate lead without an hg38-resolved row in the 65-series score file: {key}")
        ref, alt = str(a["hg38_ref"]).upper(), str(a["hg38_alt"]).upper()
        if a.get("is_indel") == "TRUE" or len(ref) != 1 or len(alt) != 1:
            continue
        chrom, pos = f"chr{r['chr']}", int(r["pos_hg38"])
        base = str(fetch(chrom, pos - 1, pos)).upper()
        if base != ref:
            raise ContractError(f"gate variant {chrom}:{pos} hg38_ref {ref} != FASTA {base}")
        ea = str(r["effect_allele"]).upper()
        if ea == alt:
            sign = 1
        elif ea == ref:
            sign = -1
        else:
            raise ContractError(f"gate effect allele {ea} not in hg38 alleles {ref}/{alt} at {chrom}:{pos}")
        out.append({"variant_uid": variant_uid(chrom, pos, ref, alt), "ensembl": key[1], "effect_allele": ea, "hg38_ref": ref, "hg38_alt": alt,
                    "orientation_sign": sign, "label": 1 if float(r["eqtl_sign"]) > 0 else 0, "variant_id_hg19": r["variant_id"],
                    "allele_note": str(a.get("note", ""))})
    return out
