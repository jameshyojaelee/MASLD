#!/usr/bin/env python
"""
parse_mirgenedb.py — Build L2b miRNA ortholog layer from MirGeneDB 2.1.

MirGeneDB 2.1 (Fromm et al. 2022 NAR) applies stricter orthology criteria than
miRBase. It uses a phylogeny-anchored naming convention:

    >Hsa-Let-7-P1b_5p          (Species-Family-Paralog_arm)
    >Hsa-Mir-122_5p
    >Hsa-Mir-15-P1a_5p

  - Species  : Hsa (Homo sapiens) / Mmu (Mus musculus)
  - Family   : Let-7, Mir-122, Mir-15, Mir-17, ...
  - Paralog  : P1a, P1b, P2a1, P3b (or empty — e.g. Mir-122 has no paralog ID)
  - Arm      : 5p or 3p

MirGeneDB pairs across species are defined by SAME family + SAME paralog ID.
This is the ortholog claim — much stricter than miRBase, because MirGeneDB
applies sequence + phylogeny + readcount support filters.

Strategy:
  1. Parse Hsa/Mmu mature FASTA from MirGeneDB.
  2. Pair Hsa ↔ Mmu by (family, paralog, arm). This is the orthology call.
  3. Map MirGeneDB hairpin name → miRBase hairpin name using cached MirGeneDB
     browse pages (the only place MirGeneDB publishes this cross-reference).
  4. From miRBase hairpin name (e.g. "hsa-mir-122"), derive HGNC/MGI candidate
     symbols (MIR122, Mir122) and match against GENCODE v49 (human) and
     GENCODE vM38 (mouse) miRNA gene metadata to land on Ensembl IDs.

Output common-schema TSV at data/external/orthologs/layers/L2b_mirgenedb.tsv.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CURATED = PROJECT / "data/external/orthologs/curated_dbs"
DEFAULT_HSA = CURATED / "mirgenedb_hsa_mat.fa"
DEFAULT_MMU = CURATED / "mirgenedb_mmu_mat.fa"
DEFAULT_OUT = PROJECT / "data/external/orthologs/layers/L2b_mirgenedb.tsv"
HUMAN_GTF_META = PROJECT / "data/gencode_v49_gene_metadata.tsv.gz"
MOUSE_GTF_META = PROJECT / "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv"
BROWSE_CACHE = CURATED / "mirgenedb_browse_cache"

# MirGeneDB header structure: Species-Family[-Paralog][-variant]_arm
# Examples seen in MirGeneDB 2.1:
#   Hsa-Mir-122_5p              (no paralog)
#   Hsa-Let-7-P1b_5p            (Paralog P1b)
#   Hsa-Let-7-P2a1_5p           (Paralog P2a1)
#   Hsa-Mir-10-P1b-v1_5p        (Paralog + variant v1)
#   Hsa-Mir-136-v1_3p           (no paralog + variant v1)
#   Hsa-Mir-219-P2-as_5p        (Paralog + antisense -as)
#   Hsa-Mir-337-as_5p           (no paralog + antisense -as)
# Family component: "Let-7" or "Mir-<number>" (digits, sometimes
#   followed by a single letter like 7a, 8b, 22a — but those don't occur
#   in current MirGeneDB hsa/mmu; family is always Mir-N with N numeric).
HEADER_RE = re.compile(
    r"^(?P<species>Hsa|Mmu)-"
    r"(?P<family>(?:Let|Mir)-[0-9]+[a-z]?)"                  # Let-7, Mir-1, Mir-122
    r"(?:-(?P<paralog>P[0-9A-Za-z]+))?"                       # optional Pxxx
    r"(?:-(?P<variant>v[0-9]+|as[0-9]*|o[0-9]+))?"            # optional -v1 / -as / -o1
    r"_(?P<arm>5p|3p)$"
)


def parse_mirgenedb_fa(fa_path: Path, species_code: str) -> pd.DataFrame:
    """Parse MirGeneDB mature FASTA. Returns DataFrame with:
    species, mature_name, hairpin_name, family, paralog, variant, arm, seq."""
    rows = []
    cur = None
    seq: list[str] = []
    unparsed = []
    with open(fa_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    rows.append((*cur, "".join(seq)))
                name = line[1:].split()[0]
                m = HEADER_RE.match(name)
                if not m:
                    unparsed.append(name)
                    cur = None
                    seq = []
                    continue
                d = m.groupdict()
                if d["species"] != species_code:
                    cur = None
                    seq = []
                    continue
                paralog = d["paralog"] or ""
                variant = d["variant"] or ""
                # Hairpin name = strip the arm suffix from mature name
                hairpin = f"{d['species']}-{d['family']}"
                if paralog:
                    hairpin += f"-{paralog}"
                if variant:
                    hairpin += f"-{variant}"
                cur = (d["species"], name, hairpin, d["family"],
                       paralog, variant, d["arm"])
                seq = []
            else:
                seq.append(line.strip())
        if cur is not None:
            rows.append((*cur, "".join(seq)))
    if unparsed:
        print(f"[mirgenedb] WARNING: could not parse {len(unparsed)} headers; "
              f"first 5: {unparsed[:5]}", file=sys.stderr)
    df = pd.DataFrame(rows, columns=["species", "mature_name", "hairpin_name",
                                     "family", "paralog", "variant", "arm",
                                     "seq"])
    print(f"[mirgenedb] {species_code}: parsed {len(df)} mature miRNAs from "
          f"{fa_path.name}", file=sys.stderr)
    return df


def build_orthologs(hsa: pd.DataFrame, mmu: pd.DataFrame) -> pd.DataFrame:
    """Pair Hsa ↔ Mmu by (family, paralog, variant, arm). Variant must
    match (Mir-N-v1 and Mir-N-v2 are distinct loci even within a species).
    Returns long-form pairs."""
    keys = ["family", "paralog", "variant", "arm"]
    h = hsa.rename(columns={
        "mature_name": "human_mature", "hairpin_name": "human_hairpin",
        "seq": "human_seq"
    })[keys + ["human_mature", "human_hairpin", "human_seq"]]
    m = mmu.rename(columns={
        "mature_name": "mouse_mature", "hairpin_name": "mouse_hairpin",
        "seq": "mouse_seq"
    })[keys + ["mouse_mature", "mouse_hairpin", "mouse_seq"]]
    pairs = h.merge(m, on=keys, how="inner")
    print(f"[mirgenedb] Hsa↔Mmu pairs by (family,paralog,variant,arm): "
          f"{len(pairs)}", file=sys.stderr)
    return pairs


def load_mirgenedb_to_mirbase_xref(cache_dir: Path) -> pd.DataFrame:
    """Parse the cached MirGeneDB browse pages to extract the published
    MirGeneDB hairpin name → miRBase hairpin name cross-reference.

    The browse table is the ONLY place MirGeneDB exposes this mapping;
    there is no downloadable mapping file.

    Subtlety: variant hairpins (Mir-10-P1b-v1) are listed in the browse
    table with their URL containing the variant suffix but the displayed
    text dropping it (e.g. <a href="/show/hsa/Mir-10-P1b-v1">Hsa-Mir-10-P1b</a>).
    We reconstruct the full hairpin name from the URL path, so that
    variant hairpins map correctly to their miRBase entries.

    Returns DataFrame with: mirgenedb_hairpin, mirbase_hairpin, mirbase_mi_id.
    """
    rows = []
    for species, fname in [("Hsa", "hsa_browse.html"),
                           ("Mmu", "mmu_browse.html")]:
        p = cache_dir / fname
        if not p.exists():
            print(f"[mirgenedb] WARNING: browse cache missing: {p}",
                  file=sys.stderr)
            continue
        with open(p) as f:
            text = f.read()
        sp_lower = species.lower()
        # Pattern: <a href="/show/{sp}/{family-paralog-variant}">{Sp}-...</a>
        # ... <a target="_blank" href="http://www.mirbase.org/hairpin/{MI...}">
        # {sp}-{let|mir}-{...}</a>
        # We capture the URL stub (after /show/{sp_lower}/) and use it to
        # reconstruct the full hairpin name including any variant suffix.
        pat = re.compile(
            rf'<a href="/show/{sp_lower}/([^"]+)">{species}-[A-Za-z0-9-]+</a>'
            r'.*?'
            rf'<a target="_blank" href="http://www\.mirbase\.org/hairpin/(MI\d+)">'
            rf'({sp_lower}-(?:let|mir)-[a-zA-Z0-9-]+)</a>',
            re.DOTALL
        )
        n_before = len(rows)
        for url_stub, miid, mb in pat.findall(text):
            # Reconstruct: Hsa-Mir-10-P1b-v1 from url_stub "Mir-10-P1b-v1"
            full = f"{species}-{url_stub}"
            rows.append((full, mb, miid))
        print(f"[mirgenedb] {species} browse xref: "
              f"{len(rows) - n_before} hairpin mappings", file=sys.stderr)
    df = pd.DataFrame(rows, columns=["mirgenedb_hairpin", "mirbase_hairpin",
                                     "mirbase_mi_id"])
    df = df.drop_duplicates()
    return df


def mirbase_hairpin_to_gencode_candidates(mb_hairpin: str) -> list[str]:
    """Convert a miRBase hairpin name (e.g. 'hsa-mir-122', 'mmu-let-7a-1')
    into candidate GENCODE gene symbols (e.g. ['MIR122','Mir122'], or
    ['MIRLET7A1','Mirlet7a1','MIRLET-7A1','Mirlet-7a1']).

    Convention used (HGNC/MGI):
      - hsa-mir-122       → MIR122          (HGNC)  / Mir122        (MGI)
      - hsa-mir-15a       → MIR15A          / Mir15a
      - hsa-mir-16-1      → MIR16-1         / Mir16-1   (cluster paralogs keep hyphen)
      - hsa-let-7a-1      → MIRLET7A1       / Mirlet7a1 (HGNC) or Mirlet7a-1 (MGI)
      - mmu-let-7c-1      → Mirlet7c-1      (MGI; keeps cluster hyphen)
      - hsa-mir-7a-2      → MIR7A2 / Mir7a-2 (MGI uses Mir7-2 — see mouse_split)
      - hsa-mir-126b      → MIR126B         / Mir126b
    """
    m = re.match(r"^[a-z]{3,4}-(.+)$", mb_hairpin.lower())
    if not m:
        return []
    rest = m.group(1)            # e.g. mir-122 / let-7a-1 / mir-16-1
    # Ordered list of candidates so we can pick the canonical-form match
    # first; the deduplication step preserves the first occurrence.
    candidates: list[str] = []

    def add(c: str) -> None:
        if c and c not in candidates:
            candidates.append(c)

    def emit(stem: str, body: str) -> None:
        """Add many capitalization+hyphenation variants for stem+body.
        Order matters: most canonical first."""
        stem_u = stem.upper()
        stem_c = stem.capitalize()
        nohyphen = body.replace("-", "")
        # Canonical/preferred forms first
        add(stem_u + nohyphen.upper())          # MIR16-1 → MIR161
        add(stem_u + body.upper())              # MIR16-1
        add(stem_c + nohyphen)                  # Mir161
        add(stem_c + body)                      # Mir16-1
        add(stem_u + "-" + body.upper())        # MIRLET-7A (rare)
        add(stem_c + "-" + body)                # Mirlet-7f
        # Strip trailing -N (cluster number) for MGI quirk:
        # mmu-mir-7a-2 → "Mir7-2" (MGI drops the 'a' since paralogs differ
        # only by cluster index).
        cluster = re.match(r"^([0-9]+)([a-z])(-[0-9]+)$", body)
        if cluster:
            no_letter = cluster.group(1) + cluster.group(3)  # 7-2
            add(stem_u + no_letter.upper())
            add(stem_c + no_letter)
        # MGI also drops a trailing single letter for single-locus families
        # (e.g. mmu-mir-146a → MGI "Mir146"). Add letter-stripped variant
        # as fallback — placed last so canonical "Mir146a" wins when it
        # exists in GENCODE.
        trail = re.match(r"^([0-9]+)([a-z])$", body)
        if trail:
            no_letter_body = trail.group(1)
            add(stem_u + no_letter_body)
            add(stem_c + no_letter_body)

    if rest.startswith("mir-"):
        body = rest[4:]
        emit("Mir", body)
    elif rest.startswith("let-"):
        body = rest[4:]
        emit("Mirlet", body)
    else:
        add(rest.upper())
        add(rest.capitalize())
    return candidates


def load_gencode_mirna(meta_path: Path, biotype_col: str,
                       symbol_col: str, id_col: str) -> pd.DataFrame:
    df = pd.read_csv(meta_path, sep=None, engine="python", dtype=str)
    df = df[df[biotype_col].astype(str) == "miRNA"].copy()
    df["symbol_norm"] = df[symbol_col].astype(str).str.upper()
    return df[["symbol_norm", symbol_col, id_col]].rename(
        columns={symbol_col: "sym", id_col: "ensembl"}
    )


def attach_gencode(pairs: pd.DataFrame, xref: pd.DataFrame) -> pd.DataFrame:
    """Annotate each MirGeneDB pair with miRBase hairpin name + GENCODE
    ensembl IDs on both sides."""
    # Build xref maps per species
    h_xref = xref[xref["mirgenedb_hairpin"].str.startswith("Hsa-")]
    m_xref = xref[xref["mirgenedb_hairpin"].str.startswith("Mmu-")]
    pairs = pairs.merge(
        h_xref.rename(columns={
            "mirgenedb_hairpin": "human_hairpin",
            "mirbase_hairpin": "human_mirbase_hairpin",
            "mirbase_mi_id": "human_mirbase_mi_id",
        }), on="human_hairpin", how="left")
    pairs = pairs.merge(
        m_xref.rename(columns={
            "mirgenedb_hairpin": "mouse_hairpin",
            "mirbase_hairpin": "mouse_mirbase_hairpin",
            "mirbase_mi_id": "mouse_mirbase_mi_id",
        }), on="mouse_hairpin", how="left")

    hsa_gc = load_gencode_mirna(HUMAN_GTF_META, "gene_biotype",
                                "gene_name", "ensembl_base")
    mmu_gc = load_gencode_mirna(MOUSE_GTF_META, "mouse_biotype",
                                "mouse_symbol_gtf", "mouse_ensembl_base")

    def match_one(mb_hairpin, gc: pd.DataFrame) -> tuple[str, str]:
        if pd.isna(mb_hairpin) or mb_hairpin in ("", "nan"):
            return ("", "")
        cands = mirbase_hairpin_to_gencode_candidates(str(mb_hairpin))
        norm = [c.upper() for c in cands]
        hit = gc[gc["symbol_norm"].isin(norm)]
        if hit.empty:
            return ("", "")
        # Prefer exact case match (HGNC over heuristic); else first.
        # Order candidates by "preferred" (HGNC uppercase first for human,
        # MGI capitalized first for mouse) — pick first matching candidate.
        for c in cands:
            row = hit[hit["sym"] == c]
            if not row.empty:
                return (row.iloc[0]["sym"], row.iloc[0]["ensembl"])
        return (hit.iloc[0]["sym"], hit.iloc[0]["ensembl"])

    h_sym, h_eid = [], []
    for hp in pairs["human_mirbase_hairpin"]:
        s, e = match_one(hp, hsa_gc)
        h_sym.append(s); h_eid.append(e)
    m_sym, m_eid = [], []
    for hp in pairs["mouse_mirbase_hairpin"]:
        s, e = match_one(hp, mmu_gc)
        m_sym.append(s); m_eid.append(e)
    pairs["human_symbol"] = h_sym
    pairs["human_ensembl"] = h_eid
    pairs["mouse_symbol"] = m_sym
    pairs["mouse_ensembl"] = m_eid
    return pairs


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hsa", type=Path, default=DEFAULT_HSA)
    p.add_argument("--mmu", type=Path, default=DEFAULT_MMU)
    p.add_argument("--browse-cache", type=Path, default=BROWSE_CACHE)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()

    print(f"[mirgenedb] parsing {args.hsa}", file=sys.stderr)
    hsa = parse_mirgenedb_fa(args.hsa, "Hsa")
    print(f"[mirgenedb] parsing {args.mmu}", file=sys.stderr)
    mmu = parse_mirgenedb_fa(args.mmu, "Mmu")

    pairs = build_orthologs(hsa, mmu)

    print(f"[mirgenedb] loading miRBase xref from {args.browse_cache}",
          file=sys.stderr)
    xref = load_mirgenedb_to_mirbase_xref(args.browse_cache)
    print(f"[mirgenedb] total xref hairpin mappings: {len(xref)}",
          file=sys.stderr)

    print("[mirgenedb] attaching GENCODE Ensembl IDs...", file=sys.stderr)
    pairs = attach_gencode(pairs, xref)

    # Stats
    both = ((pairs["human_ensembl"] != "") & (pairs["mouse_ensembl"] != "")).sum()
    h_only = ((pairs["human_ensembl"] != "") & (pairs["mouse_ensembl"] == "")).sum()
    m_only = ((pairs["human_ensembl"] == "") & (pairs["mouse_ensembl"] != "")).sum()
    neither = ((pairs["human_ensembl"] == "") & (pairs["mouse_ensembl"] == "")).sum()
    print(f"[mirgenedb]   pairs with both Ensembl IDs:  {both}", file=sys.stderr)
    print(f"[mirgenedb]   pairs with human-only:        {h_only}", file=sys.stderr)
    print(f"[mirgenedb]   pairs with mouse-only:        {m_only}", file=sys.stderr)
    print(f"[mirgenedb]   pairs with neither (symbol-only): {neither}",
          file=sys.stderr)

    # Common-schema output
    out = pd.DataFrame({
        "mouse_ensembl": pairs["mouse_ensembl"].replace("", pd.NA),
        "mouse_symbol": pairs["mouse_symbol"],
        "mouse_biotype": "miRNA",
        "human_ensembl": pairs["human_ensembl"].replace("", pd.NA),
        "human_symbol": pairs["human_symbol"],
        "human_biotype": "miRNA",
        "tier_H_biomart": 0,
        "tier_H_mirbase": 0,
        "tier_H_mirgenedb": 1,
        "tier_M_phasej_synteny": 0,
        "tier_M_lncbook": 0,
        "tier_M_blast": 0,
        "tier_M_pseudogene_parent": 0,
        "tier_L_liftover": 0,
        "confidence_tier": "H",
        "evidence_count": 1,
        "provenance_sources": "mirgenedb_2.1",
        "mirgenedb_family": pairs["family"],
        "mirgenedb_paralog": pairs["paralog"],
        "mirgenedb_variant": pairs["variant"],
        "mirgenedb_arm": pairs["arm"],
        "mirgenedb_human_hairpin": pairs["human_hairpin"],
        "mirgenedb_mouse_hairpin": pairs["mouse_hairpin"],
        "mirgenedb_human_mature": pairs["human_mature"],
        "mirgenedb_mouse_mature": pairs["mouse_mature"],
        "mirgenedb_human_mirbase_hairpin": pairs["human_mirbase_hairpin"],
        "mirgenedb_mouse_mirbase_hairpin": pairs["mouse_mirbase_hairpin"],
        "notes": "MirGeneDB 2.1 strict family+paralog+arm orthology",
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[mirgenedb] wrote {len(out)} miRNA pairs to {args.out}",
          file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
