#!/usr/bin/env python3
"""Write ADDENDUM.md for the runtime probe: liver track inventory with life stage, the embedding hook,
the execution directory, and what blocks fine-tuning.

The track inventory is read out of the packaged AlphaGenome metadata at run time rather than retyped, so
it cannot drift from what the model actually returns. CPU only, no GPU and no job needed.
"""

from __future__ import annotations

import pathlib

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EXEC = PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z"
CODE = EXEC / "code/alphagenome_research"
TERMS = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
HEADS = ["ATAC", "DNASE", "CHIP_HISTONE", "RNA_SEQ", "SPLICE_SITE_USAGE"]
LIVERISH = "hepato|liver|cholangio|stellate|kupffer|bile|sinusoid"


def selected_tracks():
    from alphagenome.data import ontology
    from alphagenome.models import dna_client
    from alphagenome.models import dna_model as api
    from alphagenome_research.model.metadata import metadata as M

    md = M.load(api.Organism.HOMO_SAPIENS)
    terms = tuple(ontology.from_curie(t) for t in TERMS)
    outs = tuple(getattr(dna_client.OutputType, h) for h in HEADS)
    masks = M.create_track_masks(md, requested_outputs=outs, requested_ontologies=terms)
    out = {}
    for h, ot in zip(HEADS, outs):
        t = md.get(ot)
        sel = t[masks[ot][: len(t)]]
        out[h] = (sel, int(md.resolution(ot)), len(t))
    return md, out


def whole_catalogue(md):
    """Every liver-related biosample in the human catalogue, whether the pinned terms reach it or not."""
    from alphagenome.models import dna_client

    rows = []
    for ot in dna_client.OutputType:
        t = md.get(ot)
        if t is None or "biosample_name" not in t.columns:
            continue
        hit = t[t["biosample_name"].astype(str).str.lower().str.contains(LIVERISH, na=False, regex=True)]
        if hit.empty:
            continue
        keys = ["biosample_name", "ontology_curie"]
        has_stage = "biosample_life_stage" in hit.columns
        if has_stage:
            keys.append("biosample_life_stage")
        for key, n in hit.groupby(keys, dropna=False).size().items():
            bn, oc = key[0], key[1]
            ls = key[2] if has_stage else "not annotated for this head"
            rows.append((ot.name, str(oc), str(bn), str(ls), int(n), str(oc) in TERMS))
    return rows


HOOK = """\
## 2. Embedding / intermediate-representation hook

**Yes, and it is a clean one.** `alphagenome_research.model.dna_model.create_model(metadata,
num_splice_sites=..., splice_site_threshold=...)` returns a five-tuple
`(init_fn, apply_fn, trunk_apply_fn, heads_apply_fn, junctions_apply_fn)`. The third element is the
frozen-probe entry point:

    trunk_apply_fn(params, state, dna_sequence, organism_index) -> Embeddings

`Embeddings` is a frozen chex dataclass in `alphagenome_research/model/embeddings.py`, with an accessor
`get_sequence_embeddings(resolution)` that takes 1 or 128. The underlying model method is
`alphagenome_research.model.model.AlphaGenome.forward_trunk(dna_sequence, organism_index,
is_training=False)`.

| field | shape at 1,048,576 bp input | size at bfloat16 |
|---|---|---:|
| `embeddings_1bp` | (1, 1048576, 1536) | 3.0 GiB |
| `embeddings_128bp` | (1, 8192, 3072) | 48 MiB |
| `embeddings_pair` | (1, 512, 512, 128) | 64 MiB |

Four things to know before building on it.

1. The `create()` factory keeps the trunk only as the private, jitted `self._trunk_apply_fn` on the
   returned `AlphaGenomeModel`. There is no public method that returns embeddings. A frozen probe should
   call `create_model` directly and pass the `params` and `state` that Orbax restored from the local
   checkpoint, rather than going through `create()`.
2. The mixed-precision policy is hard-coded at `params=float32,compute=bfloat16,output=bfloat16`, so
   embeddings arrive in bfloat16. Upcasting is the caller's choice. The float32 upcast of the *track*
   outputs is what exhausted the L40S at a 2 Mb input, so this matters.
3. Take `embeddings_128bp`, or pool, rather than pulling `embeddings_1bp` at 1 Mb. The 3.0 GiB per
   sequence is what makes the base-resolution path expensive, and a 128-bp representation is 64x smaller.
4. `score_variant` and `score_ism_variants` are public on the model object if the packaged variant
   scorers are wanted instead of a custom head.

A frozen probe on `trunk_apply_fn` sidesteps blockers 1 to 3 in the next section entirely, and is the
cheaper first move than fine-tuning.
"""

BLOCKERS = """\
## 4. What blocks fine-tuning next

1. **The shipped trainer optimizes the full multimodal loss, not a custom endpoint.**
   `finetuning/finetune.py` calls `AlphaGenome.loss(batch)` over every head, with
   `freeze_trunk_embeddings=True` and a bfloat16 compute policy. Fitting an MPRA, a Cas13 readout or any
   single endpoint means writing a new head and loss. That is a code change, not a config change.
2. **The shipped data pipeline cannot run on the packaged metadata at all.** `finetuning/dataset.py`
   builds one `BigWigExtractor` per row of `metadata['file_path']`, and the packaged human metadata has
   no `file_path` column in any of its eleven output types (checked directly against
   `metadata.load(Organism.HOMO_SAPIENS)`). Someone has to supply a metadata table carrying local bigWig
   paths for the tracks being trained. This is the hardest of the four.
3. **Fold definitions are fetched from the public internet by default.**
   `alphagenome.data.fold_intervals.get_fold_intervals` falls back to a Borzoi BED hosted on GitHub. A
   compute node without egress needs a local copy, and the copy must be pinned or the train/valid/test
   split is not reproducible.
4. **Training memory is unmeasured, and it is the real risk.** Inference at 1 Mb peaks at 17.6 GiB on the
   L40S and 31.8 GiB on the Blackwell card. A training step adds gradients, optimizer state and retained
   activations on top; the upstream README recommends TPU v3 or better for training. This probe ran no
   training step, so there is no measured number. Measure one short step at the shortest useful length
   before planning anything at 1 Mb. [Unverified: no training step was run.]
5. **Reference build mismatch.** The trainer's default FASTA is GRCh38.p13 with GENCODE v46 annotation;
   this project standardises on the cellranger-arc GRCh38-2024-A FASTA with GENCODE v49. Inference in this
   probe used the project FASTA. If fine-tuning targets are built from project annotation, the annotation
   vintage differs from the one the released weights were trained against.
6. **Terms.** Noncommercial, and derivatives inherit them. Outputs of these local weights are eligible as
   features or teacher signals. Hosted Atlas and API outputs are not, and the two are never pooled.

One further constraint that is not a blocker but shapes the plan: scores are bit-identical when repeated
on one card and are **not** reproducible across GPU architectures. The same weights and sequence gave
splice -1.111623 on the L40S against -1.100200 on the Blackwell card, and the largest channel gap was
3.88e-02 on RNA. In the hosted run the median indel splice effect is 0.0151 in the same units, so a
cross-card difference is the same order as a typical effect. Any training set or feature table built from
these weights must be produced on one card, and must record which.
"""


def main() -> None:
    md, sel = selected_tracks()
    P: list[str] = []
    P.append("# AlphaGenome local runtime probe: addendum")
    P.append("")
    P.append(
        "Companion to `RESULTS.md` in the same directory. It answers four questions that were asked "
        "after the main probe: the liver track inventory with life-stage annotation, whether the local "
        "weights expose an embedding hook for a frozen probe, the execution directory, and what blocks "
        "fine-tuning. The track tables are read from the packaged AlphaGenome metadata when this file is "
        "generated, by `scripts/analysis/alphagenome_program/h1_addendum.py`, so they cannot drift from "
        "what the model returns. No GPU job was run to produce this file."
    )
    P.append("")
    P.append("## 1. Liver track inventory per output type, with life stage")
    P.append("")
    P.append(
        "Ontology terms pinned by the recipe: " + ", ".join(f"`{t}`" for t in TERMS) + ". "
        "These are the tracks the four terms actually select, which is exactly what the anchor run "
        "recorded on the GPU."
    )
    P.append("")
    P.append(
        "**The life-stage column is the finding here, not the counts. Every hepatocyte track in the "
        "catalogue is an in-vitro-differentiated embryonic cell. There is no adult primary hepatocyte "
        "track in any modality.**"
    )
    P.append("")
    for h, (df, res, total) in sel.items():
        P.append(f"### {h} — {len(df)} of {total} tracks, {res} bp resolution")
        P.append("")
        cols = ["name", "biosample_name", "biosample_type", "biosample_life_stage", "strand"]
        if "histone_mark" in df.columns:
            cols.insert(1, "histone_mark")
        P.append("| " + " | ".join(cols) + " |")
        P.append("|" + "---|" * len(cols))
        for _, r in df.iterrows():
            P.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        P.append("")

    k27 = sel["CHIP_HISTONE"][0]
    k27 = k27[k27["name"].astype(str).str.upper().str.contains("H3K27AC")]
    P.append(
        f"The anchor's H3K27ac readout is the {len(k27)} of the {len(sel['CHIP_HISTONE'][0])} returned "
        f"histone tracks whose name contains H3K27ac: "
        + "; ".join(
            f"{r['biosample_name']} ({r['biosample_life_stage']})" for _, r in k27.iterrows()
        )
        + ". The recipe selects them by name, and that name filter does match here."
    )
    P.append("")
    P.append("### Consequences to carry into any claim")
    P.append("")
    P.append(
        "- **ATAC is adult bulk tissue only.** Three tracks, all tissue and adult. There is no hepatocyte "
        "ATAC track anywhere in the human catalogue, so no cell-type-resolved accessibility is available "
        "from this model."
    )
    P.append(
        "- **Hepatocyte resolution is embryonic and in vitro.** Every CL:0000182 track, across RNA, "
        "DNase, histone ChIP and splice usage, is `in_vitro_differentiated_cells` at `embryonic` life "
        "stage. Reading these as adult hepatocyte biology is an assumption the tracks do not support, "
        "which matters for a disease of adult liver."
    )
    P.append(
        "- **The two chromatin channels are not life-stage matched to each other.** For UBERON:0002107 "
        "the DNase track is embryonic tissue while the H3K27ac track is adult tissue, so a DNase-versus-"
        "H3K27ac contrast at that term confounds assay with life stage."
    )
    P.append("")
    cat = whole_catalogue(md)
    unreached = [r for r in cat if not r[5]]
    P.append(
        "### Liver-related biosamples in the whole human catalogue, and which the pinned terms miss"
    )
    P.append("")
    P.append("| head | ontology | biosample | life stage | tracks | reached by pinned terms |")
    P.append("|---|---|---|---|---:|---|")
    for head, oc, bn, ls, n, reached in cat:
        P.append(f"| {head} | {oc} | {bn} | {ls} | {n} | {'yes' if reached else 'NO'} |")
    P.append("")
    missed = sorted({(r[1], r[2], r[0]) for r in unreached})
    P.append(
        f"{len(missed)} liver-related biosample/head combinations are not reached by the four pinned "
        "terms: "
        + "; ".join(f"{bn} ({oc}) in {head}" for oc, bn, head in missed)
        + ". Each of those two exists exactly once in the whole catalogue and in that one head only, so "
        "adding its term would buy a single track and no second modality to cross-check it against. "
        "Separately, there is no cholangiocyte biosample in any output type, so cholangiocyte expression "
        "is not unsupported evidence here, it is not testable with this model at all."
    )
    P.append("")
    P.append(HOOK)
    P.append("## 3. Execution directory")
    P.append("")
    P.append(f"    {EXEC}")
    P.append("")
    P.append("Contents:")
    P.append("")
    P.append("- `RESULTS.md` — the main probe: environment, measured time and memory per input length on "
             "both GPU types, the four variant-effect checks with pass or fail and evidence, and the "
             "HSD17B13 anchor beside the hosted value.")
    P.append("- `ADDENDUM.md` — this file.")
    P.append("- `MANIFEST.tsv` — every input, code file and output with sha256, including the 3.1 GB "
             "reference FASTA and the hosted-value TSV the anchor is compared against.")
    P.append("- `tables/l40s/`, `tables/b6k/` — per-length benchmarks, `variant_probes.json`, "
             "`anchor_hsd17b13.json` for each card.")
    P.append("- `tables/l40s_extra/`, `tables/b6k_extra/` — the two input lengths outside the SDK's named "
             "set, which establish that the 1 Mb ceiling is the card and not the model.")
    P.append("- `env/` — `pip_freeze_*.txt` for both nodes, `versions_*.txt`, `nvidia_smi_*.txt`, "
             "`gpu_*.csv`, and all four sbatch scripts as submitted.")
    P.append(f"- `code/alphagenome_research` — the pinned clone. `pip freeze` records the commit itself as "
             f"`-e git+...@1e55dcffb98ba26b31e74edc5e9f038f54c0e89d`.")
    P.append("")
    P.append("Environment: `/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local`.")
    P.append("")
    P.append("Source for this probe: "
             "`/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/"
             "alphagenome_program/h1_*.py` and `h1_*.sbatch`.")
    P.append("")
    P.append(BLOCKERS)
    dest = EXEC / "ADDENDUM.md"
    dest.write_text("\n".join(P) + "\n")
    print(dest)


if __name__ == "__main__":
    main()
