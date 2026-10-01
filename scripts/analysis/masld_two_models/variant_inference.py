#!/usr/bin/env python3
"""Experimental fresh inference of a variant-centered liver effect.

This interface computes from GRCh38 sequence. It does not use an arbitrary
target peak: the native readout is centered on the variant, and the adapter
pools a 384-bp central embedding interval. New-variant accuracy beyond the
Currin significance-selected lead population has not been established.
Archived evaluation metrics refer to saved predictions: fresh processes failed
the 1e-4 numerical reproduction check even on the same L40S GPU.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/"scripts/analysis/alphagenome_program"))
sys.path.insert(0, str(ROOT/"scripts/analysis/alphagenome_campaign"))
import i1_common as C

BASE = ROOT/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
SEEDS = tuple(range(20260917, 20260922))
RECIPE = "adapter_r16_last5_symmetric"
VARIANT = re.compile(r"^GRCh38:(chr(?:[1-9]|1[0-9]|2[0-2])):([1-9][0-9]*):([ACGT]):([ACGT])$")


def parse_variant(value):
    matched = VARIANT.fullmatch(value)
    if matched is None:
        raise ValueError("Expected GRCh38:chr1..chr22:1-based-position:REF:ALT SNV")
    chrom, position, ref, alt = matched.groups()
    if ref == alt:
        raise ValueError("REF and ALT must differ")
    return {"chr": chrom, "pos_hg38": int(position), "ref": ref, "alt": alt,
            "key": C.key_of(chrom, position, ref, alt)}


def development_fold(chrom):
    # Only genomic metadata are parsed, including when the source has fold 0.
    source = pd.read_csv(C.C2_LABELS, sep="\t", usecols=["chr", "heldout_fold"])
    mapping = source.groupby("chr").heldout_fold.agg(lambda values: sorted(set(values)))
    if chrom not in mapping or len(mapping[chrom]) != 1:
        raise ValueError("Chromosome has no unique recorded development fold")
    fold = int(mapping[chrom][0])
    if fold not in range(1, 5):
        raise ValueError("Chromosome belongs to spent fold 0; this interface excludes it")
    return fold


def require_gpu():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("A Slurm GPU compute allocation is required")
    card = C.gpu_card()
    if card.get("card_tag") != "l40s":
        raise RuntimeError("The evaluated inference route requires L40S hardware")
    return card


def window(fasta, row, length):
    lo, hi = C.window_bounds(int(row["pos_hg38"]), length)
    check = C.validate_window(fasta, row["chr"], lo, hi, int(row["pos_hg38"]), row["ref"])
    if not check["acgt_ok"] or not check["reference_match"]:
        raise ValueError("Unsupported sequence window: " + json.dumps(check, default=str))
    return C.extract_ref_alt(fasta, row["chr"], lo, hi, int(row["pos_hg38"]), row["ref"], row["alt"]), lo, hi


class NativePredictor:
    """Single-variant, two-allele native ATAC inference with the saved readout."""
    def __init__(self, atac_only=False):
        import jax
        import pysam
        from alphagenome.models import dna_client
        self.card = require_gpu()
        tick = time.monotonic()
        self.model, self.device, self.restore_seconds = C.load_model()
        # These fields were checked in the installed AlphaGenomeModel source.
        jax.block_until_ready((self.model._params, self.model._state))
        self.fasta = pysam.FastaFile(C.FASTA_PATH)
        # Preserve the saved scorer's compiled graph by default. Dropping DNase
        # changes a static JIT argument and requires numerical reproduction.
        self.outputs = [dna_client.OutputType.ATAC]
        if not atac_only:
            self.outputs.append(dna_client.OutputType.DNASE)
        self.requested_output_names = [output.name for output in self.outputs]
        self.scorer, _ = C.atac_dnase_scorers()
        self.make_mask, self.aggregate, _, _ = C.center_mask_and_aggregation()
        self.initialization_seconds = time.monotonic()-tick

    def predict(self, row):
        import jax.numpy as jnp
        from alphagenome.data import genome
        sequences, lo, hi = window(self.fasta, row, 1048576)
        interval = genome.Interval(chromosome=row["chr"], start=lo, end=hi)
        variant = genome.Variant(chromosome=row["chr"], position=int(row["pos_hg38"]),
                                 reference_bases=row["ref"], alternate_bases=row["alt"])
        predictions = [self.model.predict_sequence(sequence=seq, requested_outputs=self.outputs,
                         ontology_terms=C.LIVER_TERMS, interval=interval).atac for seq in sequences]
        ref, alt = predictions
        if ref is None or alt is None or ref.resolution != 1 or alt.resolution != 1:
            raise ValueError("Expected 1-bp native ATAC output")
        names = ref.metadata["name"].astype(str).to_numpy()
        if not np.array_equal(names, alt.metadata["name"].astype(str).to_numpy()):
            raise ValueError("Native allele track identities differ")
        indices = []
        for name in C.AG_ATAC_TRACKS:
            hit = np.flatnonzero(names == name)
            if len(hit) != 1:
                raise ValueError("Native liver track identity differs: " + name)
            indices.append(int(hit[0]))
        mask = self.make_mask(interval.as_unstranded(), variant, width=501, resolution=1)
        if int(mask.sum()) != 501:
            raise ValueError("Native center mask has unexpected width")
        values = self.aggregate(jnp.asarray(np.asarray(ref.values, dtype=np.float32)),
                                jnp.asarray(np.asarray(alt.values, dtype=np.float32)),
                                jnp.asarray(mask), aggregation_type=self.scorer.aggregation_type)
        # Conversion synchronizes the device before returning to a timing caller.
        result = float(np.asarray(values, dtype=np.float64)[indices].mean())
        if not np.isfinite(result):
            raise ValueError("Nonfinite native effect")
        return result


class AdapterPredictor:
    """Five saved seeds, with one variant per batch and identical allele routes."""
    def __init__(self, fold):
        import jax
        import jax.numpy as jnp
        import pysam
        from alphagenome_research.model import one_hot_encoder
        from i1_extract_mpra_embeddings import load_trunk
        import model_scalar as M
        if fold not in range(1, 5):
            raise ValueError("Only development folds 1-4 are available")
        self.card = require_gpu()
        tick = time.monotonic()
        self.fold = fold
        self.trunk, self.base, self.state, self.device, self.restore_seconds = load_trunk()
        self.fasta = pysam.FastaFile(C.FASTA_PATH)
        self.encoder = one_hot_encoder.DNAOneHotEncoder()
        self.weights = []
        self.weight_sha256 = {}
        folder = BASE/f"nested_f{fold}_{21783074+fold}"
        paths = M.qv_paths(self.base, 5)
        for seed in SEEDS:
            directory = folder/f"{RECIPE}__seed{seed}"
            split = json.loads((directory/"split.json").read_text())
            if (split["seed"] != seed or split["validation_fold"] != fold or
                    split["held_fold"] != 0 or split["training_folds"] != sorted(set(range(1, 5))-{fold})):
                raise ValueError("Saved adapter split or seed differs")
            if json.loads((directory/"feasibility.json").read_text())["completed_steps"] != 5000:
                raise ValueError("Saved adapter has incomplete training")
            weight_path = directory/"weights.npz"
            metadata = json.loads(Path(str(weight_path)+".json").read_text())
            config = metadata["config"]
            if (str(metadata["checkpoint"]) != str(C.CHECKPOINT) or config["mode"] != "adapter" or
                    config["length"] != 2048 or config["pooling"] != "symmetric" or
                    config["rank"] != 16 or config["last_blocks"] != 5 or
                    [tuple(item) for item in metadata["paths"]] != paths):
                raise ValueError("Saved adapter geometry or checkpoint differs")
            self.weights.append(M.load_trainable(weight_path))
            self.weight_sha256[str(weight_path)] = hashlib.sha256(weight_path.read_bytes()).hexdigest()
        self.pool = jnp.asarray(M.pool_weights(2048, "symmetric")[None], dtype=jnp.float32)

        @jax.jit
        def effect(weights, base, state, ref, alt, pool):
            return M.sequence_effect(self.trunk, base, state, weights, paths, "adapter", ref, alt, pool)

        self.effect = effect
        jax.block_until_ready((self.base, self.state, self.weights))
        self.initialization_seconds = time.monotonic()-tick

    def predict_seeds(self, row):
        import jax.numpy as jnp
        sequences, _, _ = window(self.fasta, row, 2048)
        ref, alt = [jnp.asarray(self.encoder.encode(seq)[None], dtype=jnp.float32) for seq in sequences]
        values = np.array([float(np.asarray(self.effect(weights, self.base, self.state, ref, alt, self.pool))[0])
                           for weights in self.weights])
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite adapter effect")
        return values

    def predict(self, row):
        return float(self.predict_seeds(row).mean())


def calibration(path, fold):
    result = json.loads(path.read_text())
    weights = [entry for entry in result["weights"] if entry["held_fold"] == fold]
    if len(weights) != 1:
        raise ValueError("No unique frozen calibration for requested chromosome fold")
    if 0 in weights[0]["meta_folds"] or fold in weights[0]["meta_folds"]:
        raise ValueError("Calibration folds contain a prohibited evaluation fold")
    return weights[0]


def main(args):
    if args.out.exists():
        raise FileExistsError(args.out)
    row = parse_variant(args.variant)
    fold = development_fold(row["chr"])
    if args.assay != "liver-caQTL":
        raise ValueError("Only the liver-caQTL development prediction is supported")
    if args.model != "adapter" and args.calibration_results is None:
        raise ValueError("Frozen combination results are required for native beta calibration")
    cal = calibration(args.calibration_results, fold) if args.calibration_results else None
    values = {}
    if args.model in ("native1m", "combined"):
        native = NativePredictor()
        values["native_ATAC_log2_sum_difference"] = native.predict(row)
        values["native1m_beta"] = cal["native1m_slope"]*values["native_ATAC_log2_sum_difference"]
        if args.model == "combined":
            import gc
            import jax
            native.fasta.close()
            del native
            gc.collect()
            jax.clear_caches()
    if args.model in ("adapter", "combined"):
        adapter = AdapterPredictor(fold)
        seeds = adapter.predict_seeds(row)
        values["adapter_beta"] = float(seeds.mean())
        values["adapter_seed_predictions"] = seeds.tolist()
        values["adapter_seed_sd"] = float(seeds.std(ddof=1))
    if args.model == "combined":
        values["combined_beta"] = cal["native1m_weight"]*values["native1m_beta"]+cal["adapter_weight"]*values["adapter_beta"]
    center = row["pos_hg38"]-1
    result = {"variant": args.variant, "model": args.model, "development_fold": fold,
              "assay": "Currin liver-caQTL association-effect prediction", "allele_contrast": "ALT minus REF",
              "beta_units": "source FastQTL ALT-dosage beta", "prediction": values,
              "beta_interpretation": "Association coefficient per additional ALT allele, not an isolated causal allele-editing effect",
              "raw_native_effect_unit": "mean across the three named liver ATAC tracks of log2(1+ALT track sum) minus log2(1+REF track sum) over the 501-bp center mask",
              "native_target_interval_0based_halfopen": f"{row['chr']}:{center-250}-{center+251}",
              "adapter_pool_interval_0based_halfopen": f"{row['chr']}:{center-192}-{center+192}",
              "input_length_bp": {"native1m": 1048576, "adapter": 2048},
              "arbitrary_target_interval_supported": False, "target_gene_evidence": "unsupported",
              "new_variant_generalization": "unvalidated outside significance-selected Currin leads",
              "uncertainty": "seed spread is fit variation; it excludes runtime numerical variation and is not a calibrated predictive interval",
              "inference_status": "experimental fresh inference; archived evaluation metrics apply to saved predictions",
              "runtime_reproducibility": {
                  "status": "failed 1e-4 numerical reproduction on 32 development variants on the same L40S",
                  "within_process_repeats": "exact in the bounded diagnostic",
                  "between_process_native_max_absolute_raw_ATAC_difference": 0.013058662414550795,
                  "between_process_adapter_batch1_ensemble_max_absolute_beta_difference": 0.00571859,
                  "between_process_adapter_batch4_ensemble_max_absolute_beta_difference": 0.00347379,
                  "scope": "observed diagnostic differences, not tolerance bounds or future guarantees",
                  "source": "two-models-20260922T203900EDT/reproduction_discriminator_21849377/results.json"},
              "measured_outcomes_loaded": False, "clinical_or_MASLD_effect": "unsupported",
              "weights_state": "fixed three-development-fold research models; not an adopted release"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--assay", default="liver-caQTL")
    parser.add_argument("--model", choices=["native1m", "adapter", "combined"], default="adapter")
    parser.add_argument("--calibration-results", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
