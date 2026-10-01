#!/usr/bin/env python3
"""C3b: per-slide image features for Model C (PRESPEC_C section 5, Amendment 1).

For each C1-matched GTEx liver slide (IDC DICOM, level 0 = 0.494 um/px, 20x):
- tissue mask on the lowest-resolution level: Otsu on HSV saturation, holes filled
  so fat vacuoles count as tissue;
- fat-vacuole fraction on the 4x-downsampled level (~2 um/px): white (min RGB > 200,
  saturation < 0.1) connected regions inside tissue with equivalent diameter
  5-100 um and eccentricity < 0.9, as a fraction of tissue area;
- tiles: 224 px at level 0 on a regular grid, kept if >= 50% tissue; at most 2,000
  per slide, sampled with a seed derived from the donor ID;
- Phikon-v2 CLS embeddings (1,024-d), fp16 on GPU; mean-pooled per slide; per-tile
  embeddings saved as float16 for the secondary MIL model.
Sealed donors' outputs go to <out>/sealed/ (mode 0700) and are not read by
development fits. No RNA value is read here.
"""
import argparse
import glob
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import openslide
import pandas as pd
import torch
from scipy import ndimage
from skimage import color, filters, measure, morphology
from transformers import AutoImageProcessor, AutoModel

TILE = 224
MAX_TILES = 2000
MODEL = "owkin/phikon-v2"
REVISION = "2ae989a9c40cffaa27f0a6cb29cc94d1d6f9a5fd"


def open_series(series_dir):
    files = sorted(glob.glob(os.path.join(series_dir, "*.dcm")))
    if not files:
        raise FileNotFoundError(series_dir)
    return openslide.OpenSlide(files[0])


def tissue_mask(slide):
    lvl = slide.level_count - 1
    img = np.asarray(slide.read_region((0, 0), lvl, slide.level_dimensions[lvl]).convert("RGB"))
    sat = color.rgb2hsv(img)[..., 1]
    mask = sat > filters.threshold_otsu(sat)
    mask = morphology.remove_small_objects(mask, 64)
    return ndimage.binary_fill_holes(mask), slide.level_downsamples[lvl]


def fat_fraction(slide, mask, mask_ds, mpp0):
    lvl = 1 if slide.level_count > 2 else 0
    ds = slide.level_downsamples[lvl]
    w, h = slide.level_dimensions[lvl]
    img = np.asarray(slide.read_region((0, 0), lvl, (w, h)).convert("RGB"))
    scale = mask_ds / ds
    tissue = np.kron(mask, np.ones((int(round(scale)), int(round(scale))), bool))[:h, :w]
    tissue = np.pad(tissue, ((0, max(0, h - tissue.shape[0])), (0, max(0, w - tissue.shape[1]))))
    hsv_s = color.rgb2hsv(img)[..., 1]
    white = (img.min(axis=2) > 200) & (hsv_s < 0.1) & tissue
    lab = measure.label(white)
    um = mpp0 * ds
    keep_area = 0
    for r in measure.regionprops(lab):
        d = r.equivalent_diameter_area * um
        if 5 <= d <= 100 and r.eccentricity < 0.9:
            keep_area += r.area
    return keep_area / max(1, tissue.sum()), float(tissue.sum() * um * um / 1e6)  # fraction, tissue mm^2


def tile_coords(slide, mask, mask_ds, seed):
    w0, h0 = slide.level_dimensions[0]
    xs, ys = np.arange(0, w0 - TILE, TILE), np.arange(0, h0 - TILE, TILE)
    coords = []
    for y in ys:
        for x in xs:
            mx0, my0 = int(x / mask_ds), int(y / mask_ds)
            mx1, my1 = max(mx0 + 1, int((x + TILE) / mask_ds)), max(my0 + 1, int((y + TILE) / mask_ds))
            if mask[my0:my1, mx0:mx1].mean() >= 0.5:
                coords.append((x, y))
    coords = np.array(coords)
    rng = np.random.default_rng(seed)
    if len(coords) > MAX_TILES:
        coords = coords[np.sort(rng.choice(len(coords), MAX_TILES, replace=False))]
    return coords


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matched", required=True, help="C1 liver_slides_matched.tsv")
    ap.add_argument("--seal", required=True)
    ap.add_argument("--slides", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out); (out / "tiles").mkdir(parents=True, exist_ok=True)
    (out / "sealed" / "tiles").mkdir(parents=True, exist_ok=True)
    os.chmod(out / "sealed", 0o700)
    m = pd.read_csv(a.matched, sep="\t")
    seal = pd.read_csv(a.seal, sep="\t").set_index("SUBJID")["sealed"]
    if a.limit:
        m = m.head(a.limit)

    dev = "cuda" if torch.cuda.is_available() else "cpu"
    proc = AutoImageProcessor.from_pretrained(MODEL, revision=REVISION)
    model = AutoModel.from_pretrained(MODEL, revision=REVISION).to(dev).eval()

    rows = {False: [], True: []}
    for _, r in m.iterrows():
        subj = r["PatientID"]; sealed = bool(seal[subj])
        base = out / "sealed" if sealed else out
        slide = open_series(os.path.join(a.slides, r["SeriesInstanceUID"]))
        mpp0 = float(slide.properties.get("openslide.mpp-x", 0.4942))
        mask, mask_ds = tissue_mask(slide)
        ff, tissue_mm2 = fat_fraction(slide, mask, mask_ds, mpp0)
        seed = int(hashlib.sha256(("MASLD-C-tiles" + subj).encode()).hexdigest()[:8], 16)
        coords = tile_coords(slide, mask, mask_ds, seed)
        embs = []
        for i in range(0, len(coords), 128):
            batch = [slide.read_region((int(x), int(y)), 0, (TILE, TILE)).convert("RGB") for x, y in coords[i:i + 128]]
            with torch.inference_mode(), torch.autocast(dev, dtype=torch.float16, enabled=dev == "cuda"):
                px = proc(images=batch, return_tensors="pt")["pixel_values"].to(dev)
                embs.append(model(pixel_values=px).last_hidden_state[:, 0].float().cpu().numpy())
        e = np.concatenate(embs) if embs else np.zeros((0, 1024), np.float32)
        np.savez_compressed(base / "tiles" / f"{subj}.npz", emb=e.astype(np.float16), coords=coords)
        rows[sealed].append({"SUBJID": subj, "SeriesInstanceUID": r["SeriesInstanceUID"], "n_tiles": len(coords),
                             "fat_fraction": ff, "tissue_mm2": tissue_mm2, "mpp0": mpp0,
                             **{f"emb_{j}": v for j, v in enumerate(e.mean(0) if len(e) else np.full(1024, np.nan))}})
        print(subj, "sealed" if sealed else "dev", len(coords), round(ff, 4), flush=True)
    for sealed, base in ((False, out), (True, out / "sealed")):
        pd.DataFrame(rows[sealed]).to_csv(base / "slide_features.tsv.gz", sep="\t", index=False)
    (out / "embedding_model.json").write_text(json.dumps({"model": MODEL, "revision": REVISION, "tile_px": TILE,
                                                          "max_tiles": MAX_TILES, "level0_um_per_px": "from slide"}, indent=2))


if __name__ == "__main__":
    main()
