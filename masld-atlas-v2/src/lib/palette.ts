/**
 * Central color authority for every chart and data mark in the portal.
 *
 * Charts must import colors from HERE (not hardcode hex, not read CSS vars in
 * canvas/WebGL contexts) so the scientific color contract is enforced in one
 * place. Project doctrine baked in:
 *   - Control / healthy / non-significant is ALWAYS `#9E9E9E`.
 *   - Effect direction uses a diverging blue(down)/gray/red(up) scale, and
 *     significance is ALSO encoded by opacity (never color alone).
 *   - Ancestry uses the colorblind-safe Okabe-Ito assignment.
 *
 * Modality hex values are exact sRGB conversions of the OKLCH `--color-s*`
 * tokens in globals.css (kept in sync so SVG-via-CSS-var and canvas-via-hex
 * render identically).
 */

// ---------------------------------------------------------------------------
// Control / non-significant
// ---------------------------------------------------------------------------

/** The single canonical control / healthy / non-significant color. */
export const CONTROL = "#9E9E9E";

/** Opacity applied to a mark when it is NOT statistically significant. */
export const NONSIG_OPACITY = 0.4;
/** Opacity applied to a significant mark. */
export const SIG_OPACITY = 1;

// ---------------------------------------------------------------------------
// Diverging effect scale  (down = blue, mid = control gray, up = red)
// ---------------------------------------------------------------------------

export const DIVERGING = {
  down: "#2166AC", // downregulated / negative logFC
  mid: CONTROL,
  up: "#B2182B", // upregulated / positive logFC
} as const;

// ---------------------------------------------------------------------------
// Sequential magnitude ramp  (single-hue teal-blue, light -> deep)
// ---------------------------------------------------------------------------

export const SEQUENTIAL_RAMP = [
  "#e8f1f8",
  "#b3d1e6",
  "#7fb0d4",
  "#4a8fc2",
  "#1f6aa5",
  "#0d4a7a",
] as const;

// ---------------------------------------------------------------------------
// Perceptual heatmap ramp  (viridis anchors, low -> high)
// ---------------------------------------------------------------------------

export const HEATMAP_RAMP = [
  "#440154",
  "#414487",
  "#2a788e",
  "#22a884",
  "#7ad151",
  "#fde725",
] as const;

// ---------------------------------------------------------------------------
// Colorblind-safe ancestry palette (Okabe-Ito)
// ---------------------------------------------------------------------------

export type Ancestry = "EUR" | "AFR" | "EAS" | "AMR" | "SAS";

export const ANCESTRY_COLORS: Record<Ancestry, string> = {
  EUR: "#0072B2",
  AFR: "#E69F00",
  EAS: "#009E73",
  AMR: "#CC79A7",
  SAS: "#D55E00",
};

export const ANCESTRY_ORDER: Ancestry[] = ["EUR", "AFR", "EAS", "AMR", "SAS"];

/** Ancestry color with a gray fallback for unknown labels. */
export function ancestryColor(a: string): string {
  return ANCESTRY_COLORS[(a as Ancestry)] ?? CONTROL;
}

// ---------------------------------------------------------------------------
// Modality palette (S1–S8) — exact sRGB of the OKLCH tokens in globals.css
// ---------------------------------------------------------------------------

export type ModalityKey =
  | "s1_human"
  | "s2_genetic"
  | "s3_essential"
  | "s4_epigenomic"
  | "s5_spatial"
  | "s6_singlecell"
  | "s7_mouse"
  | "s8_proteomics";

export interface Modality {
  key: ModalityKey;
  label: string;
  /** CSS custom-property reference — use for SVG/DOM marks. */
  cssVar: string;
  /** Exact sRGB hex — use for canvas / WebGL marks. */
  hex: string;
}

export const MODALITIES: Modality[] = [
  { key: "s1_human", label: "Human RNA-seq", cssVar: "var(--color-s1-human)", hex: "#31aa40" },
  { key: "s2_genetic", label: "Genetic", cssVar: "var(--color-s2-genetic)", hex: "#5888fc" },
  { key: "s3_essential", label: "Essentiality", cssVar: "var(--color-s3-essential)", hex: "#e85a48" },
  { key: "s4_epigenomic", label: "Epigenomic", cssVar: "var(--color-s4-epigenomic)", hex: "#b16ae0" },
  { key: "s5_spatial", label: "Spatial", cssVar: "var(--color-s5-spatial)", hex: "#df6900" },
  { key: "s6_singlecell", label: "Single-cell", cssVar: "var(--color-s6-singlecell)", hex: "#00adba" },
  { key: "s7_mouse", label: "Mouse", cssVar: "var(--color-s7-mouse)", hex: "#db589e" },
  { key: "s8_proteomics", label: "Proteomics", cssVar: "var(--color-s8-proteomics)", hex: "#b59f00" },
];

export const MODALITY_HEX: Record<ModalityKey, string> = Object.fromEntries(
  MODALITIES.map((m) => [m.key, m.hex])
) as Record<ModalityKey, string>;

// ---------------------------------------------------------------------------
// Drug development-stage palette  (Approved -> Preclinical; Preclinical = gray)
// ---------------------------------------------------------------------------

export const DEV_STAGE_ORDER = [
  "Approved",
  "Phase III",
  "Phase II",
  "Phase I",
  "Preclinical",
] as const;

export const DEV_STAGE_COLORS: Record<string, string> = {
  Approved: "#0d4a7a",
  "Phase III": "#1f6aa5",
  "Phase II": "#4a8fc2",
  "Phase I": "#7fb0d4",
  // Coarse bucket used by the drug-pipeline data (approved / clinical / preclinical).
  Clinical: "#2e86ab",
  Preclinical: CONTROL,
};

/** Dev-stage color, case-insensitive, gray fallback. */
export function devStageColor(stage: string): string {
  const hit = Object.keys(DEV_STAGE_COLORS).find(
    (k) => k.toLowerCase() === stage.toLowerCase()
  );
  return hit ? DEV_STAGE_COLORS[hit] : CONTROL;
}

// ---------------------------------------------------------------------------
// Disease-stage palette  (Healthy = control gray; progression light -> dark)
// ---------------------------------------------------------------------------

export const DISEASE_STAGE_ORDER = [
  "Healthy",
  "Steatosis",
  "Steatohepatitis",
  "Cirrhosis",
] as const;

export const DISEASE_STAGE_COLORS: Record<string, string> = {
  Healthy: CONTROL,
  Steatosis: "#F6C445",
  Steatohepatitis: "#E8743B",
  Cirrhosis: "#B2182B",
};

export function diseaseStageColor(stage: string): string {
  const hit = Object.keys(DISEASE_STAGE_COLORS).find(
    (k) => k.toLowerCase() === stage.toLowerCase()
  );
  return hit ? DISEASE_STAGE_COLORS[hit] : CONTROL;
}

// ---------------------------------------------------------------------------
// Qualitative categorical palette (colorblind-aware; NO near-gray so control
// gray stays reserved). For arbitrary unordered categories — cell types,
// communities, two-class splits, etc. Cycles for large category counts.
// ---------------------------------------------------------------------------

export const CATEGORICAL = [
  "#0072B2", // blue
  "#E69F00", // orange
  "#009E73", // bluish green
  "#CC79A7", // reddish purple
  "#D55E00", // vermillion
  "#56B4E9", // sky blue
  "#F0E442", // yellow
  "#8C564B", // brown
  "#7570B3", // muted purple
  "#66A61E", // olive green
  "#A6761D", // dark goldenrod
  "#1B9E77", // teal
] as const;

/** i-th categorical color, cycling (negative-safe). */
export function categoricalColor(i: number): string {
  const n = CATEGORICAL.length;
  return CATEGORICAL[((Math.trunc(i) % n) + n) % n];
}

// ---------------------------------------------------------------------------
// Single-cell cell-type palette — the canonical Fig 3 map (publication_theme.R
// `ct_palette`). Color cells by type wherever they appear (main UMAP + landing
// mini-UMAP) so the web matches the published figures instead of an arbitrary
// index-based cycle. Unknown types fall back to the categorical cycle.
// ---------------------------------------------------------------------------

export const CELL_TYPE_COLORS: Record<string, string> = {
  Hepatocytes: "#0D47A1",
  Cholangiocytes: "#1565C0",
  "Endothelial cells": "#2E7D32",
  Fibroblasts: "#F57F17",
  Macrophages: "#C2185B",
  "Mono+mono derived cells": "#E91E63",
  "T cells": "#7B1FA2",
  "B cells": "#9C27B0",
  "Resident NK": "#00695C",
  "Circulating NK/NKT": "#00897B",
  "Plasma cells": "#5D4037",
  Neutrophils: "#FF6F00",
  cDC1s: "#AD1457",
  cDC2s: "#D81B60",
  pDCs: "#6A1B9A",
  "Mig.cDCs": "#AB47BC",
  Basophils: "#78909C",
};

/** Canonical cell-type color; categorical-cycle fallback for unmapped types. */
export function cellTypeColor(cellType: string, fallbackIndex = 0): string {
  return CELL_TYPE_COLORS[cellType] ?? categoricalColor(fallbackIndex);
}

// ---------------------------------------------------------------------------
// Network edge-layer palette (categorical layer identities). Mirrors the
// map in network-data.ts; import from here so the canvas + legend agree.
// ---------------------------------------------------------------------------

export const LAYER_COLORS: Record<string, string> = {
  ppi: "#e74c3c",
  coexpr: "#3498db",
  regulon: "#2ecc71",
  lr: "#f39c12",
  genetic: "#9b59b6",
  pathway: "#1abc9c",
  spatial: "#e67e22",
  cosmos: "#34495e",
  cerna: "#e91e63",
  xspecies: "#607d8b",
};

/** Edge-layer color, gray fallback for unknown layers. */
export function layerColor(layer: string): string {
  return LAYER_COLORS[layer] ?? CONTROL;
}

// ---------------------------------------------------------------------------
// Cross-species category palette (Not_Significant / Unclassified = control).
// ---------------------------------------------------------------------------

export const SPECIES_CATEGORY_COLORS: Record<string, string> = {
  Conserved_Core: CATEGORICAL[2], // green
  Human_Enriched: CATEGORICAL[0], // blue
  Mouse_Specific: CATEGORICAL[3], // reddish purple
  Diet_Selective: CATEGORICAL[1], // orange
  Species_Discordant: CATEGORICAL[4], // vermillion
  Not_Significant: CONTROL,
  Unclassified: CONTROL,
};

/** Cross-species category color, gray fallback. */
export function speciesCategoryColor(category: string): string {
  const hit = Object.keys(SPECIES_CATEGORY_COLORS).find(
    (k) => k.toLowerCase() === category.toLowerCase()
  );
  return hit ? SPECIES_CATEGORY_COLORS[hit] : CONTROL;
}

// ---------------------------------------------------------------------------
// Kleiner fibrosis-stage palette  (F0 = control gray; F1–F4 warm progression,
// mirroring DISEASE_STAGE_COLORS). F0 is "no fibrosis" → control gray.
// ---------------------------------------------------------------------------

export const FIBROSIS_STAGE_ORDER = ["F0", "F1", "F2", "F3", "F4"] as const;

export const FIBROSIS_STAGE_COLORS: Record<string, string> = {
  F0: CONTROL,
  F1: "#F6C445",
  F2: "#E8743B",
  F3: "#D6412B",
  F4: "#B2182B",
};

/** Fibrosis-stage color, gray fallback. */
export function fibrosisStageColor(stage: string): string {
  const key = stage.trim().toUpperCase();
  return FIBROSIS_STAGE_COLORS[key] ?? CONTROL;
}

// ---------------------------------------------------------------------------
// Interpolation helpers
// ---------------------------------------------------------------------------

const clamp01 = (t: number): number => (t < 0 ? 0 : t > 1 ? 1 : t);

function hexToRgb(hex: string): [number, number, number] {
  const h = hex.replace("#", "");
  return [
    parseInt(h.slice(0, 2), 16),
    parseInt(h.slice(2, 4), 16),
    parseInt(h.slice(4, 6), 16),
  ];
}

function rgbToHex(r: number, g: number, b: number): string {
  const c = (x: number) =>
    Math.round(Math.max(0, Math.min(255, x)))
      .toString(16)
      .padStart(2, "0");
  return `#${c(r)}${c(g)}${c(b)}`;
}

/** Linear interpolate between two hex colors. */
export function lerpHex(a: string, b: string, t: number): string {
  const [r1, g1, b1] = hexToRgb(a);
  const [r2, g2, b2] = hexToRgb(b);
  const u = clamp01(t);
  return rgbToHex(r1 + (r2 - r1) * u, g1 + (g2 - g1) * u, b1 + (b2 - b1) * u);
}

/** Sample a multi-stop ramp at t in [0, 1]. */
export function sampleRamp(ramp: readonly string[], t: number): string {
  const u = clamp01(t) * (ramp.length - 1);
  const i = Math.floor(u);
  if (i >= ramp.length - 1) return ramp[ramp.length - 1];
  return lerpHex(ramp[i], ramp[i + 1], u - i);
}

/** Sequential magnitude color for t in [0, 1]. */
export const sequentialColor = (t: number): string => sampleRamp(SEQUENTIAL_RAMP, t);

/** Perceptual heatmap color for t in [0, 1]. */
export const heatmapColor = (t: number): string => sampleRamp(HEATMAP_RAMP, t);

/**
 * Diverging effect color for a signed value (e.g. logFC).
 * `clamp` sets the magnitude that saturates to full blue/red (default 2).
 */
export function divergingColor(value: number, clamp = 2): string {
  const t = clamp01(Math.abs(value) / clamp);
  return value < 0
    ? lerpHex(DIVERGING.mid, DIVERGING.down, t)
    : lerpHex(DIVERGING.mid, DIVERGING.up, t);
}

/**
 * The primary diverging-mark encoder: direction by hue, significance by
 * opacity. Never rely on color alone — always spread both `fill` and
 * `opacity` onto the mark.
 */
export function getColor(
  logfc: number,
  sig: boolean,
  clamp = 2
): { fill: string; opacity: number } {
  return {
    fill: sig ? divergingColor(logfc, clamp) : CONTROL,
    opacity: sig ? SIG_OPACITY : NONSIG_OPACITY,
  };
}
