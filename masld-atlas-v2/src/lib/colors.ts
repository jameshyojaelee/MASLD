import type { EvidenceSource, EvidenceStrengths } from "./types";
import { GWAS_COUNT } from "./atlas-constants";

export const EVIDENCE_SOURCES: EvidenceSource[] = [
  {
    key: "s1_human",
    label: "Human Bulk RNA-seq",
    short: "Human",
    color: "var(--color-s1-human)",
    description: "Pooled (cohort-adjusted) analysis across 5 control-bearing cohorts (846 samples)",
  },
  {
    key: "s2_genetic",
    label: "Genetic Causal",
    short: "Genetic",
    color: "var(--color-s2-genetic)",
    description: `COLOC + TWAS across ${GWAS_COUNT} GWAS (5 ancestries)`,
  },
  {
    key: "s3_essential",
    label: "Essentiality",
    short: "Essential",
    color: "var(--color-s3-essential)",
    description: "DepMap CHRONOS essentiality + druggability",
  },
  {
    key: "s4_epigenomic",
    label: "Epigenomic",
    short: "Epigenomic",
    color: "var(--color-s4-epigenomic)",
    description: "SCENIC+ regulons + scATAC DA peaks",
  },
  {
    key: "s5_spatial",
    label: "Spatial",
    short: "Spatial",
    color: "var(--color-s5-spatial)",
    description: "Spatially variable genes + zonation",
  },
  {
    key: "s6_singlecell",
    label: "Single-Cell",
    short: "scRNA",
    color: "var(--color-s6-singlecell)",
    description: "Pseudobulk DE across 11 cell types + LIANA",
  },
  {
    key: "s7_mouse",
    label: "Mouse Concordance",
    short: "Mouse",
    color: "var(--color-s7-mouse)",
    description: "Cross-species validation across 5 diet models",
  },
  {
    key: "s8_proteomics",
    label: "Proteomics",
    short: "Proteins",
    color: "var(--color-s8-proteomics)",
    description: "Olink plasma + DIA-MS liver proteomics, mRNA-protein concordance",
  },
];

/** Get the CSS variable color string for a source key */
export function getSourceColor(key: keyof EvidenceStrengths): string {
  const source = EVIDENCE_SOURCES.find((s) => s.key === key);
  return source?.color ?? "var(--color-muted-foreground)";
}
