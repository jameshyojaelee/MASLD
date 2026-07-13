/** 8-axis evidence strength for a gene (0-1 normalized) */
export interface EvidenceStrengths {
  s1_human: number;      // Human Bulk RNA-seq
  s2_genetic: number;    // Genetic Causal (COLOC/TWAS)
  s3_essential: number;  // Essentiality (DepMap)
  s4_epigenomic: number; // Epigenomic (SCENIC+/ATAC)
  s5_spatial: number;    // Spatial (SVG/Moran's I)
  s6_singlecell: number; // Single-Cell (pseudobulk)
  s7_mouse: number;      // Mouse concordance
  s8_proteomics: number; // Proteomics (Olink plasma + DIA-MS liver)
}

/** Compact gene entry for search index and explorer table */
export interface GeneIndexEntry {
  symbol: string;
  ensembl_id: string;
  biotype: string | null;
  bulk_logfc: number | null;
  bulk_padj: number | null;
  is_deg: boolean;
  is_conserved_core: boolean;
  sex_class: string | null;
  zonation_class: string | null;
  ferroptosis_class: string | null;
  dgidb_druggable: boolean;
  layers_active: number;
  evidence: EvidenceStrengths;
}

/** Summary statistics for landing page */
export interface AtlasSummary {
  total_genes: number;
  total_degs: number;
  total_cohorts: number;
  total_samples: number;
  conserved_core_count: number;
  coloc_genes: number;
  drug_targets: number;
  mouse_datasets: number;
  evidence_sources: number;
}

/** A drug annotation attached to a featured gene (from `featured_genes.json`). */
export interface FeaturedGeneDrug {
  drug: string;
  stage: string;
  moa: string;
  /** Strength of atlas evidence backing the target, e.g. "Weak" | "Absent". */
  atlas_support?: string;
}

/**
 * Featured gene for landing-page cards. Shape mirrors `featured_genes.json`
 * exactly (there is no `coloc_pp4`/`drug` field — those were dead reads; drug
 * annotations live in the `drugs` array). NOTE: the emitted `evidence`
 * sub-object currently ships only s1–s7 (no `s8_proteomics`); the full
 * `EvidenceStrengths` type is retained because consumers render it through
 * `EvidenceFingerprint`, which tolerates a missing axis (`evidence[key] ?? 0`).
 */
export interface FeaturedGene {
  symbol: string;
  ensembl_id: string;
  tagline: string;
  category: string;
  bulk_logfc: number;
  bulk_padj: number;
  is_deg: boolean;
  is_conserved: boolean;
  layers_active: number;
  evidence: EvidenceStrengths;
  drugs: FeaturedGeneDrug[];
}

/** Modality metadata */
export interface EvidenceSource {
  key: keyof EvidenceStrengths;
  label: string;
  short: string;
  color: string;
  description: string;
}

/** Global store state */
export interface AppState {
  theme: "light" | "dark";
  setTheme: (theme: "light" | "dark") => void;
  sexLens: "combined" | "female" | "male";
  setSexLens: (lens: "combined" | "female" | "male") => void;
  compareGenes: string[];
  addCompareGene: (symbol: string) => void;
  removeCompareGene: (symbol: string) => void;
  clearCompareGenes: () => void;
  commandOpen: boolean;
  setCommandOpen: (open: boolean) => void;
}
