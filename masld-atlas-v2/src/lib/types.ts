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
  dream_logfc: number | null;
  dream_padj: number | null;
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

/** Featured gene for landing page cards */
export interface FeaturedGene {
  symbol: string;
  ensembl_id: string;
  tagline: string;
  evidence: EvidenceStrengths;
  dream_logfc: number;
  coloc_pp4: number | null;
  drug: string | null;
  category: string;
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
