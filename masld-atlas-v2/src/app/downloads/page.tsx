import Link from "next/link";
import { PageContainer } from "@/components/page-container";
import { PageHeader } from "@/components/page-header";
import { dataUrl } from "@/lib/data-base";
import { DEG_COUNT, DEG_GATE_LABEL, fmt } from "@/lib/atlas-constants";

// ---------------------------------------------------------------------------
// Download resources
// ---------------------------------------------------------------------------

interface DownloadResource {
  name: string;
  description: string;
  format: string;
  size: string;
  href: string | null;
  note: string | null;
}

const DOWNLOAD_RESOURCES: DownloadResource[] = [
  {
    name: "Multi-Evidence Atlas",
    description:
      "27,187 genes x 446 columns with 7 independent modalities, causal inference, and derived annotations.",
    format: "Parquet",
    size: "21 MB",
    href: dataUrl("atlas.parquet"),
    note: null,
  },
  {
    name: "Gene Index",
    description:
      "Compact gene search index with symbol and biotype for all 27,187 genes, used by the in-app search.",
    format: "JSON",
    size: "1.3 MB",
    href: dataUrl("gene_symbols.json"),
    note: null,
  },
  {
    name: "Per-Gene Profiles",
    description:
      "Long-format, symbol-keyed tables (per-cohort DE, stage trajectories, COLOC by GWAS, per-cell-type DE, drug/LINCS evidence) queried client-side via DuckDB-WASM rather than shipped as one file per gene.",
    format: "Parquet (6 tables)",
    size: "~8 MB total",
    href: null,
    note: "Query any gene from the Gene page (/gene?symbol=SYMBOL), or download the individual tables below.",
  },
  {
    name: "Pathway Gene Sets",
    description:
      "50 Hallmark gene sets from MSigDB for enrichment analysis, formatted for direct use with the atlas.",
    format: "JSON",
    size: "61 KB",
    href: dataUrl("pathway_genesets.json"),
    note: null,
  },
  {
    name: "GWAS-ATAC Variants",
    description:
      "Regulatory variant annotations: 533 credible set variants overlapping scATAC peaks with motif disruption scores.",
    format: "JSON",
    size: "121 KB",
    href: dataUrl("gwas_atac_browser.json"),
    note: null,
  },
  {
    name: "Knowledge Graph",
    description:
      "Gene-drug-TF-pathway network capturing multi-evidence relationships for network visualization and analysis.",
    format: "JSON",
    size: "61 KB",
    href: dataUrl("knowledge_graph.json"),
    note: null,
  },
];

// ---------------------------------------------------------------------------
// Modality summaries
// ---------------------------------------------------------------------------

interface EvidenceSource {
  name: string;
  abbrev: string;
  description: string;
}

const EVIDENCE_SOURCES: EvidenceSource[] = [
  {
    name: "Human Bulk RNA-seq",
    abbrev: "S1",
    description: `Pooled (cohort-adjusted) analysis across 5 control-bearing cohorts (846 samples). Per-study DE via limma-voom, followed by cross-cohort integration. ${fmt(
      DEG_COUNT
    )} DEGs by the ${DEG_GATE_LABEL}, validated by leave-one-out cross-validation.`,
  },
  {
    name: "Mouse Bulk RNA-seq",
    abbrev: "S2",
    description:
      "Integration of 5 mouse diet models (463 samples) with metafor random-effects meta-analysis and pooled (cohort-adjusted) analysis. Cross-species concordance via strict 1:1 ortholog mapping (GENCODE v49 / vM38). 723 Conserved_Core genes with concordant direction.",
  },
  {
    name: "Genetic Causal Inference",
    abbrev: "S3",
    description:
      "TWAS (elastic net, OTTERS multi-method), COLOC (SuSiE + ABF across 50 GWAS, 5 ancestries), HyPrColoc multi-trait, sc-TWAS (multi-cell-type). FinnGen + BBJ + Pan-UKBB cross-ancestry replication. 434 genes with SuSiE PP.H4 > 0.9.",
  },
  {
    name: "Essentiality (DepMap)",
    abbrev: "S4",
    description:
      "CRISPR-Cas9 essentiality screening from DepMap (Chronos scores). Identifies genes whose knockout affects hepatocyte viability, providing functional validation orthogonal to expression data.",
  },
  {
    name: "Epigenomic (ATAC / SCENIC+)",
    abbrev: "S5",
    description:
      "scATAC-seq peak accessibility, SCENIC+ gene regulatory networks (139 hepatocyte regulons), chromVAR motif activity, and cross-species promoter conservation. GWAS-ATAC integration identifies 380 regulatory variants disrupting TF motifs.",
  },
  {
    name: "Spatial Transcriptomics",
    abbrev: "S6",
    description:
      "Visium spatial transcriptomics (GSE192741) with spatially variable gene detection, hepatocyte-dominant spatial domains, and zonation classification (pericentral/periportal). 150 SVGs with BH FDR correction.",
  },
  {
    name: "Single-Cell (Pseudobulk + LIANA)",
    abbrev: "S7",
    description:
      "scRNA-seq pseudobulk DE across 5 cell types, ligand-receptor interaction analysis via LIANA, and cell-type-specific TWAS. Hepatocyte, macrophage, fibroblast, endothelial, and cholangiocyte resolution.",
  },
];

// ---------------------------------------------------------------------------
// Page
// ---------------------------------------------------------------------------

export default function DownloadsPage() {
  return (
    <PageContainer>
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <PageHeader
        title="Downloads & Documentation"
        description="Full multi-evidence atlas, gene profiles, and supplementary data. All resources are served as static files for programmatic access."
      />

      {/* ------------------------------------------------------------------ */}
      {/* Section 1: Data Downloads                                            */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Data Downloads
        </h2>
        <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-3">
          {DOWNLOAD_RESOURCES.map((res) => (
            <div
              key={res.name}
              className="flex flex-col gap-3 rounded-lg border border-border p-4 transition-colors hover:bg-muted/30"
            >
              <div>
                <h3 className="text-sm font-semibold">{res.name}</h3>
                <p className="mt-1 text-xs leading-relaxed text-muted-foreground">
                  {res.description}
                </p>
              </div>
              <div className="mt-auto flex items-center justify-between gap-2 pt-1">
                <div className="flex items-center gap-2">
                  <span className="inline-flex items-center rounded-full border border-border px-2 py-0.5 text-[10px] font-medium text-muted-foreground">
                    {res.format}
                  </span>
                  <span className="text-[10px] text-muted-foreground">
                    {res.size}
                  </span>
                </div>
                {res.href ? (
                  <a
                    href={res.href}
                    download
                    className="inline-flex items-center rounded-md bg-primary px-3 py-1.5 text-xs font-medium text-primary-foreground transition-colors hover:bg-primary/80"
                  >
                    Download
                  </a>
                ) : (
                  <span className="text-[10px] italic text-muted-foreground">
                    See note
                  </span>
                )}
              </div>
              {res.note && (
                <p className="text-[10px] leading-relaxed text-muted-foreground">
                  {res.note}
                </p>
              )}
            </div>
          ))}
        </div>
        <p className="mt-4 text-xs text-muted-foreground">
          The full multi-evidence atlas CSV (~15 MB) is available from the{" "}
          <a
            href="https://github.com/sanjana-lab/masld-atlas"
            target="_blank"
            rel="noopener noreferrer"
            className="font-medium text-primary hover:underline"
          >
            project GitHub repository
          </a>
          .
        </p>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 2: Methods & Documentation                                   */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          Methods &amp; Documentation
        </h2>

        <div className="mb-6 flex flex-wrap gap-3">
          <div className="rounded-lg border border-border px-4 py-3">
            <p className="text-sm font-semibold">Reference</p>
            <p className="mt-1 text-xs text-muted-foreground">
              Lee J et al. A multi-evidence transcriptomic atlas reveals a
              multi-step metabolic-to-inflammatory progression in MASLD.{" "}
              <span className="italic">In preparation</span> (2026).
            </p>
          </div>
          <a
            href="https://github.com/sanjana-lab/masld-atlas"
            target="_blank"
            rel="noopener noreferrer"
            className="inline-flex items-center gap-2 rounded-lg border border-border px-4 py-3 transition-colors hover:bg-muted/30"
          >
            <span className="text-sm font-semibold">GitHub Repository</span>
            <span className="text-xs text-muted-foreground">&rarr;</span>
          </a>
        </div>

        <h3 className="mb-3 text-base font-semibold">Modalities</h3>
        <div className="space-y-3">
          {EVIDENCE_SOURCES.map((src) => (
            <div
              key={src.abbrev}
              className="rounded-lg border border-border/50 px-4 py-3"
            >
              <p className="text-sm">
                <span className="font-mono text-xs font-semibold text-primary">
                  {src.abbrev}
                </span>{" "}
                <span className="font-semibold">{src.name}</span>
              </p>
              <p className="mt-1 text-xs leading-relaxed text-muted-foreground">
                {src.description}
              </p>
            </div>
          ))}
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 3: How to Cite                                               */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          How to Cite
        </h2>
        <div className="overflow-x-auto rounded-lg border border-border bg-muted/30 p-4">
          <pre className="text-xs leading-relaxed text-foreground">
{`@article{lee2026masld,
  title={A multi-evidence transcriptomic atlas reveals a multi-step
         metabolic-to-inflammatory progression in MASLD},
  author={Lee, James and others},
  journal={In preparation},
  year={2026}
}`}
          </pre>
        </div>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Section 4: API / Programmatic Access                                 */}
      {/* ------------------------------------------------------------------ */}
      <section className="mb-10">
        <h2 className="mb-4 text-xl font-semibold tracking-tight">
          API / Programmatic Access
        </h2>
        <p className="mb-4 text-sm text-muted-foreground">
          All data is served as static Parquet and JSON files from the{" "}
          <a
            href="https://huggingface.co/datasets/jameshyojaelee/masld-atlas-data"
            target="_blank"
            rel="noopener noreferrer"
            className="font-medium text-primary hover:underline"
          >
            companion Hugging Face Dataset
          </a>
          . The app itself queries per-gene evidence client-side with
          DuckDB-WASM &mdash; no per-gene endpoint is served; query the parquet
          files directly instead:
        </p>

        <div className="overflow-x-auto rounded-lg border border-border bg-muted/30 p-4">
          <p className="mb-2 text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
            Python Example (pandas)
          </p>
          <pre className="text-xs leading-relaxed text-foreground">
{`import pandas as pd

base = "https://huggingface.co/datasets/jameshyojaelee/masld-atlas-data/resolve/main"
atlas = pd.read_parquet(f"{base}/atlas_core.parquet")
row = atlas.loc[atlas["human_symbol"] == "THRB"]
print(row["coloc_best_susie_pp4"].iloc[0])  # 0.9999`}
          </pre>
        </div>

        <div className="mt-3 overflow-x-auto rounded-lg border border-border bg-muted/30 p-4">
          <p className="mb-2 text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
            JavaScript Example (DuckDB-WASM, as used in-app)
          </p>
          <pre className="text-xs leading-relaxed text-foreground">
{`import { queryParquet } from "@/lib/duck";

const rows = await queryParquet(
  "atlas_core.parquet",
  (t) => \`SELECT * FROM \${t} WHERE human_symbol = ?\`,
  ["HNF4A"]
);
console.log(rows[0].bulk_logFC);  // -0.482`}
          </pre>
        </div>

        <p className="mt-3 text-xs text-muted-foreground">
          The slim search index at{" "}
          <code className="rounded bg-muted px-1 py-0.5 text-[10px]">
            gene_symbols.json
          </code>{" "}
          contains symbol and biotype for all 27,187 genes, suitable for
          building custom search or enrichment tools.
        </p>
      </section>

      {/* ------------------------------------------------------------------ */}
      {/* Footer nav                                                           */}
      {/* ------------------------------------------------------------------ */}
      <div className="flex flex-wrap gap-3 border-t border-border pt-6">
        <Link
          href="/explore"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          &larr; Explore Genes
        </Link>
        <Link
          href="/"
          className="rounded-md border border-border px-4 py-2 text-sm font-medium transition-colors hover:bg-muted"
        >
          Home
        </Link>
      </div>
    </PageContainer>
  );
}
