import Link from "next/link";

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
      "33,943 genes x 101 columns with 7 independent modalities, causal inference, and derived annotations.",
    format: "Parquet",
    size: "4.3 MB",
    href: "/data/atlas.parquet",
    note: null,
  },
  {
    name: "Gene Index",
    description:
      "Compact gene search index with DEG status, evidence strength indicators, and effect sizes for all 33,943 genes.",
    format: "JSON",
    size: "6 MB",
    href: "/data/gene_index.json",
    note: null,
  },
  {
    name: "Per-Gene Profiles",
    description:
      "Individual evidence cards for all genes. Each file contains full multi-evidence profiles including expression, causal, drug, and epigenomic data.",
    format: "JSON (33,943 files)",
    size: "19 MB total",
    href: null,
    note: "Access individual genes at /data/genes/{SYMBOL}.json",
  },
  {
    name: "Pathway Gene Sets",
    description:
      "50 Hallmark gene sets from MSigDB for enrichment analysis, formatted for direct use with the atlas.",
    format: "JSON",
    size: "61 KB",
    href: "/data/pathway_genesets.json",
    note: null,
  },
  {
    name: "GWAS-ATAC Variants",
    description:
      "Regulatory variant annotations: 533 credible set variants overlapping scATAC peaks with motif disruption scores.",
    format: "JSON",
    size: "121 KB",
    href: "/data/gwas_atac_browser.json",
    note: null,
  },
  {
    name: "Knowledge Graph",
    description:
      "Gene-drug-TF-pathway network capturing multi-evidence relationships for network visualization and analysis.",
    format: "JSON",
    size: "61 KB",
    href: "/data/knowledge_graph.json",
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
    description:
      "Integrated mixed-model mega-analysis across 10 cohorts (1,444 samples). Per-study DE via limma-voom, followed by variance partition and cross-cohort integration. 5,484 DEGs at padj < 0.05, |logFC| > 0.3, validated by leave-one-out cross-validation (88.1% mean recovery).",
  },
  {
    name: "Mouse Bulk RNA-seq",
    abbrev: "S2",
    description:
      "Integration of 5 mouse diet models (463 samples) with metafor random-effects meta-analysis and Integrated mega-analysis. Cross-species concordance via strict 1:1 ortholog mapping (GENCODE v49 / vM38). 723 Conserved_Core genes with concordant direction.",
  },
  {
    name: "Genetic Causal Inference",
    abbrev: "S3",
    description:
      "TWAS (elastic net, OTTERS multi-method), COLOC (SuSiE + ABF across 24 European GWAS), HyPrColoc multi-trait, cTWAS, sc-TWAS (multi-cell-type), and bidirectional MR. FinnGen + BBJ cross-ancestry replication. 179 genes with PP.H4 > 0.9.",
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
    <div className="mx-auto max-w-5xl px-6 py-10">
      {/* ------------------------------------------------------------------ */}
      {/* Header                                                               */}
      {/* ------------------------------------------------------------------ */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold tracking-tight">
          Downloads &amp; Documentation
        </h1>
        <p className="mt-2 text-muted-foreground">
          Full multi-evidence atlas, gene profiles, and supplementary data.
          All resources are served as static files for programmatic access.
        </p>
      </div>

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
              Lee J et al. A multi-evidence transcriptomic atlas reveals the F2
              metabolic-to-inflammatory switch in MASLD.{" "}
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
  title={A multi-evidence transcriptomic atlas reveals the F2
         metabolic-to-inflammatory switch in MASLD},
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
          All data is served as static JSON and Parquet files. You can fetch any
          gene&apos;s evidence profile directly:
        </p>

        <div className="mb-3 rounded-lg border border-border bg-muted/30 p-4">
          <p className="mb-2 text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
            Endpoint
          </p>
          <code className="text-sm font-semibold text-primary">
            /data/genes/&#123;SYMBOL&#125;.json
          </code>
        </div>

        <div className="overflow-x-auto rounded-lg border border-border bg-muted/30 p-4">
          <p className="mb-2 text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
            Python Example
          </p>
          <pre className="text-xs leading-relaxed text-foreground">
{`import requests

gene = requests.get("https://masld-atlas.org/data/genes/THRB.json").json()
print(gene["causal"]["coloc_pp4_max"])  # 0.9999`}
          </pre>
        </div>

        <div className="mt-3 overflow-x-auto rounded-lg border border-border bg-muted/30 p-4">
          <p className="mb-2 text-[10px] font-semibold uppercase tracking-wider text-muted-foreground">
            JavaScript Example
          </p>
          <pre className="text-xs leading-relaxed text-foreground">
{`const res = await fetch("/data/genes/HNF4A.json");
const gene = await res.json();
console.log(gene.expression.bulk_logfc);  // -0.482`}
          </pre>
        </div>

        <p className="mt-3 text-xs text-muted-foreground">
          The full gene index at{" "}
          <code className="rounded bg-muted px-1 py-0.5 text-[10px]">
            /data/gene_index.json
          </code>{" "}
          contains compact records for all 33,943 genes, suitable for building
          custom search or enrichment tools.
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
    </div>
  );
}
