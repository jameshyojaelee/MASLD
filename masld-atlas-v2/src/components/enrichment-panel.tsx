"use client";

import { useEffect, useMemo, useState } from "react";
import { dataUrl } from "@/lib/data-base";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import {
  runEnrichment,
  type EnrichmentResult,
  type PathwayGeneSet,
} from "@/lib/enrichment";

interface EnrichmentPanelProps {
  genes: string[];
  universeSize: number;
  onClose: () => void;
}

function formatPval(val: number): string {
  if (val < 1e-300) return "<1e-300";
  if (val < 0.001) return val.toExponential(1);
  return val.toFixed(3);
}

const MAX_DISPLAY = 20;

export function EnrichmentPanel({
  genes,
  universeSize,
  onClose,
}: EnrichmentPanelProps) {
  const [geneSets, setGeneSets] = useState<PathwayGeneSet[]>([]);
  const [loading, setLoading] = useState(true);
  const [expandedRow, setExpandedRow] = useState<string | null>(null);
  const [showAll, setShowAll] = useState(false);

  // Load pathway gene sets once
  useEffect(() => {
    fetch(dataUrl("pathway_genesets.json"))
      .then((r) => r.json())
      .then(
        (data: {
          collections: { gene_sets: PathwayGeneSet[] }[];
        }) => {
          const sets: PathwayGeneSet[] = [];
          for (const collection of data.collections) {
            for (const gs of collection.gene_sets) {
              sets.push(gs);
            }
          }
          setGeneSets(sets);
          setLoading(false);
        }
      );
  }, []);

  const results = useMemo(() => {
    if (geneSets.length === 0 || genes.length === 0) return [];
    return runEnrichment(genes, geneSets, universeSize);
  }, [genes, geneSets, universeSize]);

  const displayResults = showAll ? results : results.slice(0, MAX_DISPLAY);

  return (
    <div className="mt-6 rounded-lg border border-border bg-muted/20 p-4">
      {/* Header */}
      <div className="mb-4 flex items-center justify-between">
        <div>
          <h2 className="text-lg font-semibold tracking-tight">
            Pathway Enrichment
          </h2>
          <p className="text-sm text-muted-foreground">
            Fisher exact test (Hallmark gene sets) &middot;{" "}
            <span className="font-medium text-foreground">
              {genes.length.toLocaleString()}
            </span>{" "}
            input genes &middot; BH-corrected
          </p>
        </div>
        <Button variant="ghost" size="sm" onClick={onClose}>
          Close
        </Button>
      </div>

      {/* Content */}
      {loading ? (
        <p className="py-8 text-center text-sm text-muted-foreground">
          Loading gene sets...
        </p>
      ) : results.length === 0 ? (
        <p className="py-8 text-center text-sm text-muted-foreground">
          No enriched pathways found (0 overlapping genes).
        </p>
      ) : (
        <>
          <div className="overflow-x-auto rounded-lg border border-border">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border bg-muted/50">
                  <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                    Pathway
                  </th>
                  <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                    Overlap
                  </th>
                  <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                    P-value
                  </th>
                  <th className="px-3 py-2 text-right text-xs font-medium text-muted-foreground">
                    Adj. P
                  </th>
                  <th className="px-3 py-2 text-left text-xs font-medium text-muted-foreground">
                    Genes
                  </th>
                </tr>
              </thead>
              <tbody>
                {displayResults.map((r) => {
                  const isExpanded = expandedRow === r.pathway_id;
                  const isSig = r.padj < 0.05;

                  return (
                    <tr
                      key={r.pathway_id}
                      className={`border-b border-border/50 transition-colors hover:bg-muted/30 ${
                        isSig ? "bg-primary/5" : ""
                      }`}
                    >
                      {/* Pathway name */}
                      <td className="px-3 py-1.5">
                        <span className="font-medium">{r.pathway_name}</span>
                        {isSig && (
                          <Badge
                            variant="default"
                            className="ml-2 text-[10px]"
                          >
                            FDR&lt;0.05
                          </Badge>
                        )}
                      </td>

                      {/* Overlap */}
                      <td className="px-3 py-1.5 text-right font-mono text-xs">
                        <span className="text-foreground">{r.overlap}</span>
                        <span className="text-muted-foreground">
                          /{r.pathway_size}
                        </span>
                      </td>

                      {/* P-value */}
                      <td className="px-3 py-1.5 text-right font-mono text-xs text-muted-foreground">
                        {formatPval(r.pvalue)}
                      </td>

                      {/* Adj P-value */}
                      <td
                        className={`px-3 py-1.5 text-right font-mono text-xs ${
                          isSig
                            ? "font-semibold text-foreground"
                            : "text-muted-foreground"
                        }`}
                      >
                        {formatPval(r.padj)}
                      </td>

                      {/* Overlap genes (click to expand) */}
                      <td className="max-w-xs px-3 py-1.5">
                        <button
                          className="text-left text-xs text-muted-foreground hover:text-foreground"
                          onClick={() =>
                            setExpandedRow(isExpanded ? null : r.pathway_id)
                          }
                        >
                          {isExpanded ? (
                            <span className="font-mono leading-relaxed">
                              {r.overlap_genes.join(", ")}
                            </span>
                          ) : (
                            <span>
                              {r.overlap_genes.slice(0, 3).join(", ")}
                              {r.overlap_genes.length > 3 && (
                                <span className="ml-1 text-primary">
                                  +{r.overlap_genes.length - 3} more
                                </span>
                              )}
                            </span>
                          )}
                        </button>
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          {/* Show more / fewer toggle */}
          {results.length > MAX_DISPLAY && (
            <div className="mt-3 text-center">
              <Button
                variant="ghost"
                size="sm"
                onClick={() => setShowAll((v) => !v)}
              >
                {showAll
                  ? "Show top 20 only"
                  : `Show all ${results.length} results`}
              </Button>
            </div>
          )}
        </>
      )}
    </div>
  );
}
