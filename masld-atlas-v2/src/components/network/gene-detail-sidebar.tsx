"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Separator } from "@/components/ui/separator";
import { EvidenceFingerprint } from "@/components/evidence-fingerprint";
import { getLayerColor } from "@/lib/network-data";
import type { CommunityData, GeneHalo, NetworkEdge, NetworkNode, EdgeLayer } from "@/lib/network-types";

interface GeneDetailSidebarProps {
  node: NetworkNode | null;
  neighbors: NetworkNode[];
  edges: NetworkEdge[];
  communities: CommunityData | null;
  onClose: () => void;
}

function formatLogFC(val: number | null): string {
  if (val == null) return "\u2014";
  return val >= 0 ? `+${val.toFixed(3)}` : val.toFixed(3);
}

function logfcColor(val: number | null): string {
  if (val == null) return "text-muted-foreground";
  if (val > 0) return "text-red-500 dark:text-red-400";
  if (val < 0) return "text-blue-500 dark:text-blue-400";
  return "text-muted-foreground";
}

/** Aggregate edges to get top neighbors sorted by max posterior */
function getTopNeighbors(
  center: string,
  edges: NetworkEdge[],
  neighbors: NetworkNode[],
  limit = 15
): Array<{ node: NetworkNode; maxPosterior: number; layers: EdgeLayer[] }> {
  const neighborMap = new Map<string, { maxPosterior: number; layers: Set<EdgeLayer> }>();

  for (const edge of edges) {
    const otherId = edge.source === center ? edge.target : edge.source;
    if (otherId === center) continue;

    const existing = neighborMap.get(otherId);
    if (existing) {
      existing.maxPosterior = Math.max(existing.maxPosterior, edge.posterior);
      existing.layers.add(edge.layer);
    } else {
      neighborMap.set(otherId, {
        maxPosterior: edge.posterior,
        layers: new Set([edge.layer]),
      });
    }
  }

  const nodeById = new Map(neighbors.map((n) => [n.id, n]));

  return Array.from(neighborMap.entries())
    .map(([id, data]) => ({
      node: nodeById.get(id)!,
      maxPosterior: data.maxPosterior,
      layers: Array.from(data.layers),
    }))
    .filter((entry) => entry.node != null)
    .sort((a, b) => b.maxPosterior - a.maxPosterior)
    .slice(0, limit);
}

function CommunityBreadcrumb({
  node,
  communities,
}: {
  node: NetworkNode;
  communities: CommunityData | null;
}) {
  const macro = communities?.macro.find((c) => c.id === node.community_macro);
  const meso = communities?.meso.find((c) => c.id === node.community_meso);

  return (
    <div className="flex flex-wrap items-center gap-1 text-xs text-muted-foreground">
      {macro && (
        <Badge
          variant="outline"
          className="text-[10px]"
          style={{ borderColor: macro.color, color: macro.color }}
        >
          {macro.label}
        </Badge>
      )}
      {macro && meso && <span className="text-muted-foreground/50">/</span>}
      {meso && (
        <Badge variant="outline" className="text-[10px]">
          {meso.label}
        </Badge>
      )}
      {node.community_label && (
        <>
          <span className="text-muted-foreground/50">/</span>
          <span className="text-xs">{node.community_label}</span>
        </>
      )}
    </div>
  );
}

export function GeneDetailSidebar({
  node,
  neighbors,
  edges,
  communities,
  onClose,
}: GeneDetailSidebarProps) {
  const [halo, setHalo] = useState<GeneHalo | null>(null);

  // Load halo data when node changes
  useEffect(() => {
    if (!node) {
      setHalo(null);
      return;
    }

    setHalo(null);
    fetch(`/data/network/halos/${node.symbol.toUpperCase()}.json`)
      .then((r) => (r.ok ? (r.json() as Promise<GeneHalo>) : null))
      .then(setHalo)
      .catch(() => setHalo(null));
  }, [node]);

  if (!node) return null;

  const topNeighbors = getTopNeighbors(node.id, edges, neighbors);

  return (
    <div className="flex h-full w-80 shrink-0 flex-col border-l border-border bg-card">
      {/* Header */}
      <div className="flex items-start justify-between p-4">
        <div className="min-w-0">
          <div className="flex items-center gap-2">
            <h2 className="truncate font-mono text-lg font-bold">{node.symbol}</h2>
            {node.biotype && (
              <Badge variant="secondary" className="shrink-0 text-[10px]">
                {node.biotype}
              </Badge>
            )}
          </div>
          <p className="mt-0.5 truncate text-xs text-muted-foreground">
            {node.ensembl_id}
          </p>
        </div>
        <button
          type="button"
          onClick={onClose}
          className="ml-2 shrink-0 rounded-md p-1 text-muted-foreground transition-colors hover:bg-muted hover:text-foreground"
          aria-label="Close sidebar"
        >
          <svg className="size-4" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
            <path strokeLinecap="round" strokeLinejoin="round" d="M6 18L18 6M6 6l12 12" />
          </svg>
        </button>
      </div>

      <Separator />

      <ScrollArea className="flex-1">
        <div className="space-y-4 p-4">
          {/* Evidence fingerprint */}
          <div>
            <p className="mb-2 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
              Evidence
            </p>
            <div className="flex items-center gap-3">
              <EvidenceFingerprint evidence={node.evidence} size={56} showTooltip />
              <div className="space-y-0.5 text-xs">
                <p>
                  <span className="text-muted-foreground">Layers:</span>{" "}
                  <span className="font-mono font-medium">{node.layers_active}/7</span>
                </p>
                <p>
                  <span className="text-muted-foreground">Degree:</span>{" "}
                  <span className="font-mono font-medium">{node.degree}</span>
                </p>
                <p>
                  <span className="text-muted-foreground">logFC:</span>{" "}
                  <span className={`font-mono font-medium ${logfcColor(node.bulk_logfc)}`}>
                    {formatLogFC(node.bulk_logfc)}
                  </span>
                </p>
              </div>
            </div>
          </div>

          {/* Tags */}
          <div className="flex flex-wrap gap-1">
            {node.is_deg && <Badge variant="default" className="text-[10px]">DEG</Badge>}
            {node.is_conserved_core && <Badge variant="secondary" className="text-[10px]">Conserved Core</Badge>}
            {node.dgidb_druggable && <Badge variant="outline" className="text-[10px]">Druggable</Badge>}
            {node.sex_class && <Badge variant="outline" className="text-[10px]">{node.sex_class}</Badge>}
            {node.progression_class && <Badge variant="outline" className="text-[10px]">{node.progression_class}</Badge>}
          </div>

          <Separator />

          {/* Community membership */}
          <div>
            <p className="mb-2 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
              Community
            </p>
            <CommunityBreadcrumb node={node} communities={communities} />
          </div>

          <Separator />

          {/* Top neighbors */}
          <div>
            <p className="mb-2 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
              Top Neighbors ({topNeighbors.length})
            </p>
            {topNeighbors.length === 0 ? (
              <p className="text-xs text-muted-foreground">No neighbors at current threshold</p>
            ) : (
              <ul className="space-y-1.5">
                {topNeighbors.map(({ node: neighbor, maxPosterior, layers }) => (
                  <li
                    key={neighbor.id}
                    className="flex items-center justify-between rounded-md px-2 py-1 text-xs transition-colors hover:bg-muted/50"
                  >
                    <div className="flex items-center gap-1.5 min-w-0">
                      <span className="truncate font-mono font-medium">{neighbor.symbol}</span>
                      <div className="flex gap-0.5">
                        {layers.map((l) => (
                          <span
                            key={l}
                            className="inline-block size-1.5 rounded-full"
                            style={{ backgroundColor: getLayerColor(l) }}
                            title={l}
                          />
                        ))}
                      </div>
                    </div>
                    <span className="shrink-0 font-mono text-muted-foreground">
                      {maxPosterior.toFixed(2)}
                    </span>
                  </li>
                ))}
              </ul>
            )}
          </div>

          <Separator />

          {/* Halo data */}
          <div>
            <p className="mb-2 text-xs font-semibold uppercase tracking-wider text-muted-foreground">
              Annotations
            </p>
            {halo == null ? (
              <p className="text-xs text-muted-foreground">Loading...</p>
            ) : (
              <div className="space-y-2">
                {halo.drugs.length > 0 && (
                  <div>
                    <p className="text-[10px] font-medium text-muted-foreground">Drugs</p>
                    <div className="mt-0.5 flex flex-wrap gap-1">
                      {halo.drugs.map((d) => (
                        <Badge key={d} variant="outline" className="text-[10px]">{d}</Badge>
                      ))}
                    </div>
                  </div>
                )}
                {halo.gwas_variants.length > 0 && (
                  <div>
                    <p className="text-[10px] font-medium text-muted-foreground">GWAS Variants</p>
                    <div className="mt-0.5 flex flex-wrap gap-1">
                      {halo.gwas_variants.map((v) => (
                        <Badge key={v} variant="outline" className="text-[10px] font-mono">{v}</Badge>
                      ))}
                    </div>
                  </div>
                )}
                {halo.pathways.length > 0 && (
                  <div>
                    <p className="text-[10px] font-medium text-muted-foreground">Pathways</p>
                    <div className="mt-0.5 flex flex-wrap gap-1">
                      {halo.pathways.map((p) => (
                        <Badge key={p} variant="secondary" className="text-[10px]">{p}</Badge>
                      ))}
                    </div>
                  </div>
                )}
                {halo.regulons.length > 0 && (
                  <div>
                    <p className="text-[10px] font-medium text-muted-foreground">Regulons</p>
                    <div className="mt-0.5 flex flex-wrap gap-1">
                      {halo.regulons.map((r) => (
                        <Badge key={r} variant="secondary" className="text-[10px]">{r}</Badge>
                      ))}
                    </div>
                  </div>
                )}
                {halo.drugs.length === 0 &&
                  halo.gwas_variants.length === 0 &&
                  halo.pathways.length === 0 &&
                  halo.regulons.length === 0 && (
                  <p className="text-xs text-muted-foreground">No annotations available</p>
                )}
              </div>
            )}
          </div>
        </div>
      </ScrollArea>

      {/* Footer — link to full gene page */}
      <Separator />
      <div className="p-4">
        <Button
          variant="outline"
          size="sm"
          className="w-full"
          render={<Link href={`/gene/${encodeURIComponent(node.symbol)}`} />}
        >
          View full gene page
          <svg className="ml-1 size-3.5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
            <path strokeLinecap="round" strokeLinejoin="round" d="M13.5 6H5.25A2.25 2.25 0 003 8.25v10.5A2.25 2.25 0 005.25 21h10.5A2.25 2.25 0 0018 18.75V10.5m-10.5 6L21 3m0 0h-5.25M21 3v5.25" />
          </svg>
        </Button>
      </div>
    </div>
  );
}
