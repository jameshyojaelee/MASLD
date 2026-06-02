"use client";

import { useCallback } from "react";
import { Badge } from "@/components/ui/badge";
import { DEFAULT_LAYERS } from "@/lib/network-data";
import type { EdgeLayer, LayerMetadata, NodeColorBy } from "@/lib/network-types";

/**
 * Per-edge-type confidence thresholds. Each threshold dispatches to a single
 * edge type via the force-graph edgeVisible() function:
 *   - `string`: applies to v4 type === "S" (STRING score)
 *   - `loco`:   applies to v4 type === "D-F2" (LOCO replication fraction)
 *   - `pp4`:    applies to v4 type === "D-COLOC" (PP.H4 min across GWAS)
 */
export interface Thresholds {
  string: number;
  loco: number;
  pp4: number;
}

export const DEFAULT_THRESHOLDS: Thresholds = {
  string: 0.7,
  loco: 0.7,
  pp4: 0.5,
};

interface ControlsPanelProps {
  layers: LayerMetadata[];
  selectedLayers: Set<EdgeLayer>;
  onToggleLayer: (layer: EdgeLayer) => void;
  thresholds: Thresholds;
  onThresholdsChange: (value: Thresholds) => void;
  colorBy: NodeColorBy;
  onColorByChange: (value: NodeColorBy) => void;
  labelTopK: number;
  onLabelTopKChange: (value: number) => void;
}

const COLOR_BY_OPTIONS: { value: NodeColorBy; label: string }[] = [
  { value: "community", label: "Community" },
  { value: "layers", label: "Evidence Layers" },
  { value: "sex_class", label: "Sex Class" },
  { value: "progression", label: "Progression" },
  { value: "logfc", label: "logFC" },
  { value: "druggability", label: "Druggability" },
];

export function ControlsPanel({
  layers,
  selectedLayers,
  onToggleLayer,
  thresholds,
  onThresholdsChange,
  colorBy,
  onColorByChange,
  labelTopK,
  onLabelTopKChange,
}: ControlsPanelProps) {
  // Hide legacy layers (regulon, pathway, spatial, cosmos) from the v4 UI.
  // Keys remain in the EdgeLayer type for backcompat — this just filters the
  // rendered toggles.
  const effectiveLayers = (layers.length > 0 ? layers : DEFAULT_LAYERS).filter(
    (l) => !l.legacy
  );

  const updateThreshold = useCallback(
    (key: keyof Thresholds, value: number) => {
      onThresholdsChange({ ...thresholds, [key]: value });
    },
    [onThresholdsChange, thresholds]
  );

  return (
    <div className="flex flex-wrap items-start gap-4 rounded-lg border border-border bg-card px-4 py-3">
      {/* Edge layer toggles */}
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-xs font-medium text-muted-foreground">Edges:</span>
        {effectiveLayers.map((layer) => {
          const active = selectedLayers.has(layer.id);
          return (
            <button
              key={layer.id}
              type="button"
              onClick={() => onToggleLayer(layer.id)}
              className={
                "inline-flex items-center gap-1.5 rounded-md border px-2 py-1 text-xs transition-colors " +
                (active
                  ? "border-border bg-muted text-foreground"
                  : "border-transparent bg-transparent text-muted-foreground opacity-50 hover:opacity-75")
              }
              title={layer.description}
            >
              <span
                className="inline-block size-2.5 shrink-0 rounded-full"
                style={{ backgroundColor: layer.color }}
              />
              {layer.label}
              {layer.edge_count > 0 && (
                <Badge variant="secondary" className="ml-0.5 text-[9px] px-1 py-0">
                  {layer.edge_count.toLocaleString()}
                </Badge>
              )}
            </button>
          );
        })}
      </div>

      {/* Separator */}
      <div className="hidden h-10 w-px bg-border sm:block" />

      {/* Edge-confidence threshold card */}
      <div className="flex flex-col gap-2 rounded-md border border-border bg-background/60 px-3 py-2">
        <span className="text-[11px] font-semibold uppercase tracking-wide text-muted-foreground">
          Edge confidence thresholds
        </span>

        {/* STRING */}
        <div className="flex items-center gap-2">
          <label
            htmlFor="string-threshold"
            className="min-w-[160px] text-xs text-foreground"
          >
            STRING &ge;{" "}
            <span className="font-mono text-xs">{thresholds.string.toFixed(2)}</span>
          </label>
          <input
            id="string-threshold"
            type="range"
            min={0.4}
            max={1.0}
            step={0.05}
            value={thresholds.string}
            onChange={(e) => updateThreshold("string", parseFloat(e.target.value))}
            className="h-1.5 w-40 cursor-pointer accent-primary"
          />
        </div>

        {/* LOCO replication */}
        <div className="flex items-center gap-2">
          <label
            htmlFor="loco-threshold"
            className="min-w-[160px] text-xs text-foreground"
          >
            LOCO replication &ge;{" "}
            <span className="font-mono text-xs">{thresholds.loco.toFixed(2)}</span>
            <span className="ml-1 text-[10px] text-muted-foreground">
              (Pre-reg Criterion 1: median 0.875)
            </span>
          </label>
          <input
            id="loco-threshold"
            type="range"
            min={0.0}
            max={1.0}
            step={0.05}
            value={thresholds.loco}
            onChange={(e) => updateThreshold("loco", parseFloat(e.target.value))}
            className="h-1.5 w-40 cursor-pointer accent-primary"
          />
        </div>

        {/* PP.H4 */}
        <div className="flex items-center gap-2">
          <label
            htmlFor="pp4-threshold"
            className="min-w-[160px] text-xs text-foreground"
          >
            PP.H4 &ge;{" "}
            <span className="font-mono text-xs">{thresholds.pp4.toFixed(2)}</span>
            <span className="ml-1 text-[10px] text-muted-foreground">
              (Pre-reg Criterion 1)
            </span>
          </label>
          <input
            id="pp4-threshold"
            type="range"
            min={0.0}
            max={1.0}
            step={0.05}
            value={thresholds.pp4}
            onChange={(e) => updateThreshold("pp4", parseFloat(e.target.value))}
            className="h-1.5 w-40 cursor-pointer accent-primary"
          />
        </div>
      </div>

      {/* Separator */}
      <div className="hidden h-10 w-px bg-border sm:block" />

      {/* Label top-K slider */}
      <div className="flex items-center gap-2">
        <label
          htmlFor="label-top-k"
          className="text-xs font-medium text-muted-foreground whitespace-nowrap"
        >
          Labels: top{" "}
          <span className="font-mono text-foreground">{labelTopK}</span>
        </label>
        <input
          id="label-top-k"
          type="range"
          min={0}
          max={100}
          step={5}
          value={labelTopK}
          onChange={(e) => onLabelTopKChange(parseInt(e.target.value, 10))}
          className="h-1.5 w-28 cursor-pointer accent-primary"
        />
      </div>

      {/* Separator */}
      <div className="hidden h-10 w-px bg-border sm:block" />

      {/* Node coloring dropdown */}
      <div className="flex items-center gap-2">
        <label
          htmlFor="color-by"
          className="text-xs font-medium text-muted-foreground whitespace-nowrap"
        >
          Color:
        </label>
        <select
          id="color-by"
          value={colorBy}
          onChange={(e) => onColorByChange(e.target.value as NodeColorBy)}
          className="h-7 rounded-md border border-input bg-transparent px-2 text-xs text-foreground outline-none focus-visible:border-ring focus-visible:ring-2 focus-visible:ring-ring/50"
        >
          {COLOR_BY_OPTIONS.map((opt) => (
            <option key={opt.value} value={opt.value}>
              {opt.label}
            </option>
          ))}
        </select>
      </div>
    </div>
  );
}
