"use client";

import { useMemo } from "react";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import { EVIDENCE_SOURCES } from "@/lib/colors";
import type { EvidenceStrengths } from "@/lib/types";

interface EvidenceFingerprintProps {
  evidence: EvidenceStrengths;
  size?: number;
  showTooltip?: boolean;
  className?: string;
}

export function EvidenceFingerprint({
  evidence,
  size = 24,
  showTooltip = true,
  className,
}: EvidenceFingerprintProps) {
  const axes = EVIDENCE_SOURCES;
  const n = axes.length;
  const cx = size / 2;
  const cy = size / 2;
  const radius = size / 2 - 1;

  const points = useMemo(() => {
    return axes.map((source, i) => {
      const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
      const val = evidence[source.key] ?? 0;
      const r = radius * val;
      return {
        x: cx + r * Math.cos(angle),
        y: cy + r * Math.sin(angle),
        val,
        source,
      };
    });
  }, [evidence, cx, cy, radius, n, axes]);

  const polygonPath = points.map((p) => `${p.x},${p.y}`).join(" ");

  const gridPath = Array.from({ length: n }, (_, i) => {
    const angle = (Math.PI * 2 * i) / n - Math.PI / 2;
    return `${cx + radius * Math.cos(angle)},${cy + radius * Math.sin(angle)}`;
  }).join(" ");

  const activeSources = points.filter((p) => p.val > 0).length;

  const svg = (
    <svg
      width={size}
      height={size}
      viewBox={`0 0 ${size} ${size}`}
      className={className}
      role="img"
      aria-label={`Evidence: ${activeSources}/${n} sources active`}
    >
      <polygon
        points={gridPath}
        fill="none"
        stroke="currentColor"
        strokeWidth={0.5}
        opacity={0.15}
      />
      <polygon
        points={polygonPath}
        fill="var(--color-primary)"
        fillOpacity={0.25}
        stroke="var(--color-primary)"
        strokeWidth={size > 40 ? 1.5 : 0.8}
        strokeOpacity={0.7}
      />
      {size > 40 &&
        points.map(
          (p, i) =>
            p.val > 0 && (
              <circle
                key={i}
                cx={p.x}
                cy={p.y}
                r={size > 80 ? 3 : 2}
                fill={p.source.color}
                stroke="var(--color-background)"
                strokeWidth={0.5}
              />
            )
        )}
    </svg>
  );

  if (!showTooltip) return svg;

  return (
    <TooltipProvider delay={200}>
      <Tooltip>
        <TooltipTrigger render={<span className="inline-flex cursor-default" />}>
          {svg}
        </TooltipTrigger>
        <TooltipContent side="right" className="max-w-xs">
          <div className="space-y-1">
            <p className="text-xs font-semibold">
              Evidence: {activeSources}/{n} sources
            </p>
            {points.map((p) => (
              <div
                key={p.source.key}
                className="flex items-center justify-between gap-4 text-xs"
              >
                <span className="text-muted-foreground">{p.source.short}</span>
                <span className="font-mono">
                  {p.val > 0 ? (p.val * 100).toFixed(0) + "%" : "—"}
                </span>
              </div>
            ))}
          </div>
        </TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}
