import type { ReactNode } from "react";
import { cn } from "@/lib/utils";
import {
  ancestryColor,
  diseaseStageColor,
  sequentialColor,
  CONTROL,
} from "@/lib/palette";

/**
 * Small pill primitives for categorical / evidence annotations. All are
 * colorblind-safe: color is carried by a leading dot, never by text alone, and
 * the label is always present.
 */

const PILL_BASE =
  "inline-flex w-fit items-center gap-1.5 whitespace-nowrap rounded-full border border-border bg-muted/40 px-2 py-0.5 text-xs font-medium text-foreground";

function Dot({ color }: { color: string }) {
  return (
    <span
      aria-hidden
      className="size-2 shrink-0 rounded-full"
      style={{ backgroundColor: color }}
    />
  );
}

interface DotPillProps {
  color: string;
  children: ReactNode;
  className?: string;
}

function DotPill({ color, children, className }: DotPillProps) {
  return (
    <span className={cn(PILL_BASE, className)}>
      <Dot color={color} />
      {children}
    </span>
  );
}

/** Ancestry pill (Okabe-Ito colorblind-safe palette). */
export function AncestryChip({
  ancestry,
  className,
}: {
  ancestry: string;
  className?: string;
}) {
  return (
    <DotPill color={ancestryColor(ancestry)} className={className}>
      {ancestry}
    </DotPill>
  );
}

/** Disease-stage pill (Healthy = control gray; progression light → dark). */
export function StageChip({
  stage,
  className,
}: {
  stage: string;
  className?: string;
}) {
  return (
    <DotPill color={diseaseStageColor(stage)} className={className}>
      {stage}
    </DotPill>
  );
}

/** Expression / decomposition program pill (caller supplies the program color). */
export function ProgramChip({
  label,
  color = CONTROL,
  className,
}: {
  label: ReactNode;
  color?: string;
  className?: string;
}) {
  return (
    <DotPill color={color} className={className}>
      {label}
    </DotPill>
  );
}

/**
 * Evidence-support badge: `count`/`total` layers, with the dot shaded by the
 * support ratio (light → deep) so stronger support reads darker. A subtle
 * decorative outer glow (chrome, NOT a data mark — the sequential dot carries
 * the data) scales with the support ratio so stronger support quietly lifts.
 */
export function SupportBadge({
  count,
  total = 8,
  label = "support",
  className,
}: {
  count: number;
  total?: number;
  label?: ReactNode;
  className?: string;
}) {
  const ratio = total > 0 ? Math.max(0, Math.min(1, count / total)) : 0;
  const glow =
    ratio > 0.05
      ? `0 0 0 1px color-mix(in oklab, var(--primary) ${Math.round(
          ratio * 22
        )}%, transparent), 0 2px 10px -4px color-mix(in oklab, var(--primary) ${Math.round(
          ratio * 42
        )}%, transparent)`
      : undefined;
  return (
    <span
      className={cn(PILL_BASE, className)}
      style={glow ? { boxShadow: glow } : undefined}
    >
      <Dot color={sequentialColor(ratio)} />
      <span className="tabular-nums">
        {count}/{total}
      </span>
      {label && <span className="text-muted-foreground">{label}</span>}
    </span>
  );
}

/**
 * Convergence tier badge (Tier 1 = strongest). Intensity encodes tier so
 * higher tiers recede visually.
 */
export function TierBadge({
  tier,
  className,
}: {
  tier: number | string;
  className?: string;
}) {
  const n =
    typeof tier === "number" ? tier : parseInt(String(tier).replace(/\D/g, ""), 10);
  const style =
    n === 1
      ? "border-primary/40 bg-primary/10 text-primary"
      : n === 2
        ? "border-border bg-muted text-foreground"
        : "border-border bg-muted/50 text-muted-foreground";
  // Tier 1 (strongest) gets a subtle decorative brand glow — chrome only; the
  // tier text + tint still carry the category.
  const glow =
    n === 1
      ? "0 0 0 1px color-mix(in oklab, var(--primary) 30%, transparent), 0 2px 12px -4px color-mix(in oklab, var(--primary) 45%, transparent)"
      : undefined;
  return (
    <span
      className={cn(
        "inline-flex w-fit items-center whitespace-nowrap rounded-full border px-2 py-0.5 text-xs font-medium tabular-nums",
        style,
        className
      )}
      style={glow ? { boxShadow: glow } : undefined}
    >
      Tier {Number.isNaN(n) ? tier : n}
    </span>
  );
}
