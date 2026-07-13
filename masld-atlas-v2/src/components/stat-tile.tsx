"use client";

import type { ReactNode } from "react";
import { ArrowDown, ArrowUp, Minus } from "lucide-react";
import { cn } from "@/lib/utils";
import { CountUp } from "@/components/motion/count-up";

export type DeltaDirection = "up" | "down" | "neutral";
export type StatTileVariant = "default" | "hero";

export interface StatDelta {
  /** Change value to display next to the arrow (e.g. "+12%" or 12). */
  value: ReactNode;
  /** Arrow + tint direction. Defaults to "neutral" (muted, no color). */
  direction?: DeltaDirection;
}

export interface StatTileProps {
  /** Small uppercase label above the value. */
  label: ReactNode;
  /**
   * The headline value (rendered in tabular numerals). Kept for callers that
   * pass a preformatted string / node (e.g. `fmt(...)`). Ignored when
   * `numericValue` is provided.
   */
  value?: ReactNode;
  /**
   * Raw number to animate up via <CountUp>. When set, takes precedence over
   * `value` and is formatted with `format`. Backward-compatible: existing
   * callers keep passing `value` and see no change.
   */
  numericValue?: number;
  /** Formatter for `numericValue` (default: rounded, locale-grouped). */
  format?: (n: number) => string;
  /** Optional secondary line below the value. */
  sublabel?: ReactNode;
  /** Optional signed delta chip (arrow + value). */
  delta?: StatDelta;
  /** Optional icon element (e.g. a lucide icon) shown top-right, muted. */
  icon?: ReactNode;
  /** Visual weight. "hero" adds a brand sheen, glow accent bar, and lift. */
  variant?: StatTileVariant;
  className?: string;
}

/**
 * Delta tint uses the decorative chrome tokens (`--color-delta-*`), NOT a
 * data-mark color — a StatTile delta is UI chrome, not a chart mark.
 */
const DELTA_STYLES: Record<DeltaDirection, string> = {
  up: "text-delta-up",
  down: "text-delta-down",
  neutral: "text-muted-foreground",
};

const DELTA_ICON = {
  up: ArrowUp,
  down: ArrowDown,
  neutral: Minus,
} as const;

/**
 * KPI tile used across atlas pages: a subtle card with a small label, a large
 * numeric value, an optional sublabel and delta chip, and an optional icon.
 * `variant="hero"` upgrades it to an immersive tile (gradient sheen + glow
 * accent bar + display-scale numeral) for landing / summary rows.
 */
export function StatTile({
  label,
  value,
  numericValue,
  format,
  sublabel,
  delta,
  icon,
  variant = "default",
  className,
}: StatTileProps) {
  const dir = delta?.direction ?? "neutral";
  const DeltaIcon = DELTA_ICON[dir];
  const isHero = variant === "hero";

  const numeral = (
    <span
      className={cn(
        "font-semibold leading-none",
        isHero
          ? "font-display-numeric text-display-sm"
          : "font-numeric text-2xl tracking-tight"
      )}
    >
      {numericValue != null ? (
        <CountUp value={numericValue} format={format} />
      ) : (
        value
      )}
    </span>
  );

  return (
    <div
      className={cn(
        "relative flex flex-col gap-1 rounded-lg border border-border bg-card px-4 py-3.5",
        isHero
          ? "hover-lift overflow-hidden border-border/60 shadow-[var(--shadow-elev-2)]"
          : "shadow-[var(--shadow-elev-1)]",
        className
      )}
      style={isHero ? { backgroundImage: "var(--gradient-surface)" } : undefined}
    >
      {isHero && (
        <span
          aria-hidden
          className="absolute inset-y-2 left-0 w-1 rounded-full"
          style={{
            background: "var(--gradient-brand)",
            boxShadow: "var(--glow-primary)",
          }}
        />
      )}
      <div className={cn("flex items-start justify-between gap-2", isHero && "pl-2")}>
        <span className="text-[11px] font-semibold uppercase tracking-wider text-muted-foreground">
          {label}
        </span>
        {icon && (
          <span className="shrink-0 text-muted-foreground [&_svg]:size-4">
            {icon}
          </span>
        )}
      </div>
      <div className={cn("flex items-baseline gap-2", isHero && "pl-2")}>
        {numeral}
        {delta && (
          <span
            className={cn(
              "inline-flex items-center gap-0.5 text-xs font-medium tabular-nums",
              DELTA_STYLES[dir]
            )}
          >
            <DeltaIcon className="size-3" />
            {delta.value}
          </span>
        )}
      </div>
      {sublabel && (
        <span className={cn("text-xs text-muted-foreground", isHero && "pl-2")}>
          {sublabel}
        </span>
      )}
    </div>
  );
}
