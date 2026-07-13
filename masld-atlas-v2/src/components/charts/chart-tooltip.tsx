"use client";

import {
  useCallback,
  useRef,
  useState,
  type ReactNode,
} from "react";
import { cn } from "@/lib/utils";

/**
 * Lightweight hover-tooltip primitive shared by every chart.
 *
 * Positions relative to a wrapper element (the chart's own `relative` div) using
 * the pointer's client coordinates, so it stays correct regardless of the
 * chart's title/caption/margins. Pair with {@link ChartTooltip} for rendering.
 */
export interface TooltipState<T> {
  /** px offset from the left edge of the wrapper. */
  left: number;
  /** px offset from the top edge of the wrapper. */
  top: number;
  data: T;
}

export function useChartTooltip<T>() {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const [tooltip, setTooltip] = useState<TooltipState<T> | null>(null);

  const show = useCallback(
    (event: { clientX: number; clientY: number }, data: T) => {
      const rect = wrapperRef.current?.getBoundingClientRect();
      if (!rect) return;
      setTooltip({
        left: event.clientX - rect.left,
        top: event.clientY - rect.top,
        data,
      });
    },
    []
  );

  const hide = useCallback(() => setTooltip(null), []);

  return { wrapperRef, tooltip, show, hide };
}

export interface ChartTooltipProps {
  left: number;
  top: number;
  className?: string;
  children: ReactNode;
}

/**
 * Themed, theme-aware tooltip card. Rendered inside the chart's `relative`
 * wrapper and floated just above/right of the cursor. `pointer-events-none` so
 * it never steals hover from the marks underneath.
 */
export function ChartTooltip({ left, top, className, children }: ChartTooltipProps) {
  return (
    <div
      className={cn(
        "pointer-events-none absolute z-30 max-w-56 rounded-md border border-border",
        "bg-popover/95 px-2 py-1 text-xs leading-tight text-popover-foreground shadow-md",
        "tabular-nums backdrop-blur-sm",
        className
      )}
      style={{
        left,
        top,
        transform: "translate(-50%, calc(-100% - 10px))",
      }}
      role="tooltip"
    >
      {children}
    </div>
  );
}
