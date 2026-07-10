import type { ReactNode } from "react";
import { Inbox } from "lucide-react";
import { cn } from "@/lib/utils";

/** A pulsing placeholder block for loading states. Size via className. */
export function SkeletonBlock({ className }: { className?: string }) {
  return (
    <div
      className={cn("h-4 w-full animate-pulse rounded-md bg-muted", className)}
    />
  );
}

export interface EmptyStateProps {
  /** Icon element (defaults to an inbox). */
  icon?: ReactNode;
  title: ReactNode;
  description?: ReactNode;
  /** Optional action (e.g. a "clear filters" button). */
  action?: ReactNode;
  className?: string;
}

/** Centered empty / no-results state inside a dashed card. */
export function EmptyState({
  icon,
  title,
  description,
  action,
  className,
}: EmptyStateProps) {
  return (
    <div
      className={cn(
        "flex min-h-[220px] flex-col items-center justify-center rounded-lg border border-dashed border-border bg-card/40 p-8 text-center",
        className
      )}
    >
      <span className="mb-3 text-muted-foreground [&_svg]:size-6">
        {icon ?? <Inbox />}
      </span>
      <p className="text-sm font-medium text-foreground">{title}</p>
      {description && (
        <p className="mt-1 max-w-md text-sm text-muted-foreground">
          {description}
        </p>
      )}
      {action && <div className="mt-4">{action}</div>}
    </div>
  );
}

export interface LegendItem {
  color: string;
  label: ReactNode;
}

/** Swatch + label legend for encoding color mappings. */
export function Legend({
  items,
  orientation = "horizontal",
  className,
}: {
  items: LegendItem[];
  orientation?: "horizontal" | "vertical";
  className?: string;
}) {
  return (
    <ul
      className={cn(
        "flex gap-x-4 gap-y-1.5 text-xs text-muted-foreground",
        orientation === "vertical" ? "flex-col" : "flex-wrap items-center",
        className
      )}
    >
      {items.map((item, i) => (
        <li key={i} className="flex items-center gap-1.5">
          <span
            aria-hidden
            className="size-2.5 shrink-0 rounded-[3px]"
            style={{ backgroundColor: item.color }}
          />
          <span>{item.label}</span>
        </li>
      ))}
    </ul>
  );
}
