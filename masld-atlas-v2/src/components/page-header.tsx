import type { ReactNode } from "react";
import { cn } from "@/lib/utils";

/**
 * Standard page title block. Matches the existing page convention
 * (`text-3xl font-bold tracking-tight` heading + muted description).
 */
export function PageHeader({
  title,
  description,
  eyebrow,
  actions,
  gradient = false,
  className,
}: {
  title: ReactNode;
  description?: ReactNode;
  /** Small uppercase kicker above the title (e.g. a section name). */
  eyebrow?: ReactNode;
  /** Right-aligned actions (buttons, toggles). */
  actions?: ReactNode;
  /** Render the title with the decorative brand gradient. */
  gradient?: boolean;
  className?: string;
}) {
  return (
    <div className={cn("mb-8 flex items-start justify-between gap-4", className)}>
      <div className="min-w-0">
        {eyebrow && (
          <p className="mb-1.5 text-[11px] font-semibold uppercase tracking-wider text-muted-foreground">
            {eyebrow}
          </p>
        )}
        <h1
          className={cn(
            "font-display text-3xl font-semibold tracking-tight text-balance",
            gradient && "text-gradient"
          )}
        >
          {title}
        </h1>
        {description && (
          <p className="mt-2 max-w-2xl text-sm text-muted-foreground">
            {description}
          </p>
        )}
      </div>
      {actions && <div className="flex shrink-0 items-center gap-2">{actions}</div>}
    </div>
  );
}

/**
 * Placeholder card for routes that are scaffolded but not yet built. Replaced
 * by the real page content as each layer lands.
 */
export function ComingSoon({ children }: { children?: ReactNode }) {
  return (
    <div className="flex min-h-[40vh] flex-col items-center justify-center rounded-lg border border-dashed border-border bg-card/40 p-10 text-center">
      <p className="text-sm font-medium text-foreground">Coming soon</p>
      <p className="mt-1 max-w-md text-sm text-muted-foreground">
        {children ?? "This view is under construction."}
      </p>
    </div>
  );
}
