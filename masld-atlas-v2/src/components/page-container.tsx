import type { ReactNode } from "react";
import { cn } from "@/lib/utils";

/**
 * Standard reading-width page container. Unifies the per-page wrappers that
 * previously drifted between `max-w-5xl` and `max-w-6xl`. Full-bleed canvas
 * tools (e.g. the network explorer) intentionally opt out and do not use this.
 */
export function PageContainer({
  children,
  className,
}: {
  children: ReactNode;
  className?: string;
}) {
  return (
    <div className={cn("mx-auto w-full max-w-6xl px-6 py-8", className)}>
      {children}
    </div>
  );
}
