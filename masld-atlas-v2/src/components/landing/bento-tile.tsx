"use client";

/**
 * Reusable bento entry tile + a viewport-gated lazy mount for its preview.
 *
 * <BentoTile> renders a card (link, button, or plain div depending on props)
 * with a title, one-line description, optional stat, and a fixed-height preview
 * slot. <LazyMount> defers the preview's fetch/canvas work until the tile is
 * within ~200px of the viewport (IntersectionObserver), showing a low-cost
 * skeleton before then — so mini-viz work never blocks the landing's first
 * paint. Visual grammar matches the rest of the portal's cards
 * (rounded-lg / border / bg-card / hover-lift).
 */

import { useEffect, useRef, useState, type ReactNode } from "react";
import { HashLink } from "@/components/hash-link";
import { cn } from "@/lib/utils";

// ---------------------------------------------------------------------------
// LazyMount
// ---------------------------------------------------------------------------

export function LazyMount({
  children,
  className,
  rootMargin = "200px",
}: {
  children: ReactNode;
  className?: string;
  rootMargin?: string;
}) {
  const ref = useRef<HTMLDivElement>(null);
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    if (typeof IntersectionObserver === "undefined") {
      setVisible(true);
      return;
    }
    const io = new IntersectionObserver(
      (entries) => {
        if (entries.some((e) => e.isIntersecting)) {
          setVisible(true);
          io.disconnect();
        }
      },
      { rootMargin }
    );
    io.observe(el);
    return () => io.disconnect();
  }, [rootMargin]);

  return (
    <div ref={ref} className={className}>
      {visible ? (
        children
      ) : (
        <div
          aria-hidden
          className="h-full w-full animate-pulse rounded-md bg-muted/40 motion-reduce:animate-none"
        />
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// BentoTile
// ---------------------------------------------------------------------------

export interface BentoTileProps {
  title: string;
  description: string;
  stat?: ReactNode;
  href?: string;
  onClick?: () => void;
  className?: string;
  /** Wrap `children` in a viewport-gated LazyMount (default true). */
  lazy?: boolean;
  children?: ReactNode;
}

export function BentoTile({
  title,
  description,
  stat,
  href,
  onClick,
  className,
  lazy = true,
  children,
}: BentoTileProps) {
  const cls = cn(
    "group flex h-full flex-col overflow-hidden rounded-lg border border-border bg-card p-4 text-left",
    "hover-lift hover:border-primary/40 focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring",
    className
  );

  const header = (
    <div className="min-w-0">
      <span className="text-sm font-semibold group-hover:text-primary">{title}</span>
      <p className="mt-0.5 truncate text-xs text-muted-foreground">{description}</p>
      {stat != null && (
        <p className="mt-1 font-numeric text-[11px] text-muted-foreground">{stat}</p>
      )}
    </div>
  );

  const body = (
    <>
      {header}
      {children != null &&
        (lazy ? (
          <LazyMount className="mt-3 min-h-0 flex-1">{children}</LazyMount>
        ) : (
          <div className="mt-3 min-h-0 flex-1">{children}</div>
        ))}
    </>
  );

  if (onClick) {
    return (
      <button type="button" onClick={onClick} className={cls}>
        {body}
      </button>
    );
  }
  if (href) {
    return (
      <HashLink href={href} className={cls}>
        {body}
      </HashLink>
    );
  }
  return <div className={cls}>{body}</div>;
}
