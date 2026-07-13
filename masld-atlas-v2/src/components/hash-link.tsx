"use client";

import { forwardRef } from "react";
import { toHash } from "@/lib/hash-router";

/**
 * Drop-in replacement for `next/link`'s default export under hash routing.
 * Renders a plain anchor whose href is the hash form of the given route, so
 * files can migrate with a single import swap:
 *   `import Link from "next/link"` ->
 *   `import { HashLink as Link } from "@/components/hash-link"`.
 * next/link-only props (prefetch/scroll/replace) are accepted and ignored so
 * existing call sites type-check unchanged.
 */
type HashLinkProps = Omit<React.ComponentPropsWithoutRef<"a">, "href"> & {
  href: string;
  prefetch?: boolean;
  scroll?: boolean;
  replace?: boolean;
};

export const HashLink = forwardRef<HTMLAnchorElement, HashLinkProps>(
  function HashLink({ href, prefetch: _p, scroll: _s, replace: _r, ...rest }, ref) {
    return <a ref={ref} href={toHash(href)} {...rest} />;
  }
);

export default HashLink;
