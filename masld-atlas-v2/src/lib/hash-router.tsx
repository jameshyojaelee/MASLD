"use client";

// ---------------------------------------------------------------------------
// Hash-based client router (static-host deep-link safe).
//
// Hugging Face Static Spaces serve exact files + `/` -> index.html ONLY (no
// `.html`-append, no directory index, no 404 fallback — empirically confirmed).
// So a pathname deep-link like `#/gene?symbol=THRB` 404s on a cold load. Routing
// the app through the URL HASH (`/#/gene?symbol=THRB`) makes every request hit
// `/` -> index.html (200); the route lives in `location.hash`, read client-side.
// Changing only the hash never triggers a network request, so in-app links are
// plain `<a href="#/...">` anchors and hard-refresh / paste / share all work.
// ---------------------------------------------------------------------------

import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useState,
} from "react";

export interface HashLocation {
  /** Plain route path, e.g. "/", "/atlas", "/gene". */
  path: string;
  /** Query string incl. leading "?", or "". */
  search: string;
}

/** Parse a raw `location.hash` ("#/gene?symbol=X") into a plain path + search. */
export function parseHash(hash: string): HashLocation {
  let h = (hash || "").replace(/^#/, "");
  if (h === "") return { path: "/", search: "" };
  if (!h.startsWith("/")) h = "/" + h;
  const qi = h.indexOf("?");
  if (qi === -1) return { path: h, search: "" };
  return { path: h.slice(0, qi) || "/", search: h.slice(qi) };
}

/** Convert an internal path ("/atlas") to its hash href ("#/atlas"). External
 *  URLs, mailto/tel, and already-hash hrefs pass through unchanged. */
export function toHash(href: string): string {
  if (!href) return "#/";
  if (/^(https?:|mailto:|tel:)/i.test(href)) return href;
  if (href.startsWith("#")) return href;
  return "#" + (href.startsWith("/") ? href : "/" + href);
}

const HashContext = createContext<HashLocation>({ path: "/", search: "" });

export function HashRouterProvider({ children }: { children: React.ReactNode }) {
  // SSR / first paint = "/" (the prerendered index). A useEffect then reads the
  // real hash, so a cold deep-load briefly shows home before the routed view.
  const [loc, setLoc] = useState<HashLocation>({ path: "/", search: "" });
  useEffect(() => {
    const update = () => setLoc(parseHash(window.location.hash));
    update();
    window.addEventListener("hashchange", update);
    return () => window.removeEventListener("hashchange", update);
  }, []);
  return <HashContext.Provider value={loc}>{children}</HashContext.Provider>;
}

export function useHashLocation(): HashLocation {
  return useContext(HashContext);
}

/** Drop-in replacement for next/navigation `useSearchParams` (reads the hash). */
export function useHashSearchParams(): URLSearchParams {
  const { search } = useHashLocation();
  return useMemo(() => new URLSearchParams(search), [search]);
}

/** Drop-in replacement for `router.push` — navigate by setting the hash. */
export function useHashNavigate() {
  return useCallback((href: string) => {
    window.location.hash = toHash(href);
    // Ensure top-of-page on route change (mirrors a normal navigation).
    if (typeof window !== "undefined") window.scrollTo({ top: 0 });
  }, []);
}

/** In-place hash update with no history entry and no scroll — the equivalent of
 *  `router.replace(url, { scroll: false })`. Use for URL/query state sync (e.g.
 *  shareable filter state) where a re-route + history push is undesirable. Does
 *  NOT emit hashchange, so it will not re-trigger the router for the same path. */
export function replaceHash(href: string): void {
  if (typeof window !== "undefined") {
    window.history.replaceState(null, "", toHash(href));
  }
}
