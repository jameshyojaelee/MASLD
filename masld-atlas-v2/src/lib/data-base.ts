/**
 * Single source of truth for the data root.
 *
 * Local dev / bundled build (env unset) resolves to `/data` (served from
 * `public/data/`). Production sets `NEXT_PUBLIC_DATA_BASE` to the Hugging Face
 * Dataset CDN (e.g. the `resolve/main` URL) at build time so the static export
 * fetches data from the CDN instead of shipping it.
 *
 * Every `/data/...` fetch/href in the app MUST go through `dataUrl()` so the
 * base is edited in exactly one place.
 */
export const DATA_BASE = process.env.NEXT_PUBLIC_DATA_BASE ?? "/data";

/**
 * Join a data-relative path onto DATA_BASE, tolerating a leading slash on the
 * argument. `dataUrl("atlas.parquet")` and `dataUrl("/atlas.parquet")` both
 * resolve to `${DATA_BASE}/atlas.parquet`.
 */
export const dataUrl = (p: string): string =>
  `${DATA_BASE}/${p.replace(/^\//, "")}`;
