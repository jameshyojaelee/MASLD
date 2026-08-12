/**
 * Fail closed while the Next.js portal still exposes the legacy rank-first UI.
 *
 * ALLOW_LEGACY_PORTAL_BUILD=true permits local provenance builds only. It is
 * not deployment or publication authorization. Plan 50/60 must remove this
 * legacy gate when the MASLD Gene Catalog portal is promoted.
 */

if (process.env.ALLOW_LEGACY_PORTAL_BUILD !== "true") {
  console.error(
    "REFUSED: this is the legacy pre-Resource portal. " +
      "Build the MASLD Gene Catalog candidate and obtain final release promotion."
  );
  process.exit(64);
}

console.warn(
  "LEGACY PROVENANCE BUILD ONLY: no publication or deployment authorization."
);
