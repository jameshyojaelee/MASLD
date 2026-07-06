// update_links.jsx
// Force-reload every LINKED (placed) file in the active Illustrator document
// from disk, preserving each item's position and scale. Run this after the
// figure PDFs have been re-synced from the HPC.
//
// Run it either way:
//   * Manually:  Illustrator > File > Scripts > Other Script...  (pick this file)
//   * Via MCP:   ask the illustrator MCP to "run the JSX at <path>"
//   * Via shell: osascript -e 'tell application "Adobe Illustrator" to do javascript POSIX file "<path>"'
//
// Notes
//   * Only LINKED placed items (PlacedItem) are touched. Embedded art is skipped.
//   * Reassigning pi.file re-imports the artwork from disk while keeping the
//     PlacedItem's transformation matrix, so panels do NOT move or resize.
//   * No alert() is used so the script never blocks a headless MCP/osascript run;
//     the summary string is returned (do javascript returns the last value).

#target illustrator

(function () {
    if (app.documents.length === 0) {
        return "ERROR: no document open.";
    }
    var doc     = app.activeDocument;
    var items   = doc.placedItems;
    var n       = items.length;
    var updated = 0, missing = 0, skipped = 0, failed = 0;
    var log     = [];

    for (var i = 0; i < n; i++) {
        var pi = items[i];
        var f;
        try { f = pi.file; }               // the linked file on disk
        catch (e) { failed++; log.push("[no-file] item " + i); continue; }

        if (f === null) { skipped++; continue; }               // embedded, skip
        if (!f.exists)  { missing++; log.push("[missing] " + f.fsName); continue; }

        try {
            pi.file = f;                    // force re-read from disk, in place
            updated++;
        } catch (e2) {
            failed++; log.push("[failed] " + f.fsName + " :: " + e2);
        }
    }

    if (doc.redraw) { try { doc.redraw(); } catch (e3) {} }

    var msg = "placed items: " + n +
              " | updated: " + updated +
              " | missing: " + missing +
              " | embedded/skipped: " + skipped +
              " | failed: " + failed;
    if (log.length) { msg += "\n" + log.join("\n"); }
    $.writeln(msg);
    return msg;
})();
