#!/usr/bin/env python3
"""Render all 5 Fig 4A v2 candidates from ONE live compute.

get_data() runs the canonical compute_counts() (atlas/DEG/proteomics read) + the
10k-draw permutation null ONCE; the result is shared across all five builders, so
the expensive null is not recomputed five times. Each builder that raises is
caught and reported (the others still render), so one broken panel never blocks
the rest. Run on a compute node in the rnaseq env (see run_fig4a_v2.sh)."""
import traceback

from fig4a_v2_common import get_data
import fig4a_v2_corefield
import fig4a_v2_physlayers
import fig4a_v2_icons
# forest / ghostflow / triangle were CUT by the user 2026-07-10 (PDFs deleted); their
# builder scripts remain on disk but are no longer auto-rendered here so a re-render
# does not resurrect the deleted candidates. corefield ("Null Lineup") is the WINNER.


def main():
    d, stats = get_data()
    builders = [
        ("corefield",  fig4a_v2_corefield.build),   # WINNER (Null Lineup)
        ("physlayers", fig4a_v2_physlayers.build),   # fallback
        ("icons",      fig4a_v2_icons.build),        # fallback (PI's assayed-vs-validated)
    ]
    ok, fail = [], []
    for name, fn in builders:
        print(f"\n----- rendering {name} -----")
        try:
            fn(d, stats)
            ok.append(name)
        except Exception as e:               # noqa: BLE001 — report, keep going
            traceback.print_exc()
            fail.append((name, repr(e)))
    print("\n[render_all] OK:", ok)
    if fail:
        print("[render_all] FAILED:", fail)
    else:
        print("[render_all] all 5 candidates rendered.")


if __name__ == "__main__":
    main()
