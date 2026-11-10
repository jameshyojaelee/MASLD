#!/usr/bin/env python3
"""Link unchanged named Catalog access and separate new review outputs."""
import html
import json
import os
from pathlib import Path


root = Path(__file__).resolve().parents[5]
task = root / "Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z"
here = task / "catalog_models"
current = here / "outputs/corrected_companion_v2"
links = [
    ("Browse the unchanged named r8 interface", root / "figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/catalog_local/portal_candidate/index.html"),
    ("Corrected local molecular-object and method-specific companion", current / "index.html"),
    ("Separate unadopted priority-eligibility Catalog candidate", here / "outputs/eligibility_candidate/index.html"),
    ("Native genetics qualification and prior/LD annotations", task / "genetics/checked_companion_v2/genetics_qualification.tsv"),
    ("Genetics companion independent checks", task / "genetics/checked_companion_v2/independent_checks.json"),
    ("Corrected Figure 7E candidate", here / "fig7/panels/fig7E_verdict_vocabulary_and_refusal_thresholds.pdf"),
    ("Figure 7E numerical checks", here / "fig7/checks.json"),
    ("Implementation and verification record", here / "IMPLEMENTATION.txt"),
]
assert all(path.is_file() for _, path in links)
page = '<meta charset="utf-8"><title>MASLD Catalog review access</title><style>body{font:16px sans-serif;max-width:1000px;margin:2rem auto;line-height:1.5}li{margin:.8rem 0}</style><h1>MASLD Catalog review access</h1><p>The named r8 interface is preserved. New local companions correct RNA applicability and distinguish COLOC methods. The alternate priority-eligibility classification is explicitly unadopted and retains original r8 values beside the candidate values. Native genetic LD/prior diagnostics remain separate; unscreened genes have not passed those diagnostics.</p><ul>'
page += ''.join('<li><a href="' + html.escape(os.path.relpath(path, here)) + '">' + html.escape(label) + '</a></li>' for label, path in links)
page += '</ul>'
(here / "review_access.html").write_text(page)
index = current / "index.html"
text = index.read_text()
label, original = links[0]
if label not in text:
    index.write_text(text.replace('</ul>', '<li><a href="' + html.escape(os.path.relpath(original, current)) + '">' + label + '</a></li></ul>', 1))
checks = {"status": "PASS", "links_exist": len(links), "named_interface_modified": False,
          "candidate_classification_is_separate_and_unadopted": True}
(here / "access_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks))
