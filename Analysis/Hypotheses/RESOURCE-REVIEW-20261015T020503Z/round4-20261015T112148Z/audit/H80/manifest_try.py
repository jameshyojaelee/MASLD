# Audit H80: try plausible manifest recipes (sha256 of sorted "path size" listing) for each pinned directory.
import os, glob, hashlib, itertools, sys
ROOTS = ["/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"]
FM = "GWAS/finemapping"
targets = {
 "eqtl_susie_polyfun": ("results/eqtl_susie_polyfun", "baa3f20a83346e7e6dcdfcada4993a1b51fcbf25ef5692f4f6a4e9a323d3ad72", "*_susie.rds"),
 "polyfun_eur": ("data/ld_ref/polyfun_eur", "5f47504b43574fc3e4fd9d8390bc4a58c3d4940304055eb367ef01a378110f86", None),
 "1kg_eur": ("data/ld_ref/1kg_eur", "4499fa48b3e090fb3ed3bd66fea9a31e66bfade80b74b1a764fee13c22bd5177", None),
 "1kg_afr_gram": ("data/ld_ref/1kg_afr_gram", "6aa5b0f4e431db2e9aece99e52220281c7a4d1456af5ad393abba7f0804e8700", None),
 "1kg_amr_gram": ("data/ld_ref/1kg_amr_gram", "6b78be707b6fd69e5bd914051a25fd03585323ce5359014c6c353d165f7aa1a6", None),
 "1kg_eas_gram": ("data/ld_ref/1kg_eas_gram", "5860814ac8758832ffa5f669b4796fdf8a25a36723a653a8191c8759e75e0df0", None),
 "1kg_sas_gram": ("data/ld_ref/1kg_sas_gram", "cf73b2ef08817bb2c987dfdd3e51035e923bc5a87befac846aeece0228b0897c", None),
}
which = sys.argv[1:] or list(targets)
for name in which:
    rel, target, pat = targets[name]
    base = f"{ROOTS[0]}/{FM}/{rel}"
    files = []
    for dp, dn, fn in os.walk(base):
        for f in fn:
            p = os.path.join(dp, f)
            if pat and not glob.fnmatch.fnmatch(f, pat):
                continue
            if os.path.isfile(p):
                files.append((p, os.path.getsize(p), int(os.path.getmtime(p))))
    print(name, "files", len(files), flush=True)
    found = None
    pathforms = {
        "abs_atlas": lambda p: p,
        "abs_lib": lambda p: p.replace(ROOTS[0], ROOTS[1]),
        "rel_root": lambda p: os.path.relpath(p, ROOTS[0]),
        "rel_fm": lambda p: os.path.relpath(p, f"{ROOTS[0]}/{FM}"),
        "rel_dir": lambda p: os.path.relpath(p, base),
        "basename": lambda p: os.path.basename(p),
    }
    for (pn, pf), sep, term, mt in itertools.product(pathforms.items(), [" ", "\t"], ["nl_trailing", "nl_join"], [False, True]):
        lines = sorted(f"{pf(p)}{sep}{s}" + (f"{sep}{m}" if mt else "") for p, s, m in files)
        txt = "\n".join(lines) + ("\n" if term == "nl_trailing" else "")
        h = hashlib.sha256(txt.encode()).hexdigest()
        if h == target:
            found = (pn, repr(sep), term, mt); break
    print(name, "MATCH" if found else "NO MATCH", found, flush=True)
