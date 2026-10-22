"""Select biopsy pairs and apply the sex rule for N2a. Reads only metadata (no allele counts).
Usage: python -I 02_select_pairs.py OUT_TSV     (refuses to overwrite)

Sets: primary = pair_type 'different_visits'; secondary = pair_type 'not_established' (reported separately).
Sex rule (person level; the two libraries are the same person by source IDs or genotype): female when at least
one library has an expression-based F call and no source in either library says M or 'both'. Libraries with an
'ambiguous' or missing call are neutral. Sources: N1 expr_sex (library_table.tsv), refit inferred_sexM
(model_design.tsv, 1 = male), GEO 'Sex' (GSE130970 only).
Libraries that N1 leaves as called mixed (GSE193066) are excluded."""
import sys, os
import pandas as pd

out = sys.argv[1]
if os.path.exists(out):
    sys.exit("refusing to overwrite " + out)
H = "Analysis/Hypotheses/"
t = pd.read_csv(H + "N4-20261007T210855Z/pair_type_table_with_n1.tsv", sep="\t")
n1lib = pd.read_csv(H + "N1-20261007T210855Z/library_table.tsv", sep="\t", usecols=["run", "cohort", "expr_sex"]).drop_duplicates("run").set_index("run")
n1call = pd.read_csv(H + "N1-20261007T210855Z/gse193066_library_calls.tsv", sep="\t", usecols=["run", "n1_status"]).set_index("run")
md = pd.read_csv("RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/model_design.tsv",
                 sep="\t", usecols=["sample_id", "inferred_sexM"]).set_index("sample_id")
geo = pd.read_csv(H + "N4-20261007T210855Z/source_metadata/GSE130970_samples_joined.tsv", sep="\t", usecols=["Run", "Sex"]).drop_duplicates("Run").set_index("Run")
man = pd.read_csv(H + "READS-20261009T141017Z/library_manifest.tsv", sep="\t").drop_duplicates("library").set_index("library")


def lib_sex(lib):
    src = {}
    if lib in n1lib.index and n1lib.loc[lib, "expr_sex"] in ("F", "M", "both"):
        src["n1_expr"] = n1lib.loc[lib, "expr_sex"]
    if lib in md.index:
        src["refit_inferred"] = "M" if int(md.loc[lib, "inferred_sexM"]) == 1 else "F"
    if lib in geo.index and geo.loc[lib, "Sex"] in ("F", "M"):
        src["geo_recorded"] = geo.loc[lib, "Sex"]
    return src


rows = []
for _, r in t.iterrows():
    libs = r["libraries"].split(",")
    stages = r["stage_codes"].split(",") if isinstance(r["stage_codes"], str) else [""] * len(libs)
    rec = dict(person_id=r["person_id"], cohort=r["cohort"], set={"different_visits": "primary", "not_established": "secondary"}.get(r["pair_type"], "other"),
               pair_type=r["pair_type"], n_libraries=len(libs), n1_eligible=r["n1_eligible"], include=False, reason="")
    if len(libs) != 2:
        rec.update(reason="three_libraries_no_prespecified_pair", libs=",".join(libs))
        rows.append(rec); continue
    rec.update(lib1=libs[0], lib2=libs[1], stage1=stages[0], stage2=stages[1])
    reasons = []
    allv, has_f_expr = set(), False
    for k, lib in enumerate(libs, 1):
        s = lib_sex(lib)
        rec["sex_sources_%d" % k] = ";".join("%s=%s" % kv for kv in s.items())
        allv |= set(s.values())
        has_f_expr |= any(s.get(x) == "F" for x in ("n1_expr", "refit_inferred"))
    if allv - {"F"}:
        reasons.append("male_or_conflicting_sex")
    elif not has_f_expr:
        reasons.append("no_expression_female_call")
    if r["cohort"] == "GSE193066":
        for k, lib in enumerate(libs, 1):
            st = n1call.loc[lib, "n1_status"] if lib in n1call.index else "missing"
            rec["n1_status_%d" % k] = st
            if st != "not_supported_as_mixed":
                reasons.append("n1_%s_lib%d" % (st, k))
    for k, lib in enumerate(libs, 1):
        if lib not in man.index:
            reasons.append("no_bam_in_manifest_lib%d" % k)
        else:
            rec["bam%d" % k] = man.loc[lib, "bam_path"]
            rec["layout%d" % k] = man.loc[lib, "layout"]
    rec["reason"] = ";".join(reasons)
    rec["include"] = not reasons
    try:
        rec["abs_stage_change"] = abs(float(stages[1]) - float(stages[0]))
    except ValueError:
        rec["abs_stage_change"] = float("nan")
    rows.append(rec)

d = pd.DataFrame(rows)
d.to_csv(out, sep="\t", index=False)
print(d.groupby(["set", "cohort", "include"]).size())
print(d[~d["include"]].groupby(["set", "cohort", "reason"]).size())
print(d[d["include"]][["set", "cohort", "person_id", "stage1", "stage2", "layout1", "layout2"]].to_string())
