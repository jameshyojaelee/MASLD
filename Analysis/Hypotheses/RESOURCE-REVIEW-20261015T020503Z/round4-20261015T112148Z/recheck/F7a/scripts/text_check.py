import json, re, sys
root="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
E=json.load(open(sys.argv[1]))
rep=json.load(open(root+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7a/exec/text_edits.json"))
repnew=[e["new_text"] for v in rep.values() for e in (v if isinstance(v,list) else v.get("text_edits",[v]))] if isinstance(rep,dict) else []
for f in E:
  for e in f["text_edits"]:
    txt=open(root+e["file"]).read()
    n=txt.count(e["old_text"]); nn=txt.count(e["new_text"])
    lineno = txt[:txt.find(e["old_text"])].count("\n")+1 if n else None
    new=e["new_text"]
    flags=[]
    if "—" in new: flags.append("EM_DASH")
    if "–" in new: flags.append("EN_DASH")
    if re.search(r"\[(Inference|Unverified|Speculation|Estimate)", new): flags.append("EVIDENCE_TAG")
    for w in ["progression","causal gene","proves","demonstrates","novel","first "]:
      if w in new.lower(): flags.append("WORD:"+w)
    inrep = any(new==r for r in repnew)
    print(f["fid"], e["file"], "given_line", e["line"], "old_count", n, "found_line", lineno, "new_already", nn, "flags", flags, "matches_repair_json", inrep)
