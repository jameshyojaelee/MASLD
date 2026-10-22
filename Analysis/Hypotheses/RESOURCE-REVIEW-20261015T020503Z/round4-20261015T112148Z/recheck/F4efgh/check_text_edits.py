import json, re, sys
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
D=R+"/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z"
d=json.load(open(D+"/repair/F4efgh/repair_output_F4efgh.json"))
tags=re.compile(r"\[(Inference|Unverified|Speculation|Estimate|Contains Unverified Content)\]")
rows=[]
for e in d["entries"]:
    print("FID",e["fid"],"new_value:",e.get("new_value"))
    for t in e.get("text_edits",[]):
        txt=open(R+"/"+t["file"]).read()
        c_old=txt.count(t["old_text"]); c_new=txt.count(t["new_text"])
        # line of first occurrence
        ln = txt[:txt.find(t["old_text"])].count("\n")+1 if c_old else None
        emd = "—" in t["new_text"]
        en = "–" in t["new_text"]
        tg = bool(tags.search(t["new_text"]))
        print(f"  {t['file']}:{t['line']} old_count={c_old} at_line={ln} new_already={c_new} emdash={emd} endash={en} tag={tg} claim={t.get('claim_level')}")
