#!/usr/bin/env python3
"""Process a mega-review team workflow result: write team report md, append to master CSV, print summary.
Usage: _process_team.py <task_output_file> <team_id> <report_basename> [<csv_team_id>]
"""
import json, sys, csv, os
from collections import Counter

out, team, report_base = sys.argv[1], sys.argv[2], sys.argv[3]
csv_team = sys.argv[4] if len(sys.argv) > 4 else team
DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/audit/2026-06-13_megareview"

w = json.load(open(out))
res = w.get("result", w)
if isinstance(res, str):
    res = json.loads(res)
finds = res.get("findings", []) or []
repro = res.get("repro", []) or []
rep = res.get("report", {}) or {}

# 1) team report
hdr = f"""# Team {team} — Mega-Review Report (2026-06-13)

*Agents: {w.get('agentCount')} · raw findings: {len(finds)} · number-reproduction checks: {len(repro)} · headline tally (not-refuted): {rep.get('tally')}*

*Read-only audit. Findings only; no code modified. Severity/verdict per `REVIEW_CONTEXT.md`.*

---

"""
open(f"{DIR}/{report_base}", "w").write(hdr + rep.get("markdown_report", "(no synthesizer report returned)"))

# 2) append to master CSV
cols = ["id","team","severity","type","adversarial_verdict","title","file","line","claim","evidence","recommendation","lit_question"]
path = f"{DIR}/00_master_findings.csv"
new = not os.path.exists(path)
# find current max index for this csv_team to avoid id collisions
start = 0
if not new:
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if row.get("team") == csv_team and row.get("id","").startswith(csv_team+"-"):
                try: start = max(start, int(row["id"].split("-")[-1]))
                except: pass
with open(path, "a", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    if new: wr.writeheader()
    for i, f in enumerate(finds, start+1):
        row = {k: f.get(k, "") for k in cols}
        row["id"] = f"{csv_team}-{i:04d}"; row["team"] = csv_team
        wr.writerow(row)
print(f"wrote {DIR}/{report_base}; appended {len(finds)} rows as {csv_team} (ids {start+1}..{start+len(finds)})")

# 3) summary
print("severity", dict(Counter(x.get("severity") for x in finds)))
print("type", dict(Counter(x.get("type") for x in finds)))
print("verdict", dict(Counter(x.get("adversarial_verdict") for x in finds)))
def sk(s): return {"Critical":0,"High":1,"Medium":2,"Low":3}.get(s,4)
conf = sorted([f for f in finds if f.get("adversarial_verdict")=="confirmed" and f.get("severity") in ("Critical","High")], key=lambda f: sk(f["severity"]))
print(f"\n=== CONFIRMED CRITICAL/HIGH ({len(conf)}) ===")
for f in conf:
    print(f"[{f['severity']}/{f['type']}] {f['title']}\n   {f.get('file','')}:{f.get('line','')}")
crit = [f for f in finds if f["severity"]=="Critical"]
print(f"\n=== ALL CRITICAL ({len(crit)}, any verdict) ===")
for f in sorted(crit, key=lambda f: f.get("adversarial_verdict","")):
    print(f"[{f.get('adversarial_verdict')}] {f['title']}  ::  {f.get('file','')}")
print("\n=== NUMBER-REPRODUCTION ISSUES (matches != yes) ===")
for c in repro:
    if c.get("matches") != "yes":
        print(f"[{c.get('matches')}] {c.get('number_name')}: claimed={str(c.get('claimed'))[:120]} | reproduced={str(c.get('reproduced'))[:160]}")
