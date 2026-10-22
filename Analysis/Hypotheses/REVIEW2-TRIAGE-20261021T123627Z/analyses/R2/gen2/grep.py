import re,sys,os
R='Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z'
files=['docs/manuscript/draft/%s.md'%x for x in 'abstract title 01_intro fig1 fig2 fig3 fig4 fig5 fig5_discussion fig7 supplementary_notes references'.split()]+['docs/manuscript/METHODS.md','docs/manuscript/05_figure_legends.md','docs/PAPER.md']
def text(f):
    p=os.path.join(R,'staged/round1',f)
    return open(p if os.path.exists(p) else f).read()
pats=sys.argv[1:]
W=int(os.environ.get('W','220'))
for f in files:
    t=text(f)
    for p in pats:
        for m in re.finditer(p,t):
            s=max(0,m.start()-W); e=min(len(t),m.end()+W)
            ln=t.count('\n',0,m.start())+1
            print(f'--- {f}:{ln} [{p}]'); print(t[s:e].replace('\n',' '))
