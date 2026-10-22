import json,sys
root='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
E=json.load(open(sys.argv[1]))
for e in E:
    t=open(root+e['file'],encoding='utf-8').read()
    n=t.count(e['old_text'])
    if n:
        pos=t.index(e['old_text']); line=t[:pos].count('\n')+1
        endline=line+e['old_text'].count('\n')
    else:
        line=endline=None
        # whitespace-insensitive search
        import re
        norm=lambda s:re.sub(r'\s+',' ',s)
        tn=norm(t); on=norm(e['old_text'])
        ws=tn.count(on)
    print(e['fid'],e['file'],'stated',e['line'],'found',n,'at',line,'-',endline, '' if n else f'ws_insensitive={ws}')
