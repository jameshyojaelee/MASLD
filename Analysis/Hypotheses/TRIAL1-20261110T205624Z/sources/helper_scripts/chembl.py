import json,sys,urllib.request,urllib.parse,time
B='https://www.ebi.ac.uk/chembl/api/data/'
def g(u): 
    time.sleep(0.2); return json.load(urllib.request.urlopen(B+u))
for name in sys.argv[1].split(','):
    try:
        r=g('molecule.json?pref_name__iexact='+urllib.parse.quote(name))['molecules'] or g('molecule/search.json?limit=5&q='+urllib.parse.quote(name))['molecules']
        r=[m for m in r if (m.get('pref_name') or '').lower()==name.lower()] or r[:1]
        if not r: print(name,'| NOT FOUND'); continue
        m=r[0]; cid=m['molecule_chembl_id']
        mech=g('mechanism.json?molecule_chembl_id='+cid)['mechanisms']
        if not mech:
            par=m.get('molecule_hierarchy',{}) or {}
            pid=par.get('parent_chembl_id')
            if pid and pid!=cid: mech=g('mechanism.json?molecule_chembl_id='+pid)['mechanisms']
        out=[]
        for x in mech:
            tg=x.get('target_chembl_id'); syms=[]
            if tg:
                t=g('target/%s.json'%tg)
                for c in t.get('target_components',[]):
                    for s in c.get('target_component_synonyms',[]):
                        if s['syn_type']=='GENE_SYMBOL': syms.append(s['component_synonym'])
            out.append('%s [%s] %s -> %s'%(x['mechanism_of_action'],x['action_type'],tg,'/'.join(syms)))
        print(name,'|',cid,'| phase',m.get('max_phase'),'|','; '.join(out) or 'NO MECHANISM')
    except Exception as e: print(name,'| ERR',e)
