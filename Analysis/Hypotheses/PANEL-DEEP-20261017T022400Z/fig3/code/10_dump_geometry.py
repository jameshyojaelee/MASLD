"""Dump vector marks (chars, rects, lines, curves) of each Figure 3 panel PDF with colours, for geometry checks.
Usage: python 10_dump_geometry.py <pdf> <outprefix>"""
import pdfplumber, sys, json, pandas as pd
pdf,outp=sys.argv[1],sys.argv[2]
with pdfplumber.open(pdf) as p:
    pg=p.pages[0]; H=pg.height
    rows=[]
    for kind in ['char','rect','line','curve']:
        for o in getattr(pg,kind+'s'):
            r=dict(kind=kind,x0=o['x0'],x1=o['x1'],top=o['top'],bottom=o['bottom'],
                   stroke=str(o.get('stroking_color')),fill=str(o.get('non_stroking_color')),
                   lw=o.get('linewidth'),text=o.get('text',''),font=o.get('fontname',''),size=o.get('size',None),
                   npts=len(o.get('pts',[]) or []), fillflag=o.get('fill'), strokeflag=o.get('stroke'))
            if kind=='curve': r['pts']=json.dumps([(round(a,3),round(b,3)) for a,b in o['pts']])
            rows.append(r)
    df=pd.DataFrame(rows); df.to_csv(outp+"_objects.tsv",sep="\t",index=False)
    print(pdf, "page", pg.width, pg.height, df.kind.value_counts().to_dict())
    # words
    w=pg.extract_words(keep_blank_chars=False,extra_attrs=['size','fontname'])
    pd.DataFrame(w).to_csv(outp+"_words.tsv",sep="\t",index=False)
    print("fonts:",df[df.kind=='char'].font.value_counts().to_dict(), "sizes:",df[df.kind=='char']['size'].round(2).value_counts().to_dict())
