import pdfplumber
p='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/main/fig6_gene_catalog/panels/fig6d_next_experiment.pdf'
with pdfplumber.open(p) as pdf:
    pg=pdf.pages[0]
    for r in sorted(pg.rects,key=lambda r:(r['top'],r['x0'])):
        print(round(r['x0'],3),round(r['x1'],3),round(r['top'],2),round(r['bottom'],2),r.get('non_stroking_color'))
    ticks=[w for w in pg.extract_words() if w['text'] in ('0','5,000','10,000','15,000')]
    for w in ticks: print('tick',w['text'],round((w['x0']+w['x1'])/2,3),round(w['top'],1))
