import pdfplumber, collections
R='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/main/fig6_gene_catalog/panels/'
for p in ['fig6a_catalog_structure','fig6b_coverage_observability','fig6c_worked_entries','fig6d_next_experiment','fig6e_portal_source_access']:
    with pdfplumber.open(R+p+'.pdf') as pdf:
        pg=pdf.pages[0]
        c=collections.Counter((ch['fontname'],round(ch['size'],2),str(ch.get('non_stroking_color'))) for ch in pg.chars)
        fills=collections.Counter(str(r.get('non_stroking_color')) for r in pg.rects)
        strokes=collections.Counter(str(l.get('stroking_color')) for l in pg.lines+pg.rects+pg.curves)
        print(p, 'size(pt)',round(pg.width,1),round(pg.height,1)); print('  chars',dict(c)); print('  rect fills',dict(fills)); print('  strokes',dict(strokes))
