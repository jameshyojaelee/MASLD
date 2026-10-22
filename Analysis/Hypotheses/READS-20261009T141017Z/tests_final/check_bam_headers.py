import csv, hashlib, sys
import pysam
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
want = {int(i) for i in sys.argv[2:]}
print("index\tcohort\tlibrary\tn_SQ\tsq_sha256_12\thas_KI270847.1\tPG_STAR_version\tGTF_or_index_hint")
for r in rows:
    if int(r["array_index"]) not in want:
        continue
    af = pysam.AlignmentFile(r["bam_path"], "rb", check_sq=False)
    h = af.header.to_dict()
    sq = [(x["SN"], x["LN"]) for x in h.get("SQ", [])]
    digest = hashlib.sha256(repr(sq).encode()).hexdigest()[:12]
    pg = h.get("PG", [{}])[0]
    co = " | ".join(c[:60] for c in h.get("CO", [])[:1])
    print(f"{r['array_index']}\t{r['cohort']}\t{r['library']}\t{len(sq)}\t{digest}\t{'KI270847.1' in dict(sq)}\t{pg.get('VN','NA')}\t{co[:80]}")
    af.close()
