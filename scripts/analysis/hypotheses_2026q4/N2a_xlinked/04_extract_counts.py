"""Count REF/ALT reads at the N2a sites, split by read-name parity, for every library of an included pair.
Usage: python -I 04_extract_counts.py RUN_ROOT N_PROCS
Reads RUN_ROOT/inputs/n2a_pairs.tsv and RUN_ROOT/inputs/sites.tsv.
Writes RUN_ROOT/counts/<library>.counts.tsv.gz (columns: chrom pos ref_p0 alt_p0 ref_p1 alt_p1 other_p0 other_p1),
RUN_ROOT/bai/<library>.bam.bai for BAMs that have no index (never written next to the source BAM),
and RUN_ROOT/counts/extraction_qc.tsv (depth totals and a MAPQ check on chr21; no allele ratios).
Libraries with a finished counts file are skipped, so a rerun never overwrites output."""
import sys, os, gzip, csv, time, zlib, collections
import multiprocessing as mp
import pysam

root, nproc = sys.argv[1], int(sys.argv[2])
os.makedirs(os.path.join(root, "counts"), exist_ok=True)
os.makedirs(os.path.join(root, "bai"), exist_ok=True)

pairs = list(csv.DictReader(open(os.path.join(root, "inputs", "n2a_pairs.tsv")), delimiter="\t"))
libs = {}
for r in pairs:
    if r["include"] == "True":
        libs[r["lib1"]] = r["bam1"]
        libs[r["lib2"]] = r["bam2"]
sites = []
with open(os.path.join(root, "inputs", "sites.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sites.append((r["chrom"], int(r["pos"]), r["ref"], r["alt"], r["cls"]))
sites.sort(key=lambda s: (s[0], s[1]))


def run_lib(item):
    lib, bam_path = item
    out = os.path.join(root, "counts", lib + ".counts.tsv.gz")
    if os.path.exists(out):
        return dict(library=lib, status="skipped_existing")
    t0 = time.time()
    f = pysam.AlignmentFile(bam_path, "rb")
    index_filename = None
    if not f.has_index():
        f.close()
        bai = os.path.join(root, "bai", lib + ".bam.bai")
        if not os.path.exists(bai):
            tmp = bai + ".tmp"
            pysam.index("-@", "2", bam_path, tmp)
            os.rename(tmp, bai)
        index_filename = bai
        f = pysam.AlignmentFile(bam_path, "rb", index_filename=index_filename)
    t_index = time.time() - t0
    # MAPQ check on a chr21 window (autosome, no X information)
    mq = collections.Counter()
    for k, rd in enumerate(f.fetch("chr21", 20000000, 20100000)):
        mq["mapq255" if rd.mapping_quality == 255 else "other"] += 1
        if k >= 20000:
            break
    depth = collections.Counter(); nrows = collections.Counter()
    tmp_out = out + ".tmp"
    with gzip.open(tmp_out, "wt") as fh:
        fh.write("chrom\tpos\tref_p0\talt_p0\tref_p1\talt_p1\tother_p0\tother_p1\n")
        for chrom, pos, ref, alt, cls in sites:
            c = collections.Counter()
            for col in f.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=20, min_mapping_quality=255,
                                max_depth=1000000, stepper="samtools", ignore_overlaps=True, ignore_orphans=True):
                if col.reference_pos != pos - 1:
                    continue
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip or pr.query_position is None:
                        continue
                    b = pr.alignment.query_sequence[pr.query_position].upper()
                    par = zlib.crc32(pr.alignment.query_name.encode()) & 1
                    c[(ref_alt_other(b, ref, alt), par)] += 1
            tot = sum(c.values())
            if tot == 0:
                continue
            fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (chrom, pos, c[("r", 0)], c[("a", 0)], c[("r", 1)], c[("a", 1)], c[("o", 0)], c[("o", 1)]))
            key = "X" if chrom == "chrX" else "auto"
            depth[key] += c[("r", 0)] + c[("a", 0)] + c[("r", 1)] + c[("a", 1)]
            nrows[key] += 1
    os.rename(tmp_out, out)
    f.close()
    return dict(library=lib, status="done", bam=bam_path, indexed_here=index_filename is not None, seconds_index=round(t_index, 1),
                seconds_total=round(time.time() - t0, 1), mapq255=mq["mapq255"], mapq_other=mq["other"],
                sites_with_reads_X=nrows["X"], sites_with_reads_auto=nrows["auto"], total_depth_X=depth["X"], total_depth_auto=depth["auto"])


def ref_alt_other(b, ref, alt):
    return "r" if b == ref else ("a" if b == alt else "o")


if __name__ == "__main__":
    print("libraries:", len(libs), "sites:", len(sites), "procs:", nproc, flush=True)
    res = []
    with mp.Pool(nproc) as pool:
        for r in pool.imap_unordered(run_lib, sorted(libs.items())):
            res.append(r)
            print(r["library"], r["status"], r.get("seconds_total", ""), flush=True)
    qc = os.path.join(root, "counts", "extraction_qc.tsv")
    cols = ["library", "status", "bam", "indexed_here", "seconds_index", "seconds_total", "mapq255", "mapq_other",
            "sites_with_reads_X", "sites_with_reads_auto", "total_depth_X", "total_depth_auto"]
    mode = "w" if not os.path.exists(qc) else "a"
    with open(qc, mode) as fh:
        if mode == "w":
            fh.write("\t".join(cols) + "\n")
        for r in res:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print("done", sum(r["status"] == "done" for r in res), "skipped", sum(r["status"] != "done" for r in res))
