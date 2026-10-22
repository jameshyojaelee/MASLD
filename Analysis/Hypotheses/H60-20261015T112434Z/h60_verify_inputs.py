# H60 step 1: recompute every input hash listed in the prespec. Exit 1 on any mismatch.
import hashlib, os, sys
RUN = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H60-20261015T112434Z"
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
PRESPEC_SHA = "420dbb486398504d2c9bcbdb59f4d4790784fa3edafa79c453942e1ae17bc8fb"
# hash of the 22 per-file hashes, newline-joined with trailing newline, chr1..chr22 (prespec Inputs table)
HOH = {"ukbb_alt": "e32b043b29693164cf32dc17cb79f6e974c6a0f4a0f77f60b11be98ee11a0e59",
       "ghodsian": "28181aef0ab2b0e999b99135a522f6a3af9a41fa8c4c4f15ea98d319ee0c8314",
       "ukbb_ast": "24357e7a47236b7ac6c0557738a94b58b89b7caf612098cdf365ef323dac256b",
       "ukbb_ggt": "a3f095583ea2e6d56c3fa722321458486cb5d2fb8d281b22567a88bb0854120d"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main(out):
    bad = []
    lines = []
    s = sha(os.path.join(RUN, "PRESPEC_H60.md"))
    lines.append(f"PRESPEC_H60.md\t{s}\t{'OK' if s == PRESPEC_SHA else 'MISMATCH'}")
    if s != PRESPEC_SHA:
        bad.append("PRESPEC")
    for lst in ["inputs_sha256.txt", "otters_p005_per_file_sha256.txt"]:
        for l in open(os.path.join(RUN, lst)):
            if not l.strip():
                continue
            exp, rel = l.split()
            got = sha(ROOT + rel)
            ok = got == exp
            lines.append(f"{rel}\t{got}\t{'OK' if ok else 'MISMATCH'}")
            if not ok:
                bad.append(rel)
    for t, exp in HOH.items():
        hs = [sha(ROOT + f"RNA-seq/results/causal_inference/otters_broadaway/{t}/chr{c}/P0.05.txt") for c in range(1, 23)]
        got = hashlib.sha256(("\n".join(hs) + "\n").encode()).hexdigest()
        ok = got == exp
        lines.append(f"OTTERS_P0.05_{t}_hash_of_hashes\t{got}\t{'OK' if ok else 'MISMATCH'}")
        if not ok:
            bad.append(t)
    # D9 input (P0.001 ALT) is not hashed in the prespec; record its hashes for provenance only.
    hs = [sha(ROOT + f"RNA-seq/results/causal_inference/otters_broadaway/ukbb_alt/chr{c}/P0.001.txt") for c in range(1, 23)]
    hoh = hashlib.sha256(("\n".join(hs) + "\n").encode()).hexdigest()
    lines.append(f"OTTERS_P0.001_ukbb_alt_hash_of_hashes_RECORDED_ONLY\t{hoh}\tRECORDED")
    open(os.path.join(out, "input_hashes_checked.tsv"), "w").write("item\tsha256\tstatus\n" + "\n".join(lines) + "\n")
    print(f"input hash check: {len(lines)} items, {len(bad)} mismatches")
    if bad:
        print("MISMATCH:", bad)
        sys.exit(1)
    print("INPUT_HASHES_OK")


if __name__ == "__main__":
    main(sys.argv[1])
