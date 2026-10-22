import json, re
RUN = '.'
s = open('scratch/main-rederive/synthesis_draft.md').read()
sec2bc = open('scratch/main-rederive/section_2b_2c.md').read()
d = json.load(open('workflow_result.json'))
v2 = {r['fid']: r['verdict'] for r in json.load(open('verify2_result.json'))['results']}

def rep(old, new, count=1):
    global s
    assert s.count(old) == count, (old[:80], s.count(old))
    s = s.replace(old, new)

# 1. Summary
rep("Sixteen review dimensions (R01 to R16) and three focused gap finders (G1 to G3) together raised 175 findings. Twenty were sent to independent verification, and all 20 were confirmed. None was refuted. The three planted decoys were all classified correctly (one true, two false). The other 112 findings remain unverified.",
    "Sixteen review dimensions (R01 to R16) and three focused gap finders (G1 to G3) together raised 175 findings. 43 were merged as duplicates, which leaves 132. 40 of the 132 were verified by blind agents in two passes. The main pass covered the 20 ranked first, and all 20 were confirmed. A second pass covered the 20 remaining high-severity findings: 13 were confirmed, 7 partly confirmed, and none refuted (Section 2B). All five planted decoys were classified correctly (one true, four false). The other 92 findings (62 medium, 30 low) remain unverified. I re-derived three confirmed findings myself and got the same values (Section 2C).")
rep("- **Claims with no analysis behind them.** Title clause 2 (\"the coordination of multicellular remodeling\") has no Results section, legend or panel. The Catalog has no regulatory DNA object, although the abstract says it links them.\n",
    "- **Claims with no analysis behind them.** Title clause 2 (\"the coordination of multicellular remodeling\") has no Results section, legend or panel. The Catalog has no regulatory DNA object, although the abstract says it links them.\n\nThe second pass adds three more:\n\n- **Results that depend on one cohort or on one analysis choice.** Without GSE130970 the Fig 3E F1 bar drops from 1,338 calls to 14. 29 of the 54 Fig 3F genes depend on the same cohort. Fig 4E's 16 fibrosis calls are 0 on the 774 units alone, and 0 under truncated Knapp-Hartung.\n- **A Catalog label that measures less than it says.** \"Genetics\" is liver-eQTL colocalization only, so MTARC1, PNPLA3, TM6SF2 and HSD17B13 read as unresolved.\n- **Wording beyond the numbers.** This includes \"resolves the architecture\", a \"direct MASLD\" coding enrichment that mostly reflects the three PDFF GWAS (NAFLD and NASH alone: ratio 3.21 [0.53, 8.25]), and four uncited prior papers.\n")
rep("The round proposes ten analyses, together 608.5 CPU hours by the plans' own estimates:",
    "The round proposes ten analyses. By the plans' own estimates they need about 252 CPU hours of expected use, and H51 alone could reach about 400 if single-thread fits are slower:")

# 2. Insert 2B/2C before section 3
rep("## 3. Refuted, unresolved and unverified findings\n", sec2bc.rstrip('\n') + "\n\n## 3. Refuted, unresolved and unverified findings\n")

# 3. Refuted list
rep("- Decoy X-03 (planted false): refuted, correct.\n",
    "- Decoy X-03 (planted false): refuted, correct.\n- Second-pass decoy G4-01 (planted false: the Fig 3E F3 bar is 6,910, not 6,190): refuted, correct. The source table gives 3,144 + 3,046 = 6,190.\n- Second-pass decoy G4-02 (planted false: r8 has 426 genetically anchored genes): refuted, correct. r8 has 400 genetically anchored and 26 concordant.\n")

# 4. Section 3.3 rewrite
i = s.index("### 3.3 Unverified findings to verify first")
j = s.index("### 3.4 Full list of unverified findings")
new33 = ("### 3.3 Unverified findings\n\n"
"92 findings remain unverified: 62 medium and 30 low. By status, 73 are new, 16 are known and still open, and 3 are known and resolved. One finder reported each one with commands, but no second agent re-derived it. Do not act on them before verification. The 43 raw findings not counted here were merged into others as duplicates (`RUN/review.json`, `merged`).\n\n"
"The 20 high-severity findings that this section listed before the second pass are now in Section 2B.\n\n"
"A further number from the portfolio critic was not verified separately. `fig5b_protein_triage_summary.tsv` gives n_full_measured 5,836, rho_full 0.34582 and direction_full 0.61361. Fig 5B prints rho 0.46 and 66% agreement on a 3,054-gene subset chosen from the superseded July COLOC universe. LEG:388 discloses the subset choice, but the full-set values are not printed.\n\n")
s = s[:i] + new33 + s[j:]
rep("### 3.4 Full list of unverified findings", "### 3.4 Full list of findings not verified in the main pass")
rep("Status key: new; known, open; known, resolved; in flight.",
    "Status key: new; known, open; known, resolved; in flight. For the 20 high-severity rows, the status cell also gives the second-pass verdict (Section 2B).")
n = 0
for fid, v in v2.items():
    if fid.startswith('G4'):
        continue
    pat = re.compile(r'^(\| ' + re.escape(fid) + r' \| high \| )([^|]*?)( \|)', re.M)
    s, k = pat.subn(lambda m: m.group(1) + m.group(2) + '; second pass: ' + ('confirmed' if v['verdict'] == 'confirmed' else 'partly confirmed') + m.group(3), s)
    assert k == 1, (fid, k)
    n += k
assert n == 20

# 5. Plans
rep("Each audited prespecification is in `RUN/scratch/audit-<ID>/` (`PRESPEC_<ID>_patched.md` or `PRESPEC_<ID>_audited.md`).",
    "Each audited prespecification is in `RUN/plans/PRESPEC_<ID>.md`. The audit working copies are in `RUN/scratch/audit-<ID>/`. CPU h is expected use, not the reservation: every plan requests the 90 h default wall time.")
rep("| 360 | Passed", "| about 4 (360 reserved) | Passed")
rep("| 116 | Passed", "| 116 (up to about 400) | Passed")
rep("Total: 608.5 CPU hours.",
    "Total expected use: about 252 CPU hours (116 + 6 + 40 + 8 + 16 + 48 + 0.5 + 4 + 8 + 6). An earlier total of 608.5 counted H58's 4 CPU x 90 h reservation as use. The H58 prespecification expects under 1 hour.")
rep("- **H58** (corrective, 360 CPU h).", "- **H58** (corrective, about 4 CPU h expected; 360 reserved).")

# 6. Section 5
i = s.index("**Pipeline counts.**"); j = s.index("**Failed the gate or validity check:**")
ng = d['hyp']['not_gated']
pc = ("**Pipeline counts.** 57 candidates were generated:\n\n"
"- 50 from lenses L01 to L09.\n- 4 from the referee lens L10, written after review findings were confirmed.\n- 3 from the portfolio critic (C-01 to C-03).\n\n"
"The code check dropped none. Triage of the 50 lens candidates kept 35, merged 14 into others and dropped 1. The dropped one is L02-05: it asks the same question as an archive candidate and has the H40 failure mode, where the estimand is nearly collinear with the nuisance. The gate cap admitted 14 of the 35, so 21 were never gated (table below). Of the 18 gated lens candidates (14 from L01 to L09 and 4 from L10), 7 passed: L01-01, L02-01, L03-04, L07-01, L10-01, L10-03 and L10-04. Six of the seven became plans H51 to H56. Of the critic's three, C-01 passed and became H60.\n\n"
"**Passed the gate, killed by the judges:** L01-01 (GGT-only vs ALT/AST-only colocalized genes in biliary vs hepatocyte compartments). All three judges scored it 2 and killed it. The MDE is 0.43 SD, but the whole 289-locus COLOC set departs from background by only 0.15 SD. GGT1 itself has an index of 0.07, and the class labels partly record which biobank measured which enzyme. It was the only candidate that tested the inherited-to-remodeling link directly.\n\n")
s = s[:i] + pc + s[j:]
rep("**Critic candidates.** The portfolio critic proposed three further candidates, gated as C-01 to C-03. C-01 (TWAS) became H60. Going by their scratch directories, C-02 concerned HMSMA spatial pairs and C-03 Catalog annotation. The input does not record their outcome.",
    "**Critic candidates.** C-01 (TWAS-predicted expression against disease direction) passed and became H60. C-02 (are bulk-coordinated program pairs co-patterned within HMSMA livers?) failed because its test is invalid for its claim: under random relabeling of the 662 pairs, the donor sign-flip test rejects at p<0.05 in 40.8% of 4,000 labelings. C-03 (Catalog per-gene columns for risk-allele direction, stage effects and liver protein) failed its concordance test on power: the MDE is 0.65 at 83 to 89 non-MHC loci, against a disclosed agreement of 0.60. Its column-building part needs no test, and it overlaps gap 6.1(1).\n\n"
"**Never gated (21).** These passed triage but were beyond the gate cap of 14. Nothing here says whether they are feasible. A follow-up gating round would cost about 21 gate agents plus validity and prior-art checks for the survivors.\n\n"
"| ID | Title |\n|---|---|\n" + "".join(f"| {c['tmp_id']} | {c['title']} |\n" for c in ng))

# 7. Section 7
rep("- Three decoys were mixed in, each checked by both agent types: X-01 (true) was confirmed, and X-02 and X-03 (false) were refuted. That is 3 of 3 correct. With three decoys this shows the verifiers can refute a false claim. It does not estimate an error rate.",
    "- Three decoys were mixed in, each checked by both agent types: X-01 (true) was confirmed, and X-02 and X-03 (false) were refuted. That is 3 of 3 correct.\n- Second pass: the 20 remaining high-severity findings each got one blind reproduce verifier. Each verifier read only its own claim (title, claim, locations) from `RUN/context/verify2_claims_blind.json`, not the finder's evidence or severity. Two false decoys (G4-01, G4-02) were given ordinary-looking IDs, and both were refuted. Over both passes, decoys were 5 of 5 correct. With five decoys this shows the verifiers can refute a false claim. It does not estimate an error rate, and single-verifier verdicts were not downgraded.\n- Main-loop re-derivation of R01-01, R03-03 and R05-02 (Hotspot part) with new code: all three match (Section 2C).")
i = s.index("**Jobs.** sbatch jobs were capped at 8 minutes"); j = s.index("`integrated_atlas.h5ad` was not loaded.")
jobs = ("**Jobs.** sbatch jobs were capped at 8 minutes on `cpu` with at most 4 CPUs and 16 GB. Each had a 4-digit numeric name, and all are listed in `RUN/job_index.tsv` (42 jobs: 32 main pass, 10 second pass). `sacct` on 2026-10-15 gave 36 COMPLETED, 3 TIMEOUT and 3 FAILED, with none still running. The non-completed jobs:\n\n"
"- 1031: TIMEOUT after the pooled models. Its per-cohort GSE126848 part was rerun on the login node (n=30, 13 s).\n"
"- 3081: TIMEOUT. The L03-04 positive-control refit printed no fit, and the gate used other evidence.\n"
"- 3831: TIMEOUT after the H51 grouped discovery refit. The permutation and SE steps did not run.\n"
"- 3051: FAILED after 59 s. The L01-04 detection and testability counts finished, but the pilot link model did not.\n"
"- 3101: FAILED after 3 s, when `fread` read the GSE225740 header as data. The L03-05 gate used login-node counts.\n"
"- 3841: FAILED after 2 min 39 s in the L10-04 caQTL pilot.\n\n"
"Other notes from the main pass:\n\n"
"- Job 1161 completed but keyed its join on the wrong GWAS, and a login-node recomputation replaced it.\n"
"- Job 1032 skipped F1 leave-one-cohort-out, the F2/F3 interaction and the complexity covariates. Second-pass job 4191 supplied the F1 leave-one-cohort-out values (Section 2B, G2-01).\n"
"- Job 3111 (H55 morphometry pilot) finished 118 of 222 slides.\n\n"
"**Agents.** 134 in the main workflow (`wf_a8d2ba5f-3d1`) and 22 in the second pass (`wf_00cf4a6c-d4e`). Two more started in a second-pass launch that was stopped within a minute because a placeholder was passed instead of the claims; they submitted no job and wrote no file. Total: 158, under the planned maximum of 160.\n\n"
"**Edits outside RUN.** One workflow agent imported a repo script in place inside a figure-candidate snapshot, which left two `.pyc` files and an empty `__pycache__` directory under `figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/scripts/figures/`. I checked that they were created after the run started, then removed them; the second pass prompts forbade in-place imports. After both passes, the sha256 manifest of `docs/`, `figures/main/` and `docs/manuscript/` (3,441 files) matched the start of the run, and `git status` was identical to the start.\n\n")
s = s[:i] + jobs + s[j:]
rep("This synthesis agent ran only directory listings and read the brief, `known_issues.tsv` and the 2026-10-07 combined review. Every other number above comes from the finder, verifier, gate, audit and critic outputs passed to it.",
    "The synthesis agent wrote Sections 1 to 7 from the finder, verifier, gate, audit and critic outputs. The main loop then added Sections 2B and 2C, the decoy and agent accounting, the job table and the never-gated list. It also corrected the candidate counts in Section 5 and the CPU total in Sections 1 and 4. Values in 2B come from `RUN/verify2_result.json`. Values in 2C come from `RUN/scratch/main-rederive/`.")

assert '—' not in s
open('scratch/main-rederive/final_doc.md', 'x').write(s)
print(len(s), s.count('\n'))
