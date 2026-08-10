# Firewall fixtures and tests

`test_firewall.py` is intentionally dependency-free and reads only blind-phase
candidate artifacts. Two synthetic CSV fixtures exercise exact HLF-versus-HLF-a
identity resolution, default-profile extraction, Entrez-aware gene headers,
duplicate-symbol failure, and explicit blanking of the pinned MEF2B/RLN2-style
cross-matrix Entrez conflicts. A dedicated availability fixture verifies that low-but-observed
TPM is covariate-complete while remaining ineligible under the frozen TPM>=1
floor, preventing missingness and biological detectability from being
conflated. The suite also tests the exact masked-column allow-list, absence of
outcome export/read markers, authenticated source hash,
frozen program weights and directions, known-hit exclusions, fail-closed
DepMap behavior, fixed draw and magnitude rules, and read-only raw-source
permissions.

It never opens a MAGeCK outcome value column. Run with:

```bash
python3 -B Analysis/Multimodal_Program_Projection/scripts/myojin_hlf/tests/test_firewall.py
```
