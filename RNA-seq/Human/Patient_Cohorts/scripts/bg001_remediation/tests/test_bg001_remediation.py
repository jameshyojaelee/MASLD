#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


sources = load_module("sources", SCRIPT_DIR / "validate_featurecounts_sources.py")
prepare = load_module("prepare", SCRIPT_DIR / "prepare_run.py")
dependency_manifest = load_module("dependency_manifest", SCRIPT_DIR / "generate_dependency_manifest.py")
figure_labels = load_module(
    "validate_rendered_figure_labels",
    SCRIPT_DIR / "validate_rendered_figure_labels.py",
)
analysis_contract = load_module(
    "verify_analysis_contract_test", SCRIPT_DIR / "verify_analysis_contract.py"
)
artifact_seal = load_module(
    "seal_analysis_artifacts_test", SCRIPT_DIR / "seal_analysis_artifacts.py"
)
safe_io = load_module("safe_io_test", SCRIPT_DIR / "safe_io.py")
analysis_invalid = load_module(
    "seal_analysis_invalid_test", SCRIPT_DIR / "seal_analysis_invalid.py"
)


class ImmutableAnalysisContractTests(unittest.TestCase):
    def valid_contract(self, root: Path):
        return {
            "schema_version": "1.2",
            "run_id": root.name,
            "candidate_root": str(root),
            "reference_root": str(root / "source_snapshot"),
            "thresholds": dict(analysis_contract.EXPECTED_THRESHOLDS),
            "arms": {arm: str(root / "arms" / arm) for arm in analysis_contract.EXPECTED_ARMS},
            "arm_modes": {
                arm: dict(mode) for arm, mode in analysis_contract.EXPECTED_ARM_MODES.items()
            },
            "analysis_model": dict(analysis_contract.EXPECTED_ANALYSIS_MODEL),
            "active_affected_cohorts": list(analysis_contract.EXPECTED_ACTIVE_COHORTS),
            "canonical_affected_cohorts": list(analysis_contract.EXPECTED_CANONICAL_COHORTS),
            "featurecounts": dict(analysis_contract.EXPECTED_FEATURECOUNTS),
            "bam_manifest_finalization_required": True,
        }

    def test_threshold_arm_mode_and_model_mutations_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "bg001-fragment-v211-gencode49-20260807T000000Z"
            contract = self.valid_contract(root)
            analysis_contract.validate_analysis_contract(contract, root)
            mutations = (
                ("thresholds", "gene_jaccard_min", 0.98),
                ("arm_modes", "F_five", {**contract["arm_modes"]["F_five"], "filter_scope": "legacy_all"}),
                ("analysis_model", "treat_lfc", 0.3),
            )
            for section, key, value in mutations:
                mutated = json.loads(json.dumps(contract))
                mutated[section][key] = value
                with self.subTest(section=section, key=key):
                    with self.assertRaisesRegex(SystemExit, "differs from frozen specification"):
                        analysis_contract.validate_analysis_contract(mutated, root)


class AnalysisArtifactSealTests(unittest.TestCase):
    def make_args(self, root: Path):
        artifact_root = root / "arms/R0"
        artifact_root.mkdir(parents=True)
        (artifact_root / "provenance").mkdir()
        (artifact_root / "result.tsv").write_text("gene\tvalue\nG1\t1\n")
        bound = root / "contract/run_contract.json"
        bound.parent.mkdir(parents=True)
        bound.write_text("{}\n")
        return SimpleNamespace(
            run_root=root,
            artifact_root=artifact_root,
            manifest=artifact_root / "provenance/artifact_manifest.tsv",
            marker=artifact_root / "ARM_COMPLETE.json",
            label="R0",
            bound_input=[bound],
        )

    def test_exact_seal_rejects_content_extra_symlink_and_bound_input_mutations(self):
        mutations = ("content", "extra", "symlink", "bound")
        for mutation in mutations:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp).resolve()
                args = self.make_args(root)
                artifact_seal.seal(args)
                artifact_seal.verify(args)
                if mutation == "content":
                    (args.artifact_root / "result.tsv").write_text("mutated\n")
                elif mutation == "extra":
                    (args.artifact_root / "extra.tsv").write_text("extra\n")
                elif mutation == "symlink":
                    (args.artifact_root / "linked.tsv").symlink_to(args.bound_input[0])
                else:
                    args.bound_input[0].write_text('{"mutated":true}\n')
                with self.assertRaises((SystemExit, RuntimeError)):
                    artifact_seal.verify(args)

    def test_seal_rejects_temporary_and_dangling_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            args = self.make_args(root)
            (args.artifact_root / ".partial.tmp").write_text("partial\n")
            with self.assertRaisesRegex(SystemExit, "temporary"):
                artifact_seal.seal(args)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            target = root / "out.tsv"
            target.symlink_to(root / "missing")
            with self.assertRaises(safe_io.SafeIOError):
                safe_io.publish_new_bytes(target, b"no overwrite\n", root)


class AnalysisInvalidStateTests(unittest.TestCase):
    def test_failed_r0_is_exclusively_sealed_and_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            for relative in analysis_invalid.COMMON_BOUND + analysis_invalid.R0_BOUND:
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("bound\n")
            report = root / "comparisons/R0_reproduction.json"
            report.write_text(json.dumps({"status": "FAIL", "failures": ["delta"]}) + "\n")
            with mock.patch("sys.argv", ["seal_analysis_invalid.py", "seal", "--run-root", str(root), "--kind", "r0"]):
                analysis_invalid.main()
            marker = json.loads((root / "ANALYSIS_INVALID.json").read_text())
            self.assertEqual(marker["status"], "INVALID")
            self.assertEqual(marker["hard_gate"], "R0_REPRODUCTION")
            self.assertEqual(set(marker["bound_inputs"]), set(analysis_invalid.COMMON_BOUND + analysis_invalid.R0_BOUND))
            with mock.patch("sys.argv", ["seal_analysis_invalid.py", "verify", "--run-root", str(root), "--kind", "r0"]):
                analysis_invalid.main()
            (root / "arms/R0/results/integration/deg_results.csv").write_text("drift\n")
            with mock.patch("sys.argv", ["seal_analysis_invalid.py", "verify", "--run-root", str(root), "--kind", "r0"]):
                with self.assertRaisesRegex(SystemExit, "drifted"):
                    analysis_invalid.main()

    def test_runner_stops_at_failed_r0_before_corrected_arms(self):
        runner = (SCRIPT_DIR / "run_four_arms.sh").read_text()
        invalid_call = 'seal_analysis_invalid.py" seal'
        self.assertIn(invalid_call, runner)
        self.assertLess(runner.index(invalid_call), runner.index("run_qc F_locked"))
        self.assertIn("ANALYSIS_INVALID.json", runner)


class ScientificGateMutationTests(unittest.TestCase):
    def test_exact_locked_equality_finite_and_expression_sensitivity_helpers(self):
        expression = f'''
          source({json.dumps(str(SCRIPT_DIR / "scientific_gate_helpers.R"))})
          library(data.table)
          a <- data.table(sample_id=c("S1","S2"), x=c(1,2), y=c(0.1,0.2))
          b <- copy(a); b[2, y := y + 1e-15]
          stopifnot(bg001_table_equal_by_key(a, a, "sample_id"))
          stopifnot(!bg001_table_equal_by_key(a, b, "sample_id"))
          meta <- data.table(sample_id=c("S1","S2"), inferred_sex=c("Female",""))
          stopifnot(bg001_missing_model_sex(meta, c("S1","S2")) == 1L)
          coef <- data.table(logFC=c(1,Inf), SE=c(0.1,0.2))
          stopifnot(bg001_nonfinite_counts(coef, c("logFC","SE"))[["logFC"]] == 1L)
          model <- data.table(delta_logFC=c(0.001,0.2), AveExpr_old=c(2,2), AveExpr_new=c(2,0))
          primary <- bg001_expression_gate_summary(model, pmin(model$AveExpr_old, model$AveExpr_new)>1, .01, .05)
          old <- bg001_expression_gate_summary(model, model$AveExpr_old>1, .01, .05)
          stopifnot(primary$overall_pass, !old$overall_pass)
          direction <- data.table(
            protected=c(TRUE, FALSE),
            claim_protected=c(FALSE, FALSE),
            logFC_old=c(0.5, 0.5),
            logFC_new=c(-0.5, -0.5)
          )
          stopifnot(bg001_protected_direction_reversals(direction) == 1L)
          expected_jobs <- c(
            G1="qc_sex_reconstruction_and_refit",
            G2="dependency_full_raw_dge_rebuild",
            G3="normalization_weight_localization",
            G4="limma_leave_one_cohort_out",
            G5="inspect_treat_threshold_and_dependencies"
          )
          observed_jobs <- vapply(names(expected_jobs), function(gate) {{
            bg001_gate_failure_specification(gate)$job
          }}, character(1))
          stopifnot(identical(unname(observed_jobs), unname(expected_jobs)))
          cat("scientific helper mutations PASS\\n")
        '''
        result = subprocess.run(
            ["micromamba", "run", "-n", "rnaseq", "Rscript", "--vanilla", "-e", expression],
            text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("mutations PASS", result.stdout)

    def test_compare_and_arm_validator_encode_required_fail_closed_seams(self):
        comparison = (SCRIPT_DIR / "compare_scientific_gates.R").read_text()
        for token in (
            "QC_MEMBERSHIP_CHANGE", "PROTECTED_TREAT_MEMBERSHIP_CHANGE",
            "PROTECTED_DIRECTION_CHANGE", "G5_FAIL_G4_PASS",
            "escalation_disposition.tsv", "coefficient_gene_order_matches_dge",
            "invalid_treat_lfc", "comparison_sample_coverage.tsv",
        ):
            self.assertIn(token, comparison)
        self.assertIn('protected_reason != "canonical_treat"', comparison)
        self.assertIn("claim_protected == TRUE", comparison)
        self.assertIn("protected_flips <- bg001_protected_direction_reversals(model)", comparison)
        direction_block = comparison[
            comparison.index("protected_direction_changes <- gene_level["):
            comparison.index("if (nrow(protected_direction_changes))", comparison.index("protected_direction_changes <- gene_level["))
        ]
        self.assertIn("protected == TRUE", direction_block)
        self.assertNotIn("claim_protected == TRUE", direction_block)
        self.assertIn('overall_status <- if (nrow(failed))', comparison)
        self.assertIn('"PENDING_TARGETED_ESCALATION"', comparison)
        self.assertIn("claim_protected = logical()", comparison)
        self.assertIn("setorder(disposition, trigger_id", comparison)
        expression_call = comparison[
            comparison.index('add_metric("G4", comparison_id, "expression_stratum_verdict_invariant"'):
            comparison.index('for (stratum_name in', comparison.index('expression_stratum_verdict_invariant'))
        ]
        self.assertIn("gating = FALSE", expression_call)
        validator = (SCRIPT_DIR / "validate_arm.R").read_text()
        for token in (
            "DGE contains a nonfinite, negative, or noninteger count",
            "Filter-stat table does not exactly cover raw genes/kept DGE genes",
            "model.matrix(~ dataset + inferred_sex + group_binary",
            "Effective count source differs from override/path/hash contract",
        ):
            self.assertIn(token, validator)
        fit = (
            SCRIPT_DIR.parents[1]
            / "analysis/integration/scripts/05h_limma_voom_qw_canonical.R"
        ).read_text()
        self.assertNotIn("if (is.null(sw)) sw <- rep(1", fit)
        self.assertIn("validate_voom_curve(v$voom.xy", fit)
        self.assertIn("any(out$treat_lfc != 0.25)", fit)
        self.assertIn("fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)", fit)
        self.assertNotIn("abs(logFC / t)", fit)


class PinnedRunnerHarnessTests(unittest.TestCase):
    def test_env_wrapper_executes_absolute_rscript_with_environment_and_arguments(self):
        runner = (SCRIPT_DIR / "run_four_arms.sh").read_text()
        start = runner.index("run_pinned_r() {")
        end = runner.index("verify_recount() {")
        functions = runner[start:end]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stub = root / "Rscript-stub"
            stub.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$TEST_TOKEN" "$@"\n')
            stub.chmod(0o750)
            harness = root / "harness.sh"
            harness.write_text(
                "#!/usr/bin/env bash\nset -euo pipefail\n"
                f"R_SCRIPT={json.dumps(str(stub))}\n"
                "verify_runtime() { :; }\n"
                f"{functions}\n"
                "run_pinned_r_env TEST_TOKEN=bound -- --vanilla script.R arg1\n"
            )
            result = subprocess.run(["bash", str(harness)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.splitlines(), ["bound", "--vanilla", "script.R", "arg1"])


class SourceParserTests(unittest.TestCase):
    def test_multiline_fragment_command(self):
        text = """featureCounts \\
          -T 8 \\
          -p --countReadPairs -B \\
          -s 2 -a ref.gtf -o out.txt sample.bam
        """
        regions = sources.command_regions(text)
        self.assertEqual(len(regions), 1)
        for token in ("-p", "--countReadPairs", "-B", "-s"):
            self.assertTrue(sources.has_token(regions[0], token))

    def test_comment_does_not_satisfy_contract(self):
        text = "# featureCounts -p --countReadPairs -B -s 2\nfeatureCounts -p -B -s 2 x.bam\n"
        regions = sources.command_regions(text)
        self.assertEqual(len(regions), 1)
        self.assertFalse(sources.has_token(regions[0], "--countReadPairs"))

    def test_embedded_false_flag_is_not_token(self):
        region = "featureCounts -p --countReadPairs=false -B -s 2 x.bam"
        self.assertFalse(sources.has_token(region, "--countReadPairs"))


class RunIdTests(unittest.TestCase):
    def test_candidate_sentinel_is_plain_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            marker = root / ".bg001_candidate_root"
            marker.write_text("test\n")
            self.assertTrue(marker.is_file())
            self.assertFalse(marker.is_symlink())


class DependencyInventoryTests(unittest.TestCase):
    def test_dependency_acceptance_state_is_fail_closed(self):
        pending = dependency_manifest.classify_verdict_state({
            "overall_status": "PENDING_TARGETED_ESCALATION",
            "completion_status": "PENDING_TARGETED_ESCALATION",
            "escalation_required": True,
        })
        self.assertFalse(pending[3])
        self.assertTrue(pending[5])
        stable = dependency_manifest.classify_verdict_state({
            "overall_status": "STABLE",
            "completion_status": "PRIMARY_GATES_COMPLETE",
            "escalation_required": False,
        })
        self.assertTrue(stable[3])
        with self.assertRaisesRegex(SystemExit, "inconsistent"):
            dependency_manifest.classify_verdict_state({
                "overall_status": "STABLE",
                "completion_status": "PRIMARY_GATES_COMPLETE",
                "escalation_required": True,
            })

    def test_script_dir_wrapper_and_custom_loader_are_resolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp)
            scripts = project / "scripts"
            scripts.mkdir()
            wrapper = scripts / "run.sh"
            child = scripts / "consumer.R"
            wrapper.write_text('SCRIPT_DIR="' + str(scripts) + '"\nRscript "$SCRIPT_DIR/consumer.R"\n')
            child.write_text(
                'deg_path <- file.path("x", "canonical_deg_results.csv")\n'
                'result <- load_project_table(deg_path)\n'
            )
            lines = prepare.noncomment_lines(wrapper)
            targets = prepare.resolve_wrapper_targets(
                project,
                "scripts/run.sh",
                lines,
                {"consumer.R": {"scripts/consumer.R"}},
            )
            self.assertEqual(targets, {"scripts/consumer.R"})
            modes, _ = prepare.role_access_modes(
                prepare.noncomment_lines(child),
                "canonical_deg",
                prepare.DEPENDENCY_PATTERNS["canonical_deg"],
            )
            self.assertIn("read", modes)

    def test_write_only_artifact_is_not_a_consumer_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "producer.R"
            path.write_text('fwrite(result, "canonical_deg_results.csv")\n')
            modes, _ = prepare.role_access_modes(
                prepare.noncomment_lines(path),
                "canonical_deg",
                prepare.DEPENDENCY_PATTERNS["canonical_deg"],
            )
            self.assertEqual(modes, {"write"})

    def test_required_derived_chain_has_detectable_input_roles(self):
        project = SCRIPT_DIR.parents[4]
        self.assertEqual(len(dependency_manifest.REQUIRED_CANONICAL_CHAIN), 9)
        for relative in dependency_manifest.REQUIRED_CANONICAL_CHAIN:
            path = project / relative
            if not path.is_file():
                path = project / "dependency_sources" / relative
            self.assertTrue(path.is_file(), relative)
            lines = prepare.noncomment_lines(path)
            detected = set()
            for role, needles in prepare.DEPENDENCY_PATTERNS.items():
                modes, _ = prepare.role_access_modes(lines, role, needles)
                if modes and modes != {"write"}:
                    detected.add(role)
            self.assertTrue(detected, f"No dependency input role detected for {relative}")


class ArmRetryTests(unittest.TestCase):
    def test_nonempty_arm_is_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            (root / ".bg001_candidate_root").write_text(root.name + "\n")
            (root / "RECOUNT_COMPLETE").write_text("x\n")
            (root / "BASELINE_FROZEN.json").write_text("{}\n")
            required = [
                "manifests/count_overrides.tsv",
                "manifests/read_count_overrides.tsv",
                "source_snapshot/config/human_datasets.yaml",
                "source_snapshot/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv",
                "frozen_sets/locked_sample_qc_report.csv",
                "frozen_sets/locked_meta_matched.rds",
                "frozen_sets/locked_merged_dge.rds",
                "frozen_sets/canonical_deg_results.csv",
                "source_snapshot/data/gencode_v49_gene_metadata.tsv.gz",
            ]
            for relative in required:
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("placeholder\n")
            for arm in ("R0", "F_locked", "F_legacy", "F_five"):
                arm_root = root / "arms" / arm
                for relative in ("qc", "results/integration", "logs", "provenance"):
                    (arm_root / relative).mkdir(parents=True, exist_ok=True)
                (arm_root / ".bg001_candidate_root").write_text(arm + "\n")
            bad = root / "arms/R0/qc/preexisting.txt"
            bad.write_text("must fail\n")
            result = subprocess.run(
                ["bash", str(SCRIPT_DIR / "run_four_arms.sh"), str(root)],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Refusing nonempty or previously attempted arm R0", result.stderr)


class ComparisonSealTests(unittest.TestCase):
    def test_comparison_manifest_uses_exact_fail_closed_sealer_and_is_verified(self):
        script = (SCRIPT_DIR / "run_four_arms.sh").read_text()
        self.assertIn('seal_analysis_artifacts.py" seal', script)
        verification = 'seal_analysis_artifacts.py" verify'
        self.assertIn(verification, script)
        self.assertIn('ANALYSIS_MARKER="$RUN_ROOT/ANALYSIS_COMPLETE.json"', script)
        sealer = (SCRIPT_DIR / "seal_analysis_artifacts.py").read_text()
        self.assertIn("Artifact tree contains temporary file", sealer)
        self.assertIn("Artifact tree file set, size, or hash differs", sealer)

    def test_runtime_comparison_and_dependency_code_use_only_frozen_source(self):
        runner = (SCRIPT_DIR / "run_four_arms.sh").read_text()
        self.assertNotIn("LIVE_PROJECT_ROOT", runner)
        self.assertIn(
            'compare_scientific_gates.R" "$SNAP" "$RUN_ROOT"',
            runner,
        )
        self.assertIn(
            'generate_dependency_manifest.py" --project-root "$SNAP" --run-root "$RUN_ROOT"',
            runner,
        )
        comparison = (SCRIPT_DIR / "compare_scientific_gates.R").read_text()
        self.assertIn(
            'expected_project_root <- normalizePath(file.path(run_root, "source_snapshot")',
            comparison,
        )

    def test_dependency_manifest_rejects_a_non_snapshot_project_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            snapshot = root / "source_snapshot"
            snapshot.mkdir()
            dependency_manifest.require_frozen_project_root(snapshot, root)
            with self.assertRaisesRegex(SystemExit, "frozen source_snapshot"):
                dependency_manifest.require_frozen_project_root(root, root)

    def test_comparison_entrypoint_rejects_a_non_snapshot_project_root(self):
        if shutil.which("micromamba") is None:
            self.skipTest("micromamba is unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            (root / "source_snapshot").mkdir()
            (root / ".bg001_candidate_root").write_text(root.name + "\n")
            result = subprocess.run(
                [
                    "micromamba",
                    "run",
                    "-n",
                    "rnaseq",
                    "Rscript",
                    str(SCRIPT_DIR / "compare_scientific_gates.R"),
                    str(root),
                    str(root),
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("frozen source_snapshot", result.stderr)


class FrozenHelperInvocationTests(unittest.TestCase):
    def test_sbatch_wrappers_do_not_require_helper_execute_bits(self):
        expected = {
            "bam_manifest_hash_array.sbatch": 'exec bash "$SCRIPT_DIR/hash_bam_manifest_task.sh"',
            "recount_array.sbatch": 'exec bash "$SCRIPT_DIR/run_featurecounts_contract.sh"',
            "analysis_arms.sbatch": 'exec bash "$SCRIPT_DIR/run_four_arms.sh"',
        }
        for filename, invocation in expected.items():
            with self.subTest(filename=filename):
                self.assertIn(invocation, (SCRIPT_DIR / filename).read_text())


class ManifestLineEndingTests(unittest.TestCase):
    def test_manifest_writers_use_lf_and_shell_defensively_strips_cr(self):
        for filename in ("build_bam_manifest.py", "finalize_bam_manifest.py"):
            with self.subTest(filename=filename):
                self.assertIn('lineterminator="\\n"', (SCRIPT_DIR / filename).read_text())
        contract = (SCRIPT_DIR / "run_featurecounts_contract.sh").read_text()
        self.assertEqual(contract.count("expected_sha=${expected_sha%$'\\r'}"), 2)


class ProtectedPanelParsingTests(unittest.TestCase):
    def test_heldout_panel_parser_requires_exact_table_header(self):
        freeze = (SCRIPT_DIR / "freeze_baseline.R").read_text()
        self.assertIn('startsWith(lines, "gene_symbol\\t")', freeze)
        self.assertIn("length(header_lines) != 1L", freeze)
        self.assertNotIn('skip = "gene_symbol"', freeze)


class ProtectedActiveFigureLabelTests(unittest.TestCase):
    # Keep this registry narrow: these are literal vectors whose entries are
    # displayed as figure labels. Dynamic/data-selected labels require frozen
    # output-side inventories and should not be inferred from all uppercase tokens.
    ACTIVE_LITERAL_LABEL_VECTORS = {
        "scripts/figures/spatial_reorganization.R": ("hep_label",),
    }

    def test_active_literal_label_vectors_are_in_inventory(self):
        project = SCRIPT_DIR.parents[4]
        inventory_path = SCRIPT_DIR / "protected_figure_genes.tsv"
        with inventory_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        inventory_symbols = {row["symbol"] for row in rows}

        for relative_path, vector_names in self.ACTIVE_LITERAL_LABEL_VECTORS.items():
            source = project / relative_path
            text = source.read_text()
            for vector_name in vector_names:
                match = re.search(
                    rf"\b{re.escape(vector_name)}\s*<-\s*c\((.*?)\)",
                    text,
                    flags=re.DOTALL,
                )
                self.assertIsNotNone(match, f"Missing literal label vector {vector_name} in {relative_path}")
                labels = set(re.findall(r'[\"\']([A-Z][A-Z0-9.-]*)[\"\']', match.group(1)))
                self.assertTrue(labels, f"No gene labels parsed from {relative_path}:{vector_name}")
                self.assertEqual(
                    labels - inventory_symbols,
                    set(),
                    f"Active figure labels missing from protected inventory: {relative_path}:{vector_name}",
                )

    def test_hmgcs1_provenance_points_to_active_label(self):
        project = SCRIPT_DIR.parents[4]
        with (SCRIPT_DIR / "protected_figure_genes.tsv").open(newline="") as handle:
            matches = [
                row for row in csv.DictReader(handle, delimiter="\t")
                if row["symbol"] == "HMGCS1"
            ]
        self.assertEqual(len(matches), 1)
        row = matches[0]
        self.assertEqual(row["source_path"], "scripts/figures/spatial_reorganization.R")
        self.assertEqual(row["source_line"], "89")
        self.assertEqual(row["baseline_status"], "testable")
        source_line = (project / row["source_path"]).read_text().splitlines()[int(row["source_line"]) - 1]
        self.assertIn('hep_label <- c(', source_line)
        self.assertIn('"HMGCS1"', source_line)


class RenderedFigureLabelInventoryTests(unittest.TestCase):
    @staticmethod
    def write_tsv(path: Path, fields, rows) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    def setUp(self):
        self.project = SCRIPT_DIR.parents[4]
        self.paths = figure_labels.default_paths(self.project)

    def test_full_active_pdf_scope_and_dynamic_inventory_validate(self):
        result = figure_labels.validate_inventory(self.project, *self.paths)
        self.assertEqual(
            {
                key: result[key]
                for key in (
                    "panels", "dynamic_panels", "static_panels", "no_gene_panels",
                    "classified_rows", "gene_label_rows", "dynamic_gene_label_rows",
                    "unique_dynamic_gene_labels", "non_gene_collision_rows",
                )
            },
            {
                "panels": 77,
                "dynamic_panels": 24,
                "static_panels": 9,
                "no_gene_panels": 44,
                "classified_rows": 485,
                "gene_label_rows": 445,
                "dynamic_gene_label_rows": 388,
                "unique_dynamic_gene_labels": 261,
                "non_gene_collision_rows": 40,
            },
        )
        self.assertRegex(result["pdftotext_version"], r"^pdftotext version ")

    def test_corrected_scope_covers_full_module_and_data_selected_figs2n(self):
        scope = {
            row["panel_id"]: row
            for row in figure_labels.read_exact_tsv(self.paths[0], figure_labels.SCOPE_FIELDS)
        }
        roster = {
            row["panel_id"]: row
            for row in figure_labels.read_exact_tsv(self.paths[1], figure_labels.ROSTER_FIELDS)
        }
        self.assertEqual(scope["figs3_module_heatmap_full"]["label_mode"], "dynamic_rendered")
        self.assertEqual(scope["figs2n_progression_coloc"]["label_mode"], "dynamic_rendered")
        self.assertEqual(
            roster["figs2n_progression_coloc"]["producer_script"],
            "scripts/figures/progression_driver_coloc_phenotype_class_heatmap.R",
        )
        producer = (self.project / roster["figs2n_progression_coloc"]["producer_script"]).read_text()
        self.assertRegex(producer, r"core\s*<-\s*g\[coloc_best_pp4\s*>\s*0[.]5\]")

    def test_generated_sidecars_are_utf8_lf_only(self):
        for path in self.paths:
            with self.subTest(path=path.name):
                payload = path.read_bytes()
                payload.decode("utf-8", errors="strict")
                self.assertNotIn(b"\r", payload)
                self.assertTrue(payload.endswith(b"\n"))

    def test_missing_rendered_token_fails_closed(self):
        rows = figure_labels.read_exact_tsv(self.paths[3], figure_labels.TOKEN_FIELDS)
        with tempfile.TemporaryDirectory() as tmp:
            mutated = Path(tmp) / "tokens.tsv"
            self.write_tsv(mutated, figure_labels.TOKEN_FIELDS, rows[:-1])
            with self.assertRaisesRegex(figure_labels.InventoryError, "missing="):
                figure_labels.validate_inventory(
                    self.project, self.paths[0], self.paths[1], self.paths[2], mutated, self.paths[4]
                )

    def test_gene_to_collision_reclassification_fails_closed(self):
        rows = figure_labels.read_exact_tsv(self.paths[3], figure_labels.TOKEN_FIELDS)
        row = next(item for item in rows if item["classification"] == "gene_label")
        row["classification"] = "non_gene_collision"
        row["reason"] = "mutated"
        with tempfile.TemporaryDirectory() as tmp:
            mutated = Path(tmp) / "tokens.tsv"
            self.write_tsv(mutated, figure_labels.TOKEN_FIELDS, rows)
            with self.assertRaisesRegex(figure_labels.InventoryError, "Unapproved token classification"):
                figure_labels.validate_inventory(
                    self.project, self.paths[0], self.paths[1], self.paths[2], mutated, self.paths[4]
                )

    def test_pdf_hash_mutation_fails_closed(self):
        rows = figure_labels.read_exact_tsv(self.paths[2], figure_labels.INVENTORY_FIELDS)
        rows[0]["pdf_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as tmp:
            mutated = Path(tmp) / "inventory.tsv"
            self.write_tsv(mutated, figure_labels.INVENTORY_FIELDS, rows)
            with self.assertRaisesRegex(figure_labels.InventoryError, "PDF/text binding changed"):
                figure_labels.validate_inventory(
                    self.project, self.paths[0], self.paths[1], mutated, self.paths[3], self.paths[4]
                )

    def test_static_glp_axis_has_exact_source_proof(self):
        rows = figure_labels.read_exact_tsv(self.paths[4], figure_labels.PROTECTED_FIELDS)
        by_symbol = {row["symbol"]: row for row in rows}
        for symbol, status in {
            "DPP4": "testable", "GCG": "NOT_TESTABLE", "GCGR": "testable",
            "GIPR": "testable", "GLP1R": "NOT_TESTABLE", "GLP2R": "testable",
        }.items():
            with self.subTest(symbol=symbol):
                row = by_symbol[symbol]
                self.assertEqual(row["source_path"], "scripts/figures/fig5e_glp1ra_incretin_axis.R")
                self.assertEqual(row["source_line"], "38")
                self.assertEqual(row["baseline_status"], status)

    def test_prepare_and_freeze_bind_rendered_inventory(self):
        prepare_text = (SCRIPT_DIR / "prepare_run.py").read_text()
        self.assertLess(
            prepare_text.index("live_figure_gate = subprocess.run"),
            prepare_text.index("run_root.mkdir(mode=0o700)"),
        )
        self.assertIn('explicit_sources.update(row["pdf_path"] for row in figure_rows)', prepare_text)
        self.assertIn("frozen_figure_gate = subprocess.run", prepare_text)
        self.assertIn('"rendered_figure_label_check_sha256"', prepare_text)

        freeze_text = (SCRIPT_DIR / "freeze_baseline.R").read_text()
        self.assertLess(
            freeze_text.index("rendered_check <- suppressWarnings(system2"),
            freeze_text.index("locked_qc_path <- copy_exact"),
        )
        self.assertIn('classification == "gene_label"', freeze_text)
        self.assertIn('label_mode == "dynamic_rendered"', freeze_text)
        self.assertIn('"active_rendered_figure_label"', freeze_text)


class ExactPerSampleMergeTests(unittest.TestCase):
    merger = (
        SCRIPT_DIR.parents[1]
        / "pipelines/shared/scripts/merge_featurecounts_exact.py"
    )

    @staticmethod
    def write_sample(directory: Path, sample: str, rows, *, paired: bool = True, assigned: bool = True):
        if paired:
            command = (
                '# Program:featureCounts v2.1.1; Command:"featureCounts" "-p" '
                '"--countReadPairs" "-B" "-s" "2" "-a" "test.gtf"\n'
            )
        else:
            command = '# Program:featureCounts v2.1.1; Command:"featureCounts" "-s" "2" "-a" "test.gtf"\n'
        path = directory / f"{sample}.counts.txt"
        with path.open("w") as handle:
            handle.write(command)
            handle.write(f"Geneid\tChr\tStart\tEnd\tStrand\tLength\t/{sample}.Aligned.sortedByCoord.out.bam\n")
            for gene, value in rows:
                handle.write(f"{gene}\tchr1\t1\t2\t+\t2\t{value}\n")
        with Path(f"{path}.summary").open("w") as handle:
            handle.write(f"Status\t/{sample}.Aligned.sortedByCoord.out.bam\n")
            status = "Assigned" if assigned else "Unassigned_NoFeatures"
            handle.write(f"{status}\t{sum(value for _, value in rows)}\n")

    def test_keyed_merge_accepts_reordered_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            per = root / "per"
            per.mkdir()
            self.write_sample(per, "S1", [("G1", 1), ("G2", 2)])
            self.write_sample(per, "S2", [("G2", 4), ("G1", 3)])
            out = root / "gene_counts.txt"
            result = subprocess.run(
                [
                    sys.executable, str(self.merger),
                    "--per-sample-dir", str(per),
                    "--sample-id", "S1", "--sample-id", "S2",
                    "--output", str(out), "--summary", f"{out}.summary",
                    "--expected-genes", "2",
                ],
                text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = out.read_text().splitlines()
            self.assertEqual(lines[1].split("\t")[-2:], ["S1", "S2"])
            self.assertEqual(lines[2].split("\t")[-2:], ["1", "3"])
            self.assertEqual(lines[3].split("\t")[-2:], ["2", "4"])
            provenance = json.loads(Path(f"{out}.merge_provenance.json").read_text())
            marker = json.loads(Path(f"{out}.merge_complete").read_text())
            self.assertEqual(provenance["status"], "PASS")
            self.assertEqual(marker["count_sha256"], hashlib.sha256(out.read_bytes()).hexdigest())
            self.assertEqual(
                marker["summary_sha256"],
                hashlib.sha256(Path(f"{out}.summary").read_bytes()).hexdigest(),
            )

    def test_annotation_drift_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            per = root / "per"
            per.mkdir()
            self.write_sample(per, "S1", [("G1", 1), ("G2", 2)])
            self.write_sample(per, "S2", [("G1", 3), ("G3", 4)])
            out = root / "gene_counts.txt"
            result = subprocess.run(
                [
                    sys.executable, str(self.merger),
                    "--per-sample-dir", str(per),
                    "--sample-id", "S1", "--sample-id", "S2",
                    "--output", str(out), "--summary", f"{out}.summary",
                    "--expected-genes", "2",
                ],
                text=True, capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(out.exists())

    def test_mixed_layout_sources_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            per = root / "per"
            per.mkdir()
            self.write_sample(per, "S1", [("G1", 1), ("G2", 2)], paired=True)
            self.write_sample(per, "S2", [("G1", 3), ("G2", 4)], paired=False)
            out = root / "gene_counts.txt"
            result = subprocess.run(
                [
                    sys.executable, str(self.merger),
                    "--per-sample-dir", str(per),
                    "--sample-id", "S1", "--sample-id", "S2",
                    "--output", str(out), "--summary", f"{out}.summary",
                    "--expected-genes", "2",
                ],
                text=True, capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("mixed single-end and paired-fragment", result.stderr)
            self.assertFalse(out.exists())

    def test_missing_assigned_summary_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            per = root / "per"
            per.mkdir()
            self.write_sample(per, "S1", [("G1", 1), ("G2", 2)], assigned=False)
            out = root / "gene_counts.txt"
            result = subprocess.run(
                [
                    sys.executable, str(self.merger),
                    "--per-sample-dir", str(per),
                    "--sample-id", "S1",
                    "--output", str(out), "--summary", f"{out}.summary",
                    "--expected-genes", "2",
                ],
                text=True, capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Assigned status", result.stderr)
            self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
