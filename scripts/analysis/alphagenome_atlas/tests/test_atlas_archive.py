"""Round-trip test for the raw Atlas response archive (offline; synthetic AnnData)."""

from __future__ import annotations

import pathlib
import sys
import tempfile
import unittest

import grpc

import anndata
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import atlas_archive as aa  # noqa: E402


def fake_results():
    obs = pd.DataFrame({"variant": ["chr1:5:A>G", "chr1:5:A>G", "chr1:9:C>T"], "gene_id": ["ENSG1", "ENSG2", "ENSG1"]}, index=["0", "1", "2"])
    var = pd.DataFrame({"name": ["t1", "t2"], "ontology_curie": ["UBERON:0002107", "EFO:1"], "biosample_name": ["liver", "HepG2"]}, index=["0", "1"])
    x = np.array([[0.1, -0.2], [0.3, 0.4], [np.nan, 1.0]], dtype=np.float32)
    gene = anndata.AnnData(X=x, obs=obs, var=var, layers={"quantiles": np.abs(x)})
    center = anndata.AnnData(X=np.array([[2.0, 3.0], [4.0, 5.0]], dtype=np.float32),
                             obs=pd.DataFrame({"variant": ["chr1:5:A>G", "chr1:9:C>T"]}, index=["0", "1"]),
                             var=var.copy())
    return {"gene_scorer": gene, "center_scorer": center}


class TestArchive(unittest.TestCase):
    def test_round_trip_preserves_scores_quantiles_and_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = pathlib.Path(tmp) / "chunk_00000"
            aa.archive_scores(fake_results(), out, request={"variants": ["chr1:5:A>G", "chr1:9:C>T"]})
            back = aa.load_archive(out)
            self.assertEqual(set(back), {"gene_scorer", "center_scorer"})
            g = back["gene_scorer"]
            np.testing.assert_array_equal(np.isnan(g.X), np.isnan(fake_results()["gene_scorer"].X))
            self.assertAlmostEqual(float(g.X[1, 1]), 0.4, places=6)
            self.assertIn("quantiles", g.layers)
            self.assertEqual(list(g.obs["variant"]), ["chr1:5:A>G", "chr1:5:A>G", "chr1:9:C>T"])
            self.assertEqual(list(g.var["ontology_curie"]), ["UBERON:0002107", "EFO:1"])
            self.assertNotIn("quantiles", back["center_scorer"].layers)
            self.assertTrue((out / "request.json").exists())

    def test_refuses_to_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = pathlib.Path(tmp) / "chunk_00000"
            aa.archive_scores(fake_results(), out, request={})
            with self.assertRaises(aa.la.ContractError):
                aa.archive_scores(fake_results(), out, request={})

    def test_long_records_one_row_per_variant_gene_track(self):
        rows = list(aa.long_records(fake_results()["gene_scorer"], scorer="gene_scorer", is_signed=True))
        self.assertEqual(len(rows), 6)
        r = rows[0]
        self.assertEqual((r["variant_uid"], r["gene_id"], r["track_name"]), ("chr1:5:A:G", "ENSG1", "t1"))
        self.assertAlmostEqual(r["raw_score"], 0.1, places=6)
        self.assertAlmostEqual(r["quantile"], 0.1, places=6)
        self.assertEqual(r["track_class"], "primary_liver")
        self.assertEqual(rows[1]["track_class"], "HepG2")


if __name__ == "__main__":
    unittest.main()


class TestQuotaRetry(unittest.TestCase):
    def test_retries_only_quota_errors_and_returns_the_value(self):
        import atlas_query as aq

        class Quota(Exception):
            pass
        calls = {"n": 0}

        def fn():
            calls["n"] += 1
            if calls["n"] < 3:
                raise Quota("per-minute quota")
            return "ok"
        slept = []
        out = aq.call_with_quota_retry(fn, is_quota=lambda e: isinstance(e, Quota), sleep=lambda s: slept.append(s))
        self.assertEqual(out, "ok"); self.assertEqual(calls["n"], 3); self.assertEqual(len(slept), 2)
        with self.assertRaises(ValueError):
            aq.call_with_quota_retry(lambda: (_ for _ in ()).throw(ValueError("other")), is_quota=lambda e: isinstance(e, Quota), sleep=lambda s: None)


class TestTransientRetry(unittest.TestCase):
    def test_unavailable_is_retried_a_bounded_number_of_times_then_raised(self):
        import atlas_query as aq

        class Unavailable(Exception):
            pass
        calls = {"n": 0}

        def flaky():
            calls["n"] += 1
            if calls["n"] < 3:
                raise Unavailable("service unavailable")
            return "ok"
        slept = []
        out = aq.call_with_quota_retry(flaky, is_quota=lambda e: False, is_transient=lambda e: isinstance(e, Unavailable), sleep=lambda s: slept.append(s))
        self.assertEqual((out, calls["n"], len(slept)), ("ok", 3, 2))

        def always():
            raise Unavailable("down")
        with self.assertRaises(Unavailable):
            aq.call_with_quota_retry(always, is_quota=lambda e: False, is_transient=lambda e: isinstance(e, Unavailable), sleep=lambda s: None, max_transient=4)


class _FakeRpcError(grpc.RpcError):
    def __init__(self, code, details):
        self._code, self._details = code, details

    def code(self):
        return self._code

    def details(self):
        return self._details


class TestMessageSizeIsNotQuota(unittest.TestCase):
    """gRPC reports an oversized RESPONSE with the same status code as a quota breach.

    Retrying an oversized request can never succeed: it re-issues the identical request and gets the identical
    response. Treating it as a quota error made two jobs sleep 65 s and retry indefinitely.
    """

    def setUp(self):
        import atlas_query as aq
        self.aq = aq

    def test_oversized_response_is_classified_as_message_size(self):
        e = _FakeRpcError(grpc.StatusCode.RESOURCE_EXHAUSTED,
                          "Stream removed (CLIENT: Received message larger than max (22031061 vs. 4194304))")
        self.assertTrue(self.aq.is_message_size_error(e))
        self.assertFalse(self.aq.is_quota_error(e))

    def test_real_quota_error_is_still_a_quota_error(self):
        e = _FakeRpcError(grpc.StatusCode.RESOURCE_EXHAUSTED, "Quota exceeded for requests per minute")
        self.assertTrue(self.aq.is_quota_error(e))
        self.assertFalse(self.aq.is_message_size_error(e))

    def test_message_size_error_propagates_instead_of_looping(self):
        e = _FakeRpcError(grpc.StatusCode.RESOURCE_EXHAUSTED,
                          "Received message larger than max (22031061 vs. 4194304)")
        calls = {"n": 0, "slept": 0}

        def boom():
            calls["n"] += 1
            raise e

        with self.assertRaises(grpc.RpcError):
            self.aq.call_with_quota_retry(boom, sleep=lambda s: calls.__setitem__("slept", calls["slept"] + 1))
        self.assertEqual(calls["n"], 1)
        self.assertEqual(calls["slept"], 0)


class TestIntervalSplitOnOversizedResponse(unittest.TestCase):
    """A window too wide for one response is halved until each piece fits, rather than retried unchanged."""

    def setUp(self):
        import atlas_query as aq
        self.aq = aq

    def test_window_is_halved_until_each_piece_fits(self):
        seen = []

        def fake_query(start, end):
            seen.append((start, end))
            if end - start > 500:
                raise _FakeRpcError(grpc.StatusCode.RESOURCE_EXHAUSTED,
                                    "Received message larger than max (9 vs. 4)")
            return {"ATAC": f"{start}-{end}"}

        out = self.aq.query_interval_split(fake_query, 1000, 3000, min_width=100)
        self.assertEqual(len(out), 4)                     # 2000 -> 1000 -> 500 each
        self.assertEqual([e - s for s, e in seen if e - s <= 500], [500, 500, 500, 500])
        self.assertEqual(seen[0], (1000, 3000))           # the full window is tried first

    def test_refuses_to_split_below_the_floor(self):
        def always_too_big(start, end):
            raise _FakeRpcError(grpc.StatusCode.RESOURCE_EXHAUSTED,
                                "Received message larger than max (9 vs. 4)")

        with self.assertRaises(grpc.RpcError):
            self.aq.query_interval_split(always_too_big, 0, 128, min_width=64)


class TestClientRaisesMessageLimit(unittest.TestCase):
    """Every script builds its client through create_client, so the raised limit must be set there."""

    def test_channel_is_built_with_a_raised_receive_limit(self):
        import atlas_query as aq
        seen = {}

        class _Chan:
            pass

        def fake_secure_channel(address, creds, options=()):
            seen["options"] = dict(options)
            return _Chan()

        orig_secure, orig_ready, orig_creds = grpc.secure_channel, grpc.channel_ready_future, grpc.ssl_channel_credentials
        grpc.secure_channel = fake_secure_channel
        grpc.channel_ready_future = lambda ch: type("F", (), {"result": lambda self, t: None})()
        grpc.ssl_channel_credentials = lambda: None
        try:
            aq.create_client_with_large_messages("fake-key")
        except Exception:
            pass  # stub construction may fail; the channel options are what this test asserts
        finally:
            grpc.secure_channel, grpc.channel_ready_future, grpc.ssl_channel_credentials = orig_secure, orig_ready, orig_creds
        self.assertGreaterEqual(seen["options"]["grpc.max_receive_message_length"], 64 * 1024 * 1024)


class TestTransientBackoff(unittest.TestCase):
    """A real Atlas outage outlasted a flat 10 x 30 s budget and killed the saturation pilot.

    The service was verified up minutes later with the identical request, so the request shape was fine and
    the retry budget was the problem. Backoff must grow, cap, and cover a materially longer outage.
    """

    def setUp(self):
        import atlas_query as aq
        self.aq = aq

    def test_backoff_grows_then_caps(self):
        d = [self.aq.transient_sleep_seconds(i) for i in range(1, 12)]
        self.assertTrue(all(a <= b for a, b in zip(d, d[1:])), f"backoff must be non-decreasing, got {d}")
        self.assertTrue(all(x <= self.aq.TRANSIENT_SLEEP_CAP for x in d))
        self.assertEqual(d[-1], self.aq.TRANSIENT_SLEEP_CAP)
        self.assertLess(d[0], d[3])

    def test_total_budget_covers_a_long_outage(self):
        total = sum(self.aq.transient_sleep_seconds(i) for i in range(1, self.aq.MAX_TRANSIENT + 1))
        self.assertGreaterEqual(total, 1800)     # at least 30 minutes of outage survived

    def test_a_transient_error_that_clears_returns_the_value(self):
        import grpc
        calls = {"n": 0}

        def flaky():
            calls["n"] += 1
            if calls["n"] < 4:
                raise _FakeRpcError(grpc.StatusCode.UNAVAILABLE, "The service is currently unavailable.")
            return "ok"

        slept = []
        out = self.aq.call_with_quota_retry(flaky, sleep=slept.append)
        self.assertEqual(out, "ok")
        self.assertEqual(len(slept), 3)
        self.assertLess(slept[0], slept[2])      # backoff actually grew


class TestOntologyArchiveFilter(unittest.TestCase):
    """Saturation archives keep only the tracks the analysis reads.

    Server-side ontology filtering breaks query_interval, so retrieval is unfiltered; filtering at archive
    time instead takes a region from 46.7 MB to 10.3 MB (4.5 TB to 0.99 TB over the 96,460-region universe)
    without touching any number, because the same filter is applied at analysis time either way.
    """

    def setUp(self):
        import atlas_query as aq
        self.aq = aq

    def _res(self):
        import anndata
        var = pd.DataFrame({"ontology_curie": ["UBERON:0002107", "CL:9999999", "EFO:0001187", "CL:8888888"],
                            "name": list("abcd")})
        X = np.arange(8, dtype=np.float32).reshape(2, 4)
        a = anndata.AnnData(X=X, var=var, obs=pd.DataFrame({"variant": ["v1", "v2"]}))
        a.layers["quantiles"] = X.copy()
        return {"ATAC": a}

    def test_keeps_only_panel_tracks_and_preserves_layers(self):
        out, counts = self.aq.filter_tracks_to_ontology(self._res(), {"UBERON:0002107", "EFO:0001187"})
        a = out["ATAC"]
        self.assertEqual(a.shape, (2, 2))
        self.assertEqual(list(a.var["ontology_curie"]), ["UBERON:0002107", "EFO:0001187"])
        self.assertIn("quantiles", a.layers)
        self.assertEqual(a.layers["quantiles"].shape, (2, 2))

    def test_records_before_and_after_counts_so_the_filter_is_auditable(self):
        _, counts = self.aq.filter_tracks_to_ontology(self._res(), {"UBERON:0002107"})
        self.assertEqual(counts["ATAC"], {"tracks_returned": 4, "tracks_archived": 1})

    def test_a_scorer_with_no_ontology_column_is_kept_whole(self):
        import anndata
        a = anndata.AnnData(X=np.ones((2, 1), dtype=np.float32),
                            var=pd.DataFrame({"name": ["avi"]}), obs=pd.DataFrame({"variant": ["v1", "v2"]}))
        out, counts = self.aq.filter_tracks_to_ontology({"AVI_SCORE": a}, {"UBERON:0002107"})
        self.assertEqual(out["AVI_SCORE"].shape, (2, 1))
        self.assertEqual(counts["AVI_SCORE"]["tracks_archived"], 1)

    def test_no_matching_track_yields_an_empty_track_axis_not_a_crash(self):
        out, counts = self.aq.filter_tracks_to_ontology(self._res(), {"CL:0000000"})
        self.assertEqual(out["ATAC"].shape[1], 0)
        self.assertEqual(counts["ATAC"]["tracks_archived"], 0)

