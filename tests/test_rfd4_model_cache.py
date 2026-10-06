"""Standard-library regression tests: python -m unittest discover -s tests."""
from copy import deepcopy
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import rfd4_model_cache as cache


def model():
    atom = {"chain": "A", "resi": "204", "resn": "ALA", "name": "N", "element": "N",
            "coord": [1.0, 2.0, 3.0], "hetero": False, "occupancy": 1.0, "charge": 0,
            "fixed": True, "sequence_fixed": False, "unindexed": True, "index_group": 3,
            "hbond": "DONOR", "rasa": "UNSPECIFIED"}
    second = {**atom, "name": "CA", "element": "C", "coord": [2.3, 2.0, 3.0], "hbond": "NONE"}
    return {"source": "/old/name.cif", "atoms": [atom, second], "bonds": [[0, 1, 1]],
            "excluded_atoms": 1, "segments": [{"chain": "A", "resi": 1, "min": 120, "max": 260}],
            "condition_fields": ["condition_coordinate"], "parse": "test-native",
            "present": {name: True for name in ("fixed", "sequence_fixed", "unindexed", "hbond", "rasa")}}


class CacheTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        self.reader = {"sha256": "r" * 64, "evidence": {"test": "native-reader"}}
        self.source_hash = "s" * 64

    def tearDown(self):
        self.temp.cleanup()

    def write(self, value=None):
        return cache.write_cached(self.source_hash, self.reader, value or model(), self.directory)

    def read(self, source="/new/name.cif", digest=None, reader=None):
        return cache.read_cached(source, digest or self.source_hash, reader or self.reader["sha256"], self.directory)

    def test_roundtrip_preserves_structure_conditions_and_updates_moved_path(self):
        self.write()
        expected = model()
        expected["source"] = "/new/name.cif"
        self.assertEqual(expected, self.read())
        self.assertEqual("/renamed/input.cif", self.read("/renamed/input.cif")["source"])

    def test_source_and_reader_changes_miss(self):
        self.write()
        self.assertIsNone(self.read(digest="new contents"))
        self.assertIsNone(self.read(reader="new parser"))

    def test_payload_corruption_is_not_used(self):
        path = self.write()
        envelope = json.loads(path.read_bytes())
        envelope["model"]["atoms"][0]["coord"][0] = 999.0
        path.write_text(json.dumps(envelope))
        self.assertIsNone(self.read())

    def test_malformed_json_recovers_as_miss(self):
        path = self.write()
        for malformed in ("[", "[]", "null", '{"schema":1}'):
            path.write_text(malformed)
            self.assertIsNone(self.read())

    def test_schema_validated_even_with_matching_checksum(self):
        path = self.write()
        envelope = json.loads(path.read_bytes())
        envelope["model"]["bonds"] = [[0, 500, 1]]
        envelope["model_sha256"] = cache.sha256(cache.canonical(envelope["model"]))
        path.write_bytes(cache.canonical(envelope))
        self.assertIsNone(self.read())

    def test_rejects_nonfinite_coordinate_and_bad_identity(self):
        broken = model()
        broken["atoms"][0]["coord"][0] = float("nan")
        with self.assertRaises(ValueError):
            self.write(broken)
        broken = model()
        broken["atoms"][1] = deepcopy(broken["atoms"][0])
        with self.assertRaises(ValueError):
            self.write(broken)

    def test_cache_does_not_save_pymol_session_indices(self):
        dirty = model()
        dirty["atoms"][0]["pymol_index"] = 900
        dirty["selections_created"] = True
        self.write(dirty)
        clean = self.read()
        self.assertNotIn("pymol_index", clean["atoms"][0])
        self.assertNotIn("selections_created", clean)
        self.assertEqual(900, dirty["atoms"][0]["pymol_index"])

    def test_atomic_json_never_leaves_pending_files(self):
        path = self.write()
        self.assertEqual(0o600, path.stat().st_mode & 0o777)
        self.assertFalse(list(self.directory.rglob(".pending-*")))

    def test_disk_hit_does_not_invoke_native_reader(self):
        source = self.directory / "test.cif"
        source.write_text("data_native_fixture")
        self.source_hash = cache.sha256(source.read_bytes())
        self.write()
        with patch.object(cache, "reader_config", return_value={}), patch.object(cache, "provenance", return_value=self.reader), \
             patch.object(cache.NativeReader, "model", side_effect=AssertionError("Native parser must not run")):
            result = cache.warm([source], directory=self.directory, jobs=2)
        self.assertEqual(1, result["hits"])
        self.assertEqual(0, result["parsed"])

    def test_missing_cache_uses_native_then_next_call_reuses(self):
        source = self.directory / "test.cif"
        source.write_text("data_native_fixture")
        with patch.object(cache, "reader_config", return_value={}), patch.object(cache, "provenance", return_value=self.reader), \
             patch.object(cache.NativeReader, "model", return_value=model()) as native:
            first = cache.warm([source], directory=self.directory)
            second = cache.warm([source], directory=self.directory)
        self.assertEqual(1, first["parsed"])
        self.assertEqual(1, second["hits"])
        self.assertEqual(1, native.call_count)

    def test_dirty_code_invalidates_digest_even_if_mtime_restored(self):
        source = self.directory / "parser.py"
        source.write_text("old\n")
        original_stat = source.stat()
        before = cache._code_digest(source)
        time.sleep(0.02)  # This filesystem coalesces metadata timestamps within a clock tick.
        source.write_text("new\n")
        os.utime(source, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
        after = cache._code_digest(source)
        self.assertNotEqual(before[1], after[1])


if __name__ == "__main__":
    unittest.main()
