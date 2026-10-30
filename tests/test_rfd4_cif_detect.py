"""Regression: saved AF3 cofactors must never enter the annotated-input filter."""
from pathlib import Path
import sys
import tempfile
import unittest
import os
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from rfd4_cif_detect import is_rfd4_cif


class DetectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.source = Path(self.temp.name) / "renamed_structure.cif"

    def tearDown(self):
        self.temp.cleanup()

    def detect(self, body):
        self.source.write_text("data_test\n" + body)
        return is_rfd4_cif(self.source)

    def test_saved_prediction_retains_components_with_unknown_occupancy(self):
        self.assertFalse(self.detect("""loop_
_atom_site.label_comp_id
_atom_site.occupancy
_atom_site.annotation_expsegmin_1_residue
_atom_site.annotation_expsegmax_1_residue
ALA 1.0 -1 -1
CEC_GS nan -1 -1
AKG_GS nan -1 -1
FE2 nan -1 -1
"""))

    def test_real_conditions_work_without_expseg_columns(self):
        for field in ("_condition_coordinate_1_atom.idx0",
                      "_condition_ChainCentroid_1_atom.val0",
                      "_atom_site.condition_sequence_1_residue",
                      "_atom_site.condition_hbond_1_atom"):
            with self.subTest(field=field):
                self.assertTrue(self.detect("loop_\n" + field + "\n0\n"))

    def test_comments_quotes_and_multiline_text_are_not_conditions(self):
        self.assertFalse(self.detect("""# _condition_coordinate_1_atom.idx0
_entry.details
;
_condition_coordinate_1_atom.idx0
_atom_site.condition_hbond_1_atom
;
_entry.id '_atom_site.condition_hbond_1_atom'
' _condition_coordinate_1_atom.idx0 '
"""))

    def test_condition_after_multiline_text_and_mmcif_extension(self):
        self.source = self.source.with_suffix(".mmcif")
        self.assertTrue(self.detect("_entry.details\n;note\n;\n"
                                    "  _ATOM_SITE.CONDITION_HBOND_1_ATOM 2\n"))

    def test_missing_files_and_other_formats_keep_original_loader(self):
        self.assertFalse(is_rfd4_cif(self.source))
        self.source = self.source.with_suffix(".pdb")
        self.assertFalse(self.detect("_atom_site.condition_hbond_1_atom 2\n"))

    def test_quoted_filename(self):
        self.detect("_atom_site.condition_hbond_1_atom 2\n")
        self.assertTrue(is_rfd4_cif('"' + str(self.source) + '"'))

    def test_inline_fields_and_quote_boundaries(self):
        self.assertTrue(self.detect("_entry.id 'a value' _condition_coordinate_1_atom.idx0 2\n"))
        self.assertFalse(self.detect("_entry.details \"a quote' _condition_coordinate_1_atom.idx0\"\n"))
        self.assertFalse(self.detect("_entry.details 'don\'t _condition_coordinate_1_atom.idx0'\n"))

    def test_marker_across_binary_chunk_boundary(self):
        self.source.write_bytes(b"data_t\n#" + b"x" * (256 * 1024 - 12) + b"\n_atom_site.condition_hbond_1_atom 2\n")
        self.assertTrue(is_rfd4_cif(self.source))

    def test_cached_negative_does_not_reopen_and_invalidates_on_replacement(self):
        self.assertFalse(self.detect("_entry.id ordinary\n"))
        with patch.object(Path, "open", side_effect=AssertionError("unchanged file reopened")):
            self.assertFalse(is_rfd4_cif(self.source))
        replacement = self.source.with_suffix('.tmp')
        replacement.write_text("data_t\n_atom_site.condition_hbond_1_atom 2\n")
        os.replace(replacement, self.source)
        self.assertTrue(is_rfd4_cif(self.source))


if __name__ == "__main__":
    unittest.main()
