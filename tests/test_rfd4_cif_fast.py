"""Portable fast-reader guards plus opt-in real native E/D/internal parity."""
from copy import deepcopy
import hashlib
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import rfd4_cif_fast as fast
import rfd4_model_cache as cache
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO


def fixture():
    d = {'data_': 'fixture', '_entity.id': ['1'], '_entity.type': ['polymer'],
         '_entity_poly.entity_id': ['1'], '_entity_poly.type': ['polypeptide(L)'],
         '_chem_comp.id': ['GLY'], '_chem_comp.type': ['PEPTIDE LINKING'],
         '_chem_comp_atom.comp_id': ['GLY'] * 4,
         '_chem_comp_atom.atom_id': ['N', 'CA', 'C', 'O'],
         '_chem_comp_atom.type_symbol': ['N', 'C', 'C', 'O'],
         '_chem_comp_atom.charge': ['0'] * 4,
         '_chem_comp_bond.comp_id': ['GLY'] * 3,
         '_chem_comp_bond.atom_id_1': ['N', 'CA', 'C'],
         '_chem_comp_bond.atom_id_2': ['CA', 'C', 'O'],
         '_chem_comp_bond.value_order': ['SING', 'SING', 'DOUB']}
    fields = {'id': list(map(str, range(8))), 'group_PDB': ['ATOM'] * 8,
              'type_symbol': ['N', 'C', 'C', 'O'] * 2, 'label_entity_id': ['1'] * 8,
              'label_alt_id': ['.'] * 8, 'pdbx_PDB_ins_code': ['.'] * 8,
              'pdbx_PDB_model_num': ['1'] * 8, 'pdbx_formal_charge': ['0'] * 8,
              'occupancy': ['1.0'] * 8,
              'annotation_expsegmin_1_residue': ['-1'] * 8,
              'annotation_expsegmax_1_residue': ['-1'] * 8,
              'condition_sequence_1_residue': ['GLY'] * 8,
              'condition_index_1_residue': ['True'] * 8,
              'condition_index-group_1_token': ['0'] * 8,
              'condition_hbond_1_atom': ['0'] * 8,
              'condition_rasa-binned_1_atom': ['0'] * 8,
              'condition_CoordinateGroup_1_atom': ['0', '7'] + ['-1'] * 6,
              'Cartn_x': list(map(str, range(8))), 'Cartn_y': ['0'] * 8, 'Cartn_z': ['0'] * 8}
    for prefix in ('label_', 'auth_'):
        fields[prefix + 'asym_id'] = ['A'] * 8
        fields[prefix + 'comp_id'] = ['GLY'] * 8
        fields[prefix + 'seq_id'] = ['1'] * 4 + ['2'] * 4
        fields[prefix + 'atom_id'] = ['N', 'CA', 'C', 'O'] * 2
    d.update({'_atom_site.' + k: v for k, v in fields.items()})
    d['_condition_coordinate_1_atom.idx0'] = list(map(str, range(8)))
    for axis in range(3):
        d['_condition_coordinate_1_atom.val' + str(axis)] = ['0'] + ['nan'] * 7
    return d


def write_dict(path, d):
    writer = MMCIFIO(); writer.set_dict(d); writer.save(str(path))


def set_centroids(d, xyz):
    count = len(d['_atom_site.id'])
    d['_condition_ChainCentroid_1_atom.idx0'] = list(map(str, range(count)))
    for axis in range(3):
        d['_condition_ChainCentroid_1_atom.val' + str(axis)] = [str(xyz[axis])] * count


class FastReaderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / 'input.cif'

    def tearDown(self):
        self.temp.cleanup()

    def parse(self, d):
        write_dict(self.path, d)
        return fast.export_model(self.path)

    def test_absolute_internal_reference_identity_and_bonds(self):
        m = self.parse(fixture())
        cache.validate_model(m)
        self.assertEqual(len(m['atoms']), 8)
        self.assertEqual(len(m['bonds']), 7)
        self.assertEqual([a['coordinate_group'] for a in m['atoms']], [0, 7] + [-1] * 6)
        self.assertTrue(all(a['sequence_fixed'] for a in m['atoms']))
        self.assertEqual(sum(a['fixed'] for a in m['atoms']), 1)
        self.assertEqual(sum(a['internal'] for a in m['atoms']), 1)

    def test_native_fallback_for_unsupported_or_malformed_fields(self):
        mutations = [
            ('_atom_site.pdbx_PDB_model_num', 0, '2'),
            ('_atom_site.label_alt_id', 0, 'A'),
            ('_atom_site.auth_seq_id', 0, '9'),
            ('_atom_site.pdbx_PDB_ins_code', 0, 'A'),
            ('_atom_site.type_symbol', 0, 'Zn'),
            ('_atom_site.group_PDB', 0, 'UNKNOWN'),
            ('_atom_site.condition_index_1_residue', 0, 'maybe'),
            ('_atom_site.condition_hbond_1_atom', 0, '9'),
            ('_atom_site.condition_CoordinateGroup_1_atom', 0, '-1'),
            ('_condition_coordinate_1_atom.val1', 0, 'nan'),
            ('_condition_coordinate_1_atom.idx0', 0, '8'),
            ('_chem_comp_bond.value_order', 0, 'AROM'),
            ('_chem_comp_atom.type_symbol', 0, 'C'),
            ('_entity.type', 0, 'non-polymer'),
            ('_atom_site.pdbx_formal_charge', 2, '1'),
        ]
        for key, index, value in mutations:
            with self.subTest(key=key):
                d = fixture(); d[key][index] = value
                with self.assertRaises(fast.UnsupportedCIF):
                    self.parse(d)

    def test_nonidentity_bookkeeping_index_requires_native_reader(self):
        d = fixture()
        d['_atom_site.atom_array_index'] = list(map(str, range(8)))
        self.assertEqual(self.parse(d), self.parse(fixture()))
        # Swap complete atom rows, retaining their old bookkeeping indices.
        # No interpretation/remapping of vector target rows is invented here.
        for key, values in d.items():
            if key.startswith('_atom_site.'):
                values[0], values[1] = values[1], values[0]
        with self.assertRaisesRegex(fast.UnsupportedCIF, 'atom_array_index'):
            self.parse(d)
        d = fixture(); d['_atom_site.atom_array_index'] = ['0'] * 8
        with self.assertRaisesRegex(fast.UnsupportedCIF, 'atom_array_index'):
            self.parse(d)

    def test_unknown_condition_and_missing_templates_fall_back(self):
        for key in ('_atom_site.condition_future_1_atom', '_atom_site.annotation_future_1_atom'):
            d = fixture(); d[key] = ['0'] * 8
            with self.assertRaises(fast.UnsupportedCIF):
                self.parse(d)
        d = fixture(); del d['_chem_comp_bond.value_order']
        with self.assertRaises(fast.UnsupportedCIF):
            self.parse(d)

    def test_duplicate_tags_and_multiple_blocks_fall_back(self):
        self.parse(fixture())
        baseline = self.path.read_text()
        for suffix in ('\n_entity.id 1\n', '\ndata_second\n_entry.id second\n'):
            self.path.write_text(baseline + suffix)
            with self.assertRaises(fast.UnsupportedCIF):
                fast.export_model(self.path)

    def test_empty_physical_selection_falls_back(self):
        d = fixture(); d['_atom_site.occupancy'] = ['0'] * 8
        with self.assertRaises(fast.UnsupportedCIF):
            self.parse(d)

    def test_atom_site_charge_is_preserved_without_template_override(self):
        d = fixture(); d['_atom_site.pdbx_formal_charge'][3] = '-1'
        self.assertEqual(self.parse(d)['atoms'][3]['charge'], -1)

    def test_centroids_preserve_zero_and_absence(self):
        self.assertEqual(self.parse(fixture())['centroids'], [])
        for xyz, expected in (([float('nan')] * 3, []),
                              ([0, 0, 0], [{'chain': 'A', 'xyz': [0, 0, 0]}]),
                              ([3.25, -8, 1], [{'chain': 'A', 'xyz': [3.25, -8, 1]}])):
            d = fixture(); set_centroids(d, xyz)
            self.assertEqual(self.parse(d)['centroids'], expected)

    def test_centroids_survive_excluding_whole_chain(self):
        d = fixture(); set_centroids(d, [0, 1, 2])
        for key in ('_atom_site.label_asym_id', '_atom_site.auth_asym_id'):
            d[key][4:] = ['B'] * 4
        d['_atom_site.occupancy'][4:] = ['0'] * 4
        parsed = self.parse(d)
        self.assertEqual({a['chain'] for a in parsed['atoms']}, {'A'})
        self.assertEqual(parsed['centroids'], [{'chain': 'A', 'xyz': [0, 1, 2]},
                                              {'chain': 'B', 'xyz': [0, 1, 2]}])

    def test_centroids_reject_partial_conflicting_and_infinite_targets(self):
        for value in ('nan', 'inf', '9'):
            d = fixture(); set_centroids(d, [0, 0, 0])
            d['_condition_ChainCentroid_1_atom.val0'][0] = value
            with self.subTest(value=value), self.assertRaises(fast.UnsupportedCIF):
                self.parse(d)
        d = fixture(); set_centroids(d, [0, 0, 0])
        for axis in range(3):
            d['_condition_ChainCentroid_1_atom.val' + str(axis)][0] = 'nan'
        with self.assertRaisesRegex(fast.UnsupportedCIF, 'broadcast'):
            self.parse(d)
        d = fixture(); set_centroids(d, [0, 0, 0])
        del d['_condition_ChainCentroid_1_atom.idx0']
        with self.assertRaisesRegex(fast.UnsupportedCIF, 'serialization'):
            self.parse(d)

    def test_centroid_placeholder_metadata_and_conflicts_are_not_filtered(self):
        self.assertEqual(fast.centroid_targets(['A', 'B'], [[1, 2, 3], [0, 0, 0]]),
                         [{'chain': 'A', 'xyz': [1, 2, 3]}, {'chain': 'B', 'xyz': [0, 0, 0]}])
        with self.assertRaisesRegex(ValueError, 'broadcast'):
            fast.centroid_targets(['A', 'A'], [[1, 2, 3], [0, 0, 0]])


@unittest.skipUnless(os.environ.get('RFD4_FAST_TEST_SOURCES'), 'Set RFD4_FAST_TEST_SOURCES for real native parity')
class NativeParityTests(unittest.TestCase):
    def test_sources_and_internal_groups_match_native(self):
        native = cache.NativeReader(cache.reader_config())
        try:
            for name in os.environ['RFD4_FAST_TEST_SOURCES'].split(os.pathsep):
                source = Path(name)
                original_hash = hashlib.sha256(source.read_bytes()).hexdigest()
                with tempfile.TemporaryDirectory() as directory:
                    d = MMCIF2Dict(str(source))
                    import math
                    fixed = {int(i) for i, x in zip(d['_condition_coordinate_1_atom.idx0'], d['_condition_coordinate_1_atom.val0']) if math.isfinite(float(x))}
                    groups = ['0' if i in fixed else '-1' for i in range(len(d['_atom_site.id']))]
                    chosen = next(i for i in range(len(groups)) if i not in fixed and float(d['_atom_site.occupancy'][i]) > 0)
                    groups[chosen] = '9'
                    d['_atom_site.condition_CoordinateGroup_1_atom'] = groups
                    mixed = Path(directory) / 'internal.cif'; write_dict(mixed, d)
                    variants = [source, mixed]
                    for label, xyz in (('zero_centroid', [0, 0, 0]), ('no_centroid', [float('nan')] * 3)):
                        variant = deepcopy(d); set_centroids(variant, xyz)
                        path = Path(directory) / (label + '.cif'); write_dict(path, variant)
                        variants.append(path)
                    variant = deepcopy(d)
                    for key in list(variant):
                        if key.startswith('_condition_ChainCentroid_1_atom.'):
                            del variant[key]
                    absent = Path(directory) / 'absent_centroid.cif'; write_dict(absent, variant)
                    variants.append(absent)
                    variant = deepcopy(d)
                    first_chain = variant['_atom_site.label_asym_id'][0]
                    variant['_atom_site.occupancy'] = ['0' if chain == first_chain else occupancy
                        for chain, occupancy in zip(variant['_atom_site.label_asym_id'], variant['_atom_site.occupancy'])]
                    excluded = Path(directory) / 'excluded_chain.cif'; write_dict(excluded, variant)
                    variants.append(excluded)
                    for path in variants:
                        with self.subTest(source=str(path)):
                            actual, expected = fast.export_model(path), native.model(path)
                            actual['bonds'].sort(); expected['bonds'].sort()
                            self.assertEqual(actual, expected)
                            cache.validate_model(actual)
                            if path == excluded:
                                self.assertNotIn(first_chain, {a['chain'] for a in actual['atoms']})
                                self.assertIn(first_chain, {target['chain'] for target in actual['centroids']})
                    for label, value in (('partial_centroid', 'nan'), ('conflicting_centroid', '1')):
                        variant = deepcopy(d); set_centroids(variant, [0, 0, 0])
                        variant['_condition_ChainCentroid_1_atom.val0'][0] = value
                        path = Path(directory) / (label + '.cif'); write_dict(path, variant)
                        with self.subTest(malformed=label):
                            with self.assertRaises(ValueError):
                                fast.export_model(path)
                            with self.assertRaises(RuntimeError):
                                native.model(path)
                self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), original_hash)
        finally:
            native.close()


if __name__ == '__main__':
    unittest.main()
