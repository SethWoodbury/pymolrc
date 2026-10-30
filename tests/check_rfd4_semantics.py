"""Headless PyMOL regression: source path via RFD4_SEMANTICS_SOURCE.

pymol -cq -r tests/check_rfd4_semantics.py
Checks scientific display semantics, native internal groups, session restoration,
ordinary raw/PML styling and camera-preserving navigation. Never edits the input.
"""
from pathlib import Path
import hashlib
import os
import tempfile
from pymol import cmd


def appearance(selection):
    result = []
    cmd.iterate(selection, 'out.append((name, color, reps))', space={'out': result})
    return result


def run_checks():
    source = Path(os.environ['RFD4_SEMANTICS_SOURCE'])
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    cmd.load(str(source), 'motif')
    data = cmd._rfd4_viewer_state['objects']['motif']
    absolute = sum(a['fixed'] for a in data['atoms'])
    internal = sum(a['internal'] for a in data['atoms'])
    assert cmd.count_atoms('motif and rep sticks') == absolute + internal
    assert cmd.count_atoms('motif_reference and rep sticks') == 0
    assert cmd.count_atoms('motif_reference and rep lines') == len(data['atoms']) - absolute - internal
    assert cmd.count_atoms('motif_reference and rep spheres') == 0
    fixed_residues = {(a['chain'], a['resi'], a['resn']) for a in data['atoms'] if a['sequence_fixed']}
    fixed_ids = [str(a['pymol_index']) for a in data['atoms'] if (a['chain'], a['resi'], a['resn']) in fixed_residues]
    cmd.select('test_identity_fixed', 'motif and index ' + '+'.join(fixed_ids))
    for selection, color in [('motif_absolute and test_identity_fixed and polymer.protein and elem C', 'cyan'),
                             ('motif_absolute and not polymer.protein and elem C', 'lig_amethyst3'),
                             ('motif_absolute and elem N', 'elemN'),
                             ('motif_absolute and elem O', 'elemO')]:
        assert cmd.count_atoms(selection), selection
        assert all(c == cmd.get_color_index(color) for _, c, _ in appearance(selection)), selection
    free_carbons = 'motif and polymer.protein and elem C and not test_identity_fixed'
    assert cmd.count_atoms(free_carbons), 'Fixture needs free identities as well as fixed ones'
    expected_tint = tuple(x + 0.30 * (1.0 - x) for x in cmd.get_color_tuple('cyan'))
    for _, color, _ in appearance(free_carbons):
        assert all(abs(x - y) < 1e-6 for x, y in zip(cmd.get_color_tuple(color), expected_tint))
    assert cmd.count_atoms('motif_absolute and not test_identity_fixed and rep sticks') > 0
    # A free identity can still have a constrained backbone; oxygen stays its
    # element color in both constrained and unconstrained atoms.
    assert all(c == cmd.get_color_index('elemO') for _, c, _ in appearance('motif and elem O'))
    # A pure annotation view remains available when element colors obscure a
    # condition on a heteroatom, especially a metal-only internal group.
    cmd.keyword['rfd4_style'][0]('motif', 'coordinates', labels=0, element_colors=0)
    assert all(color == cmd.get_color_index('cyan') for _, color, _ in appearance('motif_absolute'))
    cmd.keyword['rfd4_style'][0]('motif', 'coordinates', labels=0)
    cmd.keyword['rfd4_style'][0]('motif', 'coordinates', labels=1, context=0)
    hidden_labels = []
    cmd.iterate('motif_reference', 'out.append(label)', space={'out': hidden_labels})
    assert not any(hidden_labels), hidden_labels
    cmd.keyword['rfd4_style'][0]('motif', 'coordinates', labels=0)
    original = appearance('motif')
    cmd.keyword['style_all'][0]()
    assert appearance('motif') == original
    cmd.keyword['rfd4_style'][0]('motif', 'coordinates', labels=0, context=1)
    assert cmd.count_atoms('motif_reference and rep lines') == len(data['atoms']) - absolute - internal
    assert cmd.count_atoms('motif_reference and rep sticks') == 0
    assert all(c == cmd.get_color_index('elemO') for _, c, _ in appearance('motif_reference and elem O'))
    cmd.keyword['rfd4_style'][0]('motif', 'elements', labels=0)
    assert cmd.count_atoms('motif and rep sticks') == absolute + internal
    assert cmd.count_atoms('motif_reference and rep lines') == len(data['atoms']) - absolute - internal
    cmd.keyword['rfd4_style'][0]('motif', 'elements', labels=1, context=0)
    assert cmd.count_atoms('motif_reference and (rep lines or rep nonbonded or rep spheres)') == 0
    hidden_labels = []
    cmd.iterate('motif_reference', 'out.append(label)', space={'out': hidden_labels})
    assert not any(hidden_labels)
    with tempfile.TemporaryDirectory(prefix='rfd4-semantics-') as directory:
        temp = Path(directory)
        # A raw PML-owned object deliberately has no startup autostyle metadata.
        cmd.load_raw('HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00 10.00           O\nEND\n', 'pdb', 'raw_review')
        cmd.color('magenta', 'raw_review')
        cmd.hide('everything', 'raw_review')
        cmd.show('spheres', 'raw_review')
        raw = appearance('raw_review')
        cmd.load(str(source), 'another_motif')
        assert appearance('raw_review') == raw, 'Later loads must not restyle existing PML-owned objects'
        camera = cmd.get_view()
        cmd.keyword['structure_next'][0]()
        assert appearance('raw_review') == raw
        assert cmd.get_view() == camera
        session = temp / 'style.pse'
        cmd.save(str(session))
        before = {o: appearance(o) for o in cmd.get_object_list()}
        enabled = cmd.get_names('objects', enabled_only=1)
        cmd.delete('all')
        cmd.load(str(session))
        assert {o: appearance(o) for o in cmd.get_object_list()} == before
        assert cmd.get_names('objects', enabled_only=1) == enabled
        assert cmd.get_view() == camera
        cmd.keyword['structure_next'][0]()
        assert {o: appearance(o) for o in cmd.get_object_list()} == before
        # Native parser resolves an explicit coordinate group independently of
        # sequence constraints and absent internal-coordinate payloads.
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        from Bio.PDB.mmcifio import MMCIFIO
        cif = MMCIF2Dict(str(source))
        n = len(cif['_atom_site.id'])
        import math
        fixed_indices = {int(i) for i, x in zip(cif['_condition_coordinate_1_atom.idx0'], cif['_condition_coordinate_1_atom.val0']) if math.isfinite(float(x))}
        groups = ['0' if i in fixed_indices else '-1' for i in range(n)]
        chosen = next(i for i in range(n) if i not in fixed_indices and cif['_atom_site.occupancy'][i] == '1.0')
        groups[chosen] = '7'
        cif['_atom_site.condition_CoordinateGroup_1_atom'] = groups
        writer = MMCIFIO(); writer.set_dict(cif)
        mixed = temp / 'internal.cif'; writer.save(str(mixed))
        cmd.load(str(mixed), 'internal_test')
        assert cmd.count_atoms('internal_test_internal') == 1
        assert cmd.count_atoms('internal_test_absolute') == absolute
        cmd.keyword['rfd4_style'][0]('internal_test', 'coordinates', labels=0, element_colors=0)
        assert all(color == cmd.get_color_index('orange') for _, color, _ in appearance('internal_test_internal'))
        cmd.keyword['rfd4_style'][0]('internal_test', 'coordinates', labels=0)
        assert cmd.count_atoms('internal_test_internal and rep sticks') == 1
        # One alternate-location annotation is outside the lightweight contract.
        # It must use the native reader/cache without losing atoms or conditions.
        cif['_atom_site.label_alt_id'][chosen] = 'A'
        writer.set_dict(cif)
        fallback = temp / 'native_fallback.cif'; writer.save(str(fallback))
        cmd.load(str(fallback), 'fallback_test')
        assert cmd._rfd4_viewer_state['cache_stats']['last_route'] in {'native_parses', 'disk_hits', 'memory_hits'}
        assert cmd.count_atoms('fallback_test') == len(data['atoms']) + len(data['centroids'])
        assert cmd.count_atoms('fallback_test_absolute') == absolute
        assert cmd.count_atoms('fallback_test_internal') == 1
    # Synthetic metal-only instructions test display behavior without changing
    # any campaign source or requiring an unrelated metal-containing campaign.
    from chempy import Atom
    from chempy.models import Indexed
    model = Indexed()
    for i, element in enumerate(('Zn', 'Fe', 'Mg'), 1):
        atom = Atom(); atom.name = element.upper(); atom.symbol = element
        atom.resn = element.upper(); atom.resi = str(i); atom.chain = 'Z'
        atom.coord = [float(i * 3), 0., 0.]; atom.hetatm = 1
        model.atom.append(atom)
    cmd.load_model(model, 'metal_conditions', zoom=0)
    records = []
    for atom in cmd.get_model('metal_conditions').atom:
        record = dict(data['atoms'][0])
        record.update(name=atom.name, element=atom.symbol, resn=atom.resn, resi=atom.resi,
                      chain=atom.chain, coord=atom.coord, hetero=True, pymol_index=atom.index,
                      fixed=atom.symbol == 'Zn', internal=atom.symbol == 'Mg')
        records.append(record)
    cmd._rfd4_viewer_state['objects']['metal_conditions'] = {'atoms': records, 'centroids': [], 'centroid_indices': []}
    normal_scale = cmd.get_setting_float('sphere_scale')
    cmd.set('sphere_scale', 0.91)
    style = cmd.keyword['rfd4_style'][0]
    style('metal_conditions', 'coordinates', labels=0)
    assert cmd.count_atoms('metal_conditions and rep spheres') == 2
    assert cmd.count_atoms('metal_conditions and elem Fe and rep spheres') == 0
    assert cmd.count_atoms('metal_conditions and elem Fe and rep nonbonded') == 1
    sizes = []
    cmd.iterate('metal_conditions and (elem Zn or elem Mg)', 'out.append(s.sphere_scale)', space={'out': sizes})
    assert all(abs(x - 0.91) < 1e-6 for x in sizes), sizes
    cmd.set('sphere_scale', 0.83)
    sizes = []
    cmd.iterate('metal_conditions and (elem Zn or elem Mg)', 'out.append(s.sphere_scale)', space={'out': sizes})
    assert all(abs(x - 0.83) < 1e-6 for x in sizes), sizes
    assert all(c == cmd.get_color_index('elemZn') for _, c, _ in appearance('metal_conditions_absolute or metal_conditions_internal'))
    style('metal_conditions', 'coordinates', labels=0, element_colors=0)
    assert all(c == cmd.get_color_index('orange') for _, c, _ in appearance('metal_conditions_internal'))
    style('metal_conditions', 'elements', labels=0)
    assert cmd.count_atoms('metal_conditions and rep spheres') == 2
    style('metal_conditions', 'elements', labels=0, context=1)
    assert cmd.count_atoms('metal_conditions and rep spheres') == 2
    assert cmd.count_atoms('metal_conditions_reference and rep nonbonded') == 1
    style('metal_conditions', 'elements', labels=0, context=0)
    assert cmd.count_atoms('metal_conditions_reference and (rep nonbonded or rep lines or rep spheres)') == 0
    cmd.set('sphere_scale', normal_scale)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    print('RFD4_SEMANTICS_TEST_PASS', flush=True)


try:
    run_checks()
except Exception:
    import traceback
    traceback.print_exc()
    cmd.quit(1)
else:
    cmd.quit(0)
