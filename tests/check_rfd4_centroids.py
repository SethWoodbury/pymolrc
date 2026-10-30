"""Headless centroid-marker checks; RFD4_SEMANTICS_SOURCE supplies an annotated CIF.

pymol -cq -r tests/check_rfd4_centroids.py
No private fixture is bundled; temporary variants never change the source.
"""
import hashlib
import os
from pathlib import Path
import tempfile

from pymol import cmd


def run_checks():
    source = Path(os.environ['RFD4_SEMANTICS_SOURCE'])
    original_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    cmd.load(str(source), 'centroid_motif')
    data = cmd._rfd4_viewer_state['objects']['centroid_motif']
    assert data['centroids'], 'Use a source with centroid conditioning'
    indices = data['centroid_indices']
    model = cmd.get_model('centroid_motif')
    assert len(model.atom) == len(data['atoms']) + len(indices)
    assert len(model.bond) == len(data['bonds'])
    markers = {a.index: a for a in model.atom if a.index in indices}
    labels, colors = {}, {}
    cmd.iterate('centroid_motif_centroids', 'labels[index] = label; colors[index] = color',
                space={'labels': labels, 'colors': colors})
    assert len(markers) == len(indices) == len(data['centroids'])
    for index, target in zip(indices, data['centroids']):
        marker = markers[index]
        assert marker.name == 'ORI' and marker.symbol == 'X'
        assert labels[index] == 'Centroid ' + target['chain']
        assert all(abs(x - y) < 1e-5 for x, y in zip(marker.coord, target['xyz']))
        assert colors[index] == cmd.get_color_index('hotpink')
    for bond in model.bond:
        assert all(model.atom[i].index not in markers for i in bond.index), 'Marker must have no physical bond'
    assert cmd.count_atoms('centroid_motif_centroids and rep spheres') == len(indices)
    assert cmd.count_atoms('centroid_motif_centroids and rep sticks') == 0
    assert cmd.get_object_list() == ['centroid_motif'], 'Markers must not become separate browse targets'
    toggle = cmd.keyword['rfd4_centroids'][0]
    style = cmd.keyword['rfd4_style'][0]
    toggle('centroid_motif', 0)
    style('centroid_motif', 'elements', labels=0)
    assert cmd.count_atoms('centroid_motif_centroids and rep spheres') == 0
    toggle('centroid_motif', 1)
    style('centroid_motif', 'coordinates', labels=0, context=0)
    assert cmd.count_atoms('centroid_motif_centroids and rep spheres') == len(indices)
    cmd.translate([1, 2, 3], 'centroid_motif', camera=0)
    moved = {a.index: list(a.coord) for a in cmd.get_model('centroid_motif').atom}
    style('centroid_motif', 'coordinates', labels=0)
    assert {a.index: list(a.coord) for a in cmd.get_model('centroid_motif').atom} == moved
    for index, target in zip(indices, data['centroids']):
        assert all(abs(x - y - d) < 1e-5 for x, y, d in zip(moved[index], target['xyz'], [1, 2, 3]))
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.PDB.mmcifio import MMCIFIO
    with tempfile.TemporaryDirectory(prefix='rfd4-centroids-') as directory:
        directory = Path(directory)
        session = directory / 'markers.pse'
        cmd.save(str(session))
        cmd.delete('all')
        cmd.load(str(session))
        assert {a.index: list(a.coord) for a in cmd.get_model('centroid_motif').atom} == moved
        assert cmd.count_atoms('centroid_motif_centroids and rep spheres') == len(indices)
        cif = MMCIF2Dict(str(source))
        prefix = '_condition_ChainCentroid_1_atom.'
        for axis in range(3):
            cif[prefix + 'val' + str(axis)] = ['0'] * len(cif[prefix + 'idx0'])
        writer = MMCIFIO(); writer.set_dict(cif)
        zero = directory / 'zero.cif'; writer.save(str(zero))
        cmd.load(str(zero), 'at_origin')
        zero_data = cmd._rfd4_viewer_state['objects']['at_origin']
        assert len(zero_data['centroids']) == len(data['centroids'])
        assert all(t['xyz'] == [0, 0, 0] for t in zero_data['centroids'])
        assert cmd.count_atoms('at_origin_centroids and rep spheres') == len(indices)
        for key in list(cif):
            if key.startswith(prefix):
                del cif[key]
        writer.set_dict(cif)
        none = directory / 'none.cif'; writer.save(str(none))
        cmd.load(str(none), 'no_centroids')
        none_data = cmd._rfd4_viewer_state['objects']['no_centroids']
        assert none_data['centroids'] == [] and none_data['centroid_indices'] == []
        assert cmd.count_atoms('no_centroids') == len(none_data['atoms'])
        toggle('no_centroids', 1)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == original_hash
    print('RFD4_CENTROID_TEST_PASS', flush=True)


try:
    run_checks()
except Exception:
    import traceback
    traceback.print_exc()
    cmd.quit(1)
else:
    cmd.quit(0)
