"""Run inside PyMOL: RFD4_VIEW_TEST_DIR=/tmp/inputs pymol -cq -r this_file.

Uses temporary copies of supplied CIFs. Checks CLI inputs (if provided), native
geometry/bonds, normal PDB/CIF passthrough, repeated helper runs and batch browsing.
"""
import hashlib
import math
import os
from pathlib import Path
import shutil
import tempfile
import time

from pymol import cmd


def check_model(obj):
    model = cmd.get_model(obj)
    data = cmd._rfd4_viewer_state['objects'][obj]
    expected = {(a['chain'], a['resi'], a['resn'], a['name']): a for a in data['atoms']}
    actual = {(a.chain, a.resi, a.resn, a.name): a for a in model.atom}
    assert expected.keys() == actual.keys()
    for key, atom in actual.items():
        assert atom.coord == expected[key]['coord'], (key, atom.coord, expected[key]['coord'])
        assert all(math.isfinite(v) for v in atom.coord)
        assert atom.formal_charge == expected[key]['charge']
    native_keys = list(expected)
    actual_keys = [(a.chain, a.resi, a.resn, a.name) for a in model.atom]
    expected_bonds = {(frozenset((native_keys[i], native_keys[j])), order) for i, j, order in data['bonds']}
    actual_bonds = {(frozenset((actual_keys[b.index[0]], actual_keys[b.index[1]])), b.order) for b in model.bond}
    assert expected_bonds == actual_bonds
    assert all(math.isfinite(v) for corner in cmd.get_extent(obj) for v in corner)
    assert cmd.count_atoms('model ' + obj + ' and rep sticks') == len(model.atom)


def run_checks():
    sources = sorted(Path(os.environ['RFD4_VIEW_TEST_DIR']).glob('*.cif'))
    assert sources
    hashes = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    # With a CIF passed directly on the PyMOL CLI, this proves command-table dispatch.
    for obj in cmd.get_object_list():
        check_model(obj)
    cmd.delete('all')
    cmd.load(str(sources[0]), 'motif')
    check_model('motif')
    worker = cmd._rfd4_viewer_state['worker']
    worker_pid = worker.pid if worker is not None else None
    if os.environ.get('RFD4_EXPECT_DISK_CACHE') == '1':
        assert worker is None, 'Prewarmed inputs must not start a native reader.'
        assert cmd._rfd4_viewer_state['cache_stats']['disk_hits'] >= 1
    atom_count = cmd.count_atoms('motif')
    cmd.run(str(Path.home() / 'pymolrc/rfd4_viewer.py'))
    assert 'motif' in cmd._rfd4_viewer_state['objects']
    worker = cmd._rfd4_viewer_state['worker']
    assert (worker.pid if worker is not None else None) == worker_pid
    for mode in ('coordinates', 'sequence', 'hbond', 'burial', 'unindexed', 'elements'):
        cmd.keyword['rfd4_style'][0]('motif', mode, labels=0)
    # An existing object is preserved, never overwritten or extended.
    cmd.load(str(sources[0]), 'motif')
    assert cmd.count_atoms('motif') == atom_count
    assert len(cmd.get_object_list()) == 2
    with tempfile.TemporaryDirectory(prefix='rfd4-pymol-test-') as directory:
        temp = Path(directory)
        for suffix in ('pdb', 'cif'):
            normal = temp / ('ordinary.' + suffix)
            cmd.save(str(normal), 'motif')
            cmd.load(str(normal), 'normal_' + suffix)
            assert cmd.count_atoms('normal_' + suffix) == atom_count
            assert 'normal_' + suffix not in cmd._rfd4_viewer_state['objects']
        paths = []
        for i in range(10):
            target = temp / ('sample_%02d.cif' % i)
            shutil.copy2(sources[i % len(sources)], target)
            paths.append(str(target))
        listing = temp / 'list.txt'
        listing.write_text('\n'.join(paths) + '\n')
        cmd.delete('all')
        started = time.monotonic()
        cmd.keyword['browse_list'][0](str(listing))
        objects = cmd.get_object_list()
        assert len(objects) == 10, objects
        assert len(set(objects) & set(cmd.get_names('objects', enabled_only=1))) == 1
        for obj in objects:
            check_model(obj)
        worker = cmd._rfd4_viewer_state['worker']
        assert (worker.pid if worker is not None else None) == worker_pid
        cmd.keyword['browse_next'][0]()
        assert len(set(objects) & set(cmd.get_names('objects', enabled_only=1))) == 1
        # No long object-name typing needed for optional coloring of visible motif.
        cmd.keyword['rfd4_style'][0]('', 'hbond', labels=0)
        print('RFD4_AUTOLOAD_BATCH_SECONDS', round(time.monotonic() - started, 2), flush=True)
    assert all(hashlib.sha256(path.read_bytes()).hexdigest() == digest for path, digest in hashes.items())
    if os.environ.get('RFD4_EXPECT_DISK_CACHE') == '1':
        assert cmd._rfd4_viewer_state['cache_stats']['native_parses'] == 0
        assert cmd._rfd4_viewer_state['worker'] is None
    cmd.keyword['rfd4_cache_info'][0]()
    print('RFD4_AUTOLOAD_TEST_PASS: CLI, ordinary PDB/CIF, exact native coordinates/bonds/charges, '
          'placeholder exclusion, reload, collision protection, 10-file browse, source hashes', flush=True)


try:
    run_checks()
except Exception:
    import traceback
    traceback.print_exc()
    (Path(os.environ['RFD4_VIEW_TEST_DIR']) / 'result.txt').write_text(traceback.format_exc())
    __import__('sys').stderr.flush()
    cmd.quit(1)
else:
    (Path(os.environ['RFD4_VIEW_TEST_DIR']) / 'result.txt').write_text('PASS')
    cmd.quit(0)
