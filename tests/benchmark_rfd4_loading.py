"""Headless loading timings (exclude process/GUI startup and rendering).

RFD4_TIMING_SOURCE=/path/input.cif RFD4_TIMING_PML=/path/view.pml \
  pymol -cq -r tests/benchmark_rfd4_loading.py
Optional RFD4_TIMING_OUTPUT writes JSON; source files are never modified.
Explicit RFD4_REPO forces the native reader/cache for comparison.
"""
import hashlib
import json
import os
from pathlib import Path
import time
from pymol import cmd

source = Path(os.environ['RFD4_TIMING_SOURCE']).resolve(strict=True)
pml = Path(os.environ['RFD4_TIMING_PML']).resolve(strict=True) if os.environ.get('RFD4_TIMING_PML') else None
results = []
for mode in (('annotated', 'pml', 'combined') if pml else ('annotated',)):
    for repeat in range(3):
        cmd.delete('all')
        started = time.perf_counter()
        if mode in ('annotated', 'combined'):
            cmd.load(str(source), 'input_motif')
        if mode in ('pml', 'combined'):
            cmd.do('@' + str(pml))
        result = {'mode': mode, 'repeat': repeat, 'seconds': time.perf_counter() - started,
                  'object_count': len(cmd.get_object_list()),
                  'reader_route': cmd._rfd4_viewer_state['cache_stats']['last_route'] if mode != 'pml' else None}
        results.append(result)
        print('RFD4_LOADING_TIMING', json.dumps(result), flush=True)
payload = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
           'pml': str(pml) if pml else None, 'results': results,
           'scope': 'Headless PyMOL loading; excludes process/desktop startup and rendering.'}
if os.environ.get('RFD4_TIMING_OUTPUT'):
    Path(os.environ['RFD4_TIMING_OUTPUT']).write_text(json.dumps(payload, indent=2) + '\n')
cmd.quit()
