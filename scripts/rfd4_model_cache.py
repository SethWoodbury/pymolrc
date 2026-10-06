"""Fast, read-only RFD4 display models backed by a verified JSON cache.

This module deliberately uses only the standard library. Native scientific
imports happen only in a persistent CPU subprocess on a cache miss. The CIF is
always the authoritative input; these disposable files are visualization data.
"""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import select
import subprocess
import tempfile
import threading
import time

SCHEMA = 1
MAX_CACHE_BYTES = 128 * 1024 * 1024
_DIGESTS = {}
_DIGEST_LOCK = threading.Lock()
_DISK_DIGESTS = {}
_DIGESTS_LOADED = False
_DIGESTS_CHANGED = False


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def cache_root():
    return Path(os.environ.get("RFD4_PYMOL_CACHE", str(Path.home() / ".cache/pymol/rfd4"))).expanduser()


def reader_config(repo="", python=""):
    repo = Path(repo or os.environ.get("RFD4_REPO", "~/git/RFD4-Proteina-dev")).expanduser().resolve()
    python = Path(python or os.environ.get("RFD4_PYTHON", str(repo / ".venv/bin/python"))).expanduser().absolute()
    if not python.is_file():
        raise ValueError("RFD4 Python not found: %s; set RFD4_PYTHON." % python)
    return {"repo": str(repo), "python": str(python), "cluster": os.environ.get("RFD4_CLUSTER", "digs"),
            "exporter": str(Path(__file__).with_name("rfd4_cif_export.py")),
            "preset": "ANNOTATED_CIF / asymmetric unit / model 1 / first altloc",
            "environment": {key: value for key, value in sorted(os.environ.items()) if
                            key.startswith(("ATOMWORKS_", "CCD_", "BIOTITE_")) or
                            key in {"ALLOW_BIOTITE_CCD", "PROJECT_ROOT", "PYTHONPATH"}}}


def _signature(path):
    stat = path.stat()
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def _code_digest(path):
    global _DIGESTS_CHANGED
    stamp = _signature(path)
    key = (str(path), *stamp)
    with _DIGEST_LOCK:
        digest = _DIGESTS.get(key)
        previous = _DISK_DIGESTS.get(str(path), {})
        if not isinstance(previous, dict):
            previous = {}
        if digest is None and previous.get("stamp") == stamp:
            digest = previous.get("sha256")
            if not isinstance(digest, str) or len(digest) != 64:
                digest = None
    if digest is None:
        digest = sha256(path.read_bytes())
        with _DIGEST_LOCK:
            _DISK_DIGESTS[str(path)] = {"stamp": stamp, "sha256": digest}
            _DIGESTS_CHANGED = True
    with _DIGEST_LOCK:
        _DIGESTS[key] = digest
    return [str(path), digest]


def provenance(config):
    """Fingerprint real parser/condition code, package versions, and local config.

    Code is content-hashed (including dirty edits); within this Python process
    unchanged inode/ctime/mtime/size entries reuse their hashes. Package data and
    compiled libraries carry filesystem identity/change stamps. No git command,
    native Python startup, scientific import, or model loading is required.
    """
    global _DIGESTS_LOADED, _DIGESTS_CHANGED
    digest_file = cache_root() / "code-digests.json"
    if not _DIGESTS_LOADED:
        try:
            payload = json.loads(digest_file.read_bytes())
            if isinstance(payload, dict) and isinstance(payload.get("files"), dict) and sha256(canonical(payload["files"])) == payload.get("sha256"):
                _DISK_DIGESTS.update(payload["files"])
        except (OSError, ValueError, TypeError):
            pass
        _DIGESTS_LOADED = True
    repo, python = Path(config["repo"]), Path(config["python"])
    code = {Path(__file__), Path(config["exporter"])}
    state = {python, python.resolve()}
    roots = [repo / "src/rfproteina", repo / "lib/atomworks/src/atomworks"]
    # Editable checkout code and installed parsing dependencies are both tracked.
    site_dirs = sorted((python.parent.parent / "lib").glob("python*/site-packages"))
    for site in site_dirs:
        roots.extend(site / name for name in ("biotite", "numpy", "atomworks", "rfproteina"))
        for pattern in ("biotite-*.dist-info", "numpy-*.dist-info", "atomworks-*.dist-info", "rfproteina-*.dist-info"):
            for directory in site.glob(pattern):
                code.update(directory / name for name in ("METADATA", "RECORD", "direct_url.json") if (directory / name).is_file())
        code.update(site.glob("*atomworks*.pth"))
        code.update(site.glob("*rfproteina*.pth"))
    for root in roots:
        if not root.is_dir():
            continue
        for directory, names, filenames in os.walk(root):
            names[:] = [name for name in names if name not in {"__pycache__", ".git", "tests"}]
            for filename in filenames:
                path = Path(directory) / filename
                if path.suffix in {".py", ".pyi", ".json"}:
                    code.add(path)
                elif path.suffix in {".so", ".bcif", ".cif", ".npz", ".npy"}:
                    state.add(path)
    # Hash only: never include configuration-file contents in the cache/log.
    code.update(repo.glob(".env*"))
    code.update((repo / ".ipd").glob(".env*"))
    for candidate in (repo / "uv.lock", repo / "pyproject.toml", repo / "lib/atomworks/pyproject.toml", python.parent.parent / "pyvenv.cfg"):
        if candidate.is_file():
            code.add(candidate)
    evidence = {"schema": SCHEMA, "config": config,
                "code": [_code_digest(path) for path in sorted(code) if path.is_file()],
                "runtime_data": [[str(path), _signature(path)] for path in sorted(state)]}
    if _DIGESTS_CHANGED:
        try:
            _atomic_json(digest_file, {"files": _DISK_DIGESTS, "sha256": sha256(canonical(_DISK_DIGESTS))})
            _DIGESTS_CHANGED = False
        except OSError:
            pass  # A read-only cache must never prevent native loading.
    return {"sha256": sha256(canonical(evidence)), "evidence": evidence}


def validate_model(model):
    """Reject malformed/corrupted models before they can reach PyMOL."""
    if not isinstance(model, dict) or not isinstance(model.get("atoms"), list) or not model["atoms"]:
        raise ValueError("Cached model must contain physical atoms.")
    seen = set()
    for atom in model["atoms"]:
        if not isinstance(atom, dict):
            raise ValueError("Invalid atom record.")
        for name in ("chain", "resi", "resn", "name", "element"):
            if not isinstance(atom.get(name), str):
                raise ValueError("Invalid atom identity field: " + name)
        identity = tuple(atom[name] for name in ("chain", "resi", "resn", "name"))
        if identity in seen:
            raise ValueError("Duplicate physical atom identity.")
        seen.add(identity)
        if not isinstance(atom.get("coord"), list) or len(atom["coord"]) != 3 or not all(
                type(v) in (int, float) and math.isfinite(v) for v in atom["coord"]):
            raise ValueError("Invalid or nonfinite atom coordinates.")
        for name in ("hetero", "fixed", "sequence_fixed", "unindexed"):
            if type(atom.get(name)) is not bool:
                raise ValueError("Invalid boolean annotation: " + name)
        for name in ("charge", "index_group"):
            if type(atom.get(name)) is not int:
                raise ValueError("Invalid integer annotation: " + name)
        occupancy = atom.get("occupancy")
        if type(occupancy) not in (int, float) or not math.isfinite(occupancy) or occupancy <= 0:
            raise ValueError("Invalid occupied-atom occupancy.")
        if atom.get("hbond") not in {"NONE", "DONOR", "ACCEPTOR", "BOTH"} or atom.get("rasa") not in {"UNSPECIFIED", "BURIED", "INTERMEDIATE", "EXPOSED"}:
            raise ValueError("Unknown condition enum.")
    if not isinstance(model.get("bonds"), list):
        raise ValueError("Missing bond graph.")
    bonds = set()
    for bond in model["bonds"]:
        if not isinstance(bond, list) or len(bond) != 3 or not all(type(v) is int for v in bond):
            raise ValueError("Invalid bond record.")
        i, j, order = bond
        if not (0 <= i < len(seen) and 0 <= j < len(seen) and i != j and order in (1, 2, 3, 4)):
            raise ValueError("Invalid bond indices/order.")
        pair = tuple(sorted((i, j)))
        if pair in bonds:
            raise ValueError("Duplicate native bond.")
        bonds.add(pair)
    if type(model.get("excluded_atoms")) is not int or model["excluded_atoms"] < 0:
        raise ValueError("Invalid excluded-atom count.")
    if not isinstance(model.get("segments"), list):
        raise ValueError("Missing segment metadata.")
    for segment in model["segments"]:
        if not isinstance(segment, dict) or not isinstance(segment.get("chain"), str) or not all(type(segment.get(name)) is int for name in ("resi", "min", "max")) or not 0 <= segment["min"] <= segment["max"]:
            raise ValueError("Invalid expandable segment.")
    if not isinstance(model.get("condition_fields"), list) or not all(isinstance(v, str) for v in model["condition_fields"]):
        raise ValueError("Invalid condition-field list.")
    if not isinstance(model.get("present"), dict) or any(type(model["present"].get(name)) is not bool for name in ("fixed", "sequence_fixed", "unindexed", "hbond", "rasa")):
        raise ValueError("Invalid condition-presence metadata.")
    if not isinstance(model.get("parse"), str):
        raise ValueError("Missing parser description.")
    return model


def cache_path(source_hash, reader_hash, directory=None):
    key = sha256(canonical([SCHEMA, source_hash, reader_hash]))
    return Path(directory or cache_root()) / key[:2] / (key + ".json")


def read_cached(source, source_hash, reader_hash, directory=None):
    """Return a verified model, or None on missing, obsolete, or damaged cache."""
    path = cache_path(source_hash, reader_hash, directory)
    try:
        if path.stat().st_size > MAX_CACHE_BYTES:
            return None
        entry = json.loads(path.read_bytes())
        if not isinstance(entry, dict) or entry.get("schema") != SCHEMA or entry.get("source_sha256") != source_hash or entry.get("reader_sha256") != reader_hash:
            return None
        model = entry["model"]
        if sha256(canonical(model)) != entry["model_sha256"]:
            return None
        validate_model(model)
        model["source"] = str(source)  # Same bytes can move/rename without stale UI paths.
        return model
    except (OSError, ValueError, TypeError, KeyError):
        return None


def _atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".pending-", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(canonical(payload))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_cached(source_hash, reader, model, directory=None):
    validate_model(model)
    model = deepcopy(model)
    model.pop("source", None)
    for atom in model["atoms"]:
        atom.pop("pymol_index", None)
    model.pop("selections_created", None)
    directory = Path(directory or cache_root())
    path = cache_path(source_hash, reader["sha256"], directory)
    _atomic_json(path, {"schema": SCHEMA, "source_sha256": source_hash, "reader_sha256": reader["sha256"],
                        "model_sha256": sha256(canonical(model)), "model": model})
    evidence_file = directory / "readers" / (reader["sha256"] + ".json")
    if not evidence_file.exists():
        _atomic_json(evidence_file, reader)
    return path


class NativeReader:
    """One persistent native CPU parser; reusable from PyMOL or batch warmup."""
    def __init__(self, config):
        self.config = config
        self.worker = None

    def close(self):
        if self.worker is not None:
            worker = self.worker
            if worker.poll() is None:
                worker.terminate()
                try:
                    worker.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    worker.kill()
                    worker.wait(timeout=3)
            worker.stdin.close()
            worker.stdout.close()
            worker._rfd4_stderr.close()
            self.worker = None

    def model(self, source):
        if self.worker is None or self.worker.poll() is not None:
            self.close()
            environment = dict(os.environ)
            for name in ("PYTHONHOME", "PYTHONEXECUTABLE", "__PYVENV_LAUNCHER__"):
                environment.pop(name, None)
            environment.update(CLUSTER=self.config["cluster"], CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1",
                               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
            error_log = tempfile.TemporaryFile(mode="w+t")
            self.worker = subprocess.Popen([self.config["python"], self.config["exporter"], "--serve"],
                cwd=self.config["repo"], env=environment, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=error_log, text=True, bufsize=1)
            self.worker._rfd4_stderr = error_log
        worker = self.worker
        try:
            worker.stdin.write(json.dumps({"source": str(source)}) + "\n")
            worker.stdin.flush()
            if not select.select([worker.stdout], [], [], 180)[0]:
                raise RuntimeError("Native RFD4 reader timed out after 180 seconds.")
            line = worker.stdout.readline()
            if not line:
                worker._rfd4_stderr.seek(0)
                raise RuntimeError("Native reader exited:\n" + worker._rfd4_stderr.read()[-6000:])
            response = json.loads(line)
            if not response.get("ok"):
                raise RuntimeError("Native RFD4 parsing failed:\n" + response.get("error", "Unknown failure")[-6000:])
            return validate_model(response["model"])
        except Exception:
            self.close()
            raise


def warm(paths, *, repo="", python="", jobs=1, directory=None, progress=None):
    """Cache independent sources using up to `jobs` persistent CPU readers."""
    if type(jobs) is not int or not 1 <= jobs <= 16:
        raise ValueError("jobs must be an integer from 1 to 16.")
    paths = sorted(set(Path(path).expanduser().resolve(strict=True) for path in paths))
    started = time.monotonic()
    config = reader_config(repo, python)
    reader = provenance(config)
    queue_lock = threading.Lock()
    remaining = iter(paths)
    totals = {"files": len(paths), "hits": 0, "parsed": 0, "reader_sha256": reader["sha256"]}

    def consume():
        native = NativeReader(config)
        try:
            while True:
                with queue_lock:
                    source = next(remaining, None)
                if source is None:
                    return
                source_hash = sha256(source.read_bytes())
                hit = read_cached(source, source_hash, reader["sha256"], directory) is not None
                if not hit:
                    model = native.model(source)
                    if sha256(source.read_bytes()) != source_hash:
                        raise RuntimeError("Source changed during native parsing: " + str(source))
                    write_cached(source_hash, reader, model, directory)
                with queue_lock:
                    totals["hits" if hit else "parsed"] += 1
                    if progress:
                        progress(dict(totals))
        finally:
            native.close()

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(consume) for _ in range(min(jobs, len(paths)))]
        for future in futures:
            future.result()
    if provenance(config)["sha256"] != reader["sha256"]:
        raise RuntimeError("Native reader installation changed while warming; retry with its current fingerprint.")
    totals.update(seconds=round(time.monotonic() - started, 3), directory=str(directory or cache_root()), jobs=jobs)
    return totals


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--source-dir", type=Path, help="Nonrecursive *.cif inputs; originals are never changed.")
    parser.add_argument("--jobs", type=int, default=1, help="Persistent CPU readers (default: 1).")
    parser.add_argument("--repo", default="")
    parser.add_argument("--python", default="")
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()
    paths = args.paths + (list(args.source_dir.glob("*.cif")) if args.source_dir else [])
    if not paths:
        parser.error("Provide CIF paths or --source-dir.")

    def report(totals):
        count = totals["hits"] + totals["parsed"]
        if count % 100 == 0 or count == totals["files"]:
            print("[RFD4 cache] %s/%s (%s cached; %s newly parsed)" %
                  (count, totals["files"], totals["hits"], totals["parsed"]), flush=True)

    print(json.dumps(warm(paths, repo=args.repo, python=args.python, jobs=args.jobs,
                          directory=args.cache_dir, progress=report), sort_keys=True), flush=True)
