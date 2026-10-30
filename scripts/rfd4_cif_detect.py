"""Lightweight routing of CIFs; no scientific environment or PyMOL needed."""
from pathlib import Path
import re
from collections import OrderedDict


_CONDITION_FIELD = re.compile(
    r"^(?:_condition_[^\s.]+\.[^\s]+|_atom_site\.condition_[^\s]+)(?:\s|$)",
    re.IGNORECASE,
)
_RESULTS = OrderedDict()
_MARKER = re.compile(br"_condition_|_atom_site\.condition_", re.IGNORECASE)


def _has_marker(handle):
    """Scan bytes in C, avoiding a Python iteration over every atom row."""
    tail = b""
    while chunk := handle.read(256 * 1024):
        chunk = tail + chunk
        if _MARKER.search(chunk):
            return True
        tail = chunk[-32:]
    return False


def _declares_condition(handle):
    # CIF unquoted underscore tokens are data names. Quoted and semicolon text
    # values are opaque; unlike a line-prefix test, inline declarations work.
    in_text = False
    for raw in handle:
        line = raw.decode("utf-8", errors="replace")
        if line.startswith(";"):
            in_text = not in_text
            continue
        if in_text:
            continue
        i = 0
        while i < len(line):
            if line[i].isspace():
                i += 1
                continue
            if line[i] == "#":
                break
            if line[i] in "\"'":
                quote = line[i]
                i += 1
                while i < len(line):
                    if line[i] == quote and (i + 1 == len(line) or line[i + 1].isspace()):
                        i += 1
                        break
                    i += 1
                continue
            end = i + 1
            while end < len(line) and not line[end].isspace():
                end += 1
            if _CONDITION_FIELD.fullmatch(line[i:end]):
                return True
            i = end
    return False


def is_rfd4_cif(filename):
    """Detect explicit RFD4 condition fields, not residual bookkeeping columns.

    Saved predictions also carry annotation_expsegmin/max columns. Their
    presence alone must never send physical structures through the input
    reader, which intentionally hides scaffold and unoccupied atoms.
    """
    try:
        source = Path(str(filename).strip('"\'')).expanduser()
        if source.suffix.lower() not in {".cif", ".mmcif"}:
            return False
        stat = source.stat()
        stamp = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        key = (str(source.absolute()), stamp)
        if key in _RESULTS:
            _RESULTS.move_to_end(key)
            return _RESULTS[key]
        with source.open("rb") as handle:
            result = False
            if _has_marker(handle):
                handle.seek(0)
                result = _declares_condition(handle)
        after = source.stat()
        if stamp == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            _RESULTS[key] = result
            while len(_RESULTS) > 2048:
                _RESULTS.popitem(last=False)
        return result
    except (OSError, ValueError):
        return False
