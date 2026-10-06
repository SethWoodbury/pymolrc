"""Read an RFD4 annotated CIF with its native parser for the PyMOL viewer.

Run with the RFD4 Python environment. The output is a read-only JSON model;
neither this helper nor the PyMOL commands rewrite the input CIF.
"""
from contextlib import redirect_stdout
import json
from pathlib import Path
import sys
import select
import traceback


def export_model(filename):
    from rfproteina.utils.env import load_dotenv
    load_dotenv()
    import rfproteina.conditions  # Register native condition definitions before parsing.
    import numpy as np
    from atomworks.io import parse
    from atomworks.io.utils.standard_annotations import S_SEGMIN, S_SEGMAX
    from atomworks.ml.conditions import C_CRD, C_SEQ, C_IDX
    from biotite.structure import BondType
    from rfproteina.conditions import C_IDXGRP, C_HB, C_RAS, HBondEnum, RasaEnum
    from rfproteina.utils.parse_presets import ANNOTATED_CIF

    path = Path(filename).expanduser().resolve(strict=True)
    preset = ANNOTATED_CIF.replace(build_assembly=None, model=1, altloc="first")
    array = parse(str(path), config=preset)["asym_unit"][0]
    annotations = set(array.get_annotation_categories())
    segment = S_SEGMIN.mask(array, default="generate")
    occupancy = array.occupancy if "occupancy" in annotations else np.ones(len(array))
    # Expandable placeholders are instructions, not physical atoms at the origin.
    keep = np.isfinite(array.coord).all(axis=1) & (occupancy > 0) & ~segment
    if not keep.any():
        raise ValueError("No occupied finite physical coordinates remain after excluding placeholders.")
    masks = {
        "fixed": C_CRD.mask(array, default="generate"),
        "sequence_fixed": C_SEQ.mask(array, default="generate"),
        "unindexed": ~C_IDX.annotation(array, default="generate").astype(bool),
    }
    values = {
        "index_group": C_IDXGRP.annotation(array, default="generate"),
        "hbond": C_HB.annotation(array, default="generate"),
        "rasa": C_RAS.annotation(array, default="generate"),
    }
    atoms = []
    for i in np.flatnonzero(keep):
        atoms.append({
            "chain": str(array.chain_id[i]), "resi": str(array.res_id[i]) + str(array.ins_code[i]),
            "resn": str(array.res_name[i]), "name": str(array.atom_name[i]), "element": str(array.element[i]),
            "coord": array.coord[i].tolist(), "hetero": bool(array.hetero[i]), "occupancy": float(occupancy[i]),
            "charge": int(array.charge[i]) if "charge" in annotations else 0,
            **{key: bool(mask[i]) for key, mask in masks.items()},
            "index_group": int(values["index_group"][i]),
            "hbond": HBondEnum(int(values["hbond"][i])).name,
            "rasa": RasaEnum(int(values["rasa"][i])).name,
        })
    physical = array[keep]
    bonds = []
    for i, j, kind in physical.bonds.as_array() if physical.bonds is not None else []:
        bond_name = BondType(int(kind)).name
        order = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3, "QUADRUPLE": 4,
                 "AROMATIC_SINGLE": 1, "AROMATIC_DOUBLE": 2, "AROMATIC_TRIPLE": 3}.get(bond_name, 1)
        bonds.append([int(i), int(j), order])
    segments = []
    seen = set()
    for i in np.flatnonzero(segment):
        key = (str(array.chain_id[i]), int(array.res_id[i]))
        if key not in seen:
            seen.add(key)
            segments.append({"chain": key[0], "resi": key[1],
                             "min": int(S_SEGMIN.annotation(array)[i]), "max": int(S_SEGMAX.annotation(array)[i])})
    return {"source": str(path), "atoms": atoms, "bonds": bonds, "segments": segments,
            "excluded_atoms": int((~keep).sum()), "parse": "ANNOTATED_CIF / asymmetric unit / model 1 / first altloc",
            "condition_fields": sorted(name for name in annotations if name.startswith("condition_")),
            "present": {key: condition.full_name in annotations for key, condition in {
                "fixed": C_CRD, "sequence_fixed": C_SEQ, "unindexed": C_IDX, "hbond": C_HB, "rasa": C_RAS}.items()}}


if __name__ == "__main__":
    if sys.argv[1:] == ["--serve"]:
        # One native import per browsing session, CPU only. Exit after a quiet
        # minute (or when PyMOL closes its pipe), rather than retain idle memory.
        while select.select([sys.stdin], [], [], 60)[0]:
            line = sys.stdin.readline()
            if not line:
                break
            try:
                with redirect_stdout(sys.stderr):
                    result = export_model(json.loads(line)["source"])
                response = {"ok": True, "model": result}
            except Exception:
                response = {"ok": False, "error": traceback.format_exc()}
            print(json.dumps(response, allow_nan=False), flush=True)
        raise SystemExit(0)
    if len(sys.argv) != 2:
        raise SystemExit("Usage: RFD4_PYTHON rfd4_cif_export.py annotated.cif")
    with redirect_stdout(sys.stderr):
        result = export_model(sys.argv[1])
    print(json.dumps(result, allow_nan=False))
