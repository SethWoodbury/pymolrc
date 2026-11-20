"""Restricted, self-contained annotated-CIF display reader (no scientific imports).

Unknown schemas or chemistry raise UnsupportedCIF: the caller must use the native
reader. This is only a disposable viewer derivative, never an inference parser.
"""
from collections import defaultdict
import gzip
import io
import math
import re
from pathlib import Path


class UnsupportedCIF(ValueError):
    """The native reader is required for this input."""


_AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
_PEPTIDES = {"L-PEPTIDE LINKING", "PEPTIDE LINKING"}
_IONS = set("ZN FE MG MN CU CO NI CA NA K CD".split())
_SCALARS = {
    "condition_sequence_1_residue", "condition_index_1_residue",
    "condition_rasa-binned_1_atom", "condition_hbond_1_atom",
    "condition_pistack_1_atom", "condition_index-group_1_token",
    "condition_CoordinateGroup_1_atom",
} | {"condition_sequence-probability-" + base + "_1_residue"
     for base in _AA | {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}}
_ANNOTATIONS = {"annotation_expsegmin_1_residue", "annotation_expsegmax_1_residue",
                "annotation_src-coordinate-conditioned-region_1_atom",
                "annotation_src-component_1_residue", "annotation_src-index-group_1_atom"}
_CENTROID_FIELDS = ("condition_ChainCentroid_1_chain", "condition_ChainCentroid_1_atom")
_VECTORS = {"condition_coordinate_1_atom", "condition_InternalCoordinate_1_atom", *_CENTROID_FIELDS}


def canonical_condition_fields(names):
    """Report the same names with older and newer native condition registries."""
    return sorted({_CENTROID_FIELDS[0] if name in _CENTROID_FIELDS else name for name in names})


def centroid_targets(chains, coordinates):
    """Collect native per-chain targets before dropping any source atoms.

    All-NaN rows mean absent; zero is a valid target. A target must be finite
    and identical on every atom of its chain, including expandable placeholders.
    Never substitute the measured centroid of retained coordinates.
    """
    if len(chains) != len(coordinates):
        raise ValueError("Invalid chain-centroid annotation length")
    targets = {}
    for chain, row in zip(chains, coordinates):
        chain = str(chain)
        if len(row) != 3:
            raise ValueError("Invalid chain-centroid coordinate shape")
        xyz = tuple(float(value) for value in row)
        if all(math.isnan(value) for value in xyz):
            target = None
        elif all(math.isfinite(value) for value in xyz):
            target = xyz
        else:
            raise ValueError("Partial or infinite chain-centroid coordinates")
        if chain in targets and targets[chain] != target:
            raise ValueError("Conflicting or partially broadcast chain-centroid target for " + chain)
        targets[chain] = target
    return [{"chain": chain, "xyz": list(xyz)} for chain, xyz in targets.items() if xyz is not None]


def export_model(filename):
    """Return the native viewer model for the explicitly supported CIF subset."""
    try:
        return _export_model(filename)
    except UnsupportedCIF:
        raise
    except (ValueError, KeyError, TypeError, IndexError, OverflowError) as error:
        raise UnsupportedCIF(str(error)) from error


# Only rectangular loops of bare values qualify. Any CIF syntax with quoting,
# text, comments, control tokens or wrapped rows stays with the general parser.
_TABLE_TAG = re.compile(r"_[\w-]+\.[\w-]+$")


def _bare_table(lines):
    # Fast substring checks in C. A conservative match simply hands the
    # whole loop to Bio.PDB, including quoted values containing these words.
    raw = "\n" + "".join(lines).lower()
    syntax = ("'", '"', "#", ";", "$", " _", "\t_", "\n_", "\r_",
              "\v", "\f", "loop_", "stop_", "data_", "save_", "global_")
    return not any(token in raw for token in syntax)


def _read_columns(path):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    class UniqueFields(MMCIF2Dict):
        def __setitem__(self, key, value):
            if key in self:
                raise UnsupportedCIF("Duplicate CIF field")
            super().__setitem__(key, value)

        def _tokenize(self, handle):
            blocks = 0
            for token in super()._tokenize(handle):
                blocks += token.lower().startswith("data_")
                if blocks > 1:
                    raise UnsupportedCIF("Exactly one CIF block is required")
                yield token
            if blocks != 1:
                raise UnsupportedCIF("Exactly one CIF block is required")

    small = []
    columns = {}

    def consume(lines):
        if len(lines) > 2 and lines[0].strip() == "loop_" and _TABLE_TAG.fullmatch(lines[1].strip()):
            end = 1
            while end < len(lines) and _TABLE_TAG.fullmatch(lines[end].strip()):
                end += 1
            tags = [line.strip() for line in lines[1:end]]
            body = lines[end:]
            if body and _bare_table(body):
                rows = [row.split() for row in body]
                if all(len(row) == len(tags) for row in rows):
                    columns.update(zip(tags, map(list, zip(*rows))))
                    # Keep one row so duplicate tags and block boundaries are
                    # still validated by the same parser as every other field.
                    small.extend(lines[:end + 1])
                    return
        small.extend(lines)

    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        block = []
        in_text = False
        for line in handle:
            if line.startswith(";"):
                in_text = not in_text
            if not in_text and line.strip() == "#":
                consume(block)
                small.append(line)
                block.clear()
            else:
                block.append(line)
        consume(block)
    data = UniqueFields(io.StringIO("".join(small)))
    data.update(columns)
    return data


def _export_model(filename):
    import numpy as np

    path = Path(filename).expanduser().resolve(strict=True)
    data = _read_columns(path)
    for key in data:
        if key.startswith("_atom_site.condition_") and key.split(".", 1)[1] not in _SCALARS:
            raise UnsupportedCIF("Unknown scalar condition")
        if key.startswith("_condition_") and key.split(".", 1)[0][1:] not in _VECTORS:
            raise UnsupportedCIF("Unknown vector condition")
        if key.startswith("_atom_site.annotation_") and key.split(".", 1)[1] not in _ANNOTATIONS:
            raise UnsupportedCIF("Unknown standard annotation")

    count = len(data["_atom_site.id"])
    # Refuse incomplete loop rows rather than silently truncating with zip().
    for category in ("_atom_site.", "_chem_comp.", "_chem_comp_atom.", "_chem_comp_bond.", "_struct_conn.", "_entity.", "_entity_poly."):
        lengths = {len(value) for key, value in data.items() if key.startswith(category)}
        if len(lengths) > 1:
            raise UnsupportedCIF("Inconsistent category row counts")

    def column(name, default=None):
        key = "_atom_site." + name
        if key not in data:
            if default is None:
                raise UnsupportedCIF("Missing " + key)
            return [default] * count
        value = data[key]
        if len(value) != count:
            raise UnsupportedCIF("Inconsistent atom column length")
        return value

    def boolean(name):
        value = column(name)
        if set(value) - {"True", "False"}:
            raise UnsupportedCIF("Unknown boolean encoding")
        return np.asarray(value) == "True"

    def vector(name):
        prefix = "_" + name + "."
        if set(key for key in data if key.startswith(prefix)) != {prefix + key for key in ("idx0", "val0", "val1", "val2")}:
            raise UnsupportedCIF("Unsupported vector serialization")
        if list(map(int, data[prefix + "idx0"])) != list(range(count)):
            raise UnsupportedCIF("Vector must cover ordered atoms")
        values = np.asarray([data[prefix + "val" + str(i)] for i in range(3)], dtype=float).T
        if values.shape != (count, 3):
            raise UnsupportedCIF("Invalid vector shape")
        if not np.all(np.isfinite(values).all(axis=1) | np.isnan(values).all(axis=1)):
            raise UnsupportedCIF("Partial or infinite coordinate target")
        return values

    # idx0 is positional in the qualified native serialization. Do not infer a
    # remapping from a conflicting optional bookkeeping index: leave such files
    # to the native reader, which remains authoritative for their interpretation.
    if "_atom_site.atom_array_index" in data and list(map(int, column("atom_array_index"))) != list(range(count)):
        raise UnsupportedCIF("Nonidentity atom_array_index requires native parsing")
    if set(column("pdbx_PDB_model_num")) != {"1"} or set(column("label_alt_id")) - {".", "?"}:
        raise UnsupportedCIF("Models or alternate locations require native parsing")
    for field in ("asym_id", "seq_id", "comp_id", "atom_id"):
        if column("label_" + field) != column("auth_" + field):
            raise UnsupportedCIF("Distinct author and label identities")
    if set(column("pdbx_PDB_ins_code")) - {".", "?"}:
        raise UnsupportedCIF("Insertion codes require native parsing")
    if set(column("transformation_id", "1")) != {"1"}:
        raise UnsupportedCIF("Transformed atom instances require native parsing")
    chains, names, residues = column("label_asym_id"), column("label_atom_id"), column("label_comp_id")
    resids = list(map(int, column("label_seq_id")))
    elements = column("type_symbol")
    if {e.upper() for e in elements} - ({"C", "N", "O", "P", "S", "F", "CL", "BR", "I", "X"} | _IONS):
        raise UnsupportedCIF("Hydrogens and unfamiliar elements require native parsing")
    coords = np.asarray([column("Cartn_" + axis) for axis in "xyz"], dtype=np.float32).T
    occupancy = np.asarray(column("occupancy"), dtype=np.float32)
    # Like the native CIF reader, map missing atom-site charges to zero;
    # a reusable component template can have a different protonation state.
    charges = [0 if value in {"?", "."} else int(value) for value in column("pdbx_formal_charge")]
    hetero = column("group_PDB")
    if set(hetero) - {"ATOM", "HETATM"}:
        raise UnsupportedCIF("Unknown atom record type")
    minimum = list(map(int, column("annotation_expsegmin_1_residue")))
    maximum = list(map(int, column("annotation_expsegmax_1_residue")))
    segment = np.asarray(minimum) != -1
    if any((lo == -1 and hi != -1) or lo < -1 or (lo >= 0 and hi < lo) for lo, hi in zip(minimum, maximum)):
        raise UnsupportedCIF("Invalid expandable segment bounds")
    if any((element == "X") != bool(flag) for element, flag in zip(elements, segment)):
        raise UnsupportedCIF("Unknown placeholder encoding")
    fixed = np.isfinite(vector("condition_coordinate_1_atom")).any(axis=1)
    seq = column("condition_sequence_1_residue")
    if any(value not in {"<M>", residue} for value, residue in zip(seq, residues)):
        raise UnsupportedCIF("Nontrivial sequence condition requires native parsing")
    indexed = boolean("condition_index_1_residue")
    index_groups = list(map(int, column("condition_index-group_1_token")))
    coordinate_groups = list(map(int, column("condition_CoordinateGroup_1_atom"))) if "_atom_site.condition_CoordinateGroup_1_atom" in data else np.where(fixed, 0, -1).tolist()
    if any(group < -1 or group > 32767 or (group == 0) != bool(target)
           for group, target in zip(coordinate_groups, fixed)):
        raise UnsupportedCIF("Conflicting coordinate groups")
    hbond, rasa = list(map(int, column("condition_hbond_1_atom"))), list(map(int, column("condition_rasa-binned_1_atom")))
    if set(hbond) - {0, 1, 2, 3} or set(rasa) - {0, 1, 2, 3}:
        raise UnsupportedCIF("Unknown condition enum")
    if any(key.startswith("_condition_InternalCoordinate_1_atom.") for key in data):
        vector("condition_InternalCoordinate_1_atom")
    centroid_fields = [field for field in _CENTROID_FIELDS if any(key.startswith("_" + field + ".") for key in data)]
    if len(centroid_fields) > 1:
        raise UnsupportedCIF("Multiple serializations of chain centroids")
    centroids = centroid_targets(chains, vector(centroid_fields[0])) if centroid_fields else []
    keep = np.isfinite(coords).all(axis=1) & (occupancy > 0) & ~segment
    if not keep.any():
        raise UnsupportedCIF("No physical atoms")
    # Preserve raw identity/order; templates are local, complete, and unambiguous.
    identities = list(zip(chains, resids, residues, names))
    if len(set(identities)) != count:
        raise UnsupportedCIF("Duplicate atom identities")
    component_types = dict(zip(data["_chem_comp.id"], data["_chem_comp.type"]))
    if len(component_types) != len(data["_chem_comp.id"]):
        raise UnsupportedCIF("Duplicate component definitions")
    template_keys = list(zip(data["_chem_comp_atom.comp_id"], data["_chem_comp_atom.atom_id"]))
    template_atoms = dict(zip(template_keys, zip(data["_chem_comp_atom.type_symbol"], map(int, data["_chem_comp_atom.charge"]))))
    if len(template_atoms) != len(template_keys):
        raise UnsupportedCIF("Duplicate local atom definitions")
    entity_types = dict(zip(data["_entity.id"], data["_entity.type"]))
    polymer_types = dict(zip(data["_entity_poly.entity_id"], data["_entity_poly.type"]))
    entities = column("label_entity_id")
    peptide_components = {comp for comp, kind in component_types.items() if kind.upper() in _PEPTIDES}
    for comp in peptide_components - _AA:
        if any(template_atoms.get((comp, name), (None,))[0] != name for name in ("C", "N")):
            raise UnsupportedCIF("Nonstandard peptide attachment atoms require native parsing")
    template_sizes = defaultdict(int)
    for comp, name in template_atoms:
        template_sizes[comp] += 1
    for i in range(count):
        if not segment[i] and (residues[i], names[i]) not in template_atoms:
            raise UnsupportedCIF("Missing local atom template")
        ctype = component_types[residues[i]].upper()
        if not segment[i] and ctype != "NON-POLYMER" and residues[i] not in peptide_components:
            raise UnsupportedCIF("Unsupported polymer chemistry")
        # Atom-site charges are authoritative here (native annotated inputs
        # can intentionally differ from the reusable CCD protonation template).
        if not segment[i] and template_atoms[(residues[i], names[i])][0] != elements[i]:
            raise UnsupportedCIF("Atom element differs from local template")
        if residues[i] in _AA and (ctype not in _PEPTIDES or hetero[i] != "ATOM"):
            raise UnsupportedCIF("Nonstandard protein component/entity semantics")
        if residues[i] in peptide_components and (entity_types[entities[i]] != "polymer"
                or polymer_types.get(entities[i]) != "polypeptide(L)"):
            raise UnsupportedCIF("Unsupported peptide entity")
        if not segment[i] and ctype == "NON-POLYMER" and entity_types[entities[i]] == "polymer":
            raise UnsupportedCIF("Polymer caps require native parsing")
        if elements[i].upper() in _IONS and (ctype != "NON-POLYMER" or template_sizes[residues[i]] != 1):
            raise UnsupportedCIF("Only isolated, locally defined ions are supported")
    templates = defaultdict(list)
    template_edges = set()
    for comp, left, right, order in zip(*(data["_chem_comp_bond." + field] for field in ("comp_id", "atom_id_1", "atom_id_2", "value_order"))):
        if order.upper() not in {"SING", "DOUB", "TRIP", "QUAD"}:
            raise UnsupportedCIF("Unknown bond order")
        edge = (comp, *sorted((left, right)))
        if left == right or edge in template_edges or (comp, left) not in template_atoms or (comp, right) not in template_atoms:
            raise UnsupportedCIF("Ambiguous/incomplete local bond template")
        template_edges.add(edge)
        templates[comp].append((left, right, {"SING": 1, "DOUB": 2, "TRIP": 3, "QUAD": 4}[order.upper()]))
    original_to_kept = {int(i): j for j, i in enumerate(np.flatnonzero(keep))}
    atom_lookup = {key: i for i, key in enumerate(identities)}
    bonds = {}
    valence = defaultdict(int)

    def add_bond(i, j, order):
        if i not in original_to_kept or j not in original_to_kept:
            return
        edge = tuple(sorted((original_to_kept[i], original_to_kept[j])))
        if edge in bonds:
            if bonds[edge] != order:
                raise UnsupportedCIF("Conflicting bonds")
            return
        bonds[edge] = order
        valence[i] += order
        valence[j] += order

    residue_atoms = {}
    residue_order = []
    for i, (chain, resid, resn, name) in enumerate(identities):
        key = chain, resid, resn
        if key not in residue_atoms:
            residue_atoms[key] = {}
            residue_order.append(key)
        elif residue_order[-1] != key:
            raise UnsupportedCIF("Noncontiguous residue")
        residue_atoms[key][name] = i
    for key, atoms in residue_atoms.items():
        if any(keep[i] for i in atoms.values()) and key[2] not in templates:
            first = next(iter(atoms.values()))
            ion_or_water = elements[first].upper() in _IONS or (key[2] == "HOH" and elements[first] == "O")
            if not ion_or_water or template_sizes[key[2]] != 1:
                raise UnsupportedCIF("Missing local bond template")
        for left, right, order in templates[key[2]]:
            if left in atoms and right in atoms:
                add_bond(atoms[left], atoms[right], order)
    # An identity-symmetry disulfide between two unprotonated cysteines
    # has unambiguous connectivity. Delegate other explicit link chemistry.
    for row in range(len(data.get("_struct_conn.id", []))):
        indices = []
        for partner in (1, 2):
            prefix = "_struct_conn.ptnr" + str(partner) + "_label_"
            key = (data[prefix + "asym_id"][row], int(data[prefix + "seq_id"][row]), data[prefix + "comp_id"][row], data[prefix + "atom_id"][row])
            indices.append(atom_lookup[key])
        if all(keep[i] for i in indices):
            if (data["_struct_conn.conn_type_id"][row] != "covale"
                    or data["_struct_conn.pdbx_value_order"][row].upper() != "SING"
                    or any(residues[i] != "CYS" or names[i] != "SG" or elements[i] != "S"
                           or charges[i] != 0 or valence[i] != 1 for i in indices)
                    or indices[0] == indices[1]):
                raise UnsupportedCIF("Explicit link requires native parsing")
            for partner in (1, 2):
                for tag, allowed in (("ptnr%d_symmetry", {"1_555"}),
                                     ("ptnr%d_transformation_id", {"1"}),
                                     ("pdbx_ptnr%d_PDB_ins_code", {".", "?"}),
                                     ("pdbx_ptnr%d_label_alt_id", {".", "?"})):
                    values = data.get("_struct_conn." + tag % partner)
                    if values is not None and values[row] not in allowed:
                        raise UnsupportedCIF("Ambiguous link partner")
            add_bond(*indices, 1)
    for left, right in zip(residue_order, residue_order[1:]):
        if left[0] != right[0] or left[2] not in peptide_components or right[2] not in peptide_components:
            continue
        if right[1] <= left[1]:
            raise UnsupportedCIF("Nonincreasing polymer residue numbering")
        if right[1] == left[1] + 1 and "C" in residue_atoms[left] and "N" in residue_atoms[right]:
            if "OXT" in residue_atoms[left]:
                raise UnsupportedCIF("Polymer leaving groups require native validation")
            carbon, nitrogen = residue_atoms[left]["C"], residue_atoms[right]["N"]
            if all(i in original_to_kept for i in (carbon, nitrogen)):
                if charges[carbon] != 0 or charges[nitrogen] != 0 or valence[carbon] > 3 or valence[nitrogen] > 2:
                    raise UnsupportedCIF("Polymer link requires native valence handling")
            add_bond(carbon, nitrogen, 1)
    atoms = []
    for i in np.flatnonzero(keep):
        atoms.append({"chain": chains[i], "resi": str(resids[i]), "resn": residues[i], "name": names[i],
                      "element": elements[i], "coord": coords[i].tolist(), "hetero": hetero[i] == "HETATM",
                      "occupancy": float(occupancy[i]), "charge": charges[i], "fixed": bool(fixed[i]),
                      "internal": coordinate_groups[i] > 0, "sequence_fixed": seq[i] != "<M>",
                      "unindexed": not bool(indexed[i]), "index_group": index_groups[i],
                      "coordinate_group": coordinate_groups[i], "hbond": ("NONE", "DONOR", "ACCEPTOR", "BOTH")[hbond[i]],
                      "rasa": ("UNSPECIFIED", "BURIED", "INTERMEDIATE", "EXPOSED")[rasa[i]]})
    segments = []
    seen = set()
    for i in np.flatnonzero(segment):
        key = chains[i], resids[i]
        if key not in seen:
            seen.add(key)
            segments.append({"chain": chains[i], "resi": resids[i], "min": minimum[i], "max": maximum[i]})
    fields = {key.split(".", 1)[1] for key in data if key.startswith("_atom_site.condition_")}
    fields |= {key.split(".", 1)[0][1:] for key in data if key.startswith("_condition_")}
    return {"source": str(path), "atoms": atoms, "bonds": [[i, j, order] for (i, j), order in sorted(bonds.items())],
            "segments": segments, "centroids": centroids, "excluded_atoms": int((~keep).sum()),
            "parse": "ANNOTATED_CIF / asymmetric unit / model 1 / first altloc", "condition_fields": canonical_condition_fields(fields),
            "present": {key: value in fields for key, value in {"fixed": "condition_coordinate_1_atom",
                        "sequence_fixed": "condition_sequence_1_residue", "unindexed": "condition_index_1_residue",
                        "hbond": "condition_hbond_1_atom", "rasa": "condition_rasa-binned_1_atom"}.items()}}
