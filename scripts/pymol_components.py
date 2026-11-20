"""Color organic ligand carbons while recognizing split ATOM/HETATM residues."""


def paint_component_carbons(scope, primary, cmd, protein_color=None, assignments=None, accents=()):
    """Optional assignments map component names to stable positions in a palette."""
    carbons = "(%s) and organic and elem C" % scope
    candidates = {}
    cmd.iterate(carbons, "atoms[model, index] = (model, segi, chain, resi, resn)",
                space={"atoms": candidates})
    if not candidates:
        return

    # PyMOL's byres can separate the ATOM and HETATM portions of the same
    # modified amino acid. Match the full chemical residue identity instead.
    residues = set(candidates.values())
    names = {r[-1] for r in residues}
    polymer = "(byobject (%s)) and polymer" % carbons
    if all(n.isascii() and n.isalnum() for n in names):
        polymer += " and resn " + "+".join(sorted(names))
    in_polymer = set()
    cmd.iterate(polymer, "seen.add((model, segi, chain, resi, resn))", space={"seen": in_polymer})
    in_polymer.intersection_update(residues)
    if not in_polymer and assignments is None:
        cmd.color(primary, carbons)
        return

    ligand_names = sorted({r[-1] for r in residues - in_polymer})
    by_name = {}
    for name in ligand_names:
        color = primary
        if assignments is not None:
            slot = assignments.setdefault(name, len(assignments))
            color = accents[slot % len(accents)]
        by_name[name] = cmd.get_color_index(color)
    protein = cmd.get_color_index(protein_color) if protein_color is not None else None
    colors = {}
    for atom, residue in candidates.items():
        if residue not in in_polymer:
            colors[atom] = by_name[residue[-1]]
        elif protein is not None:
            colors[atom] = protein
    if colors:
        cmd.alter(carbons, "color = colors.get((model, index), color)", space={"colors": colors})
        cmd.recolor(carbons)
