"""Optional structure display bond filters, registered in the PyMOL namespace."""
from pymol import cmd

_AF3_CLASH_FILTER_ENABLED = False   # AF3 clash cleanup OFF by default (slow get_model/bond scan, niche); `clash_bond_filter on` to enable
_AF3_LAST_REMOVED_BONDS = []
_AF3_LIGAND_CENTER_FILTER_ENABLED = False
_AF3_LAST_REMOVED_LIGAND_CENTER_BONDS = []
_AF3_LIGAND_CENTER_ELEMENTS = {'P', 'V'}
_AF3_NON_PTM_CARBON_RESIDUES = {
    'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'TYR'
}
_AF3_AA_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEC', 'PYL'
}
_AF3_BACKBONE_ATOM_NAMES = {
    'N', 'CA', 'C', 'O', 'OXT', 'H', 'H1', 'H2', 'H3', 'HA', 'HA2', 'HA3'
}


def _atom_sel(atom):
    return 'model %s and index %d' % (atom.model, atom.index)


def _is_protein_sidechain_atom(atom):
    if atom.hetatm:
        return False
    if atom.resn not in _AF3_AA_RESIDUES:
        return False
    if atom.name in _AF3_BACKBONE_ATOM_NAMES:
        return False
    return True


def remove_clash_bonds(selection="all", _self=cmd):
    """Remove likely non-physical inter-residue sidechain bonds in proteins.
    Stores removed bonds so they can be restored."""
    global _AF3_LAST_REMOVED_BONDS
    removed = []
    objects = _self.get_object_list(selection)
    for obj in objects:
        model = _self.get_model('(%s) and (%s)' % (obj, selection))
        atoms = model.atom
        for bond in model.bond:
            a1 = atoms[bond.index[0]]
            a2 = atoms[bond.index[1]]
            if (a1.segi, a1.chain, a1.resi_number, a1.resn) == (
                a2.segi, a2.chain, a2.resi_number, a2.resn
            ):
                continue
            same_chain = (a1.segi, a1.chain) == (a2.segi, a2.chain)
            if not same_chain:
                continue
            sidechain_pair = _is_protein_sidechain_atom(a1) and _is_protein_sidechain_atom(a2)
            if not sidechain_pair:
                continue

            hydrophobic_cc = (
                a1.symbol == 'C' and a2.symbol == 'C' and
                a1.resn in _AF3_NON_PTM_CARBON_RESIDUES and
                a2.resn in _AF3_NON_PTM_CARBON_RESIDUES
            )
            c_nos_sidechain = (
                (a1.symbol == 'C' and a2.symbol in ('N', 'O', 'S')) or
                (a2.symbol == 'C' and a1.symbol in ('N', 'O', 'S'))
            )
            if not (hydrophobic_cc or c_nos_sidechain):
                continue
            sel1 = _atom_sel(a1)
            sel2 = _atom_sel(a2)
            _self.unbond(sel1, sel2)
            removed.append((a1.model, a1.index, a2.index, int(getattr(bond, 'order', 1))))
    if removed:
        _AF3_LAST_REMOVED_BONDS = removed
        print('remove_clash_bonds: removed %d suspected inter-residue sidechain bonds.' % len(removed))
    return len(removed)


def remove_ligand_center_pseudobonds(selection="all", center_elements="P,V",
                                     min_coordination=5, _self=cmd):
    # Remove neighbor-neighbor pseudo-bonds around selected ligand center atoms.
    global _AF3_LAST_REMOVED_LIGAND_CENTER_BONDS
    removed = []
    try:
        min_coord = int(min_coordination)
    except:
        min_coord = 5
    center_syms = set()
    for token in str(center_elements).replace(';', ',').split(','):
        sym = token.strip().upper()
        if sym:
            center_syms.add(sym)
    if not center_syms:
        center_syms = set(_AF3_LIGAND_CENTER_ELEMENTS)

    objects = _self.get_object_list(selection)
    for obj in objects:
        model = _self.get_model('(%s) and (%s) and hetatm' % (obj, selection))
        atoms = model.atom
        if not atoms:
            continue

        adjacency = {}
        bond_order = {}
        for bond in model.bond:
            i, j = bond.index
            adjacency.setdefault(i, set()).add(j)
            adjacency.setdefault(j, set()).add(i)
            key = (i, j) if i < j else (j, i)
            bond_order[key] = int(getattr(bond, 'order', 1))

        for center_i, center in enumerate(atoms):
            if not center.hetatm:
                continue
            if center.symbol.upper() not in center_syms:
                continue
            center_res = (center.segi, center.chain, center.resi_number, center.resn)
            neighbors = []
            for nb_i in adjacency.get(center_i, ()):
                nb = atoms[nb_i]
                if not nb.hetatm:
                    continue
                if (nb.segi, nb.chain, nb.resi_number, nb.resn) != center_res:
                    continue
                neighbors.append(nb_i)
            if len(neighbors) < min_coord:
                continue

            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    a_i = neighbors[i]
                    b_i = neighbors[j]
                    if b_i not in adjacency.get(a_i, set()):
                        continue
                    a = atoms[a_i]
                    b = atoms[b_i]
                    sel1 = _atom_sel(a)
                    sel2 = _atom_sel(b)
                    _self.unbond(sel1, sel2)
                    key = (a_i, b_i) if a_i < b_i else (b_i, a_i)
                    removed.append((a.model, a.index, b.index, bond_order.get(key, 1)))
                    adjacency[a_i].discard(b_i)
                    adjacency[b_i].discard(a_i)
    if removed:
        _AF3_LAST_REMOVED_LIGAND_CENTER_BONDS = removed
        print('remove_ligand_center_pseudobonds: removed %d ligand center-neighbor pseudo-bonds.' % len(removed))
    return len(removed)


def restore_clash_bonds(_self=cmd):
    """Restore the last set of bonds removed by remove_clash_bonds."""
    restored = 0
    for model, idx1, idx2, order in _AF3_LAST_REMOVED_BONDS:
        sel1 = 'model %s and index %d' % (model, idx1)
        sel2 = 'model %s and index %d' % (model, idx2)
        if _self.count_atoms(sel1) and _self.count_atoms(sel2):
            _self.bond(sel1, sel2)
            restored += 1
            if order != 1:
                try:
                    _self.set_bond('valence', order, sel1, sel2)
                except:
                    pass
    print('restore_clash_bonds: restored %d bonds.' % restored)
    return restored


def restore_ligand_center_bonds(_self=cmd):
    # Restore bonds removed by remove_ligand_center_pseudobonds.
    restored = 0
    for model, idx1, idx2, order in _AF3_LAST_REMOVED_LIGAND_CENTER_BONDS:
        sel1 = 'model %s and index %d' % (model, idx1)
        sel2 = 'model %s and index %d' % (model, idx2)
        if _self.count_atoms(sel1) and _self.count_atoms(sel2):
            _self.bond(sel1, sel2)
            restored += 1
            if order != 1:
                try:
                    _self.set_bond('valence', order, sel1, sel2)
                except:
                    pass
    print('restore_ligand_center_bonds: restored %d bonds.' % restored)
    return restored


def clash_bond_filter(mode="status", _self=cmd):
    # Enable or disable automatic clash-bond cleanup in color_palette.
    global _AF3_CLASH_FILTER_ENABLED
    token = str(mode).strip().lower()
    if token in ('on', '1', 'true', 'enable', 'enabled'):
        _AF3_CLASH_FILTER_ENABLED = True
    elif token in ('off', '0', 'false', 'disable', 'disabled'):
        _AF3_CLASH_FILTER_ENABLED = False
    elif token == 'toggle':
        _AF3_CLASH_FILTER_ENABLED = not _AF3_CLASH_FILTER_ENABLED
    elif token not in ('status', ''):
        print('clash_bond_filter: use on/off/toggle/status')
    state = 'on' if _AF3_CLASH_FILTER_ENABLED else 'off'
    print('clash_bond_filter: %s' % state)
    return int(_AF3_CLASH_FILTER_ENABLED)


def ligand_center_filter(mode="status", _self=cmd):
    # Enable or disable ligand center-neighbor pseudo-bond cleanup.
    global _AF3_LIGAND_CENTER_FILTER_ENABLED
    token = str(mode).strip().lower()
    if token in ('on', '1', 'true', 'enable', 'enabled'):
        _AF3_LIGAND_CENTER_FILTER_ENABLED = True
    elif token in ('off', '0', 'false', 'disable', 'disabled'):
        _AF3_LIGAND_CENTER_FILTER_ENABLED = False
    elif token == 'toggle':
        _AF3_LIGAND_CENTER_FILTER_ENABLED = not _AF3_LIGAND_CENTER_FILTER_ENABLED
    elif token not in ('status', ''):
        print('ligand_center_filter: use on/off/toggle/status')
    state = 'on' if _AF3_LIGAND_CENTER_FILTER_ENABLED else 'off'
    print('ligand_center_filter: %s' % state)
    return int(_AF3_LIGAND_CENTER_FILTER_ENABLED)


cmd.extend('remove_clash_bonds', remove_clash_bonds)
cmd.extend('restore_clash_bonds', restore_clash_bonds)
cmd.extend('clash_bond_filter', clash_bond_filter)
cmd.extend('remove_ligand_center_pseudobonds', remove_ligand_center_pseudobonds)
cmd.extend('restore_ligand_center_bonds', restore_ligand_center_bonds)
cmd.extend('ligand_center_filter', ligand_center_filter)
