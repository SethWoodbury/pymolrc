"""Run with both pymol and pymolrfd4: -cq -r tests/check_startup_modes.py."""
import math
import os
from pathlib import Path
import sys
import tempfile
import traceback

from pymol import cmd

tests_dir = Path(__script__).resolve().parent
sys.path.insert(0, str(tests_dir))
from test_rfd4_cif_fast import fixture, write_dict


def check_palette(directory):
    atoms = []

    def atom(record, name, resn, chain, resi, insertion, element):
        i = len(atoms) + 1
        atoms.append(f"{record:6}{i:5} {name:>4} {resn:3} {chain}{resi:4}{insertion}   "
                     f"{float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00 10.00          {element:>2}\n")

    for chain, residues in (("A", [(2, " "), (9, " "), (9, "A")]), ("B", [(7, " ")])):
        for resi, insertion in residues:
            for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
                atom("ATOM", name, "ALA", chain, resi, insertion, element)
    for resn, number in (("GOL", 11), ("ACT", 12)):
        atom("HETATM", "C1", resn, "L", number, " ", "C")
        atom("HETATM", "O1", resn, "L", number, " ", "O")
    atom("HETATM", "ZN1", "GOL", "L", 11, " ", "Zn")
    atom("HETATM", "CD", "LIG", "L", 13, " ", "C")
    atom("HETATM", "YB", "YB", "L", 14, " ", "Yb")
    source = Path(directory) / "palette.pdb"
    source.write_text("".join(atoms) + "END\n")
    cmd.delete("all")
    cmd.do("theme mineral, gradient=off")
    cmd.load(str(source), "pair")

    def colors(sel):
        values = []
        cmd.iterate(sel, "out.append(color)", space={"out": values})
        return set(values)

    ca = "pair and name CA"
    before = colors("pair")
    camera = cmd.get_view()
    cmd.do("theme mineral, gradient=on")
    assert cmd.get_view() == camera and colors("pair") == before
    stops = {}
    cmd.iterate(ca, "out[chain, resi] = s.cartoon_color", space={"out": stops})
    assert stops["A", "2"] == cmd.get_color_index("0x786F62")
    assert stops["A", "9A"] == cmd.get_color_index("0xA09886")
    assert stops["A", "9"] == stops["B", "7"] == cmd.get_color_index("0x8C8474")
    cmd.do("theme mineral, gradient=off")
    cleared = []
    cmd.iterate(ca, "out.append(s.cartoon_color)", space={"out": cleared})
    assert all(c in (None, -1) for c in cleared)

    # Model a modified residue with a HETATM branch sharing its ATOM residue ID.
    mixed = cmd.get_model("pair and ((chain A and resi 2) or resn ACT)")
    for record in mixed.atom:
        record.chain, record.resi, record.resn = "M", "1", "ALA"
    cmd.load_model(mixed, "modified", zoom=0)
    assert cmd.count_atoms("modified and polymer.protein") and cmd.count_atoms("modified and organic and name C1")
    cmd.do("theme mineral")
    cmd.do("ligand_colors by_type")
    assert "ALA" not in cmd._pymolrc_ligand_types
    assert colors("modified and name C1") == {cmd.get_color_index("mineral_taupe")}
    for resn, color in (("ACT", "mineral_lagoon"), ("GOL", "fern_orchid")):
        assert colors("pair and resn %s and elem C" % resn) == {cmd.get_color_index(color)}
    cmd.load(str(source), "copy")
    assert colors("copy and resn GOL and elem C") == colors("pair and resn GOL and elem C")
    assert colors("resn GOL and elem O") == {cmd.get_color_index("oxygen")}
    cmd.show_as("cartoon", "pair")
    assert cmd.count_atoms("pair and metals and rep spheres") == 2
    assert cmd.count_atoms("pair and name CD and rep spheres") == 0
    cmd.hide("spheres", "pair and metals")
    cmd.do("theme fern")
    assert cmd.count_atoms("pair and metals and rep spheres") == 0
    cmd.do("show_metals pair")
    assert cmd.count_atoms("pair and metals and rep spheres") == 2
    assert colors("pair and elem Zn") == {cmd.get_color_index("zinc")}
    cmd.do("ligand_colors uniform")
    assert colors("pair and organic and elem C") == {cmd.get_color_index("fern_orchid")}


def check():
    detail = os.environ.get("PYMOL_RFD4") == "1"
    state = cmd._rfd4_viewer_state
    assert state["autoload"] == detail
    style = cmd.keyword["rfd4_style"][0]
    theme = cmd.keyword["theme"][0]
    assert theme("list") == "mineral"
    for setting, expected in {"specular": 0, "ray_shadow": 0, "ray_trace_mode": 1,
                              "power": 0.2, "ambient": 0.4, "direct": 0.45,
                              "reflect": 1, "reflect_power": 0, "orthoscopic": 1,
                              "cartoon_fancy_helices": 1, "cartoon_fancy_sheets": 1}.items():
        assert math.isclose(cmd.get_setting_float(setting), expected, abs_tol=1e-6), setting
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "a structure.cif"
        data = fixture()
        write_dict(source, data)
        detector = style.__globals__["_RFD4_DETECT"]
        original_detector = detector.is_rfd4_cif
        if not detail:
            def unexpected_read(path):
                raise AssertionError("ordinary load inspected CIF annotations")
            detector.is_rfd4_cif = unexpected_read
        try:
            cmd.load(str(source), "sample")
        finally:
            detector.is_rfd4_cif = original_detector
        assert ("sample" in state["objects"]) == detail
        if not detail:
            assert state["cache_stats"]["last_route"] == "none"
            cmd.delete("sample")
            cmd.keyword["rfd4_load"][0](str(source), "sample")

        assert cmd.count_atoms("sample and rep cartoon") > 0
        assert cmd.count_atoms("sample_fixed and rep spheres") == 1
        assert "sample_annotations" not in cmd.get_names("objects", enabled_only=1)
        assert cmd.count_atoms("sample_internal and rep spheres") == 1
        for mask, color in (("fixed", "rfd4_gold"), ("internal", "rfd4_lilac")):
            values = []
            cmd.iterate("sample_" + mask, "out.append((s.sphere_color, s.sphere_scale))", space={"out": values})
            assert values[0][0] == cmd.get_color_index(color) and abs(values[0][1] - 0.16) < 1e-5
        style("sample", fixed_scale=0.12, fixed_color="yellow")
        cmd.hide("sticks", "sample_fixed")
        camera = cmd.get_view()
        assert theme("next") == "fern" and theme("prev") == "mineral"
        assert cmd.get_view() == camera
        assert cmd.count_atoms("sample_fixed and rep sticks") == 0
        assert state["objects"]["sample"]["style"]["fixed_color"] == "yellow"
        semantic = []
        cmd.iterate("sample and name CA", "out.append(s.cartoon_color)", space={"out": semantic})
        theme("fern", gradient="on")
        after = []
        cmd.iterate("sample and name CA", "out.append(s.cartoon_color)", space={"out": after})
        assert semantic == after and cmd.count_atoms("sample_fixed and rep sticks") == 0
        style("sample", "elements")
        assert state["objects"]["sample"]["cartoon_gradient"]
        theme("mineral", gradient="off")
        assert not state["objects"]["sample"].get("cartoon_gradient")

        # User-added group members survive deleting an annotated molecule.
        cmd.pseudoatom("mine")
        cmd.group("sample_annotations", "mine", action="add")
        cmd.delete("sample")
        assert "mine" in cmd.get_names("objects")
        assert not any(n.startswith("sample_") for n in cmd.get_names("all"))
        assert cmd.keyword["rfd4_load"][0](str(source), "sample") == "sample"

        cmd.do("@" + str(tests_dir.parent / ".pymolrc"))
        cmd.delete("all")
        cmd.keyword["rfd4_autoload"][0]("off")
        for axis in "xyz":
            data["_atom_site.Cartn_" + axis][0] = "nan"
        write_dict(source, data)
        cmd.load(str(source), "ordinary")
        assert cmd.count_atoms("ordinary") == 8 and "ordinary" not in state["objects"]
        backbone = []
        cmd.iterate("ordinary and polymer.protein and name CA", "out.append(color)", space={"out": backbone})
        assert backbone and set(backbone) == {cmd.get_color_index("mineral_stone")}
        view = cmd.get_view()
        assert all(map(math.isfinite, view))
        assert all(math.isclose(a, b, abs_tol=1e-4) for a, b in zip(view[12:15], (4, 0, 0)))
        assert not any("pymolrc_finite" in n for n in cmd.get_names("all"))
        chem = Path(directory) / "elements.pdb"
        symbols = ("N", "O", "S", "P", "H")
        chem.write_text("".join(
            f"HETATM{i:5} {sym:>4} LIG A   1    {i * 4:8.3f}{0.:8.3f}{0.:8.3f}  1.00 10.00          {sym:>2}\n"
            for i, sym in enumerate(symbols, 1)) + "END\n")
        cmd.load(str(chem), "chem")
        for palette in ("mineral", "fern", "classic", "fjord", "clay", "iris", "foliage"):
            theme(palette)
            for sym, stock in zip(symbols, ("nitrogen", "oxygen", "sulfur", "phosphorus", "hydrogen")):
                colors = []
                cmd.iterate("chem and elem " + sym, "out.append(color)", space={"out": colors})
                assert colors == [cmd.get_color_index("elem" + sym if palette == "foliage" else stock)]
        check_palette(directory)
    print("STARTUP_MODE_PASS", "detail" if detail else "ordinary", flush=True)


try:
    check()
except Exception:
    traceback.print_exc()
    cmd.quit(1)
else:
    cmd.quit(0)
