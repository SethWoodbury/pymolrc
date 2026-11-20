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


def check():
    detail = os.environ.get("PYMOL_RFD4") == "1"
    state = cmd._rfd4_viewer_state
    assert state["autoload"] == detail
    style = cmd.keyword["rfd4_style"][0]
    theme = cmd.keyword["theme"][0]
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
        assert theme("next") == "clay" and theme("prev") == "fjord"
        assert cmd.get_view() == camera
        assert cmd.count_atoms("sample_fixed and rep sticks") == 0
        assert state["objects"]["sample"]["style"]["fixed_color"] == "yellow"

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
        assert all(map(math.isfinite, cmd.get_view()))
        assert not any("pymolrc_finite" in n for n in cmd.get_names("all"))
        chem = Path(directory) / "elements.pdb"
        symbols = ("N", "O", "S", "P", "H")
        chem.write_text("".join(
            f"HETATM{i:5} {sym:>4} LIG A   1    {i * 4:8.3f}{0.:8.3f}{0.:8.3f}  1.00 10.00          {sym:>2}\n"
            for i, sym in enumerate(symbols, 1)) + "END\n")
        cmd.load(str(chem), "chem")
        for palette in ("fjord", "clay", "iris", "foliage"):
            theme(palette)
            for sym, stock in zip(symbols, ("nitrogen", "oxygen", "sulfur", "phosphorus", "hydrogen")):
                colors = []
                cmd.iterate("chem and elem " + sym, "out.append(color)", space={"out": colors})
                assert colors == [cmd.get_color_index("elem" + sym if palette == "foliage" else stock)]
    print("STARTUP_MODE_PASS", "detail" if detail else "ordinary", flush=True)


try:
    check()
except Exception:
    traceback.print_exc()
    cmd.quit(1)
else:
    cmd.quit(0)
