# ligand-from-smiles.py  --  bandicoot / Coot 0.9 extension
#
# Adds a "Cootvalent" menu with TWO ligand-placement workflows, both native in Coot:
#
#   1. "Ligand from SMILES..."           (place-at-pointer + jiggle-fit)
#      Type a SMILES, acedrg builds a restraints dictionary + 3D coords, Coot
#      loads it, places it at the view centre, and jiggle-fits it into the
#      active map.  Best when you already know WHERE the ligand goes and just
#      want it snapped into local density.
#
#   2. "Find in density (blind search)..."   (native blind ligand search)
#      Type a SMILES (+ a cluster-sigma level), acedrg builds the dictionary,
#      then Coot's native execute_ligand_search searches the WHOLE Fo-Fc
#      difference map, clusters the positive density, and drops a copy of the
#      ligand into every candidate cluster.  This is the STRONGER "where does
#      my ligand go?" path -- it finds the site for you.  ALL clusters are
#      returned / centred-on so you can pick the right one (the top-ranked
#      cluster is not always the correct pose).
#
# INSTALL:  copy this file into ~/.coot-preferences/   (Coot autoloads *.py there;
#           create the folder if it doesn't exist).  Restart bcoot.
# OR run it once from Coot's scripting window (Calculate -> Scripting... -> Python).
#
# USAGE (place-at-pointer):  have an active map (Open MTZ...), then
#           Menu Cootvalent -> Ligand from SMILES...
#           Entry accepts  "SMILES"  or  "SMILES 3LETTERNAME"  (e.g. "CCO ETH").
#           Or call directly:  ligand_from_smiles("CCO", "ETH")
#
# USAGE (blind search):  have a protein model + an MTZ open (so a Fo-Fc
#           difference map exists), then Menu Cootvalent -> Find in density...
#           Entry accepts  "SMILES"  or  "SMILES 3LETTERNAME"; second box is the
#           cluster-sigma level (default 1.0; lower it, e.g. 0.6, to catch
#           weaker / partial-occupancy density).
#           Or call directly:  ligand_blind_find("CCO", "ETH", imol_prot)
#
# Requires CCP4 (acedrg) installed.  2026-07-16.

import os, subprocess, tempfile, coot

CCP4_SETUP = "/Applications/ccp4-9/bin/ccp4.setup-sh"   # edit if CCP4 lives elsewhere


def _acedrg_from_smiles(smiles, name, workdir):
    """Run acedrg on a SMILES string -> (cif_path or None, error_or_None)."""
    out = os.path.join(workdir, "lig")
    # PYTHONNOUSERSITE=1 is belt-and-suspenders: it stops a stray python user-site
    # from shadowing CCP4's bundled numpy/gemmi and breaking servalcat/acedrg.
    cmd = ('export PYTHONNOUSERSITE=1; . "%s" >/dev/null 2>&1; '
           'acedrg -i "%s" -r %s -o "%s" > "%s/acedrg.log" 2>&1'
           % (CCP4_SETUP, smiles, name, out, workdir))
    rc = subprocess.call(["bash", "-c", cmd])
    cif = out + ".cif"
    if rc != 0 or not os.path.exists(cif):
        tail = ""
        try:
            tail = open(os.path.join(workdir, "acedrg.log")).read()[-800:]
        except Exception:
            pass
        return None, "acedrg failed (rc=%d).\n\n%s" % (rc, tail)
    return cif, None


def _build_monomer_from_smiles(smiles, name):
    """SMILES -> acedrg dictionary -> read_cif_dictionary -> get_monomer.
    Returns (imol_lig or -1, error_or_None).  Shared by both workflows."""
    workdir = tempfile.mkdtemp(prefix="pk_lig_")
    coot.add_status_bar_text("Running acedrg on %s ..." % smiles)
    cif, err = _acedrg_from_smiles(smiles, name, workdir)
    if cif is None:
        return -1, err
    coot.read_cif_dictionary(cif)
    imol = coot.get_monomer(name)     # builds ideal 3D coords from the new dict
    if imol < 0:
        return -1, ("acedrg made a dictionary but Coot could not build "
                    "monomer '%s'." % name)
    return imol, None


def ligand_from_smiles(smiles, name="LIG", imol_map=None, do_jiggle=True):
    """SMILES -> acedrg dictionary -> load -> place at view centre -> jiggle-fit
    into the active map.  Returns the new ligand molecule number, or -1.
    Safe to call from the scripting window or the socket/MCP bridge."""
    name = (name or "LIG").strip().upper()[:3] or "LIG"
    smiles = (smiles or "").strip()
    if not smiles:
        coot.info_dialog("No SMILES given.")
        return -1

    imol, err = _build_monomer_from_smiles(smiles, name)
    if imol < 0:
        coot.info_dialog("Ligand from SMILES failed:\n\n" + err)
        return -1

    # place it where you are looking. Don't use coot_utils.move_molecule_to_
    # screen_centre -- in bandicoot it routes through the undefined C name
    # is_protein_chain_p and raises. Do it with C-level calls directly.
    try:
        rcp = getattr(coot, "rotation_centre_position", None)
        mcf = getattr(coot, "molecule_centre", None) or _gui_fn("molecule_centre")
        tmb = getattr(coot, "translate_molecule_by", None) or _gui_fn("translate_molecule_by")
        if rcp and mcf and tmb:
            mc = mcf(imol)
            tmb(imol, rcp(0) - mc[0], rcp(1) - mc[1], rcp(2) - mc[2])
        else:
            print("[cootvalent] move-to-centre skipped: coord API unavailable "
                  "(ligand loaded at its built coordinates -- move it by hand).")
    except Exception as e:
        print("[cootvalent] move-to-centre skipped:", e)

    # fit into density
    if do_jiggle:
        m = imol_map if imol_map is not None else coot.imol_refinement_map()
        if m is not None and m >= 0 and coot.is_valid_map_molecule(m):
            coot.set_imol_refinement_map(m)
            coot.fit_molecule_to_map_by_random_jiggle(imol, 50, 1.0)
            coot.add_status_bar_text("Placed + jiggle-fit ligand '%s' (mol %d)"
                                     % (name, imol))
        else:
            coot.add_status_bar_text("Ligand '%s' placed (mol %d) - no active "
                                     "map to fit into" % (name, imol))

    coot.set_go_to_atom_molecule(imol)
    return imol


# ---------------- blind search: "where does my ligand go?" ----------------

def _find_difference_map():
    """Return the imol of the first valid Fo-Fc difference map, or -1.
    map_is_difference_map(i)==1 flags an Fo-Fc map."""
    try:
        n = coot.graphics_n_molecules()
    except Exception:
        n = 0
    for i in range(n):
        try:
            if coot.is_valid_map_molecule(i) and coot.map_is_difference_map(i) == 1:
                return i
        except Exception:
            pass
    return -1


def ligand_blind_find(smiles, name="LIG", imol_prot=None,
                      imol_diff_map=None, sigma=1.0):
    """Native BLIND ligand search: build the ligand from SMILES (acedrg), then
    search the WHOLE Fo-Fc difference map, cluster the positive density, and
    drop a copy of the ligand into every candidate cluster.

    This is the stronger 'where does my ligand go?' path -- it finds the site.

    Args:
      smiles        : SMILES string (optionally trailing 3-letter name is
                      handled by the GUI, not here).
      name          : 3-letter monomer code for the built ligand.
      imol_prot     : the protein model to search against.  If None, the first
                      valid model molecule is used.
      imol_diff_map : the Fo-Fc map to search.  If None, auto-detected via
                      map_is_difference_map().
      sigma         : cluster-sigma level.  REQUIRED tuning knob -- default 1.0
                      matches Coot's default, but weaker / partial-occupancy
                      density often needs a LOWER value (e.g. 0.6).

    Returns a list of dicts, one per returned cluster, ORDERED as Coot ranks
    them (rank 0 = highest scoring).  NOTE: the top-ranked cluster is NOT always
    the correct pose -- return every cluster and let the user pick.  Each dict:
        {"imol": <placed molecule number>, "rank": <0-based rank>}
    Returns [] on failure (and pops an info dialog in interactive mode).
    """
    name = (name or "LIG").strip().upper()[:3] or "LIG"
    smiles = (smiles or "").strip()
    if not smiles:
        coot.info_dialog("No SMILES given.")
        return []

    # resolve the protein model
    if imol_prot is None or imol_prot < 0:
        imol_prot = -1
        try:
            n = coot.graphics_n_molecules()
        except Exception:
            n = 0
        for i in range(n):
            try:
                if coot.is_valid_model_molecule(i):
                    imol_prot = i
                    break
            except Exception:
                pass
    if imol_prot is None or imol_prot < 0 or not coot.is_valid_model_molecule(imol_prot):
        coot.info_dialog("Blind search needs a protein model open.")
        return []

    # resolve the Fo-Fc difference map
    if imol_diff_map is None or imol_diff_map < 0:
        imol_diff_map = _find_difference_map()
    if imol_diff_map < 0 or not coot.is_valid_map_molecule(imol_diff_map):
        coot.info_dialog("Blind search needs an Fo-Fc difference map.\n"
                         "Open an MTZ (auto-read makes 2Fo-Fc + Fo-Fc) first.")
        return []

    # build the ligand from SMILES
    imol_lig, err = _build_monomer_from_smiles(smiles, name)
    if imol_lig < 0:
        coot.info_dialog("Blind search: ligand build failed:\n\n" + err)
        return []

    # wire up + run the native blind search
    try:
        sigma = float(sigma)
    except Exception:
        sigma = 1.0

    coot.set_ligand_search_protein_molecule(imol_prot)
    coot.set_ligand_search_map_molecule(imol_diff_map)
    coot.add_ligand_search_ligand_molecule(imol_lig)
    coot.set_find_ligand_do_real_space_refinement(1)
    coot.set_ligand_cluster_sigma_level(sigma)
    coot.add_status_bar_text("Blind ligand search (sigma=%.2f) ..." % sigma)

    solutions = coot.execute_ligand_search_py()

    # execute_ligand_search_py returns a list of placed molecule numbers
    # (one per accepted cluster), highest-scoring first.
    result = []
    if solutions:
        for rank, entry in enumerate(solutions):
            # each entry is typically the placed imol (an int); be defensive.
            imol = entry
            if isinstance(entry, (list, tuple)) and len(entry) > 0:
                imol = entry[0]
            try:
                imol = int(imol)
            except Exception:
                continue
            result.append({"imol": imol, "rank": rank})

    if not result:
        coot.add_status_bar_text("Blind search found no clusters at sigma=%.2f "
                                 "- try a lower sigma (e.g. 0.6)." % sigma)
        coot.info_dialog("Blind search found no clusters at sigma=%.2f.\n"
                         "Try a lower cluster-sigma (e.g. 0.6) to catch weaker "
                         "or partial-occupancy density." % sigma)
        return []

    coot.add_status_bar_text("Blind search: %d cluster(s) at sigma=%.2f "
                             "(rank 0 is not always correct - inspect all)."
                             % (len(result), sigma))
    # centre on the top-ranked solution so the user lands on something.
    try:
        coot.set_go_to_atom_molecule(result[0]["imol"])
    except Exception:
        pass
    return result


# ---------------- GUI wiring: "Cootvalent" menu ----------------
def _gui_fn(name):
    """Resolve a Coot GUI helper by name without importing coot_gui (which
    re-executes and raises in bandicoot -- see covalent-link.py for details)."""
    import sys
    srcs = [globals()]
    try:
        import __main__
        srcs.append(vars(__main__))
    except Exception:
        pass
    import builtins
    srcs.append(vars(builtins))
    cg = sys.modules.get("coot_gui")
    if cg is not None:
        srcs.append(vars(cg))
    for s in srcs:
        f = s.get(name)
        if f is not None:
            return f
    return None


def _install_menu():
    mbm = _gui_fn("coot_menubar_menu")
    asm = _gui_fn("add_simple_coot_menu_menuitem")
    gse = _gui_fn("generic_single_entry")
    gde = _gui_fn("generic_double_entry")
    if not (mbm and asm and gse and gde):
        print("[cootvalent] GUI menu API not found (no-graphics?); "
              "ligand_from_smiles()/ligand_blind_find() still callable directly.")
        return

    # ---- mode 1: place-at-pointer + jiggle-fit ----
    def _go_place(text):
        parts = (text or "").split()
        if not parts:
            return
        smi = parts[0]
        nm = parts[1] if len(parts) > 1 else "LIG"
        ligand_from_smiles(smi, nm)

    def _activate_place(*args):
        gse(
            "SMILES  (optionally:  SMILES 3-letter-name)",
            "CCO", "Build + Fit", _go_place)

    # ---- mode 2: blind search into the difference map ----
    def _go_blind(smiles_text, sigma_text):
        parts = (smiles_text or "").split()
        if not parts:
            return
        smi = parts[0]
        nm = parts[1] if len(parts) > 1 else "LIG"
        try:
            sigma = float((sigma_text or "1.0").strip())
        except Exception:
            sigma = 1.0
        sols = ligand_blind_find(smi, nm, None, None, sigma)
        # summarise the clusters for the user (all of them; pick, don't auto-accept).
        if sols:
            lines = ["Blind search placed %d cluster(s) at sigma=%.2f."
                     % (len(sols), sigma),
                     "The top rank is NOT always correct -- inspect each:", ""]
            for s in sols:
                lines.append("  rank %d  ->  molecule %d" % (s["rank"], s["imol"]))
            coot.info_dialog("\n".join(lines))

    def _activate_blind(*args):
        # generic_double_entry(label1, label2, default1, default2,
        #                      check_button_label, check_fn, go_label, go_fn)
        # go_fn is called with (entry1_text, entry2_text).
        gde(
            "SMILES  (optionally:  SMILES 3-letter-name)",
            "Cluster sigma level (default 1.0; lower for weak density)",
            "CCO", "1.0",
            False, False,          # no check button
            "Find in density", _go_blind)

    try:
        menu = mbm("Cootvalent")
        if menu is None:
            try:
                import coot_python  # noqa: F401
                tail = "RESTART bcoot (menu installs at startup)."
            except Exception:
                tail = ("this Coot build has no coot_python; use the console: "
                        "ligand_from_smiles(smiles) / ligand_blind_find(smiles).")
            print("[cootvalent] smiles: menu not installed -- " + tail)
            return
        asm(menu, "Ligand from SMILES...", _activate_place)
        asm(menu, "Find in density (blind search)...", _activate_blind)
        print("[cootvalent] OK: 'Ligand from SMILES...' + 'Find in density "
              "(blind search)...' added to the Cootvalent menu.")
    except Exception:
        import traceback
        print("[cootvalent] smiles MENU INSTALL FAILED:\n" + traceback.format_exc())


_install_menu()
