# cootvalent-keys.py  --  bandicoot / Coot key bindings for the covalent workflow
#
# Bandicoot has no coot_python, a stubbed coot_toolbar_button, and invisible-stub
# GTK dialogs -- so a menu/toolbar/input-dialog cannot be added from Python here.
# Key bindings DO work (C-level dispatch, no coot_python needed), so this file
# wires the zero-input covalent actions to Ctrl+<letter>:
#
#   Ctrl+L  declare covalent link (auto-detect the placed warhead + declare)
#   Ctrl+P  propagate covalent ligand to every NCS copy of the Cys
#   Ctrl+R  full: propagate + refine via refmac.sh (-W waters)
#
# All three are INPUT-FREE: they auto-detect the ligand heavy atom sitting
# 1.2-2.6 A from a CYS SG, so the ligand must already be built/placed near the
# Cys (see "Ligand from SMILES" / ligand_from_smiles() in the console for that).
#
# Ctrl+R needs an MTZ + ligand dictionary. It guesses them from the model's
# directory (<LIGCODE>.cif or LIG.cif; the sole/newest *.mtz). Override once
# from the console if the guess is wrong:
#     cootvalent_set_refine(mtz="/path/data.mtz", lig_dict="/path/LIG.cif")
#
# Load order does not matter: the callbacks resolve declare/propagate/detect at
# key-press time, not at import.  ASCII-only.

import os, glob, coot

DECLARABLE = ("F1", "F2", "F3", "CAA")   # classes the declare side can build
KEY_REFINE_MTZ = None                    # Ctrl+R overrides (None => auto-guess)
KEY_REFINE_DICT = None


def cootvalent_set_refine(mtz=None, lig_dict=None):
    """Set the MTZ / ligand-dict that Ctrl+R uses (persists for the session)."""
    global KEY_REFINE_MTZ, KEY_REFINE_DICT
    if mtz is not None:
        KEY_REFINE_MTZ = mtz
    if lig_dict is not None:
        KEY_REFINE_DICT = lig_dict
    print("[cootvalent] Ctrl+R will use  mtz=%s  dict=%s"
          % (KEY_REFINE_MTZ, KEY_REFINE_DICT))


def _resolve(name):
    """Find a name in the shared startup namespace (globals/__main__/builtins)."""
    g = globals()
    if name in g and g[name] is not None:
        return g[name]
    try:
        import __main__
        if hasattr(__main__, name):
            return getattr(__main__, name)
    except Exception:
        pass
    import builtins
    if hasattr(builtins, name):
        return getattr(builtins, name)
    return None


def _imol():
    try:
        return coot.first_coords_imol()
    except Exception:
        return 0


def _status(msg):
    print("[cootvalent] " + msg)
    try:
        coot.add_status_bar_text("[cootvalent] " + msg)
    except Exception:
        pass


def _norm_family(f):
    return "CAA" if f == "F3" else f


def _model_molecules():
    """All valid coordinate molecules, lowest first (falls back to first_coords)."""
    n = getattr(coot, "graphics_n_molecules", None)
    ok = getattr(coot, "is_valid_model_molecule", None)
    if n is None or ok is None:
        return [_imol()]
    try:
        return [i for i in range(n()) if ok(i)]
    except Exception:
        return [_imol()]


def _detect():
    """Scan every loaded model for a covalent ligand; return the detection dict
    (with '_imol' set) for the first model that has one, else None."""
    fn = _resolve("detect_covalent_link")
    if fn is None:
        _status("covalent-detect.py not loaded"); return None
    cands = _model_molecules()
    # try the active molecule first, then the rest
    act = _imol()
    if act in cands:
        cands = [act] + [c for c in cands if c != act]
    for imol in cands:
        try:
            det = fn(imol, None)
        except Exception:
            det = None
        if det and det.get("family"):
            det["_imol"] = imol
            return det
    _status("no covalent warhead found near a CYS SG in any loaded model "
            "(build/place the ligand at the Cys first)")
    return None


def _declarable(det):
    fam = det["family"]
    if fam not in DECLARABLE:
        _status("detected %s (%s) -- not declarable yet; only F1/F2/CAA are built"
                % (fam, det.get("family_name", "")))
        return None
    return _norm_family(fam)


def cv_key_declare():
    det = _detect()
    if det is None:
        return
    imol = det["_imol"]
    fam = _declarable(det)
    if fam is None:
        return
    declare = _resolve("declare_covalent_link")
    if declare is None:
        _status("covalent-link.py not loaded"); return
    try:
        declare(imol, det["sg_cid"], det["warhead_cid"], family=fam, do_refine=False)
        _status("declared %s: %s <-> %s" % (fam, det["sg_cid"], det["warhead_cid"]))
    except Exception as e:
        _status("declare failed: %s" % e)


def _model_dir(imol):
    nm = ""
    try:
        nm = coot.molecule_name(imol)
    except Exception:
        pass
    d = os.path.dirname(nm) if nm else ""
    return d if (d and os.path.isdir(d)) else os.getcwd()


def _guess_refine_inputs(imol, det):
    d = _model_dir(imol)
    dic = KEY_REFINE_DICT
    if not dic:
        for cand in (det.get("ligand_comp", "LIG") + ".cif", "LIG.cif"):
            p = os.path.join(d, cand)
            if os.path.exists(p):
                dic = p; break
    mtz = KEY_REFINE_MTZ
    if not mtz:
        mtzs = sorted(glob.glob(os.path.join(d, "*.mtz")), key=os.path.getmtime)
        if mtzs:
            mtz = mtzs[-1]   # newest
            if len(mtzs) > 1:
                _status("multiple MTZs in %s; using newest (%s). Override with "
                        "cootvalent_set_refine(mtz=...)." % (d, os.path.basename(mtz)))
    return mtz, dic


def cv_key_propagate(do_refine=False):
    det = _detect()
    if det is None:
        return
    imol = det["_imol"]
    fam = _declarable(det)
    if fam is None:
        return
    prop = _resolve("propagate_covalent_ncs")
    if prop is None:
        _status("covalent-link.py not loaded"); return
    kwargs = dict(family=fam, do_refine=False)
    if do_refine:
        mtz, dic = _guess_refine_inputs(imol, det)
        if not mtz or not dic:
            _status("Ctrl+R needs an MTZ + ligand dict; couldn't guess both. Set "
                    "them: cootvalent_set_refine(mtz='...', lig_dict='...'), retry.")
            return
        kwargs.update(do_refine=True, use_wrapper=True, add_waters=True,
                      mtz=mtz, lig_dict=dic)
        _status("full refine using mtz=%s dict=%s"
                % (os.path.basename(mtz), os.path.basename(dic)))
    try:
        r = prop(imol, det["warhead_cid"], det["cys_resno"], **kwargs)
        _status("propagated to %d covalent copy(ies)%s%s"
                % (r.get("n_copies", 0),
                   " + refined" if do_refine else "",
                   ("  -> " + r["refined_pdb"]) if r.get("refined_pdb") else ""))
    except Exception as e:
        _status("propagate failed: %s" % e)


def _install_keys():
    akb = _resolve("add_key_binding")
    if akb is None:
        print("[cootvalent] add_key_binding unavailable (no-graphics?); "
              "keys not bound. Console: cv_key_declare()/cv_key_propagate().")
        return
    try:
        akb("Cootvalent: declare covalent link (auto)", "Control_l", cv_key_declare)
        akb("Cootvalent: propagate covalent ligand to NCS copies", "Control_p",
            lambda: cv_key_propagate(False))
        akb("Cootvalent: full (propagate + refine + waters)", "Control_r",
            lambda: cv_key_propagate(True))
        print("[cootvalent] keys bound over the graphics window: "
              "Ctrl+L declare, Ctrl+P propagate, Ctrl+R full refine.")
    except Exception:
        import traceback
        print("[cootvalent] key install failed:\n" + traceback.format_exc())


_install_keys()
