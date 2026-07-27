# cootvalent-keys.py  --  bandicoot / Coot key bindings for the covalent workflow
#
# Bandicoot has no coot_python, a stubbed coot_toolbar_button, and invisible-stub
# GTK dialogs -- so a menu/toolbar/input-dialog cannot be added from Python here.
# Key bindings DO work (C-level dispatch, no coot_python needed), so this file
# wires the zero-input covalent actions to Ctrl+<letter>:
#
#   Ctrl+L  declare covalent link (auto-detect the placed warhead + declare)
#   Ctrl+P  propagate covalent ligand to every NCS copy of the Cys
#   Ctrl+U  full: propagate + refine via refmac.sh (-W waters)
# (Ctrl+R is deliberately avoided -- Coot's native rock-view owns it; the other
#  native Ctrl combos are W/F/B/C/M.)
#
# All three are INPUT-FREE: they auto-detect the ligand heavy atom sitting
# 1.2-2.6 A from a CYS SG, so the ligand must already be built/placed near the
# Cys (see "Ligand from SMILES" / ligand_from_smiles() in the console for that).
#
# Ctrl+U needs an MTZ + ligand dictionary. It guesses them from the model's
# directory (<LIGCODE>.cif or LIG.cif; the sole/newest *.mtz). Override once
# from the console if the guess is wrong:
#     cootvalent_set_refine(mtz="/path/data.mtz", lig_dict="/path/LIG.cif")
#
# Load order does not matter: the callbacks resolve declare/propagate/detect at
# key-press time, not at import.  ASCII-only.

import os, glob, subprocess, coot

CCP4_SETUP = "/Applications/ccp4-9/bin/ccp4.setup-sh"
CCP4_PYTHON = "/Applications/ccp4-9/bin/ccp4-python"

# Per-family warhead reduction: transform the UNREACTED SMILES to the covalent
# PRODUCT form (minus the S, which the link adds), so the built ligand has the
# correct reacted geometry (sp2/sp3) and bond order for fitting -- rather than
# the linear alkyne/planar alkene you'd otherwise place and then have to bend.
_REACT_RXN = {
    # butynamide / ynamide: warhead C#C conjugated to amide -> vinyl C=C
    "F2": "[C:1]#[C:2][C:3](=[O:4])[N:5]>>[C:1]=[C:2][C:3](=[O:4])[N:5]",
    # acrylamide: warhead C=C conjugated to amide -> saturated C-C
    "F1": "[C:1]=[C:2][C:3](=[O:4])[N:5]>>[C:1][C:2][C:3](=[O:4])[N:5]",
}

_REACT_WORKER = r'''
import sys, json
from rdkit import Chem
from rdkit.Chem import AllChem
job = json.load(sys.stdin)
m = Chem.MolFromSmiles(job["smiles"])
if m is None:
    print(json.dumps({"error": "unparseable SMILES"})); sys.exit(0)
r = AllChem.ReactionFromSmarts(job["rxn"])
prods = r.RunReactants((m,))
if not prods:
    print(json.dumps({"error": "warhead pattern not matched"})); sys.exit(0)
p = prods[0][0]
try:
    Chem.SanitizeMol(p)
except Exception as e:
    print(json.dumps({"error": "sanitize: %s" % e})); sys.exit(0)
print(json.dumps({"smiles": Chem.MolToSmiles(p)}))
'''


def _reacted_smiles(smiles, family):
    """Return the covalent-product SMILES for a warhead family, or (None, reason).
    Runs RDKit under ccp4-python (Coot's own Python lacks RDKit)."""
    rxn = _REACT_RXN.get(family)
    if rxn is None:
        return None, "no reaction defined for family %s" % family
    import json as _json, tempfile
    src = os.path.join(tempfile.gettempdir(), "cv_react_worker.py")
    open(src, "w").write(_REACT_WORKER)
    cmd = '. "%s" >/dev/null 2>&1; exec "%s" "%s"' % (CCP4_SETUP, CCP4_PYTHON, src)
    try:
        p = subprocess.Popen(["bash", "-c", cmd], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate(_json.dumps({"smiles": smiles, "rxn": rxn}).encode(),
                                 timeout=120)
        res = _json.loads(out.decode("utf-8", "ignore").strip().splitlines()[-1])
    except Exception as e:
        return None, "reaction subprocess failed: %s" % e
    if res.get("error"):
        return None, res["error"]
    return res.get("smiles"), None

DECLARABLE = ("F1", "F2", "F3", "CAA")   # classes the declare side can build
KEY_REFINE_MTZ = None                    # Ctrl+U overrides (None => auto-guess)
KEY_REFINE_DICT = None


def cootvalent_set_refine(mtz=None, lig_dict=None):
    """Set the MTZ / ligand-dict that Ctrl+U uses (persists for the session)."""
    global KEY_REFINE_MTZ, KEY_REFINE_DICT
    if mtz is not None:
        KEY_REFINE_MTZ = mtz
    if lig_dict is not None:
        KEY_REFINE_DICT = lig_dict
    print("[cootvalent] Ctrl+U will use  mtz=%s  dict=%s"
          % (KEY_REFINE_MTZ, KEY_REFINE_DICT))


def _resolve(name):
    """Find a Coot/plugin function by name across every namespace it might live
    in: this module, __main__, builtins, and the already-loaded coot modules
    (coot_utils holds bare helpers like merge_molecules; use sys.modules so we
    never re-import/re-exec them)."""
    import sys
    g = globals()
    if name in g and g[name] is not None:
        return g[name]
    srcs = []
    try:
        import __main__
        srcs.append(vars(__main__))
    except Exception:
        pass
    import builtins
    srcs.append(vars(builtins))
    for modname in ("coot_utils", "coot_load_modules", "coot_gui", "coot"):
        mod = sys.modules.get(modname)
        if mod is not None:
            srcs.append(vars(mod))
    for s in srcs:
        f = s.get(name)
        if f is not None:
            return f
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


def _first_map():
    """Return a valid map molecule number, or -1."""
    n = getattr(coot, "graphics_n_molecules", None)
    ok = getattr(coot, "is_valid_map_molecule", None)
    if n and ok:
        try:
            maps = [i for i in range(n()) if ok(i)]
            if maps:
                return maps[-1]   # most recently made
        except Exception:
            pass
    return -1


def cv_build_at_cys(smiles, chain, resno, name="LIG", family="F2"):
    """Build a ligand from SMILES near a given Cys, as a SEPARATE molecule.

    By default builds the REACTED (covalent-product) form of the warhead for the
    given family (F2 alkyne->vinyl, F1 acrylamide->saturated), so the placed
    ligand already has the correct sp2/sp3 geometry and bond order for fitting
    -- you're not fitting a linear alkyne and hoping refinement bends it. The S
    is added later by the link (which deletes the warhead's placeholder H). Pass
    family=None to build the unreacted SMILES verbatim.

    Deliberately does NOT merge into the protein: a lone ligand molecule is far
    easier to Rotate/Translate + Real Space Refine into the density, and merging
    before fitting piles the ligand centroid onto the SG (unphysical clashes).

    Assumes the protein + map are already loaded. Sets the refinement map,
    centres on <chain>/<resno> SG (a starting point -- the ligand lands there),
    and builds. Then, in the GUI: drag/orient the ligand so its warhead sits
    ~1.8 A from the SG in density. Finally:  cv_merge_ligand()  ->  cv_declare().
    Returns the ligand molecule number.
    """
    # transform the SMILES to the reacted product form for this warhead family
    if family:
        rsmi, err = _reacted_smiles(smiles, family)
        if rsmi is None:
            _status("could not build reacted form (%s); building the SMILES as "
                    "given. Reason: %s" % (family, err))
        else:
            _status("reacted %s warhead: %s" % (family, rsmi))
            smiles = rsmi
    prot = _imol()
    # make sure a refinement map is set (so the build can jiggle-fit)
    sirm = getattr(coot, "set_imol_refinement_map", None)
    irm = getattr(coot, "imol_refinement_map", None)
    if sirm:
        cur = irm() if irm else -1
        if cur is None or cur < 0:
            m = _first_map()
            if m >= 0:
                sirm(m)
    # centre on the target Cys SG (starting location for the build)
    try:
        coot.set_go_to_atom_molecule(prot)
        coot.set_go_to_atom_chain_residue_atom_name(chain, resno, " SG ")
    except Exception as e:
        _status("could not centre on %s/%d SG: %s" % (chain, resno, e))
    build = _resolve("ligand_from_smiles")
    if build is None:
        _status("ligand-from-smiles.py not loaded"); return None
    lig = build(smiles, name)
    if lig is None or lig < 0:
        _status("ligand build failed"); return None
    _status("built %s as molecule %d near %s/%d (NOT merged). Drag/orient it in "
            "the density so the warhead is ~1.8 A from the SG, then "
            "cv_merge_ligand() and cv_declare()." % (name, lig, chain, resno))
    return lig


def cv_clear_ligands(comp="LIG"):
    """Remove ALL copies of a ligand comp id: delete matching residues from
    protein molecules and close any lone ligand-only molecules. Use to reset
    after piled-up/overlapping build attempts, then rebuild cleanly."""
    import tempfile, os as _os
    delres = _resolve("delete_residue")
    closem = _resolve("close_molecule")
    n_del = 0
    n_closed = 0
    for imol in list(_model_molecules()):
        pdb = _os.path.join(tempfile.gettempdir(), "cv_clear_%d.pdb" % imol)
        try:
            coot.write_pdb_file(imol, pdb)
        except Exception:
            continue
        has_protein = False
        lig_res = set()
        for line in open(pdb):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            rn = line[17:20].strip(); ch = line[21]
            try:
                rs = int(line[22:26])
            except ValueError:
                continue
            if rn in ("ALA","GLY","SER","CYS","VAL","LEU","ILE","PRO","THR",
                      "MET","PHE","TYR","TRP","ASP","GLU","ASN","GLN","HIS",
                      "LYS","ARG"):
                has_protein = True
            if rn == comp:
                lig_res.add((ch, rs))
        if not lig_res:
            continue
        if has_protein:
            if delres is None:
                _status("delete_residue unavailable; remove %s from mol %d by hand"
                        % (comp, imol)); continue
            for (ch, rs) in sorted(lig_res):
                try:
                    delres(imol, ch, rs, ""); n_del += 1
                except Exception:
                    pass
        else:
            if closem is None:
                _status("close_molecule unavailable; close mol %d by hand" % imol)
                continue
            try:
                closem(imol); n_closed += 1
            except Exception:
                pass
    _status("cleared %s: deleted %d residue(s) from protein(s), closed %d "
            "ligand-only molecule(s)." % (comp, n_del, n_closed))


def cv_merge_ligand(lig_imol=None, protein_imol=None):
    """Merge the built ligand molecule into the protein molecule so cv_declare()/
    cv_propagate() can see them together. With no args, protein = first coords
    molecule and ligand = the newest other model molecule. Call this if a build
    left the ligand as a separate molecule."""
    prot = protein_imol if protein_imol is not None else _imol()
    if lig_imol is None:
        others = [m for m in _model_molecules() if m != prot]
        if not others:
            _status("no separate ligand molecule to merge (only protein loaded)")
            return
        lig_imol = max(others)
    merge = _resolve("merge_molecules")
    if merge is None:
        _status("merge_molecules not found in any namespace; merge by hand: "
                "Edit > Merge Molecules (ligand mol %d into protein mol %d)."
                % (lig_imol, prot))
        return
    try:
        merge([lig_imol], prot)
        _status("merged ligand mol %d into protein mol %d. Now cv_declare()."
                % (lig_imol, prot))
    except Exception as e:
        _status("merge failed: %s" % e)


def cv_warhead_dist(cutoff=6.0):
    """Report the closest ligand-CARBON <-> CYS-SG distances ACROSS all loaded
    molecules (works whether or not the ligand is merged yet), so you can check
    the fit before merging. The nearest carbon is the warhead pick; a real
    S-Cbond is ~1.8 A. Also notes if the winning pair spans two molecules
    (merge with cv_merge_ligand() before cv_declare())."""
    import tempfile, os as _os, math as _math
    _skip = {"HOH", "WAT", "DOD"}
    sgs, cs = [], []   # pooled across every molecule; each tagged with imol
    for imol in _model_molecules():
        pdb = _os.path.join(tempfile.gettempdir(), "cv_probe_%d.pdb" % imol)
        try:
            coot.write_pdb_file(imol, pdb)
        except Exception:
            continue
        for line in open(pdb):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            an = line[12:16].strip(); rn = line[17:20].strip()
            ch = line[21]; el = (line[76:78].strip() or an[:1])
            try:
                rs = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            if rn == "CYS" and an == "SG":
                sgs.append((imol, ch, rs, xyz))
            elif line.startswith("HETATM") and rn not in _skip and el == "C":
                cs.append((imol, rn, ch, rs, an, xyz))
    if not sgs or not cs:
        _status("need both a Cys SG and a ligand carbon loaded (found %d SG, "
                "%d ligand C)." % (len(sgs), len(cs)))
        return
    pairs = []
    for (limol, rn, lch, lrs, an, lxyz) in cs:
        for (simol, sch, srs, sxyz) in sgs:
            d = _math.sqrt(sum((lxyz[i]-sxyz[i])**2 for i in range(3)))
            if d <= cutoff:
                pairs.append((d, limol, rn, lch, lrs, an, simol, sch, srs))
    pairs.sort()
    if not pairs:
        _status("no ligand carbon within %.1f A of any Cys SG -- move the "
                "ligand closer to the SG, then re-run." % cutoff)
        return
    _status("closest ligand-carbon <-> Cys-SG (nearest = warhead; aim ~1.8 A):")
    for i, (d, limol, rn, lch, lrs, an, simol, sch, srs) in enumerate(pairs[:6]):
        tag = "  <-- WARHEAD" if i == 0 else ""
        span = "" if limol == simol else "  [mols %d/%d]" % (limol, simol)
        print("   %s %s/%d %-4s  <->  CYS %s/%d SG   %.2f A%s%s"
              % (rn, lch, lrs, an, sch, srs, d, tag, span))
    d0, limol0, _, _, _, _, simol0, _, _ = pairs[0]
    if limol0 != simol0:
        _status("warhead pick spans molecules %d (ligand) and %d (protein) -- "
                "run cv_merge_ligand() before cv_declare()." % (limol0, simol0))


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
            _status("Ctrl+U needs an MTZ + ligand dict; couldn't guess both. Set "
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


# ---------------------------------------------------------------------------
# Console entry points (the reliable interface on bandicoot). These are bare
# names in the shared startup namespace, so call them directly in the console.
# ---------------------------------------------------------------------------
def cv_declare():
    """Auto-detect the placed warhead and declare the covalent link."""
    return cv_key_declare()


def cv_propagate():
    """Declare + propagate the covalent ligand to every NCS copy."""
    return cv_key_propagate(False)


def cv_full():
    """Propagate + refine via refmac.sh (-W waters)."""
    return cv_key_propagate(True)


# Keys are DELIBERATELY NOT bound automatically. This bandicoot build has
# compiled-in native accelerators that intercept keys before Python (observed:
# L = go-to-ligand, Ctrl+R = rock, Ctrl+D = DELETE), some destructive, and they
# aren't discoverable from Python -- so auto-binding risks silent shadowing or
# worse. Use the console functions above, or opt in explicitly to a key you
# have confirmed is free with cootvalent_bind_keys(...).
_KNOWN_NATIVE = {
    "l", "L", "Control_l", "Control_r", "Control_d", "Control_w",
    "Control_f", "Control_b", "Control_c", "Control_m", "Control_s",
    "Control_q", "Control_z",
}


def cootvalent_bind_keys(declare=None, propagate=None, full=None):
    """Opt in to key bindings, e.g.
        cootvalent_bind_keys(propagate="Control_p", full="Control_u")
    Only binds the keys you pass. Warns on keys known to be grabbed natively.
    Test a candidate first (worst case a free key just does nothing):
        add_key_binding("t", "Control_o", lambda: print("fired"))
    """
    akb = _resolve("add_key_binding")
    if akb is None:
        print("[cootvalent] add_key_binding unavailable; use the console "
              "functions cv_declare() / cv_propagate() / cv_full().")
        return
    plan = [(declare, "declare", cv_key_declare),
            (propagate, "propagate", lambda: cv_key_propagate(False)),
            (full, "full refine", lambda: cv_key_propagate(True))]
    for key, label, thunk in plan:
        if not key:
            continue
        if key in _KNOWN_NATIVE:
            print("[cootvalent] refusing %r for %s -- it's grabbed by Coot "
                  "natively (may be destructive). Pick another key." % (key, label))
            continue
        try:
            akb("Cootvalent: " + label, key, thunk)
            print("[cootvalent] bound %s -> %s" % (key, label))
        except Exception as e:
            print("[cootvalent] could not bind %r: %s" % (key, e))


print("[cootvalent] loaded. Console (safe): cv_build_at_cys(smiles, chain, "
      "resno) then cv_declare() / cv_propagate() / cv_full().")
print("[cootvalent] keys are OFF by default on this build (native accelerators "
      "can shadow/destroy). Opt in with cootvalent_bind_keys(propagate='Control_o', ...).")
