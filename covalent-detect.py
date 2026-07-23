# covalent-detect.py  --  bandicoot / Coot 0.9 extension
#
# AUTO-DETECTOR for PyKeko's Cys-covalent warhead workflow.  Removes the need
# to hand-specify the warhead family + the two link atoms when declaring a
# Cys-covalent link.  Companion to covalent-link.py (declare + refmac side).
#
# Adds  Menu "Cootvalent" -> "Auto-detect + declare covalent link"  :  scan the
# loaded model, find a Cys-SG that bonds a ligand warhead carbon, classify the
# warhead family (F1-F6), and hand the result to declare_covalent_link().
# Declines (info_dialog / returns None) when no F1-F6 warhead is present.
#
# ------------------------------------------------------------------ design ---
# The reactive-group classification needs bond ORDERS, which a PDB does not
# carry (CONECT is connectivity only, and coot.n_atoms() is 0 headless).  The
# authoritative source of bond order is the ligand's monomer DICTIONARY
# (_chem_comp_bond.type), which is exactly what Coot/refmac already use for
# restraints.  So the detector:
#
#   (1) dumps the current model to PDB (real coords + atom names);
#   (2) loads the ligand's monomer dict (CLIBD_MON / Coot lib / caller-supplied)
#       and builds an RDKit mol with correct bond orders in a CCP4 ccp4-python
#       subprocess (Coot's own Python has no RDKit);
#   (3) finds the warhead beta-carbon Cb = the ligand HEAVY atom nearest a CYS
#       SG (1.2-2.6 A, i.e. within covalent-bonding distance);
#   (4) classifies the family from graph topology around Cb PLUS the MEASURED
#       Cb-Ca bond length from the deposited coords (single ~1.5 A -> F1,
#       double ~1.33 A -> F2), and applies veto rules for warheads that are
#       NOT in the F1-F6 taxonomy (e.g. alpha-cyanoacrylamide / cyanoacetamide);
#   (5) returns short-form CIDs //chain/resno/atom for SG and the warhead + the
#       family + a variant tag, or None (declines).
#
# Why not raw registry SMARTS alone: with an S grafted onto the dict graph, the
# bare F1/F2/F3 amide SMARTS collide -- the free-drug 8E8 dict is stored reduced
# (no C=C) while YY3 keeps its acrylamide C=C, so F1 osimertinib would mis-hit
# F2, and cyanoacetamide (4C9) sails through F1.  Measured geometry + a nitrile
# veto resolve both.  The registry (cov-links/index.json) still drives the
# family -> link_id / target-distance mapping downstream in covalent-link.py.
#
# ASCII-only.  2026-07-16.

import os, re, sys, json, math, tempfile, subprocess, coot

CCP4_SETUP = "/Applications/ccp4-9/bin/ccp4.setup-sh"

# Monomer-library roots to search for a ligand dict, in priority order.
_MON_ROOTS = [
    os.environ.get("CLIBD_MON"),
    "/Applications/ccp4-9/lib/data/monomers",
    os.path.join(os.environ.get("COOT_PREFIX", ""), "share/coot/lib/data/monomers")
        if os.environ.get("COOT_PREFIX") else None,
    "/Applications/bandicoot-0.1.4.8/share/coot/lib/data/monomers",
]

# Covalent S-C bonding-distance window (A) used to locate the warhead carbon.
_SC_MIN, _SC_MAX = 1.2, 2.6
# Bond-length thresholds (A) separating single vs double Cb-Ca.
_DOUBLE_MAX = 1.42     # <= this is treated as a double bond (F2 vinyl / ynamide)
_TRIPLE_MAX = 1.27     # (informational; triple never survives in bound state)


# ---------------------------------------------------------------------------
# PDB coordinate reader (ground truth; coot.n_atoms() is 0 headless).
# ---------------------------------------------------------------------------
def _read_pdb_atoms(pdb_path):
    """Return list of dicts: rec, resn, chain, resno, atom, elem, xyz, altloc."""
    out = []
    for line in open(pdb_path):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        try:
            resno = int(line[22:26])
        except ValueError:
            continue
        atom = line[12:16].strip()
        elem = (line[76:78].strip() or atom[0:1]).strip()
        out.append({
            "rec": line[:6].strip(),
            "resn": line[17:20].strip(),
            "chain": line[21],
            "resno": resno,
            "atom": atom,
            "elem": elem,
            "altloc": line[16],
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        })
    return out


# NOTE: prefixed _pkd_ to avoid clobbering covalent-link.py's own _dist (which
# takes 5-tuples) when both extensions autoload into a shared global namespace.
def _pkd_dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def _canon_local(e):
    """Canonical titlecase element symbol (C, O, N, Cl, Br ...)."""
    e = (e or "").strip()
    if len(e) >= 2 and e[0:2].upper() in ("CL", "BR", "SE"):
        return e[0].upper() + e[1].lower()
    return e[0:1].upper()


# ---------------------------------------------------------------------------
# Locate a ligand's monomer dict.
# ---------------------------------------------------------------------------
def _find_ligand_dict(comp, extra=None):
    """Return a path to comp's monomer .cif, or None."""
    if extra and os.path.exists(extra):
        return extra
    lc = comp[0].lower() if comp else ""
    for root in _MON_ROOTS:
        if not root or not os.path.isdir(root):
            continue
        cand = os.path.join(root, lc, comp + ".cif")
        if os.path.exists(cand):
            return cand
        cand2 = os.path.join(root, lc, comp.upper() + ".cif")
        if os.path.exists(cand2):
            return cand2
    return None


# ---------------------------------------------------------------------------
# Ask Coot to write the ligand's own dictionary (covers ligands not in any
# library, e.g. XQQ), best-effort.  Returns a path or None.
# ---------------------------------------------------------------------------
def _coot_write_ligand_dict(comp, workdir):
    """Try Coot's dictionary-export entry points for one comp id."""
    out = os.path.join(workdir, comp + "_coot.cif")
    # cif_file_for_comp_id_py returns the path of the dict Coot has loaded for
    # this comp id (present when a ligand dict was auto-read with the model).
    fn = getattr(coot, "cif_file_for_comp_id_py", None)
    if fn is not None:
        try:
            p = fn(comp)
            if p and os.path.exists(p) and os.path.getsize(p) > 0:
                return p
        except Exception:
            pass
    fn = getattr(coot, "write_dictionary_from_residue", None)
    if fn is not None:
        try:
            fn(comp, out)
            if os.path.exists(out) and os.path.getsize(out) > 0:
                return out
        except Exception:
            pass
    return None


# ---------------------------------------------------------------------------
# The classification worker.  Runs in ccp4-python (RDKit) as a subprocess.
# It is emitted to a temp .py file and reads a JSON job on stdin.
# ---------------------------------------------------------------------------
_WORKER_SRC = r'''
import sys, json, math
from rdkit import Chem
from rdkit.Chem import RWMol

ORDER = {"SINGLE": Chem.BondType.SINGLE, "DOUBLE": Chem.BondType.DOUBLE,
         "TRIPLE": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC,
         "AROMATIC": Chem.BondType.AROMATIC, "DELOC": Chem.BondType.AROMATIC}


def parse_dict(path, comp):
    """Return ([(atom_id, element)], [(a1, a2, order)]) from a monomer dict."""
    atoms = []
    bonds = []
    lines = open(path).read().splitlines()

    def loop_cols(start):
        c = []
        j = start
        while j < len(lines) and lines[j].strip().startswith("_"):
            c.append(lines[j].strip())
            j += 1
        return c, j

    for i, l in enumerate(lines):
        s = l.strip()
        if s.startswith("_chem_comp_atom.comp_id"):
            cols, j = loop_cols(i)
            idx = {c.split(".", 1)[1]: k for k, c in enumerate(cols)}
            k = j
            while k < len(lines):
                r = lines[k].split()
                st = lines[k].strip()
                if not r or st.startswith(("_", "loop_", "data_", "#")):
                    break
                if r[0] != comp:
                    break
                atoms.append((r[idx["atom_id"]], r[idx["type_symbol"]]))
                k += 1
        if s.startswith("_chem_comp_bond.comp_id"):
            cols, j = loop_cols(i)
            idx = {c.split(".", 1)[1]: k for k, c in enumerate(cols)}
            k = j
            while k < len(lines):
                r = lines[k].split()
                st = lines[k].strip()
                if not r or st.startswith(("_", "loop_", "data_", "#")):
                    break
                if r[0] != comp:
                    break
                bonds.append((r[idx["atom_id_1"]], r[idx["atom_id_2"]],
                              r[idx["type"]].upper()))
                k += 1
    return atoms, bonds


# Covalent radii (A) for geometry-based connectivity perception (fallback when
# no monomer dict is available, e.g. a ligand absent from every library).
_COV = {"C": 0.77, "N": 0.70, "O": 0.66, "S": 1.04, "CL": 0.99, "BR": 1.14,
        "F": 0.57, "P": 1.07, "H": 0.31, "I": 1.33, "SE": 1.20}


def _dist3(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def _order_from_len(e1, e2, d):
    """Crude bond order from length for the C/N/O bonds classification needs."""
    key = frozenset((e1, e2))
    if key == frozenset(("C", "C")):
        return "TRIPLE" if d < 1.27 else ("DOUBLE" if d < 1.42 else "SINGLE")
    if key == frozenset(("C", "O")):
        return "DOUBLE" if d < 1.28 else "SINGLE"
    if key == frozenset(("C", "N")):
        return "TRIPLE" if d < 1.21 else ("DOUBLE" if d < 1.35 else "SINGLE")
    return "SINGLE"


def perceive_from_geometry(geom_atoms):
    """geom_atoms: [[atom_id, element, [x,y,z]], ...] (heavy atoms).
    Returns ([(atom_id, element)], [(a1, a2, order)]) by distance perception."""
    atoms = [(a[0], a[1].upper()) for a in geom_atoms]
    bonds = []
    n = len(geom_atoms)
    for i in range(n):
        for j in range(i + 1, n):
            e1 = geom_atoms[i][1].upper()
            e2 = geom_atoms[j][1].upper()
            d = _dist3(geom_atoms[i][2], geom_atoms[j][2])
            cut = _COV.get(e1, 0.77) + _COV.get(e2, 0.77) + 0.45
            if 0.4 < d < cut:
                bonds.append((geom_atoms[i][0], geom_atoms[j][0],
                              _order_from_len(e1, e2, d)))
    return atoms, bonds


def build_mol(atoms, bonds):
    """RWMol with heavy atoms only; return (mol, atomid->idx, idx->atomid)."""
    m = RWMol()
    amap = {}
    for aid, el in atoms:
        if el == "H" or el == "D":
            continue
        sym = el[0:1].upper() + el[1:].lower() if len(el) > 1 else el.upper()
        amap[aid] = m.AddAtom(Chem.Atom(sym))
    for a1, a2, o in bonds:
        if a1 in amap and a2 in amap:
            m.AddBond(amap[a1], amap[a2], ORDER.get(o, Chem.BondType.SINGLE))
    mol = m.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(mol,
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except Exception:
            pass
    inv = {v: k for k, v in amap.items()}
    return mol, amap, inv


def neigh(bonds, a):
    out = []
    for a1, a2, o in bonds:
        if a1 == a:
            out.append((a2, o))
        elif a2 == a:
            out.append((a1, o))
    return out


def _canon_elem(e):
    """Canonical titlecase element symbol (C, O, N, Cl, Br ...)."""
    e = (e or "").strip()
    if len(e) >= 2 and e[0:2].upper() in ("CL", "BR", "SE"):
        return e[0].upper() + e[1].lower()
    return e[0:1].upper()


def elem_of(atoms, aid):
    for a, e in atoms:
        if a == aid:
            return _canon_elem(e)
    return "?"


def is_carbonyl_c(atoms, bonds, c):
    """True if carbon c has a =O (double-bond oxygen)."""
    if elem_of(atoms, c) != "C":
        return False
    for nb, o in neigh(bonds, c):
        if elem_of(atoms, nb) == "O" and o == "DOUBLE":
            return True
    return False


def amide_c(atoms, bonds, c):
    """True if c is a carbonyl carbon that also bears an N (i.e. an amide)."""
    if not is_carbonyl_c(atoms, bonds, c):
        return False
    for nb, o in neigh(bonds, c):
        if elem_of(atoms, nb) == "N":
            return True
    return False


def has_nitrile(atoms, bonds, c):
    """True if carbon c bears a nitrile substituent (c-C#N)."""
    for nb, o in neigh(bonds, c):
        if elem_of(atoms, nb) == "C":
            for nb2, o2 in neigh(bonds, nb):
                if nb2 != c and elem_of(atoms, nb2) == "N" and o2 == "TRIPLE":
                    return True
    return False


def classify(job):
    """job: dict(dict_path|geom_atoms, comp, warhead_atom, cb_ca_dist, ...)
    Returns dict(family, variant, ca, reason) or dict(family=None, reason=...)."""
    source = "dict"
    if job.get("dict_path"):
        atoms, bonds = parse_dict(job["dict_path"], job["comp"])
        if not atoms or not bonds:
            # dict was unreadable/empty; fall back to geometry if provided.
            if job.get("geom_atoms"):
                atoms, bonds = perceive_from_geometry(job["geom_atoms"])
                source = "geometry"
            else:
                return {"family": None, "reason": "empty/unreadable ligand dict"}
    elif job.get("geom_atoms"):
        atoms, bonds = perceive_from_geometry(job["geom_atoms"])
        source = "geometry"
    else:
        return {"family": None, "reason": "no dict and no geometry supplied"}
    if not atoms or not bonds:
        return {"family": None, "reason": "no atoms/bonds after %s parse" % source}
    cb = job["warhead_atom"]
    if cb not in [a for a, e in atoms]:
        return {"family": None,
                "reason": "warhead atom %s not in dict" % cb}
    if elem_of(atoms, cb) != "C":
        return {"family": None, "reason": "warhead atom is not carbon"}

    heavy_nb = [(n, o) for n, o in neigh(bonds, cb)
                if elem_of(atoms, n) != "H"]
    carbon_nb = [(n, o) for n, o in heavy_nb if elem_of(atoms, n) == "C"]
    oxy_nb = [(n, o) for n, o in heavy_nb if elem_of(atoms, n) == "O"]

    # ---- F3 (chloroacetamide SN2): Cb bonded DIRECTLY to an amide carbonyl ----
    # (one carbon between S and C=O).  Pre-form still has the Cl on Cb.
    direct_amide = [n for n, o in carbon_nb if amide_c(atoms, bonds, n)]
    has_cl = any(elem_of(atoms, n) == "Cl" for n, o in heavy_nb)
    if direct_amide and not has_nitrile(atoms, bonds, cb):
        # make sure this isn't actually the F1 case where Cb->Ca->C=O; here Cb
        # itself neighbours the carbonyl, so it is F3 topology.
        return {"family": "F3", "variant": "post" if not has_cl else "chloride",
                "ca": None, "reason": "Cb directly bonded to amide carbonyl"}

    # ---- F4 (epoxide) / F6 (reversible carbonyl) via OH/O on Cb ----
    # F6 post: Cb carries an -OH (hemithioketal) and NO amide reachable.
    # F4 post: Cb-Ca-OH (hydroxyl on the ADJACENT carbon).
    cb_oh = any(elem_of(atoms, n) == "O" and o == "SINGLE" for n, o in oxy_nb)

    # ---- Identify Ca = the carbon neighbour of Cb that leads to an amide ----
    ca = None
    ca_is_amide_bearing = False
    for n, o in carbon_nb:
        for n2, o2 in neigh(bonds, n):
            if n2 != cb and amide_c(atoms, bonds, n2):
                ca = n
                ca_is_amide_bearing = True
                break
        if ca:
            break

    # ---- VETO: alpha-cyanoacrylamide / cyanoacetamide (edge case, NOT F1-F6) --
    # The 4C9 warhead: Cb-Ca(C#N)-C(=O)N.  Ca (or Cb) carries a nitrile.
    if ca and (has_nitrile(atoms, bonds, ca) or has_nitrile(atoms, bonds, cb)):
        return {"family": None,
                "reason": "alpha-cyano(acrylamide) warhead -- not in F1-F6"}

    # ---- F4 epoxide: Ca carries an -OH, no amide involved ----
    if ca is None:
        # look for a hydroxyl-bearing adjacent carbon (epoxide-opened form)
        for n, o in carbon_nb:
            for n2, o2 in neigh(bonds, n):
                if elem_of(atoms, n2) == "O" and o2 == "SINGLE":
                    return {"family": "F4", "variant": "post", "ca": n,
                            "reason": "Cb-Ca-OH beta-hydroxy thioether"}
        # F6 reversible carbonyl: Cb has an OH (hemithioketal) or =O and no N
        if cb_oh:
            return {"family": "F6", "variant": "post", "ca": None,
                    "reason": "Cb hemithioketal (-OH), reversible carbonyl"}
        # ketone-specific: Cb is a carbonyl carbon flanked by two carbons and
        # NO nitrogen (a true ketone, not an amide) -> F6 pre.
        if is_carbonyl_c(atoms, bonds, cb) and len(carbon_nb) >= 2 \
                and not any(elem_of(atoms, n) == "N" for n, o in heavy_nb):
            return {"family": "F6", "variant": "carbonyl", "ca": None,
                    "reason": "Cb ketone carbonyl (no N) -- reversible"}
        return {"family": None,
                "reason": "no amide/warhead reachable from Cb"}

    # ---- F5 maleimide: Cb and Ca in a 5-ring with TWO carbonyls sharing an N --
    # detect ring closure: Cb has a second carbonyl-C neighbour whose amide-N is
    # shared with the Ca-side carbonyl.
    ring_carbonyls = [n for n, o in carbon_nb if is_carbonyl_c(atoms, bonds, n)]
    if ca_is_amide_bearing and ring_carbonyls:
        return {"family": "F5", "variant": "post", "ca": ca,
                "reason": "maleimide 3-thiosuccinimide ring"}

    # ---- F1 vs F2: measured Cb-Ca bond length decides ----
    d = job.get("cb_ca_dist")
    if d is not None and d <= job.get("double_max", 1.42):
        return {"family": "F2", "variant": "post", "ca": ca,
                "reason": "vinyl thioether (Cb=Ca %.3f A, sp2)" % d}
    return {"family": "F1", "variant": "post", "ca": ca,
            "reason": "saturated beta-thioether (Cb-Ca %.3f A, sp3)"
                      % (d if d is not None else -1)}


def main():
    job = json.load(sys.stdin)
    try:
        res = classify(job)
    except Exception as e:
        res = {"family": None, "reason": "worker error: %s" % e}
    if "source" not in res:
        res["source"] = "dict" if job.get("dict_path") else "geometry"
    json.dump(res, sys.stdout)


main()
'''


def _run_worker(job):
    """Run the RDKit classification worker under ccp4-python; return its dict."""
    wd = job.get("_workdir") or tempfile.mkdtemp(prefix="pk_detect_")
    wpath = os.path.join(wd, "detect_worker.py")
    with open(wpath, "w") as fh:
        fh.write(_WORKER_SRC)
    cmd = ('. "%s" >/dev/null 2>&1; exec ccp4-python "%s"'
           % (CCP4_SETUP, wpath))
    p = subprocess.Popen(["bash", "-c", cmd], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate(json.dumps(job).encode("utf-8"), timeout=120)
    txt = out.decode("utf-8", "ignore").strip()
    if not txt:
        return {"family": None,
                "reason": "worker produced no output: %s"
                          % err.decode("utf-8", "ignore")[-400:]}
    try:
        return json.loads(txt)
    except Exception:
        return {"family": None, "reason": "worker output not JSON: %s" % txt[-400:]}


# ---------------------------------------------------------------------------
# Family display names (kept in sync with the registry).
# ---------------------------------------------------------------------------
_FAMILY_NAMES = {
    "F1": "acrylamide -> saturated beta-thioether",
    "F2": "alpha,beta-ynamide -> vinyl thioether",
    "F3": "alpha-chloro/cyano-acetamide SN2 thioether",
    "F4": "epoxide -> beta-hydroxy thioether",
    "F5": "maleimide -> 3-thiosuccinimide",
    "F6": "reversible ketone/aldehyde hemithioketal",
}


def _short_cid(chain, resno, atom):
    return "//%s/%d/%s" % (chain, resno, atom)


# ---------------------------------------------------------------------------
# MAIN CALLABLE.
# ---------------------------------------------------------------------------
def detect_covalent_link(imol_protein, ligand_comp_or_cid, lig_dict=None,
                         workdir=None, verbose=True):
    """Auto-detect a Cys-covalent warhead link on molecule imol_protein.

    ligand_comp_or_cid : the target ligand.  Accepts a 3-letter comp id
        ("8E8"), a residue CID ("//A/701" or "//A/701(8E8)"), or None/"" to
        scan every non-standard ligand in the model.
    lig_dict           : optional path to the ligand's monomer dict (used when
                         the ligand isn't in any library, e.g. XQQ from acedrg).

    Returns dict(sg_cid, warhead_cid, family, variant, ...) on a confident
    F1-F6 detection, or None (declines).
    """
    wd = workdir or tempfile.mkdtemp(prefix="pk_detect_")
    src_pdb = os.path.join(wd, "model.pdb")
    try:
        coot.write_pdb_file(imol_protein, src_pdb)
    except Exception as e:
        if verbose:
            print("[cootvalent] detect: write_pdb_file failed:", e)
        return None
    atoms = _read_pdb_atoms(src_pdb)
    if not atoms:
        if verbose:
            print("[cootvalent] detect: model has no atoms")
        return None

    # Candidate ligand residues (rec=HETATM, not water/ion) matching the request.
    want_comp, want_chain, want_resno = _parse_ligand_request(ligand_comp_or_cid)
    ligres = _ligand_residues(atoms, want_comp, want_chain, want_resno)
    if not ligres:
        if verbose:
            print("[cootvalent] detect: no candidate ligand residue for '%s'"
                  % ligand_comp_or_cid)
        return None

    # CYS SG atoms in the model.
    sgs = [a for a in atoms if a["resn"] == "CYS" and a["atom"] == "SG"]
    if not sgs:
        if verbose:
            print("[cootvalent] detect: no CYS SG atoms in model")
        return None

    for (comp, chain, resno) in ligres:
        res_atoms = [a for a in atoms if a["resn"] == comp
                     and a["chain"] == chain and a["resno"] == resno
                     and a["altloc"] in (" ", "A")]
        heavy = [a for a in res_atoms if a["elem"] not in ("H", "D")
                 and not a["atom"].startswith("H")]

        # find the warhead atom = ligand heavy atom nearest a CYS SG in-window.
        best = None
        for a in heavy:
            for s in sgs:
                d = _pkd_dist(a["xyz"], s["xyz"])
                if _SC_MIN < d < _SC_MAX:
                    if best is None or d < best[0]:
                        best = (d, a, s)
        if best is None:
            continue
        sc_dist, wh_atom, sg_atom = best

        # get the ligand dict (bond orders).  Preferred source of truth; when
        # absent (ligand in no library, e.g. XQQ), the worker falls back to
        # geometry perception from the deposited coordinates we pass below.
        dpath = _find_ligand_dict(comp, extra=lig_dict)
        if dpath is None:
            dpath = _coot_write_ligand_dict(comp, wd)
        if dpath is None and verbose:
            print("[cootvalent] detect: no monomer dict for %s -- using geometry "
                  "perception fallback" % comp)

        # measure Cb-Ca candidate distances (Ca = each heavy carbon neighbour
        # of the warhead atom within ~1.7 A).  Pass the SHORTEST C-C as cb_ca.
        cb_ca = None
        for a in heavy:
            if a is wh_atom:
                continue
            if _canon_local(a["elem"]) != "C":
                continue
            d = _pkd_dist(a["xyz"], wh_atom["xyz"])
            if 1.1 < d < 1.75:
                if cb_ca is None or d < cb_ca:
                    cb_ca = d

        geom_atoms = [[a["atom"], a["elem"], list(a["xyz"])] for a in heavy]
        job = {
            "dict_path": dpath, "comp": comp, "warhead_atom": wh_atom["atom"],
            "cb_ca_dist": cb_ca, "double_max": _DOUBLE_MAX, "_workdir": wd,
            "geom_atoms": geom_atoms,
        }
        res = _run_worker(job)
        fam = res.get("family")
        if fam is None:
            if verbose:
                print("[cootvalent] detect: %s/%d(%s) atom %s -> DECLINE (%s)"
                      % (chain, resno, comp, wh_atom["atom"],
                         res.get("reason")))
            continue

        result = {
            "sg_cid": _short_cid(sg_atom["chain"], sg_atom["resno"], "SG"),
            "warhead_cid": _short_cid(chain, resno, wh_atom["atom"]),
            "family": fam,
            "variant": res.get("variant"),
            "family_name": _FAMILY_NAMES.get(fam, ""),
            "ligand_comp": comp,
            "cys_resno": sg_atom["resno"], "cys_chain": sg_atom["chain"],
            "warhead_atom": wh_atom["atom"],
            "sg_warhead_dist": round(sc_dist, 3),
            "cb_ca_dist": round(cb_ca, 3) if cb_ca else None,
            "reason": res.get("reason"),
            "graph_source": res.get("source"),
            "dict_path": dpath,
        }
        if verbose:
            print("[cootvalent] detect: %s  SG(%s/%d) <-> %s(%s/%d/%s)  family=%s (%s)"
                  % (comp, sg_atom["chain"], sg_atom["resno"], comp, chain,
                     resno, wh_atom["atom"], fam, res.get("reason")))
        return result

    if verbose:
        print("[cootvalent] detect: no F1-F6 warhead detected on any candidate ligand")
    return None


def _parse_ligand_request(req):
    """Return (comp_or_None, chain_or_None, resno_or_None)."""
    if not req:
        return None, None, None
    req = req.strip()
    # residue CID //A/701 or //A/701(8E8)
    m = re.match(r"^/*/?([A-Za-z0-9])/(-?\d+)(?:\(([^)]+)\))?$", req)
    if m:
        return (m.group(3) or None), m.group(1), int(m.group(2))
    # bare comp id
    if re.match(r"^[A-Za-z0-9]{1,3}$", req):
        return req.upper(), None, None
    return None, None, None


# Residue names never treated as covalent-ligand candidates.
_SKIP_RESN = set("""HOH WAT NA CL K MG CA ZN MN FE CU CO NI CD BR IOD
    SO4 PO4 GOL EDO PEG DMS ACT NAG BMA MAN FUC GAL BGC EPE MES TRS
    CIT FMT ACY DOD OXY""".split())
# standard amino acids + nucleotides
_SKIP_RESN |= set("""ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET
    PHE PRO SER THR TRP TYR VAL MSE SEC PYL A C G U T DA DC DG DT DU""".split())


def _ligand_residues(atoms, want_comp, want_chain, want_resno):
    """Distinct (comp, chain, resno) ligand residues, filtered by request."""
    seen = []
    order = []
    for a in atoms:
        if a["rec"] != "HETATM":
            continue
        if a["resn"] in _SKIP_RESN:
            continue
        key = (a["resn"], a["chain"], a["resno"])
        if key in seen:
            continue
        if want_comp and a["resn"] != want_comp.upper():
            continue
        if want_chain and a["chain"] != want_chain:
            continue
        if want_resno is not None and a["resno"] != want_resno:
            continue
        seen.append(key)
        order.append(key)
    return order


# ---------------------------------------------------------------------------
# GUI wiring + bridge to declare_covalent_link (from covalent-link.py).
# ---------------------------------------------------------------------------
def auto_detect_and_declare(imol=None, ligand=None, mtz=None, lig_dict=None,
                            do_declare=True):
    """Detect then hand off to declare_covalent_link.  Returns the detect dict
    (augmented with declare-result keys) or None."""
    if imol is None:
        try:
            imol = coot.first_coords_imol()
        except Exception:
            imol = 0
    det = detect_covalent_link(imol, ligand, lig_dict=lig_dict)
    if det is None:
        return None
    if do_declare:
        declare = globals().get("declare_covalent_link")
        if declare is None:
            print("[cootvalent] detect: declare_covalent_link not loaded "
                  "(load covalent-link.py first) -- returning detection only")
        else:
            try:
                dres = declare(imol, det["sg_cid"], det["warhead_cid"],
                               family=det["family"], mtz=mtz, lig_dict=lig_dict,
                               do_refine=bool(mtz))
                det["declare_result"] = dres
            except Exception as e:
                print("[cootvalent] detect: declare_covalent_link failed:", e)
                det["declare_error"] = str(e)
    return det


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


def _pkd_install_menu():
    mbm = _gui_fn("coot_menubar_menu")
    asm = _gui_fn("add_simple_coot_menu_menuitem")
    gse = _gui_fn("generic_single_entry")
    if not (mbm and asm and gse):
        print("[cootvalent] detect: GUI menu API not found (no-graphics?); "
              "detect_covalent_link() still callable directly.")
        return

    def _go(ligand_text):
        try:
            imol = coot.first_coords_imol()
        except Exception:
            imol = 0
        lig = (ligand_text or "").strip() or None
        det = detect_covalent_link(imol, lig)
        if det is None:
            coot.info_dialog("No covalent warhead detected.\n\n"
                             "No F1-F6 Cys-warhead was found near a CYS SG on "
                             "the selected ligand (or Cootvalent declined an "
                             "unsupported chemistry). See the scripting console "
                             "for the reason.")
            return
        declare = globals().get("declare_covalent_link")
        msg = ("Detected covalent link:\n\n"
               "  family:  %s  (%s)\n"
               "  Cys SG:  %s\n"
               "  warhead: %s\n"
               "  S-C dist: %.2f A\n\n"
               % (det["family"], det["family_name"], det["sg_cid"],
                  det["warhead_cid"], det["sg_warhead_dist"]))
        if declare is None:
            coot.info_dialog(msg + "(covalent-link.py not loaded, so the link "
                             "was NOT declared. Load it to enable declare.)")
            return
        try:
            imol = coot.first_coords_imol()
        except Exception:
            imol = 0
        try:
            declare(imol, det["sg_cid"], det["warhead_cid"],
                    family=det["family"], mtz=None, do_refine=True)
            coot.info_dialog(msg + "Link declared. See the scripting console "
                             "for the augmented PDB + link CIF paths and the "
                             "refmac command.")
        except Exception as e:
            coot.info_dialog(msg + "Detect OK, but declare failed:\n\n%s" % e)

    def _activate(*args):
        gse(
            "Ligand comp id or CID (blank = scan all ligands), "
            "e.g. 8E8 or //A/701",
            "", "Auto-detect + declare", _go)

    try:
        menu = mbm("Cootvalent")
        if menu is None:
            try:
                import coot_python  # noqa: F401
                tail = "RESTART bcoot (menu installs at startup)."
            except Exception:
                tail = ("this Coot build has no coot_python; use the console: "
                        "auto_detect_and_declare(...).")
            print("[cootvalent] detect: menu not installed -- " + tail)
            return
        asm(menu, "Auto-detect + declare covalent link", _activate)
        print("[cootvalent] OK: 'Auto-detect + declare covalent link' added to "
              "the Cootvalent menu.")
    except Exception:
        import traceback
        print("[cootvalent] detect MENU INSTALL FAILED:\n" + traceback.format_exc())


_pkd_install_menu()
