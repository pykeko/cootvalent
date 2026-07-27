# covalent-link.py  --  bandicoot / Coot 0.9 extension
#
# Adds  Menu "Cootvalent" -> "Declare covalent link..."  :  declare a Cys-warhead
# covalent bond the PyKeko way, native in classic Coot.
#
# It mirrors PyKeko's crown-jewel covalent-ligand mechanism, ported off the
# Moorhen/JS mmCIF-surgery path onto Coot + external refmac5:
#
#   (a) writes an AUGMENTED MODEL  -- the current model saved to PDB, with a
#       LINK record (+ LINKR link-id suffix after col 80) for the Cys-SG <-> warhead
#       bond, plus a mmCIF twin carrying a _struct_conn row.  This is the load-
#       bearing artifact refmac matches on.
#   (b) writes a refmac-ready LINK-CIF restraint dictionary for the chosen
#       warhead family (data_link_list + data_mod_list catalog blocks +
#       link_id columns in every loop -- refmac silently drops unlabeled loops).
#   (c) provides a REFMAC5 refine action: spawns external refmac5 with the
#       link CIF (+ the ligand's own monomer dict) as LIBIN, the augmented
#       PDB as XYZIN, an MTZ as HKLIN, and loads the refined model back.
#
# Families supported in this prototype (the covalent-ligand-plan taxonomy):
#   F1  CYS-ACR   acrylamide -> saturated sp3 beta-thioether   (S-Cb 1.81 A)
#   F2  CYS-YNA   ynamide/butynamide -> sp2 vinyl thioether    (S-Cb 1.78 A)
#   CAA CYS-CAA   alpha-cyano/chloro-acetamide SN2 thioether   (S-Cb 1.81 A)
#
# INSTALL:  copy this file into ~/.coot-preferences/   (Coot autoloads *.py there;
#           create the folder if it doesn't exist).  Restart bcoot.
# OR run it once from Coot's scripting window (Calculate -> Scripting... -> Python).
#
# USAGE (GUI):  Menu Cootvalent -> Declare covalent link...  Enter the two atom CIDs
#           (Cys SG and the ligand warhead carbon), e.g.
#             SG CID:       //A/481(CYS)/SG
#             warhead CID:  //A/701(8E8)/CAA
#           pick a family, Declare.  Then use the printed refmac command, or
#           call declare_covalent_link(..., mtz=...) to refine in one shot.
#
# USAGE (scripting / socket / MCP):
#           declare_covalent_link(imol, "//A/481(CYS)/SG", "//A/701(8E8)/CAA",
#                                 family="F1", mtz="/path/data.mtz",
#                                 lig_dict="/path/8E8.cif")
#
# Requires CCP4 (refmac5) for the refine step.  ASCII-only.  2026-07-16.

import os, re, subprocess, tempfile, math, json, coot

CCP4_SETUP = "/Applications/ccp4-9/bin/ccp4.setup-sh"
CCP4_PYTHON = "/Applications/ccp4-9/bin/ccp4-python"       # has gemmi (NCS step)
REFMAC_WRAPPER = os.path.expanduser("~/xtal/refmac.sh")   # legacy default


def _gui_fn(name):
    """Resolve a Coot GUI helper (coot_menubar_menu, add_simple_coot_menu_menuitem,
    generic_single_entry, generic_double_entry) by name.

    Do NOT `import coot_gui`: in bandicoot that re-executes coot_gui.py, which
    references startup-only injected globals (e.g. set_found_coot_python_gui)
    and raises NameError. Coot instead exposes these helpers as names in its
    main namespace / builtins, which is where its own extensions find them.
    """
    import sys
    srcs = [globals()]
    try:
        import __main__
        srcs.append(vars(__main__))
    except Exception:
        pass
    import builtins
    srcs.append(vars(builtins))
    cg = sys.modules.get("coot_gui")   # if already loaded cleanly at startup
    if cg is not None:
        srcs.append(vars(cg))
    for s in srcs:
        f = s.get(name)
        if f is not None:
            return f
    return None


def _find_refmac_wrapper():
    """Locate the refmac.sh wrapper: env override, ~/bin, ~/xtal, then PATH."""
    cands = [os.environ.get("COOTVALENT_REFMAC_SH"),
             os.path.expanduser("~/bin/refmac.sh"),
             REFMAC_WRAPPER]
    for c in cands:
        if c and os.path.exists(c):
            return c
    for d in os.environ.get("PATH", "").split(os.pathsep):
        p = os.path.join(d, "refmac.sh")
        if os.path.exists(p):
            return p
    return None

# ---------------------------------------------------------------------------
# Warhead family registry.  Each family carries the canonical S-Cbeta target
# distance/esd + a link-CIF builder.  Kept minimal for the prototype.
# ---------------------------------------------------------------------------
FAMILIES = {
    "F1": {
        "link_id": "CYS-ACR",
        "name": "Cys-S to acrylamide post-Michael adduct (sat. beta-thioether)",
        "target": 1.81, "esd": 0.02, "kind": "sp3",
    },
    "F2": {
        "link_id": "CYS-YNA",
        "name": "Cys-S to alpha,beta-ynamide post-Michael adduct (vinyl thioether)",
        "target": 1.78, "esd": 0.02, "kind": "sp2",
    },
    "CAA": {
        "link_id": "CYS-CAA",
        "name": "Cys-S to alpha-cyano/chloro-acetamide SN2 thioether",
        "target": 1.81, "esd": 0.02, "kind": "sp3",
    },
}
# The auto-detector (covalent-detect.py) labels the SN2 chloro/cyano-acetamide
# class "F3"; alias it to the CAA link builder so detect->declare works for it.
FAMILIES["F3"] = FAMILIES["CAA"]


# ---------------------------------------------------------------------------
# CID parsing.  Accept Coot/Moorhen-style "//CHAIN/RESNO(RESNAME)/ATOM"
# short forms; also tolerate "/1/A/481/SG" long forms.
# ---------------------------------------------------------------------------
def _parse_cid(cid):
    """Return dict(chain, resno, resname, atom) from a CID string. resname may be ''."""
    cid = (cid or "").strip()
    # atom = last path element; may be plain "SG"
    parts = [p for p in cid.split("/")]
    # drop leading empties from //
    parts = [p for p in parts if p != ""]
    # a long form starts with a model number (all digits) - drop it
    if parts and parts[0].isdigit():
        parts = parts[1:]
    if len(parts) < 3:
        raise ValueError("cannot parse CID '%s' (need chain/res/atom)" % cid)
    chain = parts[0]
    resfield = parts[1]
    atom = parts[2]
    resname = ""
    m = re.match(r"^(-?\d+)\s*\(([^)]+)\)\s*$", resfield)
    if m:
        resno = int(m.group(1)); resname = m.group(2).strip()
    else:
        resno = int(re.sub(r"[^\-\d]", "", resfield))
    return {"chain": chain, "resno": resno, "resname": resname, "atom": atom}


# ---------------------------------------------------------------------------
# Read atom coords/names for a residue out of a written PDB (ground truth;
# coot.n_atoms() is 0 headless).  Returns list of (atomname, element, x,y,z, altloc).
# ---------------------------------------------------------------------------
def _residue_atoms(pdb_path, chain, resno, resname=""):
    out = []
    for line in open(pdb_path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != chain:
            continue
        try:
            rs = int(line[22:26])
        except ValueError:
            continue
        if rs != resno:
            continue
        rn = line[17:20].strip()
        if resname and rn != resname:
            continue
        an = line[12:16].strip()
        el = line[76:78].strip() or an[0]
        out.append((an, el, float(line[30:38]), float(line[38:46]),
                    float(line[46:54]), line[16]))
    return out


def _find_atom(atoms, name):
    for a in atoms:
        if a[0] == name and a[5] in (" ", "A"):
            return a
    for a in atoms:
        if a[0] == name:
            return a
    return None


def _dist(a, b):
    return math.sqrt((a[2]-b[2])**2 + (a[3]-b[3])**2 + (a[4]-b[4])**2)


def _neighbours(atoms, center_name, exclude_h=True, cutoff=1.95):
    c = _find_atom(atoms, center_name)
    if c is None:
        return []
    res = []
    for a in atoms:
        if a[0] == center_name:
            continue
        if exclude_h and (a[1] == "H" or a[0].startswith("H")):
            continue
        d = _dist(c, a)
        if 0.4 < d < cutoff:
            res.append((a[0], d))
    res.sort(key=lambda x: x[1])
    return res


# ---------------------------------------------------------------------------
# Build a refmac-ready link CIF for one Cys-warhead bond.
#
# Emits: data_link_list + data_mod_list catalog blocks (with the mod ids),
# a data_link_<ID> block with link_id columns in every loop_, and two
# data_mod_<ID>-modN blocks (delete Cys HG; delete/retype warhead H).
#
# ca = Calpha atom (the warhead C's carbonyl-side neighbour)
# cg = the "third" heavy substituent on Cbeta (methyl for F2 butynamide,
#      or the second chain carbon for CAA); may be None.
# ---------------------------------------------------------------------------
def _build_link_cif(family, lig_comp, cb, ca=None, cg=None,
                    cys_hg="HG", cb_h_to_delete=None, pre_form=False,
                    ca_h_to_add=None):
    """pre_form=True means the ligand's own dict still draws the warhead in the
    unreacted (alkene/alkyne) form, i.e. Ca=Cb is a DOUBLE bond in the dict.
    The mod2 then reconciles it to the post-Michael product:
       F1 (acrylamide): Ca=Cb double -> single, retype both sp2->sp3, +H on Ca
       F2 (ynamide):    Ca#Cb triple -> double, retype sp->sp2, +H on Ca
    Without this, refmac sees the link's SG-Cb bond AND the dict's Ca=Cb double,
    over-valences Cb, and the local geometry is pulled off target."""
    F = FAMILIES[family]
    link_id = F["link_id"]
    mod1 = link_id + "-mod1"
    mod2 = link_id + "-mod2"
    tgt = "%.2f" % F["target"]
    esd = "%.2f" % F["esd"]

    L = []
    # ---- catalog: data_link_list ----
    L.append("data_link_list")
    L.append("loop_")
    for c in ["_chem_link.id", "_chem_link.comp_id_1", "_chem_link.mod_id_1",
              "_chem_link.group_comp_1", "_chem_link.comp_id_2",
              "_chem_link.mod_id_2", "_chem_link.group_comp_2", "_chem_link.name"]:
        L.append(c)
    L.append('%s CYS %s L-peptide %s %s non-polymer "%s"'
             % (link_id, mod1, lig_comp, mod2, F["name"]))
    L.append("")
    # ---- catalog: data_mod_list ----
    L.append("data_mod_list")
    L.append("loop_")
    for c in ["_chem_mod.id", "_chem_mod.name", "_chem_mod.comp_id",
              "_chem_mod.group_id"]:
        L.append(c)
    L.append('%s "%s-side1" CYS L-peptide' % (mod1, link_id))
    L.append('%s "%s-side2" %s non-polymer' % (mod2, link_id, lig_comp))
    L.append("")
    # ---- per-link block ----
    L.append("data_link_%s" % link_id)
    L.append("_chem_link.id           %s" % link_id)
    L.append('_chem_link.name         "%s"' % F["name"])
    L.append("_chem_link.comp_id_1    CYS")
    L.append("_chem_link.mod_id_1     %s" % mod1)
    L.append("_chem_link.group_comp_1 L-peptide")
    L.append("_chem_link.comp_id_2    %s" % lig_comp)
    L.append("_chem_link.mod_id_2     %s" % mod2)
    L.append("_chem_link.group_comp_2 non-polymer")
    L.append("")
    # bonds (link_id column first)
    L.append("loop_")
    for c in ["_chem_link_bond.link_id", "_chem_link_bond.atom_1_comp_id",
              "_chem_link_bond.atom_id_1", "_chem_link_bond.atom_2_comp_id",
              "_chem_link_bond.atom_id_2", "_chem_link_bond.type",
              "_chem_link_bond.value_dist", "_chem_link_bond.value_dist_esd"]:
        L.append(c)
    L.append("%s 1 SG    2 %-4s single %s %s" % (link_id, cb, tgt, esd))
    if F["kind"] == "sp2" and ca and pre_form:
        # vinyl thioether from a PRE-form (alkyne) monomer: the link establishes
        # the sp2 C=C. If the monomer is already post-form (reacted, C=C in its
        # own dict) we must NOT redeclare it here -- that double-defines the bond.
        L.append("%s 2 %-4s 2 %-4s double 1.34 0.02" % (link_id, cb, ca))
    L.append("")
    # angles
    L.append("loop_")
    for c in ["_chem_link_angle.link_id", "_chem_link_angle.atom_1_comp_id",
              "_chem_link_angle.atom_id_1", "_chem_link_angle.atom_2_comp_id",
              "_chem_link_angle.atom_id_2", "_chem_link_angle.atom_3_comp_id",
              "_chem_link_angle.atom_id_3", "_chem_link_angle.value_angle",
              "_chem_link_angle.value_angle_esd"]:
        L.append(c)
    if F["kind"] == "sp2":
        L.append("%s 1 CB    1 SG    2 %-4s 104.2 3.0" % (link_id, cb))
        if ca:
            L.append("%s 1 SG    2 %-4s 2 %-4s 120.7 1.5" % (link_id, cb, ca))
        if cg:
            L.append("%s 1 SG    2 %-4s 2 %-4s 120.3 1.5" % (link_id, cb, cg))
    else:
        L.append("%s 1 CB    1 SG    2 %-4s 100.0 3.0" % (link_id, cb))
        if ca:
            L.append("%s 1 SG    2 %-4s 2 %-4s 113.0 2.0" % (link_id, cb, ca))
        if cg:
            L.append("%s 1 SG    2 %-4s 2 %-4s 113.0 2.0" % (link_id, cb, cg))
    L.append("")
    # torsions (soft hinges)
    L.append("loop_")
    for c in ["_chem_link_tor.link_id", "_chem_link_tor.id",
              "_chem_link_tor.atom_1_comp_id", "_chem_link_tor.atom_id_1",
              "_chem_link_tor.atom_2_comp_id", "_chem_link_tor.atom_id_2",
              "_chem_link_tor.atom_3_comp_id", "_chem_link_tor.atom_id_3",
              "_chem_link_tor.atom_4_comp_id", "_chem_link_tor.atom_id_4",
              "_chem_link_tor.value_angle", "_chem_link_tor.value_angle_esd",
              "_chem_link_tor.period"]:
        L.append(c)
    if ca:
        if F["kind"] == "sp2":
            L.append("%s hinge 1 CA 1 CB 1 SG 2 %-4s 180.0 20.0 3" % (link_id, cb))
        else:
            L.append("%s hinge 1 CA 1 CB 1 SG 2 %-4s 180.0 20.0 3" % (link_id, cb))
    L.append("")
    # planarity (sp2 vinyl thioether only)
    if F["kind"] == "sp2" and ca and cg:
        L.append("loop_")
        for c in ["_chem_link_plane.link_id", "_chem_link_plane.plane_id",
                  "_chem_link_plane.atom_comp_id", "_chem_link_plane.atom_id",
                  "_chem_link_plane.dist_esd"]:
            L.append(c)
        L.append("%s PLN_VINYL 1 SG   0.02" % link_id)
        L.append("%s PLN_VINYL 2 %-4s 0.02" % (link_id, cb))
        L.append("%s PLN_VINYL 2 %-4s 0.02" % (link_id, ca))
        L.append("%s PLN_VINYL 2 %-4s 0.02" % (link_id, cg))
        L.append("")

    # ---- mod1: delete Cys HG ----
    L.append("data_mod_%s" % mod1)
    L.append("loop_")
    for c in ["_chem_mod_atom.mod_id", "_chem_mod_atom.function",
              "_chem_mod_atom.atom_id", "_chem_mod_atom.new_atom_id",
              "_chem_mod_atom.new_type_symbol", "_chem_mod_atom.new_type_energy",
              "_chem_mod_atom.new_partial_charge"]:
        L.append(c)
    L.append("%s delete %s . . . ." % (mod1, cys_hg))
    L.append("")
    # ---- mod2: reconcile the ligand's warhead to the post-Michael product ----
    # Two things can happen on the ligand side:
    #   (a) delete the warhead H that SG replaces (post-form dicts), and/or
    #   (b) pre_form reconciliation: change Ca=Cb order + retype sp2->sp3 (F1)
    #       or Ca#Cb -> Ca=Cb + retype sp->sp2 (F2), and add an H on Ca.
    L.append("data_mod_%s" % mod2)
    atom_rows = []
    bond_rows = []
    if cb_h_to_delete:
        atom_rows.append("%s delete %s . . . ." % (mod2, cb_h_to_delete))
    if pre_form and ca:
        if F["kind"] == "sp3":
            # F1 acrylamide: Ca=Cb double -> single, both sp2->sp3
            bond_rows.append("%s change %-4s %-4s single 1.54 0.02"
                             % (mod2, ca, cb))
            atom_rows.append("%s change %-4s . C CT ." % (mod2, ca))
            atom_rows.append("%s change %-4s . C CT ." % (mod2, cb))
        else:
            # F2 ynamide: Ca#Cb triple -> double, both sp->sp2
            bond_rows.append("%s change %-4s %-4s double 1.34 0.02"
                             % (mod2, ca, cb))
            atom_rows.append("%s change %-4s . C C2 ." % (mod2, ca))
            atom_rows.append("%s change %-4s . C C2 ." % (mod2, cb))
        if ca_h_to_add:
            atom_rows.append("%s add %-4s %-4s H HCH1 ."
                             % (mod2, ca_h_to_add, ca_h_to_add))
    if bond_rows:
        L.append("loop_")
        for c in ["_chem_mod_bond.mod_id", "_chem_mod_bond.function",
                  "_chem_mod_bond.atom_id_1", "_chem_mod_bond.atom_id_2",
                  "_chem_mod_bond.new_type", "_chem_mod_bond.new_value_dist",
                  "_chem_mod_bond.new_value_dist_esd"]:
            L.append(c)
        L.extend(bond_rows)
    if atom_rows:
        L.append("loop_")
        for c in ["_chem_mod_atom.mod_id", "_chem_mod_atom.function",
                  "_chem_mod_atom.atom_id", "_chem_mod_atom.new_atom_id",
                  "_chem_mod_atom.new_type_symbol",
                  "_chem_mod_atom.new_type_energy",
                  "_chem_mod_atom.new_partial_charge"]:
            L.append(c)
        L.extend(atom_rows)
    if not bond_rows and not atom_rows:
        L.append("# no atom change on the ligand side")
    L.append("")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# PDB LINK / LINKR record.  Strict column placement (PDB v3.3); the link-id
# suffix after col 80 is what refmac matches templates by (LINKR extension).
# ---------------------------------------------------------------------------
def _link_record(sg, wh, link_id):
    """sg, wh: parsed-CID dicts with chain/resno/resname/atom."""
    def field(atom, res):
        return (atom, res["resname"] or "LIG", res["chain"], res["resno"])
    a1, r1, c1, n1 = field(sg["atom"], sg)
    a2, r2, c2, n2 = field(wh["atom"], wh)
    buf = [" "] * 81
    def put(s, start):  # 1-based columns like the PDB spec
        for i, ch in enumerate(s):
            buf[start - 1 + i] = ch
    # Atom-name field (cols 13-16) uses the PDB element-justification rule:
    # a name <=3 chars starts at col 14 (leading space in col 13), matching the
    # deposited LINK format cross-checked against 5P9J.pdb / 8FD9.pdb.  This is
    # the off-by-one the covalent-link memory warns about -- refmac is lenient
    # but we emit the canonical placement anyway.
    def name_field(a):
        return a if len(a) >= 4 else (" %-3s" % a)
    put("LINK", 1)
    put(name_field(a1), 13)         # name1 cols 13-16
    put(("%-3s" % r1), 18)          # resName1 18-20
    put(c1[:1], 22)                 # chainID1 22
    put(("%4d" % n1), 23)           # resSeq1 23-26
    put(name_field(a2), 43)         # name2 43-46
    put(("%-3s" % r2), 48)          # resName2 48-50
    put(c2[:1], 52)                 # chainID2 52
    put(("%4d" % n2), 53)           # resSeq2 53-56
    put("1555", 60)                 # sym1 60-65
    put("1555", 67)                 # sym2 67-72
    line = "".join(buf).rstrip()
    # LINKR link-id suffix after col 80 (whitespace-delimited)
    line = ("%-80s %s" % (line, link_id))
    return line


def _insert_link_record(pdb_in, pdb_out, link_line):
    """Write a copy of pdb_in with link_line inserted (after existing LINKs,
    else before the first ATOM/HETATM)."""
    lines = open(pdb_in).read().splitlines()
    out = []
    inserted = False
    # insert after last LINK if any, else before first coordinate record
    last_link = -1
    for i, l in enumerate(lines):
        if l.startswith("LINK"):
            last_link = i
    if last_link >= 0:
        for i, l in enumerate(lines):
            out.append(l)
            if i == last_link:
                out.append(link_line); inserted = True
    else:
        for l in lines:
            if not inserted and l.startswith(("ATOM", "HETATM")):
                out.append(link_line); inserted = True
            out.append(l)
    if not inserted:
        out.insert(0, link_line)
    open(pdb_out, "w").write("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# Ligand-dict inspection: is the warhead drawn in its pre-reaction form, and
# which H's sit on Cbeta / Calpha?  Drives the mod2 reconciliation.
# ---------------------------------------------------------------------------
def _monomer_lib_path(comp):
    """Best-effort path to a ligand's dict in the CCP4 monomer library."""
    root = os.environ.get("CLIBD_MON")
    if not root:
        # common CCP4-9 location
        cand = "/Applications/ccp4-9/lib/data/monomers/"
        root = cand if os.path.isdir(cand) else None
    if not root:
        return None
    p = os.path.join(root, comp[0].lower(), comp + ".cif")
    return p if os.path.exists(p) else None


def _analyse_ligand_dict(dict_path, cb, ca, family):
    """Return (pre_form, cb_h_to_delete, ca_h_to_add).

    pre_form      : True if the Ca-Cb bond in the dict is double (F1) or
                    triple (F2) -- i.e. the unreacted warhead is still drawn.
    cb_h_to_delete: an H bonded to Cb that SG will replace (post-form only;
                    for pre-form the retype handles valence so we don't delete).
    ca_h_to_add   : a synthetic H id to add on Ca (pre-form only).
    """
    if not dict_path or not os.path.exists(dict_path):
        return False, None, None
    bonds = []      # (a1, a2, order)
    try:
        in_bond = False
        cols = []
        for line in open(dict_path):
            s = line.strip()
            if s.startswith("_chem_comp_bond."):
                in_bond = True
                cols.append(s.split(".", 1)[1])
                continue
            if in_bond:
                if s.startswith("_") or s == "loop_":
                    continue
                if s.startswith("data_") or s.startswith("#") or s == "":
                    if bonds:
                        in_bond = False
                    continue
                parts = s.split()
                # rows look like: COMP a1 a2 order ...
                if len(parts) >= 4:
                    a1, a2, order = parts[1], parts[2], parts[3]
                    bonds.append((a1, a2, order.upper()))
    except Exception:
        return False, None, None

    def order_between(x, y):
        for a1, a2, o in bonds:
            if {a1, a2} == {x, y}:
                return o
        return None

    # "pre_form" = the Cb-Ca bond is still drawn in its UNREACTED order, which is
    # family-specific: F2 ynamide is unreacted at TRIPLE (reacts to DOUBLE), F1
    # acrylamide is unreacted at DOUBLE (reacts to SINGLE). A blanket "any
    # double/triple => pre_form" would wrongly flag a REACTED F2 vinyl (DOUBLE)
    # as unreacted and then double-define the bond.
    pre_form = False
    if ca:
        o = order_between(cb, ca)
        kind = FAMILIES.get(family, {}).get("kind")
        if kind == "sp2":            # F2: unreacted == triple
            pre_form = (o == "TRIPLE")
        elif kind == "sp3":          # F1 / CAA: unreacted == double
            pre_form = o in ("DOUBLE", "AROM", "AROMATIC")
        else:                        # unknown: fall back to the old heuristic
            pre_form = o in ("DOUBLE", "TRIPLE", "AROM", "AROMATIC")

    # Hs on Cb / Ca
    def hs_on(c):
        out = []
        for a1, a2, o in bonds:
            if a1 == c and (a2.startswith("H")):
                out.append(a2)
            elif a2 == c and (a1.startswith("H")):
                out.append(a1)
        return out
    cb_hs = hs_on(cb)
    cb_h_del = cb_hs[0] if (cb_hs and not pre_form) else None
    # for a terminal pre-form warhead (Cb has 2 H) SG replaces one; still no
    # explicit delete needed because retype re-derives Cb valence, but if the
    # post-form dict simply has spare H's we delete one.
    ca_h_add = None
    if pre_form and ca:
        # invent an H id that doesn't collide
        base_h = "H%sX" % ca.replace("C", "")
        existing = set(x for c in (cb, ca) for x in hs_on(c))
        ca_h_add = base_h if base_h not in existing else base_h + "1"
    return pre_form, cb_h_del, ca_h_add


# ---------------------------------------------------------------------------
# Pre-form library ligand workaround.
#
# If a ligand is present in the CCP4 monomer library in its UNREACTED
# (pre-Michael, alkene/alkyne) form -- e.g. YY3 (osimertinib) ships with
# C8=C9 drawn double -- refmac prefers that library entry over any post-form
# dict you pass in LIBIN (same comp id => library wins), so the internal
# Ca=Cb double bond fights the SG-Cb link and the geometry lands ~0.1 A long.
#
# The reliable fix is to rename that ligand copy to a NON-library comp id in
# the model, and hand refmac a post-Michael (saturated sp3 / vinyl sp2) dict
# under the new id.  MAKE NEWLIGAND NOEXIT (set in the keyword file) then lets
# refmac accept the "new" ligand.  This is the native analogue of PyKeko's
# JS-side mod2 transform, which rewrites the ligand chem_comp before handoff.
#
# rename_ligand_copy() rewrites one residue's comp id in a PDB (+ its LINK
# records) so declare_covalent_link(... warhead_cid pointing at the new id,
# lig_dict=<post-form dict under the new id> ...) then Just Works.
# ---------------------------------------------------------------------------
def rename_ligand_copy(pdb_in, pdb_out, chain, resno, old_comp, new_comp):
    """Rewrite residue (chain,resno,old_comp) -> new_comp in a PDB + LINKs."""
    new3 = ("%-3s" % new_comp)[:3]
    out = []
    for line in open(pdb_in):
        if line.startswith(("ATOM", "HETATM")):
            if (line[21] == chain and line[17:20].strip() == old_comp):
                try:
                    if int(line[22:26]) == resno:
                        line = line[:17] + new3 + line[20:]
                except ValueError:
                    pass
        elif line.startswith("LINK"):
            tag = "%s %s%4d" % (old_comp, chain, resno)
            if tag in line:
                line = line.replace(old_comp, new3, 1)
        out.append(line)
    open(pdb_out, "w").write("".join(out))


def rename_comp_in_dict(dict_in, dict_out, old_comp, new_comp):
    """Rewrite every occurrence of a comp id in a monomer dict."""
    txt = open(dict_in).read()
    txt = txt.replace("comp_%s" % old_comp, "comp_%s" % new_comp)
    # data rows begin with the comp id token
    txt = re.sub(r"(?m)^(\s*)%s\b" % re.escape(old_comp),
                 r"\g<1>%s" % new_comp, txt)
    open(dict_out, "w").write(txt)


# ---------------------------------------------------------------------------
# Main callable.
# ---------------------------------------------------------------------------
def declare_covalent_link(imol, sg_cid, warhead_cid, family="F1",
                          mtz=None, lig_dict=None, workdir=None,
                          ncyc=10, do_refine=True,
                          use_wrapper=False, add_waters=False):
    """Declare a Cys-warhead covalent link on molecule imol.

    sg_cid, warhead_cid : CID strings //CHAIN/RESNO(RESNAME)/ATOM
    family              : "F1" (acrylamide), "F2" (ynamide), "CAA" (cyano/chloro)
    mtz                 : if given (and do_refine), run refmac5 and load the result
    lig_dict            : path to the ligand's own monomer dict (.cif). If None,
                          refmac falls back to the CCP4 monomer lib.
    Returns dict with paths (augmented_pdb, link_cif, refined_pdb, refined_mtz,
    sg_c_refined) or raises on hard failure.  Safe from console/socket/MCP.
    """
    if family not in FAMILIES:
        raise ValueError("unknown family '%s' (use %s)"
                         % (family, "/".join(FAMILIES)))
    sg = _parse_cid(sg_cid)
    wh = _parse_cid(warhead_cid)
    wd = workdir or tempfile.mkdtemp(prefix="pk_cov_")
    base = os.path.join(wd, "model")

    # 1) dump the current model to PDB (ground truth for atom names/coords)
    src_pdb = base + "_src.pdb"
    coot.write_pdb_file(imol, src_pdb)

    # 2) resolve residue names + warhead neighbourhood from coords
    sg_atoms = _residue_atoms(src_pdb, sg["chain"], sg["resno"], sg["resname"])
    wh_atoms = _residue_atoms(src_pdb, wh["chain"], wh["resno"], wh["resname"])
    if not sg_atoms:
        raise RuntimeError("no atoms for Cys %s/%d" % (sg["chain"], sg["resno"]))
    if not wh_atoms:
        raise RuntimeError("no atoms for ligand %s/%d" % (wh["chain"], wh["resno"]))
    if not sg["resname"]:
        # fill from file
        sg["resname"] = sg_atoms[0][0] and _res_name(src_pdb, sg)
    if not wh["resname"]:
        wh["resname"] = _res_name(src_pdb, wh)
    lig_comp = wh["resname"]

    # warhead-carbon neighbours (heavy) -> pick Calpha (carbonyl side) + Cgamma
    nbrs = _neighbours(wh_atoms, wh["atom"])
    carbons = [n for n in nbrs if n[0][0] == "C" or n[0][0] == "c"]
    ca = carbons[0][0] if len(carbons) >= 1 else None
    cg = carbons[1][0] if len(carbons) >= 2 else None

    # measure deposited SG-warhead distance (informational)
    sgA = _find_atom(sg_atoms, sg["atom"])
    whA = _find_atom(wh_atoms, wh["atom"])
    dep_d = _dist(sgA, whA) if (sgA and whA) else None

    # 3) inspect the ligand's own dict to decide whether it is pre-form
    #    (Ca=Cb still drawn double/triple) and which warhead H's to delete/add.
    dict_path = lig_dict or _monomer_lib_path(lig_comp)
    pre_form, cb_h_del, ca_h_add = _analyse_ligand_dict(
        dict_path, wh["atom"], ca, family)

    # 4) build the link CIF
    link = FAMILIES[family]["link_id"]
    link_cif = base + "_link.cif"
    cif_text = _build_link_cif(family, lig_comp, wh["atom"], ca=ca, cg=cg,
                               pre_form=pre_form, cb_h_to_delete=cb_h_del,
                               ca_h_to_add=ca_h_add)
    open(link_cif, "w").write(cif_text)

    # 5) augmented model with LINK/LINKR record
    aug_pdb = base + "_augmented.pdb"
    link_line = _link_record(sg, wh, link)
    _insert_link_record(src_pdb, aug_pdb, link_line)

    # 6) mmCIF twin with a _struct_conn row (parity with PyKeko surgery output)
    aug_mmcif = base + "_augmented.mmcif"
    _write_struct_conn_mmcif(imol, aug_mmcif, sg, wh, link, base)

    result = {
        "augmented_pdb": aug_pdb, "augmented_mmcif": aug_mmcif,
        "link_cif": link_cif, "link_id": link, "family": family,
        "lig_comp": lig_comp, "warhead_atom": wh["atom"],
        "ca": ca, "cg": cg, "deposited_sg_c": dep_d,
        "refined_pdb": None, "refined_mtz": None, "sg_c_refined": None,
    }

    print("[cootvalent] declared %s: SG(%s/%d) <-> %s(%s/%d)  family=%s target=%.2f A"
          % (link, sg["chain"], sg["resno"], wh["atom"], wh["chain"],
             wh["resno"], family, FAMILIES[family]["target"]))
    if dep_d is not None:
        print("[cootvalent]   current model SG-%s distance: %.3f A" % (wh["atom"], dep_d))
    print("[cootvalent]   augmented PDB: %s" % aug_pdb)
    print("[cootvalent]   link CIF     : %s" % link_cif)

    # 6) refine with external refmac5 (inline) or the refmac.sh wrapper (+waters)
    if do_refine and mtz:
        if use_wrapper:
            rp, rm, sgc = refmac_refine_via_wrapper(
                aug_pdb, mtz, link_cif, lig_dict, sg, wh, wd,
                ncyc=ncyc, add_waters=add_waters)
        else:
            rp, rm, sgc = refmac_refine(aug_pdb, mtz, link_cif, lig_dict,
                                        sg, wh, wd, ncyc=ncyc)
        result["refined_pdb"] = rp
        result["refined_mtz"] = rm
        result["sg_c_refined"] = sgc
        if rp and coot.is_valid_model_molecule(coot.read_pdb(rp)) is not None:
            pass
    elif do_refine and not mtz:
        cmd = _refmac_cmd(aug_pdb, "<YOUR.mtz>", link_cif, lig_dict, base, ncyc)
        print("[cootvalent]   to refine, supply an MTZ. Example refmac command:")
        print("           " + cmd)

    return result


def _res_name(pdb, cid):
    for line in open(pdb):
        if line.startswith(("ATOM", "HETATM")) and line[21] == cid["chain"]:
            try:
                if int(line[22:26]) == cid["resno"]:
                    return line[17:20].strip()
            except ValueError:
                pass
    return "LIG"


def _write_struct_conn_mmcif(imol, path, sg, wh, link, base):
    """Write an mmCIF of the model + a _struct_conn covalent row (PyKeko parity).
    Best-effort: if coot can't write mmCIF, skip silently."""
    try:
        coot.write_cif_file(imol, path)
    except Exception:
        try:
            coot.mmcif_file_name_to_molecule  # noqa
        except Exception:
            return
        return
    block = [
        "", "loop_",
        "_struct_conn.id", "_struct_conn.conn_type_id",
        "_struct_conn.ptnr1_auth_asym_id", "_struct_conn.ptnr1_auth_seq_id",
        "_struct_conn.ptnr1_label_comp_id", "_struct_conn.ptnr1_label_atom_id",
        "_struct_conn.ptnr2_auth_asym_id", "_struct_conn.ptnr2_auth_seq_id",
        "_struct_conn.ptnr2_label_comp_id", "_struct_conn.ptnr2_label_atom_id",
        "_struct_conn.ccp4_link_id",
        "covale1 covale %s %d %s %s %s %d %s %s %s"
        % (sg["chain"], sg["resno"], sg["resname"], sg["atom"],
           wh["chain"], wh["resno"], wh["resname"], wh["atom"], link),
        "",
    ]
    try:
        with open(path, "a") as fh:
            fh.write("\n".join(block) + "\n")
    except Exception:
        pass


def _refmac_cmd(pdb, mtz, link_cif, lig_dict, base, ncyc):
    if os.path.exists(REFMAC_WRAPPER):
        extra = (' "%s"' % lig_dict) if lig_dict else ""
        return ('%s "%s" "%s" "%s"%s -c %d -o "%s_refined"'
                % (REFMAC_WRAPPER, pdb, mtz, link_cif, extra, ncyc, base))
    return "refmac5 XYZIN %s HKLIN %s LIBIN %s ..." % (pdb, mtz, link_cif)


def refmac_refine(pdb, mtz, link_cif, lig_dict, sg, wh, wd, ncyc=10):
    """Spawn external refmac5.  Merge link CIF + ligand dict into one LIBIN,
    run, measure the refined SG-warhead distance.  Returns (pdb, mtz, dist)."""
    # merge link CIF + ligand dict into a single LIBIN (refmac reads one LIBIN)
    libin = os.path.join(wd, "libin.cif")
    parts = [open(link_cif).read()]
    if lig_dict and os.path.exists(lig_dict):
        parts.append(open(lig_dict).read())
    open(libin, "w").write("\n\n".join(parts) + "\n")

    out_base = os.path.join(wd, "refined")
    out_pdb = out_base + ".pdb"
    out_mtz = out_base + ".mtz"
    keywords = os.path.join(wd, "refmac.keys")
    # detect free-R label
    labin = _guess_labin(mtz)
    # MAKE NEWLIGAND NOEXIT: don't bail when a renamed/non-library ligand id is
    #   seen (we sometimes rename a library ligand so our post-Michael dict wins
    #   over CLIBD's pre-form entry -- see the covalent-link README notes).
    # MAKE CHECK NONE: trust our supplied dict geometry.
    open(keywords, "w").write(
        "LABIN %s\nMAKE HYDR NO\nMAKE CHECK NONE\nMAKE NEWLIGAND NOEXIT\n"
        "NCYC %d\nEND\n" % (labin, ncyc))
    cmd = ('export PYTHONNOUSERSITE=1; . "%s" >/dev/null 2>&1; '
           'refmac5 XYZIN "%s" HKLIN "%s" XYZOUT "%s" HKLOUT "%s" '
           'LIBIN "%s" < "%s" > "%s/refmac.log" 2>&1'
           % (CCP4_SETUP, pdb, mtz, out_pdb, out_mtz, libin, keywords, wd))
    coot.add_status_bar_text("Running refmac5 (%d cycles)..." % ncyc)
    rc = subprocess.call(["bash", "-c", cmd])
    if rc != 0 or not os.path.exists(out_pdb):
        tail = ""
        try:
            tail = open(os.path.join(wd, "refmac.log")).read()[-1200:]
        except Exception:
            pass
        print("[cootvalent] refmac5 FAILED (rc=%d)\n%s" % (rc, tail))
        return None, None, None

    # measure refined SG-warhead distance
    sga = _residue_atoms(out_pdb, sg["chain"], sg["resno"], sg["resname"])
    wha = _residue_atoms(out_pdb, wh["chain"], wh["resno"], wh["resname"])
    a = _find_atom(sga, sg["atom"]); b = _find_atom(wha, wh["atom"])
    d = _dist(a, b) if (a and b) else None
    print("[cootvalent] refmac5 done -> %s" % out_pdb)
    if d is not None:
        print("[cootvalent]   refined SG-%s distance: %.3f A" % (wh["atom"], d))
    # load refined model back into Coot
    try:
        coot.read_pdb(out_pdb)
        coot.make_and_draw_map(out_mtz, "FWT", "PHWT", "", 0, 0)
    except Exception as e:
        print("[cootvalent]   (load-back skipped:", e, ")")
    return out_pdb, out_mtz, d


def _guess_labin(mtz):
    """Read MTZ column labels via mtzdump-free header scan; return a LABIN string."""
    # Try to read column names from the mtz binary header cheaply.
    try:
        data = open(mtz, "rb").read()
        txt = data.decode("latin-1", "ignore")
        cols = set(re.findall(r"COLUMN\s+(\S+)", txt))
        def pick(cands, default):
            for c in cands:
                if c in cols:
                    return c
            return default
        fp = pick(["FP", "F", "F_obs", "FOBS", "F-obs"], "FP")
        sig = pick(["SIGFP", "SIGF", "SIGF_obs", "SIGFOBS"], "SIGFP")
        free = pick(["FREE", "FREER", "FreeR_flag", "FREER_flag", "R-free-flags"], "FREER")
        return "FP=%s SIGFP=%s FREE=%s" % (fp, sig, free)
    except Exception:
        return "FP=FP SIGFP=SIGFP FREE=FREER"


# ---------------------------------------------------------------------------
# Refine via the refmac.sh wrapper (adds automatic water picking with -W and
# the wrapper's tested keyword set) instead of the inline refmac5 call.
# Falls back to the inline path if the wrapper isn't found.
# ---------------------------------------------------------------------------
def refmac_refine_via_wrapper(pdb, mtz, link_cif, lig_dict, sg, wh, wd,
                              ncyc=10, add_waters=False, wrapper=None):
    """Merge link CIF + ligand dict into one LIBIN and drive refmac.sh."""
    wrapper = wrapper or _find_refmac_wrapper()
    if not wrapper:
        print("[cootvalent] refmac.sh not found (set COOTVALENT_REFMAC_SH); "
              "using inline refmac5.")
        return refmac_refine(pdb, mtz, link_cif, lig_dict, sg, wh, wd, ncyc=ncyc)

    libin = os.path.join(wd, "libin.cif")
    parts = [open(link_cif).read()]
    if lig_dict and os.path.exists(lig_dict):
        parts.append(open(lig_dict).read())
    open(libin, "w").write("\n\n".join(parts) + "\n")

    out_base = os.path.join(wd, "refined_cov")
    out_pdb, out_mtz = out_base + ".pdb", out_base + ".mtz"
    log = os.path.join(wd, "refmac_wrapper.log")
    cmd = [wrapper, pdb, mtz, libin, "-o", out_base, "-c", str(ncyc),
           "-L", _guess_labin(mtz)]
    if add_waters:
        cmd.append("-W")
    # MAKE NEWLIGAND NOEXIT: accept a renamed/non-library ligand id (post-Michael
    # dict) without bailing -- same reason as the inline path.
    cmd += ["--", "MAKE NEWLIGAND NOEXIT"]
    coot.add_status_bar_text("Running refmac.sh (%d cyc%s)..."
                             % (ncyc, ", +waters" if add_waters else ""))
    with open(log, "w") as lf:
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0 or not os.path.exists(out_pdb):
        tail = ""
        try:
            tail = open(log).read()[-1500:]
        except Exception:
            pass
        print("[cootvalent] refmac.sh FAILED (rc=%d)\n%s" % (rc, tail))
        return None, None, None

    sga = _residue_atoms(out_pdb, sg["chain"], sg["resno"], sg["resname"])
    wha = _residue_atoms(out_pdb, wh["chain"], wh["resno"], wh["resname"])
    a = _find_atom(sga, sg["atom"]); b = _find_atom(wha, wh["atom"])
    d = _dist(a, b) if (a and b) else None
    print("[cootvalent] refmac.sh done -> %s" % out_pdb)
    if d is not None:
        print("[cootvalent]   refined SG-%s distance: %.3f A" % (wh["atom"], d))
    try:
        coot.read_pdb(out_pdb)
        coot.make_and_draw_map(out_mtz, "FWT", "PHWT", "", 0, 0)
    except Exception as e:
        print("[cootvalent]   (load-back skipped:", e, ")")
    return out_pdb, out_mtz, d


# ---------------------------------------------------------------------------
# NCS propagation: given one built+placed covalent ligand (the reference copy),
# place symmetry-equivalent copies at the same Cys in every other protein chain
# of the ASU, by superposing the reference chain onto each target chain and
# applying that transform to the ligand.  Deterministic geometry (gemmi via
# ccp4-python) -- but each copy should still be eyeballed against its density.
# ---------------------------------------------------------------------------
def _parse_res_cid(cid):
    """Parse a residue CID //CHAIN/RESNO(RESNAME)[/ATOM]; atom optional."""
    parts = [p for p in (cid or "").split("/") if p != ""]
    if parts and parts[0].isdigit():
        parts = parts[1:]
    if len(parts) < 2:
        raise ValueError("cannot parse residue CID '%s'" % cid)
    d = {"chain": parts[0], "resname": "", "atom": None}
    m = re.match(r"^(-?\d+)\s*\(([^)]+)\)\s*$", parts[1])
    if m:
        d["resno"] = int(m.group(1)); d["resname"] = m.group(2).strip()
    else:
        d["resno"] = int(re.sub(r"[^\-\d]", "", parts[1]))
    if len(parts) >= 3:
        d["atom"] = parts[2]
    return d


def _cys_sg_by_chain(pdb, resno, cys_atom="SG"):
    """Return {chain: (x,y,z)} for every CYS <resno> SG in the model."""
    out = {}
    for line in open(pdb):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != "CYS":
            continue
        try:
            if int(line[22:26]) != resno:
                continue
        except ValueError:
            continue
        if line[12:16].strip() != cys_atom:
            continue
        out[line[21]] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out


def _resno_used(pdb, chain, resno):
    for line in open(pdb):
        if line.startswith(("ATOM", "HETATM")) and line[21] == chain:
            try:
                if int(line[22:26]) == resno:
                    return True
            except ValueError:
                pass
    return False


# gemmi worker (runs under ccp4-python): superpose ref chain -> each target
# chain and write the transformed ligand copies into one PDB.
_NCS_WORKER_SRC = r'''
import sys, json, gemmi
job = json.load(sys.stdin)
st = gemmi.read_structure(job["model"])
model = st[0]

def get_chain(name):
    for ch in model:
        if ch.name == name:
            return ch
    return None

def ca_positions(ch):
    d = {}
    for res in ch:
        for at in res:
            if at.name == "CA":
                d[res.seqid.num] = at.pos
                break
    return d

def find_res(ch, resno, resname):
    for res in ch:
        if res.seqid.num == resno and (not resname or res.name == resname):
            return res
    return None

ref = get_chain(job["ref_chain"])
ligch = get_chain(job["lig_chain"])
if ref is None or ligch is None:
    print(json.dumps({"error": "ref/lig chain not found"})); sys.exit(0)
lig = find_res(ligch, job["lig_resno"], job.get("lig_resname", ""))
if lig is None:
    print(json.dumps({"error": "ligand residue not found"})); sys.exit(0)

refca = ca_positions(ref)
center = job["cys_resno"]; window = job.get("window", 0)
out_st = gemmi.Structure()
out_st.cell = st.cell
out_st.spacegroup_hm = st.spacegroup_hm
om = gemmi.Model("1")
report = []

for tgt in job["targets"]:
    tch = get_chain(tgt["tgt_chain"])
    if tch is None:
        report.append({"tgt_chain": tgt["tgt_chain"], "error": "chain missing"}); continue
    tgtca = ca_positions(tch)
    def pairs(win):
        fx = []; mv = []
        for num, pos in refca.items():
            if win and abs(num - center) > win:
                continue
            if num in tgtca:
                mv.append(pos); fx.append(tgtca[num])
        return fx, mv
    fx, mv = pairs(window)
    if len(fx) < 3:
        fx, mv = pairs(0)   # fall back to all matched CA
    if len(fx) < 3:
        report.append({"tgt_chain": tgt["tgt_chain"],
                       "error": "too few matched CA (%d)" % len(fx)}); continue
    sup = gemmi.superpose_positions(fx, mv)   # maps mv (ref) -> fx (target)
    T = sup.transform
    ch = gemmi.Chain(tgt["new_chain"])
    nres = gemmi.Residue()
    nres.name = lig.name
    nres.seqid = gemmi.SeqId(tgt["new_resno"], " ")
    nres.het_flag = "H"
    for at in lig:
        v = T.apply(at.pos)
        na = gemmi.Atom()
        na.name = at.name; na.element = at.element
        na.pos = gemmi.Position(v.x, v.y, v.z)
        na.occ = at.occ; na.b_iso = at.b_iso; na.altloc = at.altloc
        nres.add_atom(na)
    ch.add_residue(nres)
    om.add_chain(ch)
    report.append({"tgt_chain": tgt["tgt_chain"], "new_chain": tgt["new_chain"],
                   "new_resno": tgt["new_resno"], "rmsd": round(sup.rmsd, 3),
                   "n_ca": len(fx)})

out_st.add_model(om)
out_st.write_pdb(job["out_pdb"])
print(json.dumps({"copies": report}))
'''


def _run_ncs_worker(job):
    """Run the gemmi superposition worker under ccp4-python; return its JSON."""
    src = os.path.join(job["_wd"], "ncs_worker.py")
    open(src, "w").write(_NCS_WORKER_SRC)
    cmd = ('. "%s" >/dev/null 2>&1; exec "%s" "%s"'
           % (CCP4_SETUP, CCP4_PYTHON, src))
    p = subprocess.Popen(["bash", "-c", cmd], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate(json.dumps(job).encode("utf-8"), timeout=120)
    try:
        return json.loads(out.decode("utf-8", "ignore").strip().splitlines()[-1])
    except Exception:
        return {"error": "worker failed: %s" % err.decode("utf-8", "ignore")[-500:]}


def _hetatm_lines(pdb):
    return [l for l in open(pdb).read().splitlines()
            if l.startswith(("ATOM", "HETATM"))]


def propagate_covalent_ncs(imol, lig_cid, cys_resno, family="F2",
                           target_chains=None, ref_chain=None, window=15,
                           warhead_atom=None, cys_atom="SG",
                           mtz=None, lig_dict=None, workdir=None, ncyc=10,
                           do_refine=False, add_waters=False, use_wrapper=True):
    """Propagate a built covalent ligand to the same Cys in every ASU chain.

    lig_cid   : residue CID of the reference (already placed) ligand,
                e.g. //A/601(LIG)
    cys_resno : the reacting Cys residue number (same in every chain)
    Returns a dict with the augmented PDB, per-copy report, and (if refined)
    the refined model paths.
    """
    if family not in FAMILIES:
        raise ValueError("unknown family '%s'" % family)
    wd = workdir or tempfile.mkdtemp(prefix="pk_ncs_")
    src_pdb = os.path.join(wd, "model_src.pdb")
    coot.write_pdb_file(imol, src_pdb)

    lig = _parse_res_cid(lig_cid)
    if not lig["resname"]:
        lig["resname"] = _res_name(src_pdb, lig)
    lig_atoms = _residue_atoms(src_pdb, lig["chain"], lig["resno"], lig["resname"])
    if not lig_atoms:
        raise RuntimeError("reference ligand %s/%d not found"
                           % (lig["chain"], lig["resno"]))

    sg_by_chain = _cys_sg_by_chain(src_pdb, cys_resno, cys_atom)
    if not sg_by_chain:
        raise RuntimeError("no CYS %d /%s found in any chain" % (cys_resno, cys_atom))

    # warhead atom: given, else the ligand heavy atom closest to any SG
    if warhead_atom is None:
        # warhead = the ligand CARBON closest to any Cys SG (Cys-S adducts bond
        # through carbon; no distance window -- just the nearest carbon).
        best = None
        for a in lig_atoms:
            elem = (a[1] or a[0][:1]).strip()
            if elem not in ("C", "c"):
                continue
            for sgxyz in sg_by_chain.values():
                d = math.sqrt(sum((a[2 + i] - sgxyz[i]) ** 2 for i in range(3)))
                if best is None or d < best[1]:
                    best = (a[0], d)
        if best is None:
            raise RuntimeError("could not infer warhead carbon; pass warhead_atom=")
        warhead_atom = best[0]
        print("[cootvalent] inferred warhead carbon = %s (%.2f A from an SG)"
              % (warhead_atom, best[1]))

    # reference chain = CYS chain whose SG is closest to the reference warhead
    wa = _find_atom(lig_atoms, warhead_atom)
    if ref_chain is None:
        ref_chain = min(sg_by_chain,
                        key=lambda c: math.sqrt(sum((wa[2 + i] - sg_by_chain[c][i]) ** 2
                                                    for i in range(3))))
    # targets = every other CYS-bearing chain (or the user's subset)
    tgts = target_chains or [c for c in sorted(sg_by_chain) if c != ref_chain]
    tgts = [c for c in tgts if c != ref_chain and c in sg_by_chain]
    if not tgts:
        print("[cootvalent] no other CYS %d chains to propagate to "
              "(chains seen: %s)" % (cys_resno, ",".join(sorted(sg_by_chain))))

    # assign a free residue number per target copy (reuse lig resno if free)
    targets = []
    for c in tgts:
        rn = lig["resno"]
        while _resno_used(src_pdb, c, rn):
            rn += 1
        targets.append({"tgt_chain": c, "new_chain": c, "new_resno": rn})

    copies_pdb = os.path.join(wd, "copies.pdb")
    report = {"copies": []}
    if targets:
        job = {"_wd": wd, "model": src_pdb, "ref_chain": ref_chain,
               "lig_chain": lig["chain"], "lig_resno": lig["resno"],
               "lig_resname": lig["resname"], "cys_resno": cys_resno,
               "window": window, "targets": targets, "out_pdb": copies_pdb}
        report = _run_ncs_worker(job)
        if "error" in report:
            raise RuntimeError("NCS worker: %s" % report["error"])

    # augmented model: src + propagated ligand copies + LINK records for ALL
    # copies (reference + propagated) so the model is self-consistent.
    link_id = FAMILIES[family]["link_id"]
    src_lines = open(src_pdb).read().splitlines()
    add_lines = _hetatm_lines(copies_pdb) if (targets and os.path.exists(copies_pdb)) else []

    link_lines = []
    made = []
    # reference copy
    made.append((ref_chain, lig["chain"], lig["resno"]))
    for cp in report.get("copies", []):
        if "new_chain" in cp:
            made.append((cp["tgt_chain"], cp["new_chain"], cp["new_resno"]))
    for cys_ch, lc, lr in made:
        sg = {"chain": cys_ch, "resno": cys_resno, "resname": "CYS", "atom": cys_atom}
        wh = {"chain": lc, "resno": lr, "resname": lig["resname"], "atom": warhead_atom}
        link_lines.append(_link_record(sg, wh, link_id))

    aug_pdb = os.path.join(wd, "model_ncs_augmented.pdb")
    out = []
    inserted = False
    for l in src_lines:
        if not inserted and l.startswith(("ATOM", "HETATM")):
            out.extend(link_lines); inserted = True
        # drop END / stale MASTER, and any pre-existing LINK/LINKR (we re-emit
        # links for every copy above, so keeping the old ones would duplicate
        # the reference link); TER/CONECT/anisou are kept.
        if l.startswith(("END", "MASTER", "LINK")):
            continue
        out.append(l)
    if not inserted:
        out = link_lines + out
    out.extend(add_lines)
    out.append("END")
    open(aug_pdb, "w").write("\n".join(out) + "\n")

    print("[cootvalent] NCS propagation: reference Cys chain %s, ligand %s/%d"
          % (ref_chain, lig["chain"], lig["resno"]))
    for cp in report.get("copies", []):
        if "new_chain" in cp:
            print("[cootvalent]   -> chain %s: ligand %s/%d  (CA superpose rmsd %.3f A, %d CA)"
                  % (cp["tgt_chain"], cp["new_chain"], cp["new_resno"],
                     cp.get("rmsd", -1), cp.get("n_ca", 0)))
        else:
            print("[cootvalent]   -> chain %s: SKIPPED (%s)"
                  % (cp.get("tgt_chain", "?"), cp.get("error", "?")))
    print("[cootvalent]   augmented PDB: %s   (%d covalent copies)"
          % (aug_pdb, len(made)))

    result = {"augmented_pdb": aug_pdb, "link_id": link_id, "family": family,
              "warhead_atom": warhead_atom, "ref_chain": ref_chain,
              "copies": report.get("copies", []), "n_copies": len(made),
              "refined_pdb": None, "refined_mtz": None}

    try:
        coot.read_pdb(aug_pdb)
    except Exception as e:
        print("[cootvalent]   (load-back skipped:", e, ")")

    if do_refine and mtz:
        # one link CIF + one ligand monomer dict cover every copy (same comp id).
        # Build the link geometry the same way declare_covalent_link does.
        nbrs = _neighbours(lig_atoms, warhead_atom)
        carbons = [n for n in nbrs if n[0][0] in "Cc"]
        ca = carbons[0][0] if carbons else None
        cg = carbons[1][0] if len(carbons) >= 2 else None
        dict_path = lig_dict or _monomer_lib_path(lig["resname"])
        pre_form, cb_h_del, ca_h_add = _analyse_ligand_dict(
            dict_path, warhead_atom, ca, family)
        link_cif = os.path.join(wd, "link.cif")
        open(link_cif, "w").write(
            _build_link_cif(family, lig["resname"], warhead_atom, ca=ca, cg=cg,
                            pre_form=pre_form, cb_h_to_delete=cb_h_del,
                            ca_h_to_add=ca_h_add))
        sg0 = {"chain": ref_chain, "resno": cys_resno, "resname": "CYS", "atom": cys_atom}
        wh0 = {"chain": lig["chain"], "resno": lig["resno"],
               "resname": lig["resname"], "atom": warhead_atom}
        if use_wrapper:
            rp, rm, _ = refmac_refine_via_wrapper(
                aug_pdb, mtz, link_cif, lig_dict, sg0, wh0, wd,
                ncyc=ncyc, add_waters=add_waters)
        else:
            rp, rm, _ = refmac_refine(
                aug_pdb, mtz, link_cif, lig_dict, sg0, wh0, wd, ncyc=ncyc)
        result["refined_pdb"] = rp
        result["refined_mtz"] = rm
    return result


# ---------------------------------------------------------------------------
# GUI wiring: "Cootvalent" menu -> "Declare covalent link..."
# ---------------------------------------------------------------------------
def _install_menu():
    mbm = _gui_fn("coot_menubar_menu")
    asm = _gui_fn("add_simple_coot_menu_menuitem")
    gde = _gui_fn("generic_double_entry")
    if not (mbm and asm and gde):
        print("[cootvalent] GUI menu API not found (no-graphics?); "
              "declare_covalent_link()/propagate_covalent_ncs() still callable "
              "directly from the console.")
        return

    def _go(sg_cid, warhead_cid):
        # family is inferred from a trailing token, else defaults to F1.
        fam = "F1"
        wc = warhead_cid
        m = re.search(r"\s+(F1|F2|CAA)\s*$", warhead_cid or "", re.I)
        if m:
            fam = m.group(1).upper()
            wc = warhead_cid[:m.start()].strip()
        try:
            imol = coot.first_coords_imol()
        except Exception:
            imol = 0
        try:
            declare_covalent_link(imol, sg_cid.strip(), wc.strip(),
                                  family=fam, mtz=None, do_refine=True)
            coot.info_dialog("Covalent link declared (family %s).\n"
                             "See the scripting console for the augmented PDB "
                             "and link CIF paths + the refmac command." % fam)
        except Exception as e:
            coot.info_dialog("Declare covalent link failed:\n\n%s" % e)

    def _activate(*args):
        gde(
            "Cys SG CID   (e.g. //A/481(CYS)/SG)",
            "warhead C CID  (append F1/F2/CAA to choose family)",
            "//A/481(CYS)/SG", "//A/701(LIG)/CAA F1",
            False, False,
            "Declare", _go)

    # ---- NCS propagation menu item ----
    def _go_prop(lig_cid, spec):
        # spec = "<cys_resno> [F1|F2|CAA]", e.g. "547 F2"
        toks = (spec or "").split()
        if not toks:
            coot.info_dialog("Enter the Cys residue number (e.g. 547 F2)."); return
        try:
            cys_resno = int(re.sub(r"[^\-\d]", "", toks[0]))
        except ValueError:
            coot.info_dialog("Could not read Cys residue number from '%s'." % spec); return
        fam = "F2"
        for t in toks[1:]:
            if t.upper() in FAMILIES:
                fam = t.upper()
        try:
            imol = coot.first_coords_imol()
        except Exception:
            imol = 0
        try:
            r = propagate_covalent_ncs(imol, lig_cid.strip(), cys_resno,
                                       family=fam, do_refine=False)
            coot.info_dialog("Propagated to %d covalent copy(ies) at Cys %d "
                             "(family %s).\nAugmented model loaded; check each "
                             "copy in density, then refine.\nSee console for "
                             "per-copy superposition RMSDs.\n\n%s"
                             % (r["n_copies"], cys_resno, fam, r["augmented_pdb"]))
        except Exception as e:
            coot.info_dialog("NCS propagation failed:\n\n%s" % e)

    def _activate_prop(*args):
        gde(
            "Reference ligand CID   (e.g. //A/601(LIG))",
            "Cys resno [+ family]   (e.g. 547 F2)",
            "//A/601(LIG)", "547 F2",
            False, False,
            "Propagate", _go_prop)

    # Create the menu LOUDLY -- a swallowed exception here is why a menu can
    # silently fail to appear; print the real traceback so it's diagnosable.
    try:
        menu = mbm("Cootvalent")
        if menu is None:
            # coot_menubar_menu() returns None when coot_python.main_menubar()
            # is unavailable. If coot_python is simply not ready yet (classic
            # Coot 0.9 during startup) a restart fixes it; if the build has no
            # coot_python at all (some bandicoot builds) no Python menu is
            # possible and the console API is the interface.
            try:
                import coot_python  # noqa: F401
                hint = ("menubar not ready -- RESTART bcoot; the menu installs "
                        "at startup.")
            except Exception:
                hint = ("this Coot build has no coot_python, so a Python menu "
                        "cannot be added. Use the console API: "
                        "declare_covalent_link(...) / propagate_covalent_ncs(...).")
            print("[cootvalent] menu not installed -- " + hint)
            return
        asm(menu, "Declare covalent link...", _activate)
        asm(menu, "Propagate covalent ligand to NCS copies...", _activate_prop)
        print("[cootvalent] OK: 'Cootvalent' menu ready "
              "(Declare covalent link... + Propagate covalent ligand to NCS copies...).")
    except Exception:
        import traceback
        print("[cootvalent] MENU INSTALL FAILED -- the Cootvalent menu will not "
              "appear. Traceback:\n" + traceback.format_exc())


_install_menu()
