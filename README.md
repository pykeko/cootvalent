# Cootvalent

**Covalent-ligand and SMILES-to-fitted-ligand workflows for native Coot 0.9**
(tested on [bandicoot](https://github.com/fraser-lab/bandicoot), the macOS-native Coot
0.9 fork — should also work on any classic Coot 0.9.x with a Python console).

Pure-Python Coot extensions. No fork, no C++, no rebuild. Drop them in
`~/.coot-preferences/` and a **Cootvalent** menu appears in Coot.

The covalent workflow is the point: declare a Cys–warhead bond, generate a
refmac-ready link CIF + LINK/LINKR records, and refine to a canonical S–C
distance — the thing classic Coot has no built-in path for. The SMILES→ligand
and blind-search items are the supporting cast (you need a built ligand before
you can bond it).

Derived from the covalent-ligand mechanism developed in
[PyKeko](https://github.com/pykeko) (a branded Moorhen/Coot-WASM fork), ported
here off the web stack to run in the native desktop app.

## Menu items (Cootvalent →)

| Item | What it does |
|---|---|
| **Ligand from SMILES…** | SMILES → acedrg dict → place at pointer + jiggle-fit into density |
| **Find in density (blind search)…** | native Fo−Fc cluster search over the whole map, σ knob, returns all clusters to pick |
| **Declare covalent link…** | Cys-SG ↔ warhead-carbon: emits refmac-ready link CIF + LINK/LINKR + `_struct_conn`, spawns external refmac5 |
| **Auto-detect + declare covalent link** | classifies the warhead family (F1–F6) and finds the two link atoms automatically, then declares |
| **Propagate covalent ligand to NCS copies…** | takes one built+placed covalent ligand and copies it to the same Cys in every other protein chain of the ASU, by superposing the reference chain onto each target chain (gemmi via `ccp4-python`) and applying that transform; writes a LINK record per copy and loads the augmented model for you to check each in density |

### Refinement: inline refmac5 or the `refmac.sh` wrapper

The covalent refine step runs either the inline `refmac5` call (default) or an
external [`refmac.sh`](https://github.com/pykeko) wrapper. Pass
`use_wrapper=True` (and optionally `add_waters=True`) to
`declare_covalent_link()` / `propagate_covalent_ncs()` to route through the
wrapper, which adds automatic water picking (`findwaters`, via `-W`) after the
covalent refinement — so the link geometry and ordered solvent come out of one
command. The wrapper is located via `$COOTVALENT_REFMAC_SH`, then
`~/bin/refmac.sh`, then `PATH`.

### Console interface (bandicoot: where menus/keys can't be trusted)

Some bandicoot builds ship **without `coot_python`** (and with a stubbed
`coot_toolbar_button` and invisible-stub GTK dialogs), so a menu/toolbar/input
dialog can't be created from Python. They also have **compiled-in native key
accelerators** that intercept keys before Python and are sometimes destructive
(observed: `L` = go-to-ligand, `Ctrl+R` = rock, `Ctrl+D` = delete). So
`cootvalent-keys.py` **does not bind any keys automatically**; the reliable,
collision-free interface is three console functions:

| Call | Action |
|---|---|
| `cv_declare()` | auto-detect the placed warhead (ligand atom 1.2–2.6 Å from a CYS SG) and declare the covalent link |
| `cv_propagate()` | declare + propagate the ligand to the same Cys in every NCS copy |
| `cv_full()` | propagate **and** refine via `refmac.sh` (`-W` waters) |

These auto-detect, so the ligand must already be built/placed at the Cys and
**merged into the protein molecule** (`ligand_from_smiles()` builds it;
`merge_molecules([lig], protein_imol)` merges it). `cv_full()` guesses the MTZ +
ligand dict from the model's directory; override with
`cootvalent_set_refine(mtz="…", lig_dict="…")`.

If you have confirmed a key is free on your build, opt in explicitly:
`cootvalent_bind_keys(propagate="Control_o", full="Control_u")` (it refuses keys
known to be grabbed natively). On classic Coot (with `coot_python`) the
"Cootvalent" menu is added as well.

### MCP bridge — drive a live bcoot from an agent

`cootvalent-mcp-bridge.py` (in-Coot) + `cootvalent_mcp_server.py` (MCP server)
let an MCP client (e.g. Claude Code) drive your **live GUI bcoot** directly,
instead of pasting console lines. It uses a small JSON **file queue** polled on
Coot's GLib main loop (`coot.bandicoot_python_timeout_add`) — non-blocking and
GUI-safe — rather than a socket that would need a pump loop and freeze the GUI.

Setup:
1. In the bcoot Python console (once per session): `cootvalent_mcp_start()`
2. Register the server with your MCP client. Example (`mcp.json.example`):
   ```
   claude mcp add cootvalent /opt/anaconda3/bin/python \
       /Users/you/bin/cootvalent_mcp_server.py \
       -e COOTVALENT_MCP_DIR=/Users/you/.cootvalent_mcp
   ```
   (the server needs an interpreter with the `mcp` SDK; the queue dir must match
   the bridge's `$COOTVALENT_MCP_DIR`, default `~/.cootvalent_mcp`).

Tools: `coot_eval`, `cv_build_at_cys`, `cv_warhead_dist`, `cv_merge_ligand`,
`cv_declare`, `cv_propagate`, `cv_full`, `cv_clear_ligands`, `load_pdb`,
`read_dict`, `close_all`. Note `coot_eval` runs arbitrary Python in your session
— only register this for local, single-user use.

`bandicoot_backend.py` is a separate, older **headless** driver
(`BandicootBackend`) that spawns its own no-graphics Coot over a socket; it is
*not* installed and is kept for anyone wiring headless Coot into an MCP plane.

## Warhead families (Cys-S covalent)

| Family | Chemistry | Detect | Declare / refine | Reacted build |
|---|---|---|---|---|
| F1 | acrylamide → saturated β-thioether | ✅ | ✅ (validated) | ✅ C=C→C–C |
| F2 | α,β-ynamide/butynamide → vinyl thioether | ✅ | ✅ (validated) | ✅ C≡C→C=C |
| F3 | α-halo-acetamide SN2 thioether (incl. **α-chlorofluoroacetamide, CFA**) | ✅ | ✅ (validation pending real structure) | ✅ dehalogenate (F retained) |
| F4 | epoxide → β-hydroxy thioether | ✅ | ❌ | — |
| F5 | maleimide → 3-thiosuccinimide | ✅ | ❌ | — |
| F6 | reversible ketone/aldehyde hemithioketal | ✅ | ❌ | — |

**Reacted build**: `cv_build_at_cys(smiles, chain, resno, family=...)` transforms
the unreacted SMILES to the covalent-product form for the family *before* acedrg
(F2 alkyne→vinyl, F1 acrylamide→saturated, F3 α-halo-acetamide→dehalogenated),
so you fit the correct sp²/sp³ geometry with the right bond order — the S is
added later by the link, which deletes the warhead's placeholder H. Pass
`family=None` to build the SMILES verbatim.

## Install

```bash
git clone https://github.com/pykeko/cootvalent.git
cd cootvalent
./install.sh          # copies the 3 extensions into ~/.coot-preferences/
```

Then restart Coot/bcoot. Uninstall by deleting the three files from
`~/.coot-preferences/` (install.sh prints the exact command).

## Requirements

- **Coot 0.9.x** with a Python console (bandicoot recommended on macOS Tahoe).
- **CCP4** for `acedrg` (SMILES→dict) and `refmac5` (covalent refinement). The
  extensions source a CCP4 setup script; the default path is
  `/Applications/ccp4-9/bin/ccp4.setup-sh` — edit `CCP4_SETUP` at the top of
  each file if yours differs.
- The auto-detector uses **RDKit via CCP4's `ccp4-python`** (bundled with CCP4-9);
  Coot's own Python does not need RDKit.

## Validated

On real crystallographic data (2026-07-16):

- **Ligand blind-find** — 1NHZ 0.29 Å, 2P54 0.41 Å RMSD to deposited pose.
- **Covalent declare + refmac** — 5P9J/6JX0 (F1), 8FD9 (F2), 4YHF (edge) all
  refine to canonical S–C within ~0.1 Å.
- **Auto-detect** — F1/F2 correct on 5P9J/6JX0/8FD9; 4YHF (α-cyanoacetamide,
  outside the F1–F6 taxonomy) correctly declines.

Full logs in [`overnight/`](overnight/).

## Caveats

- **F4/F5/F6 detector paths** are implemented from the taxonomy but not yet
  exercised on real data (no epoxide/maleimide/ketone test set to hand).
- **NCS propagation is a geometric starting point, not a verdict.** It places a
  copy at every equivalent Cys assuming the pocket is NCS-conserved; occupancy
  and even presence vary per copy (partial reaction is common), so check each
  propagated ligand against its own density before refining. It uses `gemmi`
  from `ccp4-python` (superposes on Cα within ±15 residues of the Cys, falling
  back to all matched Cα).
- **Pre-form library ligands** (a monomer dict that still carries the *unreacted*
  warhead, e.g. an alkene where the bound form is saturated) can make refmac
  prefer the wrong geometry; `covalent-link.py` renames the ligand to a
  non-library id and supplies the saturated dict to sidestep this.
- Paths (`CCP4_SETUP`, monomer-library roots, the bandicoot install dir) are
  macOS/CCP4-9 defaults — adjust the constants at the top of each file for a
  different layout.

## License

GPL-3.0-or-later — these extensions call Coot's GPL-3 Python API and are a
derivative work of Coot. See [LICENSE](LICENSE).
