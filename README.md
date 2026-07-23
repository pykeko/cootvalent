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

`bandicoot_backend.py` is separate — an **external MCP driver** (`BandicootBackend`)
that drives headless Coot over a socket (`eval(code)→value`). It is *not* a Coot
extension and `install.sh` does not install it; it lives here for anyone wiring
native Coot into an MCP control plane.

## Warhead families (Cys-S covalent)

| Family | Chemistry |
|---|---|
| F1 | acrylamide → saturated β-thioether |
| F2 | α,β-ynamide → vinyl thioether |
| F3 | α-chloro/cyano-acetamide SN2 thioether |
| F4 | epoxide → β-hydroxy thioether |
| F5 | maleimide → 3-thiosuccinimide |
| F6 | reversible ketone/aldehyde hemithioketal |

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
