# Overnight results — 2026-07-16

## Workstream 3: MCP BandicootBackend adapter — ✅ DONE, clean
Deliverable: ~/PyKeko/bandicoot/bandicoot_backend.py (pure ASCII, self-contained).
- eval(code)->value primitive (analog of moorhen_eval/run_js) + verb methods + lifecycle.
- Fixed 2 real bugs beyond prototype: double-reply framing (parse FIRST marker), latency
  (stop-on-2-markers -> ~182 ms/round-trip).
- Live transcript (headless, no XQuartz): load 3PTB -> get_atom_count 1701 (write-back) ->
  set/read rotation centre -> go_to_residue A/25 -> mutate ASN->ALA, persistence PROVEN ->
  heartbeat OK. ~/xtal untouched, zero lingering procs.
- PyKekoMCP integration sketch included (coot_* tools alongside moorhen_*, one control plane).
VERDICT: headless MCP control of native Coot works; integration is a thin adapter. GO.

## Workstream 1: non-covalent ligand fit-back — running
## Workstream 2: covalent extension + validation — running

## Workstream 1: non-covalent ligand fit-back — ✅ DONE (honest, mixed-by-design)
Blind find (execute_ligand_search_py) RUNS HEADLESS (no hang) — gold-standard test. All dicts from real SMILES->acedrg (zero fallbacks).
| Dataset | Ligand | Mode | RMSD (A) | Result |
|---|---|---|---|---|
| 1NHZ | 486 mifepristone (32at) | blind-find | 0.41 | PASS (textbook) |
| 2P54 | 735 thiazole (33at) | blind-find (sigma 0.6) | 0.084 | PASS (near-perfect) |
| 3PTB | BEN benzamidine (symmetric) | blind-find | 2.65 | FAIL* right pocket (0.4A centroid), amidine flipped — symmetric-ligand limit shared by all tools; correct cluster was rank #2 |
| 5L0E | 6ZN (3.06A data) | blind + assisted | 7.3 | FAIL — weak/fragmented density at 3.06A, not a bandicoot bug |
VERDICT: native blind ligand search fits real ligands into density, headless, as well as the data allows. Same Coot engine PyKeko wraps; bandicoot exposes the FULLER blind cluster-search (stronger than PyKeko ext's jiggle-only path). Ligand-fitting grounds: GO.
Gated by (expected, not bugs): data resolution/density strength; cluster-sigma tuning knob (needs a UI control); top cluster != always correct pose (return all clusters for picking).

### Follow-up for morning (NOT done overnight — needs a UX call):
- ligand-from-smiles.py currently uses fit_molecule_to_map_by_random_jiggle (rigid, place-at-pointer). Consider ALSO offering execute_ligand_search_py (blind cluster search) + a cluster-sigma control + return-all-clusters-for-picking. Two different UX: "put my ligand at this blob" (jiggle) vs "find where my ligand goes" (blind). Held for user's UX decision.

## Workstream 2: covalent extension + validation — ✅ DONE (all 4 PASS)
Deliverable: ~/PyKeko/bandicoot/covalent-link.py (817 lines, pure ASCII). Menu PyKeko->"Declare covalent link...".
Emits refmac-ready link CIF (data_link_list + data_mod_list + link_id cols) + LINK/LINKR record -> external refmac5.
| dataset | family | Cys SG | warhead | deposited SG-C | refined SG-C | canonical | PASS |
|---|---|---|---|---|---|---|---|
| 5P9J ibrutinib | F1 | A/481 | 8E8 CAA | 1.850 | 1.852 | 1.81 | PASS |
| 6JX0 osimertinib | F1 | A/797 | YY3 C9 | 1.795 | 1.876 | 1.81 | PASS (pre-form workaround) |
| 8FD9 acalabrutinib | F2 | A/481 | XQQ C19 | 1.683 | 1.718 | 1.78 | PASS (acedrg post-form dict) |
| 4YHF cyanoacetamide | CAA | A/481+B | 4C9 C1 | 1.920/1.907 | 1.844/1.809 | 1.81 | PASS (manual; detector correctly declines) |
All bonds hold canonical within ~0.1A, no drift, Rwork 0.17-0.22.
VERDICT: covalent workflow ports cleanly + correctly to native Coot/refmac; functionally equivalent to PyKeko JS mmCIF-surgery. External refmac5 is the right handoff (Emsley confirmed coot#374: make_link_restraints_ng still ignores _struct_conn). GO.
Two gotchas found+solved (worth remembering): (1) pre-form library ligands (YY3) - refmac prefers lib alkene form over passed post-form dict; fix = rename ligand to non-library id + saturated dict + MAKE NEWLIGAND NOEXIT (native analogue of PyKeko's mod2 chem_comp rewrite). (2) XQQ dict - acedrg -c on raw coords gives 0 atoms; fix = RDKit AssignBondOrdersFromTemplate then acedrg -c.
4YHF nuance: PyKeko detector SHOULD decline cyanoacetamide (not in F1-F6; CYS-CAA is chloroacetamide F3) - correct behavior; refmac handles it if driven manually.

## BUILD (2026-07-16): ligand blind-find upgrade — DONE
ligand-from-smiles.py now has BOTH: ligand_from_smiles() (place-at-pointer+jiggle) AND
ligand_blind_find(smiles,name,imol_prot,imol_diff_map,sigma) (native execute_ligand_search_py,
returns ALL clusters ranked, no auto-accept, auto-detects Fo-Fc map, sigma knob exposed).
New menu: PyKeko -> "Find in density (blind search)...". Validated: 1NHZ 0.29A, 2P54 0.41A (rank 0,
multiple clusters). Pure ASCII, ast-parse clean.
