# Bandicoot overnight prototype+validation — 2026-07-16

Three parallel workstreams (all headless, all read-only on ~/xtal):

1. NON-COVALENT ligand fit-back validation (lig_test): 3PTB/BEN, 1NHZ/486, 2P54/735, 5L0E/?
   For each: load apo_refined (ligand removed -> +density) + map, generate ligand from
   acedrg, find/fit into the positive density, RMSD vs ground-truth pose in the full PDB.
2. COVALENT workstream: build a "Declare covalent link" Coot extension AND validate
   declare->refmac->geometry on cov_test: 5P9J/8E8(F1), 6JX0/YY3(F1), 8FD9/XQQ(F2), 4YHF/4C9(edge).
3. MCP BandicootBackend adapter: turn the proven socket driver into a real PyKekoMCP backend + test.

Environment fixed earlier: CCP4 acedrg works (orphaned ~/Library/Python/3.9 disabled).
Ligand extension exists: ~/PyKeko/bandicoot/ligand-from-smiles.py
Socket driver: scratchpad server_side.py / coot_side.py
