# Bandicoot overnight validation — morning report (2026-07-16)

Three parallel workstreams ran overnight against your real `~/xtal` data. **All three came back positive.** Details below; per-workstream logs in `RESULTS.md`; deliverables in `~/PyKeko/bandicoot/`.

## Bottom line

The three PyKeko-critical workflows — **MCP control, ligand fitting, and the covalent workflow** — all **work natively on bandicoot**, validated on real crystallographic data. Combined with the earlier assessment (~80% of PyKeko is already native Coot 0.9; you've discounted the GPL-3 and web-feature concerns), **bandicoot is a technically viable platform for your work.** The remaining reservation is not technical — it's the **bus-factor-1 / frozen-upstream maintenance risk**, which is a risk to accept with eyes open, not a blocker.

These are **feasibility prototypes**, not production. Two real covalent gotchas surfaced (and were solved) that show the port isn't free — it needs the same careful chemistry handling PyKeko already does.

## 1. MCP control — `bandicoot_backend.py` ✅

A clean, reusable `BandicootBackend` class: spawns headless `coot --no-graphics`, drives it over the socket with `eval(code)->value` (the direct analog of `moorhen_eval`/`run_js`) + verb methods + lifecycle. Fixed two bugs beyond the prototype (double-reply framing, latency → ~182 ms). Live transcript proved headless control: load 3PTB → atom count → set/read rotation centre → **mutate ASN→ALA with live persistence** → heartbeat. Includes a PyKekoMCP integration sketch (`coot_*` tools alongside `moorhen_*` = one control plane over both engines).
**Verdict: the hybrid "drive native Coot from PyKekoMCP" vision is real and cheap.**

## 2. Non-covalent ligand fitting ✅ (honest, data-limited)

Native blind `execute_ligand_search_py` **runs headless** (no hang) — the gold-standard test. All dicts from the real SMILES→acedrg pipeline.

| Dataset | Ligand | RMSD (Å) | Result |
|---|---|---|---|
| 1NHZ | mifepristone (32 at) | **0.41** | PASS — textbook |
| 2P54 | thiazole drug (33 at) | **0.084** | PASS — near-perfect (σ=0.6) |
| 3PTB | benzamidine (symmetric) | 2.65 | FAIL\* — right pocket (0.4 Å centroid), amidine flipped |
| 5L0E | 6ZN (3.06 Å data) | 7.3 | FAIL — weak/fragmented density |

\*The two "fails" are crystallographic limits **every** fitting tool shares (symmetric-ligand orientation; genuinely weak 3.06 Å density), not bandicoot bugs. **Verdict: native fitting is as good as the data allows — same Coot engine PyKeko wraps, and bandicoot exposes the fuller blind cluster-search.** Needs a cluster-σ UI knob and "return all clusters for picking" (top cluster isn't always the right pose).

## 3. Covalent workflow — `covalent-link.py` ✅ (all 4 PASS)

Extension that declares a Cys–warhead link (refmac-ready link CIF + LINK/LINKR record → external refmac5). Validated across all four families:

| Dataset | Family | Deposited SG–C (Å) | Refined (Å) | Canonical | Result |
|---|---|---|---|---|---|
| 5P9J ibrutinib | F1 | 1.850 | **1.852** | 1.81 | PASS |
| 6JX0 osimertinib | F1 | 1.795 | **1.876** | 1.81 | PASS |
| 8FD9 acalabrutinib | F2 | 1.683 | **1.718** | 1.78 | PASS |
| 4YHF cyanoacetamide | edge | 1.920/1.907 | **1.844/1.809** | 1.81 | PASS |

Every covalent bond holds at canonical distance (~0.1 Å), no drift, Rwork 0.17–0.22. **Functionally equivalent to PyKeko's JS mmCIF-surgery** — both hand refmac the same LINK + link CIF. External refmac5 is confirmed the right handoff (Emsley, coot#374: Coot still ignores `_struct_conn`). **Verdict: the crown-jewel workflow ports cleanly and correctly.**

**Two gotchas found + solved (remember these):**
1. **Pre-form library ligands** (YY3): refmac prefers the lib's *unreacted alkene* over your post-form dict → bond lands long. Fix: rename the ligand to a non-library id + supply saturated dict + `MAKE NEWLIGAND NOEXIT` (the native analogue of PyKeko's mod2 chem_comp rewrite — which sidesteps this trap in the WASM layer).
2. **XQQ dict**: `acedrg -c` on raw coords → 0 atoms. Fix: RDKit `AssignBondOrdersFromTemplate` (bound-form SMILES) → then `acedrg -c`.

The **4YHF edge case** is instructive: α-cyanoacetamide is genuinely *not* in the F1–F6 registry, so PyKeko's detector should (correctly) decline it — the extension mirrors that. refmac handles it fine when driven manually.

## What "committing to bandicoot" would still require (honest scope)

Feasibility is proven; a real move is not free. Remaining work: the PyMOL selection/color layer (deferred), packaging the extensions into a cohesive suite + keybindings, hardening the covalent chemistry for all families (the pre-form/dict gotchas), a cluster-σ + cluster-picking ligand UI, and — the standing reservation — accepting sole maintenance of a single-maintainer, frozen-upstream, GTK2 stack.

## Deliverables produced overnight
- `~/PyKeko/bandicoot/bandicoot_backend.py` — MCP adapter (+ integration sketch)
- `~/PyKeko/bandicoot/covalent-link.py` — covalent-link extension (817 lines)
- `~/PyKeko/bandicoot/ligand-from-smiles.py` — ligand extension (from before; consider adding blind-find)
- `~/PyKeko/bandicoot/overnight/RESULTS.md` — full per-workstream logs

## Suggested next decision
The evaluation has now answered every feasibility question with a "yes." The call in front of you is no longer *"can we?"* but *"do we want to own this stack?"* — a judgment about maintenance risk and where you want to spend effort, not a technical unknown. Good place to step back and decide direction before building more.
