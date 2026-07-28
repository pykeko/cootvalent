#!/usr/bin/env python3
"""Cootvalent MCP server -- drive a live bcoot session over a file queue.

Pairs with cootvalent-mcp-bridge.py running inside bcoot (arm it once with
cootvalent_mcp_start()). This server (run under an interpreter that has the mcp
SDK, e.g. /opt/anaconda3/bin/python) writes Python snippets to the queue and
reads back the result, exposing them as MCP tools so an agent can drive Coot
directly instead of dictating console lines.

Queue dir: $COOTVALENT_MCP_DIR (default ~/.cootvalent_mcp), same as the bridge.

Tools:
  coot_eval(code)               run arbitrary Python in the Coot namespace
  cv_build_at_cys(smiles,chain,resno)  build ligand at a Cys (separate molecule)
  cv_warhead_dist()             report closest ligand-carbon <-> Cys-SG
  cv_merge_ligand()             merge the built ligand into the protein
  cv_declare()                  auto-detect + declare the covalent link
  cv_propagate()                declare + propagate to all NCS copies
  cv_full()                     propagate + refine (refmac.sh -W)
  cv_clear_ligands(comp)        remove ligand copies
  load_pdb(path) / read_dict(path) / close_all()
"""

import os
import json
import time

from mcp.server.fastmcp import FastMCP

MCP_DIR = os.environ.get("COOTVALENT_MCP_DIR",
                         os.path.expanduser("~/.cootvalent_mcp"))
_REQ = os.path.join(MCP_DIR, "request.json")
_RESP = os.path.join(MCP_DIR, "response.json")
_TIMEOUT = float(os.environ.get("COOTVALENT_MCP_TIMEOUT", "600"))  # s (refine is slow)

mcp = FastMCP("cootvalent")


def _next_id():
    """Monotonic id from the last response (or 0), +1."""
    last = 0
    try:
        with open(_RESP) as fh:
            last = int(json.load(fh).get("id", 0))
    except Exception:
        pass
    # also consider a pending request id
    try:
        with open(_REQ) as fh:
            last = max(last, int(json.load(fh).get("id", 0)))
    except Exception:
        pass
    return last + 1


_NOT_ARMED = (
    "The cootvalent bridge is not responding, so no live bcoot is connected.\n"
    "In your bcoot Python console, run:\n"
    "    cootvalent_mcp_start()\n"
    "(first  exec(open('~/.coot-preferences/cootvalent-mcp-bridge.py').read())  "
    "if it's undefined), then retry. Queue dir: %s" % MCP_DIR
)


def _roundtrip(code, timeout):
    """Send one request, wait for its response, or None if no reply in time."""
    os.makedirs(MCP_DIR, exist_ok=True)
    rid = _next_id()
    tmp = _REQ + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"id": rid, "code": code}, fh)
    os.replace(tmp, _REQ)
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with open(_RESP) as fh:
                resp = json.load(fh)
            if int(resp.get("id", -1)) == rid:
                return resp
        except Exception:
            pass
        time.sleep(0.2)
    return None


def _armed(probe_timeout=2.5):
    """Fast heartbeat: is a live bcoot polling the queue right now?"""
    r = _roundtrip("1", probe_timeout)   # cheap expression the bridge evals
    return bool(r and r.get("ok"))


def _send(code, timeout=None):
    """Write a request, wait for the matching response, return a text summary.
    Fails fast with an actionable message if the bridge isn't armed, so tools
    don't hang for the full timeout against a dead/absent session."""
    if not _armed():
        return _NOT_ARMED
    resp = _roundtrip(code, timeout or _TIMEOUT)
    if resp is None:
        return ("No response after %ss. The command may still be running in "
                "bcoot, or the session stopped. Check bcoot." % (timeout or _TIMEOUT))
    parts = []
    if resp.get("stdout"):
        parts.append(resp["stdout"].rstrip())
    if resp.get("result"):
        parts.append("=> " + resp["result"])
    if not resp.get("ok"):
        parts.append("ERROR:\n" + resp.get("error", ""))
    return "\n".join(parts) if parts else "(no output)"


@mcp.tool()
def coot_eval(code: str) -> str:
    """Run arbitrary Python in the live bcoot namespace and return stdout/value."""
    return _send(code)


@mcp.tool()
def cv_build_at_cys(smiles: str, chain: str, resno: int) -> str:
    """Build a ligand from SMILES near a Cys (as a separate molecule to fit)."""
    return _send("cv_build_at_cys(%r, %r, %d)" % (smiles, chain, int(resno)))


@mcp.tool()
def cv_warhead_dist() -> str:
    """Report closest ligand-carbon <-> Cys-SG distances across all molecules."""
    return _send("cv_warhead_dist()")


@mcp.tool()
def cv_merge_ligand() -> str:
    """Merge the built ligand molecule into the protein molecule."""
    return _send("cv_merge_ligand()")


@mcp.tool()
def cv_declare() -> str:
    """Auto-detect the placed warhead and declare the covalent link."""
    return _send("cv_declare()")


@mcp.tool()
def cv_propagate() -> str:
    """Declare + propagate the covalent ligand to every NCS copy."""
    return _send("cv_propagate()")


@mcp.tool()
def cv_full() -> str:
    """Propagate + refine via refmac.sh (-W waters). May take minutes."""
    return _send("cv_full()")


@mcp.tool()
def cv_clear_ligands(comp: str = "LIG") -> str:
    """Remove all copies of a ligand comp id (reset after messy attempts)."""
    return _send("cv_clear_ligands(%r)" % comp)


@mcp.tool()
def load_pdb(path: str) -> str:
    """Read a coordinate file into bcoot; returns the molecule number."""
    return _send("coot.read_pdb(%r)" % path)


@mcp.tool()
def read_dict(path: str) -> str:
    """Read a restraint/monomer CIF dictionary into bcoot."""
    return _send("coot.read_cif_dictionary(%r)" % path)


@mcp.tool()
def close_all() -> str:
    """Close every loaded molecule (models and maps)."""
    return _send("[coot.close_molecule(i) for i in range(coot.graphics_n_molecules())]")


if __name__ == "__main__":
    mcp.run()
