# cootvalent-mcp-bridge.py  --  in-Coot side of the Cootvalent MCP bridge.
#
# Lets an external MCP server drive THIS live (GUI) bcoot session without the
# copy-paste dance. Instead of a socket (which needs a pump loop that would
# freeze the GUI), it uses a tiny file queue polled on the GLib main loop via
# coot.bandicoot_python_timeout_add -- non-blocking and GUI-safe (the same timer
# bandicoot uses for live model/map updates).
#
# Protocol (JSON files in ~/.cootvalent_mcp/, override with COOTVALENT_MCP_DIR):
#   request.json  = {"id": <int>, "code": "<python source>"}
#   response.json = {"id": <int>, "ok": <bool>, "result": <str>, "stdout": <str>,
#                    "error": <str>}
# The bridge runs each new id exactly once, evaluating an expression (returns its
# repr) or executing statements, in Coot's own namespace (coot + all cv_* funcs).
#
# USAGE (once per GUI session, in the bcoot Python console):
#     cootvalent_mcp_start()
# Stop with cootvalent_mcp_stop().  ASCII-only.

import os, io, json, time, traceback, contextlib

MCP_DIR = os.environ.get("COOTVALENT_MCP_DIR",
                         os.path.expanduser("~/.cootvalent_mcp"))
_REQ = "request.json"
_RESP = "response.json"

_mcp_state = {"running": False, "last_id": 0, "interval_ms": 300}


def _mcp_dir():
    try:
        os.makedirs(MCP_DIR, exist_ok=True)
    except Exception:
        pass
    return MCP_DIR


def _run_one(code):
    """Exec/eval one snippet in the shared namespace; return (ok, result, out, err)."""
    g = globals()
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            try:
                # expression? -> eval and repr the value
                compiled = compile(code, "<mcp>", "eval")
                val = eval(compiled, g)
                res = "" if val is None else repr(val)
            except SyntaxError:
                exec(compile(code, "<mcp>", "exec"), g)
                res = ""
        return True, res, buf.getvalue(), ""
    except Exception:
        return False, "", buf.getvalue(), traceback.format_exc()


def _mcp_poll():
    """Main-loop timer callback: run any new request, write its response.
    Returns True to keep the timer alive (GLib convention)."""
    if not _mcp_state["running"]:
        return False
    d = _mcp_dir()
    reqp = os.path.join(d, _REQ)
    try:
        if os.path.exists(reqp):
            with open(reqp) as fh:
                req = json.load(fh)
            rid = int(req.get("id", 0))
            if rid > _mcp_state["last_id"]:
                _mcp_state["last_id"] = rid
                ok, res, out, err = _run_one(req.get("code", ""))
                resp = {"id": rid, "ok": ok, "result": res,
                        "stdout": out, "error": err}
                tmp = os.path.join(d, _RESP + ".tmp")
                with open(tmp, "w") as fh:
                    json.dump(resp, fh)
                os.replace(tmp, os.path.join(d, _RESP))
    except Exception:
        # never let a bad request kill the timer
        pass
    return True


def cootvalent_mcp_start(interval_ms=300):
    """Arm the MCP bridge on this session's main loop (call once)."""
    import coot
    if _mcp_state["running"]:
        print("[cootvalent-mcp] already running (dir: %s)" % _mcp_dir())
        return
    # reset the queue so a stale request isn't replayed
    d = _mcp_dir()
    for f in (_REQ, _RESP):
        try:
            os.remove(os.path.join(d, f))
        except OSError:
            pass
    _mcp_state.update(running=True, last_id=0, interval_ms=int(interval_ms))
    add = getattr(coot, "bandicoot_python_timeout_add", None)
    if add is None:
        _mcp_state["running"] = False
        print("[cootvalent-mcp] no main-loop timer (coot.bandicoot_python_timeout_add "
              "missing); cannot run the bridge on this build.")
        return
    add(int(interval_ms), _mcp_poll)
    print("[cootvalent-mcp] bridge armed. dir=%s  interval=%dms" % (d, interval_ms))
    print("[cootvalent-mcp] point the MCP server at COOTVALENT_MCP_DIR=%s" % d)


def cootvalent_mcp_stop():
    """Stop the bridge (the timer callback returns False and is removed)."""
    _mcp_state["running"] = False
    print("[cootvalent-mcp] bridge stopped.")
