#!/usr/bin/env bash
# Install the Cootvalent extension suite into ~/.coot-preferences/
# (Coot autoloads every *.py there at startup.) Additive + reversible.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/.coot-preferences"
mkdir -p "$DEST"
for f in ligand-from-smiles.py covalent-link.py covalent-detect.py cootvalent-keys.py cootvalent-mcp-bridge.py; do
  cp "$SRC/$f" "$DEST/$f" && echo "installed  $f  ->  $DEST/"
done
echo
echo "Done. Restart bcoot (or coot)."
echo "  - Classic Coot (with coot_python): a 'Cootvalent' menu appears."
echo "  - bandicoot (no coot_python): use the console functions"
echo "      cv_declare()  cv_propagate()  cv_full()"
echo "    (keys are OFF by default -- native accelerators can shadow/destroy;"
echo "     opt in with cootvalent_bind_keys(...) once you confirm a free key)."
echo "(bandicoot_backend.py is a separate EXTERNAL MCP driver, not a Coot extension -- not installed here.)"
echo "Uninstall:  rm $DEST/{ligand-from-smiles,covalent-link,covalent-detect,cootvalent-keys}.py"
