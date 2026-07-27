#!/usr/bin/env bash
# Install the Cootvalent extension suite into ~/.coot-preferences/
# (Coot autoloads every *.py there at startup.) Additive + reversible.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
DEST="$HOME/.coot-preferences"
mkdir -p "$DEST"
for f in ligand-from-smiles.py covalent-link.py covalent-detect.py cootvalent-keys.py; do
  cp "$SRC/$f" "$DEST/$f" && echo "installed  $f  ->  $DEST/"
done
echo
echo "Done. Restart bcoot (or coot)."
echo "  - Classic Coot (with coot_python): a 'Cootvalent' menu appears."
echo "  - bandicoot (no coot_python): the menu can't be added, but key bindings"
echo "    work -- Ctrl+L declare, Ctrl+P propagate, Ctrl+R full refine (over the"
echo "    graphics window); the functions are also callable from the console."
echo "(bandicoot_backend.py is a separate EXTERNAL MCP driver, not a Coot extension -- not installed here.)"
echo "Uninstall:  rm $DEST/{ligand-from-smiles,covalent-link,covalent-detect,cootvalent-keys}.py"
