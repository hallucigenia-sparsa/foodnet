"""PyInstaller entry for the Windows build: see foodnet.app.

Not called foodnet.py: a script of that name would itself be the module `foodnet` inside the program, and
`from foodnet.app import run` would then fail ("'foodnet' is not a package"), as CI showed on the rename.
The program is still foodnet.exe, from build.py's --name."""
import sys

from foodnet.app import run

if __name__ == "__main__":
    sys.exit(run())
