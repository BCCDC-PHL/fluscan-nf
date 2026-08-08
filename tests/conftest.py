"""Shared pytest setup.

`bin/` goes on `sys.path` so the entry points under test can resolve their `from lib.…`
imports. The scripts themselves must not do this: Nextflow puts `bin/` on `PATH` and Python
puts a script's own directory on `sys.path`, so a production script needs no path
manipulation — only pytest, which imports them by file location, does.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(REPO_ROOT / "bin"))
