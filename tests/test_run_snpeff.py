"""Unit tests for bin/run_snpeff.py.

The behaviour under test is genome selection. The genome key is conventionally
influenza_<subtype>, but which keys exist is a property of the reference bundle: a
BLAST hit to a subtype the bundle has no database for used to reach snpEff as a
database name that does not exist, which fails deep inside the task with snpEff's own
error rather than naming the subtype.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "bin" / "run_snpeff.py"

BUNDLE_GENOMES = ["influenza_h1n1", "influenza_h3n2", "influenza_h5n1", "influenza_h9n2"]


def _load_module():
    # bin/ is on PATH in a task, and run_snpeff.py imports `from lib.tools import ...`
    # the way Python resolves a script's own directory. Reproduce that here.
    sys.path.insert(0, str(REPO_ROOT / "bin"))
    spec = importlib.util.spec_from_file_location("run_snpeff", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rse = _load_module()


@pytest.mark.parametrize("subtype,expected", [
    ("H5N1", "influenza_h5n1"),
    ("h5n1", "influenza_h5n1"),
    ("H1N1", "influenza_h1n1"),
])
def test_selects_the_bundle_genome_for_a_supported_subtype(subtype, expected):
    assert rse.select_genome(subtype, BUNDLE_GENOMES) == expected


def test_unsupported_subtype_names_the_subtype_and_the_alternatives():
    with pytest.raises(SystemExit) as excinfo:
        rse.select_genome("H7N9", BUNDLE_GENOMES)

    message = str(excinfo.value)
    assert "H7N9" in message
    assert "influenza_h7n9" in message
    # Listing what the bundle does have is the point: the fix is either a different
    # bundle or a rebuilt one, and the message has to say which are available.
    assert "influenza_h5n1" in message


def test_empty_genome_list_is_reported_rather_than_matching_nothing():
    with pytest.raises(SystemExit) as excinfo:
        rse.select_genome("H5N1", [])
    assert "none" in str(excinfo.value)
