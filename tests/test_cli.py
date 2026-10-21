"""Smoke tests for the resdel command line interface.

These do not run any GROMACS commands; they check that the package imports
cleanly and that every subcommand is wired up and responds to --help.
"""

import subprocess
import sys

import pytest


def run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "resdel.cli.main", *args],
        capture_output=True,
        text=True,
    )


def test_package_imports():
    import resdel

    assert resdel.__version__


@pytest.mark.parametrize("command", ["prepare", "transform", "analyze"])
def test_subcommand_help(command):
    """Every advertised subcommand exists and documents itself."""
    result = run_cli(command, "--help")
    assert result.returncode == 0, result.stderr
    assert command in result.stdout


def test_missing_config_fails_cleanly(tmp_path):
    """A missing config file is an error, not a traceback."""
    result = run_cli("prepare", str(tmp_path / "does_not_exist.yaml"))
    assert result.returncode != 0
