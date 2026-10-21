"""Shared pytest configuration for the resdel test suite."""

import shutil

import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_gromacs: test needs a working `gmx` executable on PATH",
    )


def pytest_collection_modifyitems(config, items):
    """Skip GROMACS-dependent tests when no `gmx` binary is available."""
    if shutil.which("gmx") is not None:
        return

    skip_gmx = pytest.mark.skip(reason="no `gmx` executable on PATH")
    for item in items:
        if "requires_gromacs" in item.keywords:
            item.add_marker(skip_gmx)
