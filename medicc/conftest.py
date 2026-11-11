"""Test setup for FastME, which MEDICC2 uses to infer tree topologies.

FastME is taken from the environment variable FASTME (path to the executable) or from PATH.
If FASTME is set, a temporary directory with a link named 'fastme' is put in front of PATH, so
that the command line and the API find it under its default name. Tests marked with
@pytest.mark.fastme are skipped if FastME is not found.
"""
import os
import shutil
import tempfile

import pytest

_link_dir = None


def find_test_fastme():
    """Returns the FastME executable for the tests (FASTME or PATH), or None."""
    return shutil.which(os.environ.get('FASTME', 'fastme'))


def pytest_configure(config):
    global _link_dir
    config.addinivalue_line('markers', 'fastme: the test needs the FastME executable '
                                       '(environment variable FASTME or fastme on PATH)')
    fastme = find_test_fastme()
    if 'FASTME' in os.environ and fastme is not None:
        _link_dir = tempfile.mkdtemp(prefix='medicc-test-fastme-')
        os.symlink(os.path.abspath(fastme), os.path.join(_link_dir, 'fastme'))
        os.environ['PATH'] = _link_dir + os.pathsep + os.environ.get('PATH', '')


def pytest_unconfigure(config):
    if _link_dir is not None:
        shutil.rmtree(_link_dir, ignore_errors=True)


def pytest_collection_modifyitems(config, items):
    if find_test_fastme() is not None:
        return
    skip = pytest.mark.skip(reason='FastME not found: set FASTME to the fastme executable or put '
                                   'fastme on PATH (conda install -c bioconda fastme)')
    for item in items:
        if item.get_closest_marker('fastme') is not None:
            item.add_marker(skip)
