"""PGO configuration and compiler guards; real compilation is checked separately."""
from pathlib import Path
import runpy
import sys
from types import SimpleNamespace

import Cython.Build
import pytest
import setuptools
from setuptools.command.build_ext import build_ext

ROOT = Path(__file__).resolve().parents[1]


def configure(monkeypatch, tmp_path, platform, mode):
    captured = {}
    monkeypatch.setattr(sys, 'platform', platform)
    monkeypatch.setenv('MEDICC_PGO', mode)
    monkeypatch.setenv('MEDICC_PGO_DIR', str(tmp_path))
    monkeypatch.setattr(Cython.Build, 'cythonize', lambda extensions, **kwargs: extensions)
    monkeypatch.setattr(setuptools, 'setup', lambda **kwargs: captured.update(kwargs))
    namespace = runpy.run_path(str(ROOT / 'setup.py'))
    return namespace, captured


@pytest.mark.parametrize('platform', ['darwin', 'linux'])
@pytest.mark.parametrize('mode', ['', 'generate', 'use'])
def test_pgo_flags_only_apply_to_native_engine(monkeypatch, tmp_path, platform, mode):
    (tmp_path / ('default.profdata' if platform == 'darwin' else 'training.gcda')).touch()
    _, setup = configure(monkeypatch, tmp_path, platform, mode)
    native = setup['ext_modules'][-1]
    flags = native.extra_compile_args
    assert native.name == 'medicc.native._engine'
    for extension in setup['ext_modules'][:-1]:
        assert not any('profile' in flag for flag in extension.extra_compile_args + extension.extra_link_args)
    if mode == 'generate':
        assert f'-fprofile-generate={tmp_path}' in flags
        assert '-fprofile-update=atomic' in flags
    elif mode == 'use':
        target = tmp_path / 'default.profdata' if platform == 'darwin' else tmp_path
        assert f'-fprofile-use={target}' in flags
    else:
        assert not any('profile' in flag for flag in flags + native.extra_link_args)


@pytest.mark.parametrize('platform', ['darwin', 'linux'])
def test_pgo_use_requires_profiles(monkeypatch, tmp_path, platform):
    with pytest.raises(RuntimeError, match='No .*profile'):
        configure(monkeypatch, tmp_path, platform, 'use')


@pytest.mark.parametrize('platform,macros,accepted', [
    ('darwin', '#define __clang__ 1\n#define __GNUC__ 4', True),
    ('darwin', '#define __GNUC__ 12', False),
    ('linux', '#define __GNUC__ 12', True),
    ('linux', '#define __clang__ 1\n#define __GNUC__ 4', False),
])
def test_pgo_compiler_guard(monkeypatch, tmp_path, platform, macros, accepted):
    namespace, _ = configure(monkeypatch, tmp_path, platform, 'generate')
    monkeypatch.setattr(namespace['subprocess'], 'run',
                        lambda *args, **kwargs: SimpleNamespace(stdout=macros))
    called = []
    monkeypatch.setattr(build_ext, 'build_extensions', lambda self: called.append(True))
    command = namespace['NativeBuildExt'](setuptools.Distribution())
    command.compiler = SimpleNamespace(compiler_cxx=['test-c++'])
    if accepted:
        command.build_extensions()
        assert called == [True]
    else:
        with pytest.raises(RuntimeError, match='PGO'):
            command.build_extensions()
        assert called == []


@pytest.mark.parametrize('platform,mode', [('darwin', 'invalid'), ('win32', 'generate')])
def test_invalid_pgo_configuration(monkeypatch, tmp_path, platform, mode):
    with pytest.raises(RuntimeError, match='PGO'):
        configure(monkeypatch, tmp_path, platform, mode)
