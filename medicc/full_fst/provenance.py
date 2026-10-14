"""Local code verification and runtime provenance; never import sibling repositories."""

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_code():
    for name in ["code_manifest.json", "build_manifest.json"]:
        path = ROOT / name
        if not path.is_file():
            raise RuntimeError(f"Missing {name}; run scripts/build.py")
        for relative, expected in json.loads(path.read_text()).items():
            if sha256(ROOT / relative) != expected:
                raise RuntimeError(
                    f"Local artifact changed: {relative}; review changes and rerun scripts/build.py"
                )


def record_dependencies(out):
    modules = {}
    for name in [
        "medicc",
        "fstlib",
        "medicc.full_fst.model",
        "medicc.full_fst.pair_scorer",
        "medicc.full_fst.scheduler",
        "medicc.full_fst.ancestors",
        "medicc.full_fst.core",
        "medicc.full_fst._bindings",
    ]:
        module = sys.modules.get(name)
        if module is not None and getattr(module, "__file__", None):
            path = Path(module.__file__).resolve()
            if ROOT not in path.parents:
                raise RuntimeError(
                    f"Unexpected external implementation: {name}: {path}"
                )
            modules[name] = str(path.relative_to(ROOT))
    versions = {}
    for name in [
        "numpy",
        "pandas",
        "scipy",
        "biopython",
        "matplotlib",
        "pyyaml",
        "joblib",
    ]:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    (Path(out) / "runtime.json").write_text(
        json.dumps(
            dict(
                python=sys.version,
                platform=platform.platform(),
                modules=modules,
                dependencies=versions,
                code_manifest=sha256(ROOT / "code_manifest.json"),
                build_manifest=sha256(ROOT / "build_manifest.json"),
            ),
            indent=2,
        )
        + "\n"
    )
