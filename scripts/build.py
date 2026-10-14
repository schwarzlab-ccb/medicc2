#!/usr/bin/env python3
"""Build all native engines locally; optionally rebuild the vendored OpenFST wrappers."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
PREFIX = Path(os.getenv('CONDA_PREFIX', sys.prefix))

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--record-only', action='store_true', help='Refresh reviewed code hashes without recompiling')
    ap.add_argument('--legacy', action='store_true', help='Also rebuild FSTlib and MEDICC Cython wrappers against this Python/OpenFST environment')
    args = ap.parse_args()
    if args.record_only and args.legacy: ap.error('--record-only and --legacy cannot be combined')
    for name, source, library in [('ancestors', 'ancestors.cpp', 'libdirect_wgd2.so'),
                                  ('direct_pair', 'pairwise.cpp', 'libpair_direct.so'),
                                  ('reference_pair', 'pairwise.cpp', 'libpair_direct.so')]:
        if args.record_only: continue
        folder = ROOT/'medicc/cext/full_fst'/name
        build = folder/'build'; build.mkdir(exist_ok=True)
        table = folder/'assets/ancestral_parts.txt'
        (build/'ancestral_parts.txt').write_bytes(table.read_bytes())
        temporary = build/(library+'.tmp')
        command = [os.getenv('CXX', 'g++'), '-O3', '-std=c++17', '-pthread', '-shared', '-fPIC',
                   '-DANCESTRAL_MAX_WGDS=2', str(folder/'src'/source), '-o', str(temporary)]
        if name == 'ancestors':
            command[command.index('-o'):command.index('-o')] = [
                str(ROOT/'medicc/cext/full_fst/common/src/openfst_language.cpp'),
                f'-I{PREFIX}/include', f'-L{PREFIX}/lib', f'-Wl,-rpath,{PREFIX}/lib', '-lfst']
        print(' '.join(command), flush=True)
        subprocess.run(command, check=True)
        temporary.replace(build/library)
    if args.legacy:
        subprocess.run([sys.executable, str(ROOT/'scripts/build_legacy.py'), 'build_ext',
                        '--build-lib', str(ROOT), '--build-temp', str(ROOT/'build/legacy')], cwd=ROOT, check=True)
    for relative in ['medicc', 'fstlib/cext']:
        if not list((ROOT/relative).glob('*.so')):
            raise RuntimeError(f'Missing Python extensions in {relative}; build with --legacy')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    code = {}
    binaries = {}
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() or any(v in p.parts for v in ['__pycache__', '.git']): continue
        rel = p.relative_to(ROOT)
        if rel.parts[0] == 'build' or any(part.endswith('.egg-info') for part in rel.parts): continue
        if p.suffix == '.so' or ('build' in rel.parts and p.name == 'ancestral_parts.txt'):
            binaries[str(rel)] = digest(p)
        elif 'build' not in rel.parts and p.suffix in {'.py', '.pyx', '.pxd', '.cpp', '.hpp', '.h', '.yaml', '.sh', '.sbatch', '.txt'}:
            code[str(rel)] = digest(p)
    (ROOT/'code_manifest.json').write_text(json.dumps(code, indent=2)+'\n')
    (ROOT/'build_manifest.json').write_text(json.dumps(binaries, indent=2)+'\n')
    print('Recorded local code/build hashes.' if args.record_only else 'Built three native engines and recorded local code/build hashes.', flush=True)

if __name__ == '__main__': main()
