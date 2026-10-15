#!/usr/bin/env python3
"""Check feasible bounds and strict cutoffs against literal fixed-WGD relations.

The test-only C ABI is compiled in a temporary directory. No test hooks are
added to the shipped library, and no external dataset or baseline checkout is
needed. Run after scripts/build.py in the normal MEDICC build environment.
"""
import ctypes as C
import itertools
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from oracle import get_parts, Primitive, relation, minplus


def main():
    parts, symbols = get_parts()
    count = 0
    with tempfile.TemporaryDirectory(prefix='medicc-pair-bounds-') as folder:
        source = Path(folder)/'probe.cpp'
        library = Path(folder)/'probe.so'
        source.write_text('#include "' + str(ROOT/'medicc/cext/full_fst/direct_pair/src/pairwise.cpp') + '"\n' + r'''
extern "C" int check_score(void *handle, const char *x, const char *y, int branch, int room) {
    auto *model = static_cast<AM *>(handle);
    model->scalar_cache.clear();
    double stats[6] = {};
    factorized::Chunk chunk(model, x, y);
    return chunk.score(model, branch, room, stats);
}
extern "C" int check_witness(void *handle, const char *x, const char *y, int branch) {
    auto *model = static_cast<AM *>(handle);
    factorized::Chunk chunk(model, x, y);
    chunk.prepare(model);
    int left = branch / 3, right = branch % 3;
    int upper = pair_direct::chromosome_witness(x, y, left, right);
    return pair_direct::greedy_witness(*chunk.left[left], *chunk.right[right], chunk.mask, upper);
}
''')
        subprocess.run([os.getenv('CXX', 'g++'), '-O3', '-std=c++17', '-pthread',
                        '-shared', '-fPIC', '-DANCESTRAL_MAX_WGDS=2',
                        str(source), '-o', str(library)], check=True)
        lib = C.CDLL(str(library))
        lib.am_load.argtypes = [C.c_char_p]; lib.am_load.restype = C.c_void_p
        lib.am_free.argtypes = [C.c_void_p]
        lib.check_score.argtypes = [C.c_void_p, C.c_char_p, C.c_char_p, C.c_int, C.c_int]
        lib.check_score.restype = C.c_int
        lib.check_witness.argtypes = [C.c_void_p, C.c_char_p, C.c_char_p, C.c_int]
        lib.check_witness.restype = C.c_int
        handle = lib.am_load(str(ROOT/'medicc/cext/full_fst/direct_pair/build/ancestral_parts.txt').encode())
        assert handle
        try:
            rng = random.Random(20261015)
            for length in [1, 2]:
                words = [''.join(w) for w in itertools.product('012345678', repeat=length)]
                ids = {word: i for i, word in enumerate(words)}
                a = minplus(relation(Primitive(parts['L_LOH']), words, symbols),
                            relation(Primitive(parts['G']), words, symbols))
                doubling = np.full_like(a, np.inf)
                for i, word in enumerate(words):
                    doubled = ''.join(str(min(8, 2*int(c))) for c in word)
                    doubling[i, ids[doubled]] = 0
                branches = [a]
                for _ in range(2):
                    branches.append(minplus(minplus(branches[-1], doubling), a))
                pairs = (list(itertools.product(range(9), repeat=2)) if length == 1
                         else [(rng.randrange(len(words)), rng.randrange(len(words))) for _ in range(256)])
                for left in range(3):
                    for right in range(3):
                        exact = minplus(branches[left].T, branches[right])
                        branch = 3*left + right
                        for i, j in pairs:
                            x, y = words[i].encode(), words[j].encode()
                            value = int(exact[i, j])
                            assert lib.check_witness(handle, x, y, branch) >= value, (x, y, branch, 'witness')
                            count += 1
                            for room in sorted({0, max(0, value-1), value, value+1}):
                                result = lib.check_score(handle, x, y, branch, room)
                                expected = value if value < room else 100000000
                                assert result == expected, (x, y, branch, room, result, expected)
                                count += 1
        finally:
            lib.am_free(handle)
    print(f'PASS: {count} literal fixed-WGD witness/cutoff checks', flush=True)


if __name__ == '__main__':
    main()
