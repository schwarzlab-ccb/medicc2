#!/usr/bin/env bash
# Build and run the C++ lesson tests (tests/migration/test_*.cpp), in step order.
#
#   ./test_fstengine.sh            all C++ tests
#   ./test_fstengine.sh 05         only tests whose name contains "05" (e.g. test_05_types.cpp)
#   ./test_fstengine.sh 05 07      several filters
#
# Needs the medicc-migration env active (OpenFST headers and libfst in $CONDA_PREFIX).
# Binaries go to build/migration/. Exit code is 0 only if every test built and passed.
set -u
cd "$(dirname "$0")"

if [[ -z "${CONDA_PREFIX:-}" || ! -f "$CONDA_PREFIX/include/fst/fst.h" ]]; then
    echo "OpenFST not found: run 'mamba activate medicc-migration' first" >&2
    exit 2
fi
CXX="${CXX:-c++}"
mkdir -p build/migration

# Keep CLion's compile database in sync with new headers/tests (local helper; skipped if absent).
if [[ -f tools/clion_compile_commands.py ]]; then
    "$CONDA_PREFIX/bin/python" tools/clion_compile_commands.py > /dev/null \
        || echo "warning: could not update compile_commands.json" >&2
fi

tests=()
for f in $(ls tests/migration/test_*.cpp 2>/dev/null | sort -V); do
    if [[ $# -eq 0 ]]; then
        tests+=("$f")
    else
        for pattern in "$@"; do
            [[ "$(basename "$f")" == *"$pattern"* ]] && { tests+=("$f"); break; }
        done
    fi
done
if [[ ${#tests[@]} -eq 0 ]]; then
    echo "no C++ tests found${1:+ matching: $*}" >&2
    exit 2
fi

passed=0
failed=()
for f in "${tests[@]}"; do
    name="$(basename "$f" .cpp)"
    printf '%-40s ' "$name"
    if ! log="$("$CXX" -std=c++17 -I. -I"$CONDA_PREFIX/include" "$f" \
                -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lfst \
                -o "build/migration/$name" 2>&1)"; then
        echo "BUILD FAILED"
        echo "$log" | head -40 | sed 's/^/    /'
        failed+=("$name")
        continue
    fi
    if out="$("build/migration/$name" 2>&1)"; then
        echo "PASS"
        [[ -n "$out" ]] && echo "$out" | sed 's/^/    /'
        passed=$((passed + 1))
    else
        echo "FAIL (exit $?)"
        [[ -n "$out" ]] && echo "$out" | sed 's/^/    /'
        failed+=("$name")
    fi
done

echo
echo "$passed passed, ${#failed[@]} failed"
[[ ${#failed[@]} -eq 0 ]] || { printf '  failed: %s\n' "${failed[@]}"; exit 1; }
