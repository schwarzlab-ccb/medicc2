// REFERENCE for migration Step 10 (not part of REF; untracked). It is REF's native.cpp reduced
// to what Steps 8-11 need. Left out until later steps:
//   am_start / am_final / am_arc / am_import (message inspection for test_fstengine.py, Step 41)
//   fstengine_history (Step 33), fstengine_finish_row / fstengine_clear (score cache, Steps 24/31)
//   loss / gain / doubling / a_even in Model (history and bounds); only A and AD are read here.
//
// fstengine: the C interface used by medicc/fstengine/model.py (ctypes).
// The model is P = A + ADA + ADADA with A = L_LOH o G and the capped doubling D(c) = min(2c, 8);
// each WGD costs 1. All headers are included here once: one translation unit, one OpenFST build.
// C functions never throw (an exception must not cross into C / Python): on an error they store
// the text in `error` (thread_local, fst_operations.hpp) and return nullptr or -1; am_error()
// then returns that text.
//
// Build by hand (Step 12 moves this into setup.py):
//   "${CXX:-c++}" -std=c++17 -O2 -shared -fPIC -I"$CONDA_PREFIX/include" medicc/fstengine/src/fstengine.cpp \
//     -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" -lfst -o medicc/fstengine/_fstengine.so
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "types.hpp"          // Model, Message, Vec, Arc, W, INF
#include "fst_operations.hpp" // error, output, linear
#include "search.hpp"
#include "messages.hpp"       // join, resolve
#include "scoring.hpp"        // score

// Read the model parts file of Step 4. Blocks: "NAME n start narcs", one line of n final
// weights (INF = not final), then narcs lines "source ilabel olabel weight target".
// File labels are 0..8 and X = 9; the engine adds 1 because label 0 is epsilon in OpenFST.
// Only A and AD are kept; other blocks (L_LOH, G, ...) are read and skipped.
static Model load_model(const char *path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error(std::string("Cannot open model file ") + path);
    Model m;
    std::string name;
    int n, start, narcs;
    while (in >> name >> n >> start >> narcs) {
        if (n <= 0 || start < 0 || start >= n || narcs < 0)
            throw std::runtime_error("Bad part header: " + name);
        Vec f;
        for (int i = 0; i < n; ++i)
            f.AddState();
        f.SetStart(start);
        for (int i = 0, w; i < n; ++i) {
            if (!(in >> w))
                throw std::runtime_error("Bad final weight in " + name);
            if (w < INF)
                f.SetFinal(i, W(w));
        }
        for (int i = 0, s, a, b, w, t; i < narcs; ++i) {
            if (!(in >> s >> a >> b >> w >> t))
                throw std::runtime_error("Bad arc in " + name);
            f.AddArc(s, Arc(a + 1, b + 1, W(w), t));
        }
        fst::ArcSort(&f, fst::ILabelCompare<Arc>());
        if (name == "A")
            m.a = std::move(f);
        else if (name == "AD")
            m.ad = std::move(f);
    }
    if (m.a.NumStates() == 0 || m.ad.NumStates() == 0)
        throw std::runtime_error("Model file has no A or AD");
    return m;
}

extern "C" {
const char *am_error() {
    return error.c_str();
}

// Load the model; no_wgd != 0 restricts it to A (no WGD at all).
void *fstengine_load(const char *path, int no_wgd) {
    try {
        auto m = std::make_unique<Model>(load_model(path));
        m->max_wgds = no_wgd ? 0 : 2;
        return m.release();
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}
void am_free(void *m) {
    delete static_cast<Model *>(m);
}

// A message that accepts only `s`, at cost 0 (a leaf of the tree).
void *am_string(const char *s) {
    try {
        return new Message{linear(s), 0, int(std::string(s).size())};
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}
void am_drop(void *p) {
    delete static_cast<Message *>(p);
}
double am_offset(void *p) {
    return static_cast<Message *>(p)->offset;
}

void *fstengine_join(void *m, void *l, void *r) {
    try {
        auto &model = *static_cast<Model *>(m);
        return new Message(join(model, *static_cast<Message *>(l), *static_cast<Message *>(r)));
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}

// Returns the best child profile; its cost goes to *cost. The string lives in `output`
// (thread_local) until the next call on this thread, so Python copies it at once.
const char *fstengine_resolve(void *m, const char *parent, void *child, double *cost) {
    try {
        auto &model = *static_cast<Model *>(m);
        auto r = resolve(model, parent, *static_cast<Message *>(child));
        if (!std::isfinite(r.second))
            throw std::runtime_error("No reachable child");
        output = r.first;
        *cost = r.second;
        return output.c_str();
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}

// directed != 0: T(x, y); otherwise K(x, y). inf if impossible, -1 on an error.
double fstengine_score(void *m, const char *x, const char *y, int directed) {
    try {
        auto &model = *static_cast<Model *>(m);
        return score(model, x, y, directed != 0);
    } catch (const std::exception &e) {
        error = e.what();
        return -1;
    }
}
}
