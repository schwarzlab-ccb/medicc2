// Native OpenFST copy-number engine: C interface used by medicc/native/model.py.
// The model is P = A + ADA + ADADA with A = L_LOH o G and the capped doubling
// D(c) = min(2c, 8); each WGD costs 1. The headers below are internal to this
// single translation unit (one OpenFST compilation) and in dependency order.
// Functions return nullptr / a negative value on error; am_error() then holds
// the message of the calling thread.
#include "types.hpp"
#include "fst_operations.hpp"
#include "astar.hpp"
#include "pair_filter.hpp"
#include "search.hpp"
#include "bounds.hpp"
#include "messages.hpp"
#include "scoring.hpp"
#include "history.hpp"

extern "C" {
const char *am_error() {
    return error.c_str();
}
static void *load_model(const char *path) {
    try {
        auto m = std::make_unique<Model>();
        std::ifstream in(path);
        std::string name;
        int n, st, na;
        while (in >> name >> n >> st >> na) {
            if (n <= 0 || st < 0 || st >= n || na < 0)
                throw std::runtime_error("Bad primitive header");
            Vec f;
            for (int i = 0; i < n; ++i)
                f.AddState();
            f.SetStart(st);
            for (int i = 0, w; i < n; ++i) {
                if (!(in >> w))
                    throw std::runtime_error("Bad final");
                if (w < INF)
                    f.SetFinal(i, W(w));
            }
            for (int i = 0, s, a, b, w, t; i < na; ++i) {
                if (!(in >> s >> a >> b >> w >> t))
                    throw std::runtime_error("Bad arc");
                f.AddArc(s, Arc(a + 1, b + 1, W(w), t));
            }
            fst::ArcSort(&f, fst::ILabelCompare<Arc>());
            if (name == "L_LOH")
                m->loss = f;
            if (name == "G")
                m->gain = f;
            if (name == "A")
                m->a = std::move(f);
            if (name == "AD")
                m->ad = std::move(f);
        }
        if (m->a.NumStates() == 0 || m->ad.NumStates() == 0)
            throw std::runtime_error("Missing A/AD");
        m->doubling.AddState();
        m->doubling.SetStart(0);
        m->doubling.SetFinal(0, W::One());
        for (int c = 0; c <= 9; ++c) {
            int target = c == 9 ? 9 : c == 0 ? 0 : std::min(8, c * 2);
            m->doubling.AddArc(0, Arc(c + 1, target + 1, W::One(), 0));
        }
        for (int q = 0; q < m->a.NumStates(); ++q)
            m->a_even.AddState();
        m->a_even.SetStart(m->a.Start());
        for (int q = 0; q < m->a.NumStates(); ++q) {
            m->a_even.SetFinal(q, m->a.Final(q));
            for (fst::ArcIterator<Vec> it(m->a, q); !it.Done(); it.Next()) {
                auto a = it.Value();
                if (a.ilabel == 10 || (a.ilabel - 1) % 2 == 0)
                    m->a_even.AddArc(q, a);
            }
        }
        m->a4 = quotient(m->a, false, 4);
        m->ad4 = quotient(m->ad, true, 4);
        m->a2 = quotient(m->a, false, 2);
        m->ad2 = quotient(m->ad, true, 2);
        return m.release();
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}
// Load the model primitives; no_wgd != 0 restricts the model to A (no WGD).
void *nf_load(const char *path, int no_wgd) {
    auto *model = static_cast<Model *>(load_model(path));
    if (model)
        model->max_wgds = no_wgd ? 0 : 2;
    return model;
}
void am_free(void *m) {
    delete static_cast<Model *>(m);
}
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
// Message inspection and import, used by the tests to enumerate and build
// weighted ancestor languages.
double am_offset(void *p) {
    return static_cast<Message *>(p)->offset;
}
int am_start(void *p) {
    return static_cast<Message *>(p)->fst.Start();
}
int am_final(void *p, int s) {
    double v = static_cast<Message *>(p)->fst.Final(s).Value();
    return std::isfinite(v) ? int(v) : INF;
}
int am_arc(void *p, int s, int c, int *w) {
    const auto &f = static_cast<Message *>(p)->fst;
    for (fst::ArcIterator<Vec> it(f, s); !it.Done(); it.Next())
        if (it.Value().ilabel == c + 1) {
            *w = int(it.Value().weight.Value());
            return it.Value().nextstate;
        }
    return -1;
}
void *am_import(int len, int n, int start, int offset, const int *rows) {
    try {
        auto m = std::make_unique<Message>();
        m->offset = offset;
        m->length = len;
        for (int s = 0; s < n; ++s)
            m->fst.AddState();
        if (start >= 0)
            m->fst.SetStart(start);
        for (int s = 0; s < n; ++s) {
            const int *r = rows + 21 * s;
            if (r[0] < INF)
                m->fst.SetFinal(s, W(r[0]));
            for (int c = 0; c < 10; ++c)
                if (r[1 + 2 * c] >= 0)
                    m->fst.AddArc(s, Arc(c + 1, c + 1, W(r[2 + 2 * c]), r[1 + 2 * c]));
        }
        fst::Connect(&m->fst);
        if (m->fst.Start() >= 0) {
            std::vector<W> d;
            fst::ShortestDistance(m->fst, &d, true);
            double total = d[m->fst.Start()].Value();
            if (!std::isfinite(total) || total + m->offset >= INF)
                throw std::runtime_error("Imported cost exceeds range");
            fst::Push(&m->fst, fst::REWEIGHT_TO_INITIAL, fst::kShortestDelta, true);
            m->offset += int(total);
            check(m->fst);
        }
        return m.release();
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}
void *nf_join(void *m, void *l, void *r) {
    try {
        auto &model = *static_cast<Model *>(m);
        return new Message(join(model, *static_cast<Message *>(l), *static_cast<Message *>(r)));
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}
const char *nf_resolve(void *m, const char *p, void *c, double *cost) {
    try {
        auto &model = *static_cast<Model *>(m);
        auto r = resolve(model, p, *static_cast<Message *>(c));
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
double nf_score(void *m, const char *x, const char *y, int dir) {
    try {
        auto &model = *static_cast<Model *>(m);
        return score(model, x, y, dir);
    } catch (const std::exception &e) {
        error = e.what();
        return -1;
    }
}
const char *nf_history(void *p, const char *x, const char *y) {
    try {
        auto &model = *static_cast<Model *>(p);
        output = history(model, x, y);
        return output.c_str();
    } catch (const std::exception &e) {
        error = e.what();
        return nullptr;
    }
}
// Release the per-handle score caches (between rows of the distance matrix).
void nf_clear(void *p) {
    auto &m = *static_cast<Model *>(p);
    m.chromosome_scores.clear();
    m.suffix_cache.clear();
}
}
