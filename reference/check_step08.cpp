// Checks the Step 8 reference (or your own messages.hpp) against the oracle's two-bin answers:
//   join cost == K and resolve cost == T on all 81x81 pairs, join ancestor sets == Oracle.optimum().
// Build and run from NEW's root, with the env active:
//   g++ -std=c++17 -O2 -I. -Imedicc/fstengine/src -isystem $CONDA_PREFIX/include \
//     ../medicc-full-fst-lloh-wgd2-capped-x2/reference/check_step08.cpp \
//     -L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -lfst -o build/check_step08 && build/check_step08
// To check YOUR messages.hpp instead, include "medicc/fstengine/src/messages.hpp" below instead of the reference.
#define REFDIR "../medicc-full-fst-lloh-wgd2-capped-x2/reference/"
#include <cassert>
#include <iostream>
#include <set>
#include <chrono>
#include <sstream>
#include "medicc/fstengine/src/search.hpp"
#include "../medicc-full-fst-lloh-wgd2-capped-x2/reference/step08_messages.hpp"  // or your messages.hpp

static Model load(const char *path) {
    Model m;
    std::ifstream in(path);
    std::string name; int n, st, na;
    while (in >> name >> n >> st >> na) {
        Vec f;
        for (int i = 0; i < n; ++i) f.AddState();
        f.SetStart(st);
        for (int i = 0, w; i < n; ++i) { in >> w; if (w < INF) f.SetFinal(i, W(w)); }
        for (int i = 0, s, a, b, w, t; i < na; ++i) { in >> s >> a >> b >> w >> t; f.AddArc(s, Arc(a + 1, b + 1, W(w), t)); }
        fst::ArcSort(&f, fst::ILabelCompare<Arc>());
        if (name == "A") m.a = std::move(f);
        if (name == "AD") m.ad = std::move(f);
    }
    assert(m.a.NumStates() > 0 && m.ad.NumStates() > 0);
    return m;
}
static Message word(const std::string &s) { return {linear(s), 0, int(s.size())}; }
static std::set<std::string> words(const Vec &f) {          // all words of an acyclic acceptor
    std::set<std::string> out;
    std::vector<std::pair<int, std::string>> stack{{f.Start(), ""}};
    while (!stack.empty()) {
        auto [s, w] = stack.back(); stack.pop_back();
        if (f.Final(s) != W::Zero()) out.insert(w);
        for (fst::ArcIterator<Vec> it(f, s); !it.Done(); it.Next())
            stack.push_back({it.Value().nextstate, w + symbol(it.Value().ilabel)});
    }
    return out;
}
int main(int argc, char **argv) {
    Model m = load("medicc/objects/ancestral_parts.txt");
    // ties and correlations (same cases as the oracle's Step 2 test)
    auto j1 = join(m, word("1"), word("2"));
    assert(j1.offset == 1 && words(j1.fst) == std::set<std::string>({"1", "2"}));
    auto j2 = join(m, word("11"), word("22"));
    assert(j2.offset == 1 && words(j2.fst) == std::set<std::string>({"11", "22"}));
    auto j3 = join(m, word("1X1"), word("2X2"));
    std::cout << "join(1X1, 2X2): offset " << j3.offset << ", " << words(j3.fst).size() << " ancestors\n";
    // offsets carry: join of two joins
    auto j4 = join(m, j1, word("1"));
    assert(j4.offset == 1 && words(j4.fst) == std::set<std::string>({"1"}));
    // resolve: one best child, cost includes the child's offset
    auto [c, cost] = resolve(m, "1", j1);
    assert(c == "1" && cost == 1); // 1 -> 1 costs 0, plus the message's offset 1
    std::cout << "resolve(1, {1,2}+1) = " << c << " cost " << cost << "\n";
    assert(resolve(m, "0", word("1")).second == infinity);

    std::ifstream g(REFDIR "step08_golden.txt");
    std::string x, y, Ts, Ks; int n = 0, bad = 0;
    auto t0 = std::chrono::steady_clock::now();
    while (g >> x >> y >> Ts >> Ks) {
        double T = std::stod(Ts), K = std::stod(Ks);
        double kj = infinity;
        try { kj = join(m, word(x), word(y)).offset; } catch (const std::runtime_error &) {}
        double tr = resolve(m, x, word(y)).second;
        if (kj != K || tr != T) { if (bad++ < 5) std::cout << "DIFF " << x << " " << y << " K " << K << "/" << kj << " T " << T << "/" << tr << "\n"; }
        ++n;
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << n << " pairs: join == K and resolve == T for all but " << bad << " (" << secs << " s)\n";
    // full ancestor sets vs the oracle's optimum() (ties across and within WGD branches)
    std::ifstream op(REFDIR "step08_optimum.txt");
    std::string line; int sets = 0, setbad = 0;
    while (std::getline(op, line)) {
        std::istringstream ls(line); std::string a, b, list; ls >> a >> b >> list;
        std::set<std::string> want; std::istringstream ws(list); std::string w;
        while (std::getline(ws, w, ',')) want.insert(w);
        auto got = words(join(m, word(a), word(b)).fst);
        if (got != want && setbad++ < 5) std::cout << "SET DIFF " << a << " " << b << ": " << got.size() << " vs " << want.size() << "\n";
        ++sets;
    }
    std::cout << sets << " pairs: join ancestor set == oracle optimum for all but " << setbad << "\n";
    return bad != 0 || setbad != 0;
}
