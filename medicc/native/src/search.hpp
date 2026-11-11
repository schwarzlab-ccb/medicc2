#pragma once

// Bounded A* over delayed OpenFST graphs. OpenFST owns heap operations and
// relaxation; this adapter supplies only an incumbent stopping rule and records
// settled IDs for finite tight-path export (all optimal paths are kept). Also
// canonical() (determinize/minimize) and choose() (one shortest path).
// Nonnegative costs imply no state with distance > incumbent can be optimal.
// Strict > is essential: equality states contain additional optimal words.
struct Priority {
    const std::vector<W> *distance;
    const Heuristic *heuristic;
    double value(int s) const {
        return (*distance)[s].Value() + (heuristic ? (*heuristic)(s) : 0);
    }
    bool operator()(int a, int b) const {
        double x = value(a), y = value(b);
        if (x != y)
            return x < y;
        // Prefer deeper/cheaper-to-go states on f ties; all ties are still kept.
        return (*distance)[a].Value() > (*distance)[b].Value();
    }
};
class BoundQueue : public fst::ShortestFirstQueue<int, Priority> {
    using Base = fst::ShortestFirstQueue<int, Priority>;
    const Fst &graph;
    const std::vector<W> &distance;
    const Heuristic *heuristic;
    bool ties;
    ComposeBudget *budget;
    double priority() const {
        int s = Base::Head();
        return distance[s].Value() + (heuristic ? (*heuristic)(s) : 0);
    }

  public:
    double best;
    bool found = false;
    std::vector<int> settled;
    BoundQueue(const Fst &f, const std::vector<W> &d, double cap, bool all, const Heuristic *h,
               ComposeBudget *context = nullptr)
        : Base(Priority{&d, h}), graph(f), distance(d), heuristic(h), ties(all), budget(context),
          best(cap) {}
    bool Empty() const override {
        if (Base::Empty())
            return true;
        double f = priority();
        return !std::isfinite(f) || f > best || (!ties && found && f >= best);
    }
    void Dequeue() override {
        int s = Base::Head();
        double value = distance[s].Value() + graph.Final(s).Value();
        if (std::isfinite(value) && value <= best) {
            best = value;
            found = true;
        }
        if (budget) {
            budget->cap = best;
            budget->active_state = s;
            budget->active_g = distance[s].Value();
        }
        settled.push_back(s);
        Base::Dequeue();
    }
};
struct Search {
    double cost = infinity;
    Vec tight;
};
static Search search(const Fst &graph, double cap, bool retain,
                     const Heuristic *heuristic = nullptr, ComposeBudget *budget = nullptr,
                     std::vector<int> *original_ids = nullptr) {
    Search out;
    if (cap < 0)
        return out;
    std::vector<W> distance;
    BoundQueue queue(graph, distance, cap, retain, heuristic, budget);
    fst::ShortestDistanceOptions<Arc, BoundQueue, fst::AnyArcFilter<Arc>> opts(&queue, {});
    fst::ShortestDistance(graph, &distance, opts);
    check(graph);
    if (distance.size() == 1 && !distance[0].Member())
        throw std::runtime_error("OpenFST shortest distance failed");
    if (!queue.found)
        return out;
    out.cost = queue.best;
    if (!retain)
        return out;
    std::sort(queue.settled.begin(), queue.settled.end());
    queue.settled.erase(std::unique(queue.settled.begin(), queue.settled.end()),
                        queue.settled.end());
    // Preserve all tight arcs, including ties that did not trigger a relaxation.
    // This is a finite extraction of visited OpenFST IDs, not another composer.
    std::vector<int> remap(distance.size(), -1);
    for (int s : queue.settled)
        if (distance[s].Value() <= out.cost)
            remap[s] = out.tight.AddState();
    int start = graph.Start();
    if (start < 0 || remap[start] < 0)
        throw std::runtime_error("Missing start");
    out.tight.SetStart(remap[start]);
    for (int s : queue.settled) {
        if (remap[s] < 0)
            continue;
        if (distance[s].Value() + graph.Final(s).Value() == out.cost)
            out.tight.SetFinal(remap[s], W::One());
        for (fst::ArcIterator<Fst> it(graph, s); !it.Done(); it.Next()) {
            const auto &a = it.Value();
            if (a.nextstate < int(remap.size()) && remap[a.nextstate] >= 0 &&
                distance[s].Value() + a.weight.Value() == distance[a.nextstate].Value())
                out.tight.AddArc(remap[s], Arc(a.ilabel, a.olabel, W::One(), remap[a.nextstate]));
        }
    }
    if (original_ids) {
        original_ids->resize(out.tight.NumStates());
        for (int s : queue.settled)
            if (remap[s] >= 0)
                (*original_ids)[remap[s]] = s;
    } else
        fst::Connect(&out.tight);
    check(out.tight);
    return out;
}
static Vec canonical(Vec f) {
    Vec out;
    fst::RmEpsilon(&f);
    fst::Determinize(f, &out);
    check(out);
    fst::Minimize(&out);
    check(out);
    fst::ArcSort(&out, fst::ILabelCompare<Arc>());
    return out;
}
static std::pair<std::string, double> choose(const Vec &f, bool output_label = false) {
    Vec path;
    fst::ShortestPath(f, &path);
    check(path);
    std::string word;
    double cost = 0;
    int s = path.Start();
    if (s < 0)
        return {word, infinity};
    while (true) {
        fst::ArcIterator<Vec> it(path, s);
        if (it.Done()) {
            cost += path.Final(s).Value();
            break;
        }
        auto a = it.Value();
        int l = output_label ? a.olabel : a.ilabel;
        if (l)
            word += symbol(l);
        cost += a.weight.Value();
        s = a.nextstate;
    }
    return {word, cost};
}
