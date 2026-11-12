#pragma once

// Bounded shortest-path search (Dijkstra) over delayed OpenFST graphs. OpenFST's
// ShortestDistance does the heap and the relaxation; our queue adds two things:
// a stopping rule (stop once the next state costs more than the best answer so far)
// and a record of every settled state. From those, search() copies out the "tight"
// graph: every path whose cost equals the optimum, so ties are kept.
// Costs are nonnegative, so a state that already costs more than the best answer
// cannot lead to a better one. The test is strictly ">": a state with equal cost
// can still lead to another, equally good answer.
#include "fst_operations.hpp"

// Queue order: smaller distance first. 'distance' is the vector ShortestDistance fills.
struct Priority {
    const std::vector<W> *distance;
    bool operator()(int a, int b) const {
        return (*distance)[a].Value() < (*distance)[b].Value();
    }
};

class BoundQueue : public fst::ShortestFirstQueue<int, Priority> {
    using Base = fst::ShortestFirstQueue<int, Priority>;
    const Fst &graph;
    const std::vector<W> &distance;
    bool ties; // true: keep every optimal path; false: the cost alone is enough

public:
    double best; // starts at the cap, then the best complete answer found
    bool found = false; // has any final state been reached within the cap?
    std::vector<int> settled;

    // constructor in C++
    BoundQueue(const Fst &f, const std::vector<W> &d, double cap, bool all)
        : Base(Priority{&d}), graph(f), distance(d), ties(all), best(cap) {}

    // OpenFST stops when the queue says it is empty. Say so early when the cheapest
    // waiting state cannot beat (or, without ties, cannot equal) the best answer.
    bool Empty() const override {
        if (Base::Empty())
            return true;
        double next = distance[Base::Head()].Value();
        return !std::isfinite(next) || next > best || (!ties && found && next >= best);
    }
    // Called when OpenFST settles the head state: update the best answer, record it
    void Dequeue() override {
        int s = Base::Head();
        double value = distance[s].Value() + graph.Final(s).Value();
        if (std::isfinite(value) && value <= best) {
            best = value;
            found = true;
        }
        settled.push_back(s);
        Base::Dequeue();
    }
};

struct Search {
    double cost = infinity; // infinity: nothing within the cap
    Vec tight; // all optimal paths (only if retain = true)
};

// Cheapest path through 'graph' that costs at most 'cap'. With 'retain', also return
// the tight graph: every optimal path, with costs reset to 0 (the cost is in 'cost').
static Search search(const Fst &graph, double cap, bool retain) {
    Search out;
    if (cap < 0)
        return out;
    std::vector<W> distance;
    BoundQueue queue(graph, distance, cap, retain);
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

    // Dijkstra settles each state once, but OpenFST re-enqueues a state whose
    // distance improves after it was dequeued (possible with A*, Step 27).
    // Keep each settled state once.
    // std::sort(queue.settled.begin(), queue.settled.end());
    // queue.settled.erase(std::unique(queue.settled.begin(), queue.settled.end()),
    //                     queue.settled.end());

    // Reconstruct the equally good FSA
    // Copy out the settled states that are not too expensive ...
    std::vector<int> remap(distance.size(), -1);
    for (int s : queue.settled)
        if (distance[s].Value() <= out.cost)
            remap[s] = out.tight.AddState();
    int start = graph.Start();
    if (start < 0 || remap[start] < 0)
        throw std::runtime_error("Missing start");
    out.tight.SetStart(remap[start]);
    // ... and only the arcs and finals that lie on an optimal path:
    // d(u) + w ==  d(v) for an arc, d(u) + final(u) == best for a final state
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
    fst::Connect(&out.tight);
    check(out.tight);
    return out;
}

// Clean up an acceptor: remove epsilons, determinize, minimize, sort arcs
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

// One cheapest path of 'f', spelled out as a string, and its cost.
static std::pair<std::string, double> choose(const Vec &f, bool output_label = false) {
    Vec path;
    fst::ShortestPath(f, &path, 1);
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