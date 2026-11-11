#pragma once

// Ancestor messages. join() splits two child messages into independent blocks
// (chromosome boundaries where both sides have a unique state), searches every
// WGD branch combination in bound order and keeps the union of all minimum-cost
// ancestor languages. resolve() picks one optimal child profile for a fixed
// parent profile (down pass).

static Message slice(const Message &m, const Layers &ls, int begin, int end) {
    Message out;
    out.length = end - begin;
    std::vector<int> ids(m.fst.NumStates(), -1);
    for (int d = begin; d <= end; ++d)
        for (int s : ls.states[d])
            ids[s] = out.fst.AddState();
    out.fst.SetStart(ids[ls.states[begin].at(0)]);
    for (int d = begin; d <= end; ++d)
        for (int s : ls.states[d]) {
            if (d == end) {
                out.fst.SetFinal(ids[s], end == m.length ? m.fst.Final(s) : W::One());
                continue;
            }
            for (fst::ArcIterator<Vec> it(m.fst, s); !it.Done(); it.Next()) {
                auto a = it.Value();
                a.nextstate = ids[a.nextstate];
                out.fst.AddArc(ids[s], a);
            }
        }
    return out;
}
static std::vector<std::pair<Message, Message>> split(const Message &l, const Message &r) {
    if (l.length != r.length)
        throw std::runtime_error("Length mismatch");
    auto a = layers(l), b = layers(r);
    std::vector<std::pair<Message, Message>> parts;
    int begin = 0;
    for (int end = 1; end <= l.length; ++end) {
        // A unique state on each side after X proves independence; otherwise
        // retain correlations by leaving the full block intact.
        if (end != l.length && !(a.states[end].size() == 1 && b.states[end].size() == 1 &&
                                 a.after_x[end] && b.after_x[end]))
            continue;
        parts.emplace_back(slice(l, a, begin, end), slice(r, b, begin, end));
        begin = end;
    }
    if (l.length == 0)
        parts.emplace_back(l, r);
    return parts;
}
static Message join(Model &m, const Message &l, const Message &r) {
    auto parts = split(l, r);
    auto bounds = message_bounds(m, parts, false);
    std::vector<int> order(parts.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return parts[a].first.length < parts[b].first.length; });
    double best = infinity;
    Vec combined;
    auto branches = branch_order(m, bounds, false, best);
    for (auto [lower, branch] : branches) {
        if (lower > best)
            continue; // Equality must retain all optimal words.
        int k = branch / 3, j = branch % 3;
        double total = k + j, remaining = lower - total;
        std::vector<Vec> blocks(parts.size());
        bool feasible = true;
        for (int i : order) {
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            double room = best - total - remaining;
            auto graph = make_graph(m, parts[i].first, parts[i].second, k, j, false);
            auto result = search(*graph, room, true);
            if (!std::isfinite(result.cost)) {
                feasible = false;
                break;
            }
            total += result.cost;
            blocks[i] = canonical(std::move(result.tight));
        }
        if (!feasible)
            continue;
        Vec language = std::move(blocks[0]);
        for (size_t i = 1; i < blocks.size(); ++i)
            fst::Concat(&language, blocks[i]);
        if (total < best) {
            combined.DeleteStates();
            best = total;
        }
        if (total == best)
            fst::Union(&combined, language);
    }
    if (!std::isfinite(best) || combined.Start() < 0 || best + l.offset + r.offset >= INF)
        throw std::runtime_error("No finite common ancestor or cost exceeds range");
    return {canonical(std::move(combined)), int(best) + l.offset + r.offset, l.length};
}
static std::pair<std::string, double> resolve(Model &m, const std::string &parent,
                                              const Message &child) {
    Message p{linear(parent), 0, int(parent.size())};
    auto parts = split(p, child);
    auto bounds = message_bounds(m, parts, true);
    std::string bestword;
    double best = infinity;
    bool found = false;
    auto branches = branch_order(m, bounds, true, best);
    for (auto [lower, branch] : branches) {
        if (lower > best || (found && lower == best))
            continue;
        double total = branch, remaining = lower - total;
        std::string word;
        bool feasible = true;
        for (size_t i = 0; i < parts.size(); ++i) {
            remaining -= bounds[i].lower[branch];
            if (total + bounds[i].lower[branch] + remaining > best) {
                feasible = false;
                break;
            }
            auto graph = make_graph(m, parts[i].first, parts[i].second, 0, branch, true);
            auto result =
                search(*graph, best - total - remaining, true);
            if (!std::isfinite(result.cost)) {
                feasible = false;
                break;
            }
            total += result.cost;
            word += choose(result.tight, true).first;
        }
        if (feasible && (!found || total < best)) {
            best = total;
            bestword = std::move(word);
            found = true;
        }
    }
    if (!found)
        return {"", infinity};
    return {bestword, best + child.offset};
}
