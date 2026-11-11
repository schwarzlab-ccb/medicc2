#pragma once
// Pair graphs contain no epsilon arcs: all primitives and linear profiles
// consume one symbol on each side. Use OpenFST's exact-label composition filter.
using PairMatcher = fst::SortedMatcher<Fst>;
using PairFilter = fst::NullComposeFilter<PairMatcher>;
using PairStateTable = fst::GenericComposeStateTable<Arc, PairFilter::FilterState>;

// Scalar pair scores through OpenFST's composition-filter extension. Each CN
// string is composed with A (or the even-input part of A below a doubling) once
// and cached with its forward/backward distances. A composition filter drops
// every arc whose prefix + arc + suffix bound exceeds the current incumbent,
// so the pair product is expanded only within the current search budget.
struct SuffixData {
    Vec fst;
    std::vector<W> forward, backward;
};
struct ComposeBudget {
    double cap = infinity, active_g = 0;
    int active_state = -1;
    PairStateTable *top_table = nullptr;
};
using Potential = std::function<double(int)>;
static double distance_at(const std::vector<W> &values, int state) {
    return state >= 0 && state < int(values.size()) ? values[state].Value() : infinity;
}
class CostComposeFilter : public PairFilter {
    std::shared_ptr<ComposeBudget> budget_;
    Potential left_h_, right_h_, left_p_, right_p_;
    bool top_ = false;
    double allowance_ = 0, prefix_ = 0;

  public:
    CostComposeFilter(const Fst &left, const Fst &right, PairMatcher *lm = nullptr,
                      PairMatcher *rm = nullptr)
        : PairFilter(left, right, lm, rm) {}
    CostComposeFilter(const CostComposeFilter &other, bool safe = false)
        : PairFilter(other, safe), budget_(other.budget_), left_h_(other.left_h_),
          right_h_(other.right_h_), left_p_(other.left_p_), right_p_(other.right_p_),
          top_(other.top_), allowance_(other.allowance_) {}
    void Configure(std::shared_ptr<ComposeBudget> budget, Potential lh, Potential rh, Potential lp,
                   Potential rp, bool top, double allowance) {
        budget_ = std::move(budget);
        left_h_ = std::move(lh);
        right_h_ = std::move(rh);
        left_p_ = std::move(lp);
        right_p_ = std::move(rp);
        top_ = top;
        allowance_ = allowance;
    }
    void SetState(int left, int right, const FilterState &fs) {
        PairFilter::SetState(left, right, fs);
        if (!budget_)
            return;
        prefix_ = left_p_(left) + right_p_(right);
        if (top_ && budget_->active_state >= 0) {
            auto tuple = budget_->top_table->Tuple(budget_->active_state);
            if (tuple.StateId1() == left && tuple.StateId2() == right &&
                tuple.GetFilterState() == fs)
                prefix_ = std::max(prefix_, budget_->active_g);
        }
    }
    FilterState FilterArc(Arc *left, Arc *right) const {
        auto state = PairFilter::FilterArc(left, right);
        if (state == FilterState::NoState() || !budget_)
            return state;
        double tail = left_h_(left->nextstate) + right_h_(right->nextstate);
        double bound = prefix_ + left->weight.Value() + right->weight.Value() + tail + allowance_;
        if (!std::isfinite(bound) || bound > budget_->cap)
            return FilterState::NoState();
        return state;
    }
    uint64_t Properties(uint64_t props) const {
        // Filtering removes arcs; do not advertise accessibility/completeness.
        return PairFilter::Properties(props) &
               (fst::kAcceptor | fst::kIDeterministic | fst::kODeterministic | fst::kNoEpsilons |
                fst::kNoIEpsilons | fst::kNoOEpsilons | fst::kILabelSorted | fst::kOLabelSorted |
                fst::kAcyclic | fst::kInitialAcyclic | fst::kTopSorted | fst::kUnweighted |
                fst::kError);
    }
};
struct SuffixSide {
    Ptr fst;
    Potential suffix, prefix;
};
static Ptr filtered_compose(const Fst &left, const Fst &right, PairStateTable **table,
                            std::shared_ptr<ComposeBudget> budget, Potential lh, Potential rh,
                            Potential lp, Potential rp, bool top = false, double allowance = 0) {
    auto sorted_right = (right.Properties(fst::kILabelSorted, false) & fst::kILabelSorted)
                            ? Ptr(right.Copy())
                            : sorted(right);
    fst::ComposeFstOptions<Arc, PairMatcher, CostComposeFilter, PairStateTable> options(
        fst::CacheOptions(true, 4 << 20));
    options.filter = new CostComposeFilter(left, *sorted_right);
    options.filter->Configure(std::move(budget), std::move(lh), std::move(rh), std::move(lp),
                              std::move(rp), top, allowance);
    options.state_table = new PairStateTable(left, *sorted_right);
    *table = options.state_table;
    return std::make_unique<fst::ComposeFst<Arc>>(left, *sorted_right, options);
}
static std::shared_ptr<const SuffixData> suffix_data(Model &model, const std::string &word,
                                                     bool even) {
    std::string key = (even ? "E" : "A") + word;
    auto found = model.suffix_cache.find(key);
    if (found != model.suffix_cache.end())
        return found->second;
    auto data = std::make_shared<SuffixData>();
    auto child = linear(word);
    data->fst = Vec(*compose(even ? model.a_even : model.a, child));
    fst::Connect(&data->fst);
    if (data->fst.Start() >= 0) {
        fst::ShortestDistance(data->fst, &data->backward, true);
        fst::ShortestDistance(data->fst, &data->forward);
    }
    fst::ArcSort(&data->fst, fst::ILabelCompare<Arc>());
    check(data->fst);
    model.suffix_cache.emplace(key, data);
    return data;
}
static SuffixSide suffix_side(Model &model, const std::string &word, int k,
                              std::shared_ptr<ComposeBudget> budget, double allowance) {
    auto data = suffix_data(model, word, k > 0);
    SuffixSide side{Ptr(new Vec(data->fst)),
                    [data](int s) { return distance_at(data->backward, s); },
                    [data](int s) { return distance_at(data->forward, s); }};
    Potential zero = [](int) { return 0.0; };
    for (int i = 0; i < k; ++i) {
        PairStateTable *table;
        auto graph = filtered_compose(model.ad, *side.fst, &table, budget, zero, side.suffix, zero,
                                      side.prefix, false, allowance);
        auto suffix = [table, h = side.suffix](int s) { return h(table->Tuple(s).StateId2()); };
        auto prefix = [table, p = side.prefix](int s) { return p(table->Tuple(s).StateId2()); };
        side = {std::move(graph), std::move(suffix), std::move(prefix)};
    }
    side.fst = std::make_unique<fst::ProjectFst<Arc>>(*side.fst, fst::ProjectType::INPUT);
    return side;
}
struct PairGraph {
    Ptr exact;
    std::shared_ptr<Heuristic> heuristic;
    std::shared_ptr<ComposeBudget> budget;
};
static PairGraph make_pair_graph(Model &model, const std::string &x, const std::string &y, int k,
                                 int j, double cap) {
    auto budget = std::make_shared<ComposeBudget>();
    auto ld = suffix_data(model, x, k > 0), rd = suffix_data(model, y, j > 0);
    budget->cap = cap;
    double lmin = distance_at(ld->backward, ld->fst.Start()),
           rmin = distance_at(rd->backward, rd->fst.Start());
    auto left = suffix_side(model, x, k, budget, rmin),
         right = suffix_side(model, y, j, budget, lmin);
    PairStateTable *table;
    auto exact = filtered_compose(*left.fst, *right.fst, &table, budget, left.suffix, right.suffix,
                                  left.prefix, right.prefix, true);
    budget->top_table = table;
    auto heuristic = std::make_shared<Heuristic>();
    heuristic->evaluate = [table, l = left.suffix, r = right.suffix](int state) {
        auto tuple = table->Tuple(state);
        return l(tuple.StateId1()) + r(tuple.StateId2());
    };
    return {std::move(exact), std::move(heuristic), std::move(budget)};
}
