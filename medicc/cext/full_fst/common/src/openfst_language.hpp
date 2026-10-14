#pragma once
#include "automata_types.hpp"
#include <cstddef>
namespace library_language {
DFA determinized(const std::vector<direct_nfa::Tight> &input, int start,
                 int len, int offset, std::size_t &states);
DFA united(const std::vector<DFA> &inputs);
}
