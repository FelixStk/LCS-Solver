#ifndef LCS_SOLVER_ALGORITHMS_LCS_OUTPUTTYPE_H_
#define LCS_SOLVER_ALGORITHMS_LCS_OUTPUTTYPE_H_

#include <tuple>
#include <cstddef>
#include <concepts>
#include <vector>
#include "structures/Embedding.h"

namespace lcs_solver::algorithms::lcs {

using PointVector = std::vector<std::pair<size_t, size_t>>;
using Embedding = ::lcs_solver::structures::Embedding;
using String = std::string;

template <typename T>
concept OutputType = std::same_as<T, PointVector> || std::same_as<T, Embedding> || std::same_as<T, String>;

}
#endif //LCS_SOLVER_ALGORITHMS_LCS_OUTPUTTYPE_H_
