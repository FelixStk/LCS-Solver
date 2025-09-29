#ifndef LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_RT_H_
#define LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_RT_H_

#include <sstream>

#include "algorithms/BaseAlgorithm.h"
#include "algorithms/LCS/OutputType.h"
#include "constraints/BaseConstraint.h"
#include "constraints/ConstraintType.h"
#include "structures/Matrix.h"
#include "structures/RangeTree2D.h"
#include "util/CommonTypes.hpp"

namespace lcs_solver::algorithms::lcs {

template<OutputType T>
class Reconstructor_RT {
 public:
  using uint = std::size_t;
  using Pair = std::pair<uint, uint>;
  using KeyPointVec = std::vector<std::vector<Pair>>;

  using AlgorithmPtr = ::lcs_solver::algorithms::BaseAlgorithm *;
  using BaseConstraint = ::lcs_solver::constraints::BaseConstraint;
  using ConstraintType = ::lcs_solver::constraints::ConstraintType;
  using ConstraintMap = ::lcs_solver::algorithms::BaseAlgorithm::ConstraintMap;
  using Tree = ::lcs_solver::structures::RangeTree2D<uint, uint>;
  using TreeVec = std::vector<Tree>;
  using StringView = ::lcs_solver::util::StringView;
  using StringViewVec = std::vector<StringView>;

  //== Declaration of Iterator =================================================
  class Iterator {
   public:
    using RangeIter = Tree::Result::Iterator;

    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type *;
    using reference = value_type &;

    // Iterator(RangeIter start, RangeIter end, const TreeVec &t);
    // Iterator &operator++();
    // bool operator!=(const Iterator &other) const;
    // bool operator==(const Iterator &other) const;
    // [[nodiscard]] value_type operator*() const;
    // [[nodiscard]] std::string DebugString() const;

   private:
    // void advance();
    // RangeIter current, end;
  };//== end of Iterator ======================================================

  explicit Reconstructor_RT(const AlgorithmPtr &algo);

  Iterator cbegin() const;
  Iterator cend() const;
  auto begin() const { return cbegin(); }
  auto end() const { return cend(); }
  [[nodiscard]] bool empty() const;
  [[nodiscard]] std::string DebugString() const;

  std::tuple<uint, uint, uint, uint> getWindow(uint i, uint j, uint gap) const;

 private:
  Tree::Result queryStartPoints() const {
    const Pair r0 = std::make_pair(0, svv[0].size());
    const Pair r1 = std::make_pair(0, svv[1].size());
    const uint llcs = kpv[0].size();
    return tree[llcs - 1].query(r0, r1);
  }

  std::vector<StringView> svv;
  KeyPointVec kpv;
  TreeVec tree;
  AlgorithmPtr algo;
  Tree::Result startPoints;
};

#include "algorithms/LCS/Reconstructor_RT.tpp"

}// namespace lcs_solver::algorithms::lcs
#endif//LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_RT_H_