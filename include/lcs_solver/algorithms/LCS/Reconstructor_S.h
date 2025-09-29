#ifndef LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_S_H_
#define LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_S_H_

#include <stack>

#include "algorithms/BaseAlgorithm.h"
#include "algorithms/LCS/OutputType.h"
#include "structures/Matrix.h"
#include "constraints/ConstraintType.h"
#include "constraints/BaseConstraint.h"
#include "structures/RangeTree2D.h"

namespace lcs_solver::algorithms::lcs {

template <OutputType T>
class Reconstructor_S {
 public:
  using uint = std::size_t;
  using Pair = std::pair<uint, uint>;
  using Row = std::vector<uint>;
  using Matrix = std::vector<Row>;

  using AlgorithmPtr = ::lcs_solver::algorithms::BaseAlgorithm *;
  using BaseConstraint = ::lcs_solver::constraints::BaseConstraint;
  using ConstraintType = ::lcs_solver::constraints::ConstraintType;
  using ConstraintMap = ::lcs_solver::algorithms::BaseAlgorithm::ConstraintMap;
  using Tree = ::lcs_solver::structures::RangeTree2D<uint, uint>;
  using TreeArray = std::vector<Tree>;

  //== Declaration of Iterator =================================================
  class Iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    Iterator(const TreeArray& t);
    Iterator &operator++();
    bool operator!=(const Iterator &other) const;
    bool operator==(const Iterator &other) const;
    [[nodiscard]] value_type operator*() const;
    [[nodiscard]] std::string DebugString() const;
   private:
    void advance();
    std::stack<Pair> stack;
  }; //== end of Iterator ======================================================

  Reconstructor_S(
      const AlgorithmPtr p
  );

  std::string DebugString() const;

  Iterator cbegin() const;
  Iterator cend() const;
  auto begin() const { return cbegin(); }
  auto end() const { return cend(); }
  [[nodiscard]] bool empty() const;

 private:
  Matrix mat;
  AlgorithmPtr p;
};

#include "algorithms/LCS/Reconstructor_S.tpp"

} // lcs_solver::algorithms::lcs

#endif //LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_S_H_
