#ifndef LCS_SOLVER_ALGORITHMS_LCS_LCS2_RT_H_
#define LCS_SOLVER_ALGORITHMS_LCS_LCS2_RT_H_

#include <memory>
#include <vector>

#include "algorithms/BaseAlgorithm.h"            // Required for BaseAlgorithm
#include "algorithms/solutions/BaseCollector.h"  // Required for Iterator
#include "structures/RangeTree2D.h"              // Required due to unique_ptr

// Forward declarations to reduce includes
namespace lcs_solver::algorithms::llcs {
class LLCS2_Algorithm;
}

namespace lcs_solver::algorithms::solutions {
class BaseCollector;
class BaseIterator;
class AnyIterator;
}  // namespace lcs_solver::algorithms::solutions

namespace lcs_solver::algorithms::lcs {

/*******************************************************************************
 * @brief LCS Reconstruction with the help of RangeTrees
 ******************************************************************************/
class LCS2_RT final : public BaseAlgorithm {
 public:
  class RangeTreeSolution;

  static constexpr const char* name = "LCS2_RT";

  LCS2_RT(const StrPtrVector& spv, const ConstraintMap& map, AlgoType type);

  // BaseAlgorithm Implementation
  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

  [[nodiscard]] const BaseAlgorithm* getBase() const;

 private:
  std::shared_ptr<llcs::LLCS2_Algorithm> algo_;
};

/*******************************************************************************
 * @brief Solution type of LCS2_RT.query()
 * @details Implements the BaseCollector interface.
 ******************************************************************************/
class LCS2_RT::RangeTreeSolution final : public solutions::BaseCollector {
 public:
  using AnyIterator = solutions::AnyIterator;
  using BaseIterator = solutions::BaseIterator;

  class Iterator;

  explicit RangeTreeSolution(const std::shared_ptr<llcs::LLCS2_Algorithm>& ptr);

  [[nodiscard]] AnyIterator begin() const override;
  [[nodiscard]] AnyIterator end() const override;
  [[nodiscard]] BaseSolution* clone() const override;

 private:
  const std::shared_ptr<llcs::LLCS2_Algorithm> algo_;
};

/*******************************************************************************
 * @brief Iterator used for LCS2_RT::RangeTreeSolution objects
 ******************************************************************************/
class LCS2_RT::RangeTreeSolution::Iterator final : public BaseIterator {
 public:
  using RangeTree2D = structures::RangeTree2D<uint, uint>;
  using Pair = std::pair<uint, uint>;
  Iterator();
  explicit Iterator(const std::weak_ptr<llcs::LLCS2_Algorithm>& algo, bool end);

  void advance() override;
  [[nodiscard]] std::unique_ptr<BaseIterator> clone() const override;

 private:
  [[nodiscard]] std::vector<std::shared_ptr<RangeTree2D>> genTrees(
      const std::shared_ptr<llcs::LLCS2_Algorithm>& ptr) const;
  void setEnd();
  void reversePre();
  void reversePost();

  std::shared_ptr<llcs::LLCS2_Algorithm> algo_;
  uint llcs;
  const bool aligned_;
  bool do_reversing_ = false;
  std::vector<std::shared_ptr<RangeTree2D>> tree_;// tree_[i] contains a point q = (x,y) if LLCS(s1[1:x],s2[1:y]) == i and
                                                  //          LLCS2_Algorithm::isExtensible(Pair p, Pair q, uint gap) where
                                                  //          p is some point in tree_[i-1]. These points are tracked while
                                                  //          a dp-LLCS2_Algorithm is running. (see LLCS2_Algorithm::track)
  std::vector<RangeTree2D::Result> result_;       // result_[i] stores responses when we query tree_[i]
  std::vector<RangeTree2D::Result::Iterator> it_; // it_[i] iterates through result_[i]
  std::vector<RangeTree2D::Result::Iterator> end_;// end_[i] is set to result_[i].end()
  std::unique_ptr<solutions::Points> current_;    // stack and vector like adaptor to store an embedding
  uint lvl;                                       // index in tree_ where we try to advance to the next lcs
};

}  // namespace lcs_solver::algorithms::lcs
#endif  // LCS_SOLVER_ALGORITHMS_LCS_LCS2_RT_H_