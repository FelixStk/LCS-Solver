#ifndef LCS_SOLVER_ALGORITHMS_LCS_LCS2_STD_H_
#define LCS_SOLVER_ALGORITHMS_LCS_LCS2_STD_H_

#include <memory>
#include <stack>
#include <unordered_set>

#include "algorithms/BaseAlgorithm.h"          // Required for BaseAlgorithm
#include "algorithms/solutions/BaseCollector.h"// Required for Iterator

// Forward declaration to reduce includes
namespace lcs_solver::algorithms::llcs {
class LLCS2_Algorithm;
}

namespace lcs_solver::algorithms::solutions {
class BaseCollector;
class BaseIterator;
class AnyIterator;
}// namespace lcs_solver::algorithms::solutions

namespace lcs_solver::algorithms::lcs {
class LCS2_STD_S final : public BaseAlgorithm {
 public:
  class StackSolution final : public solutions::BaseCollector {
   public:
    using AnyIterator = solutions::AnyIterator;
    using BaseIterator = solutions::BaseIterator;

    class Iterator final : public BaseIterator {
     public:
      using Pair = std::pair<uint, uint>;
      struct PairHash {
        std::size_t operator()(const Pair &p) const {
          return std::hash<uint>{}(p.first) ^ (std::hash<uint>{}(p.second) << 1);
        }
      };
      using State = std::pair<Pair, uint>;
      using Stack = std::stack<State>;
      Iterator();
      explicit Iterator(const std::weak_ptr<llcs::LLCS2_Algorithm> &algo, bool end);

      void advance() override;
      [[nodiscard]] std::unique_ptr<BaseIterator> clone() const override;

     private:
      std::shared_ptr<llcs::LLCS2_Algorithm> algo_;
      std::unique_ptr<solutions::Points> current_;
      std::vector<std::vector<uint>> dp_;
      uint llcs_{};
      Stack stack_;
      std::vector<std::unordered_set<Pair, PairHash>> processed_;
    };//== end of Iterator =====================================================

    explicit StackSolution(const std::shared_ptr<llcs::LLCS2_Algorithm> &ptr);

    [[nodiscard]] AnyIterator begin() const override;
    [[nodiscard]] AnyIterator end() const override;
    [[nodiscard]] BaseSolution *clone() const override;

   private:
    const std::shared_ptr<llcs::LLCS2_Algorithm> algo_;
  };//== end of StackSolution ==============================================

  static constexpr const char *name = "LCS2_STD_S";

  LCS2_STD_S(const StrPtrVector &spv, const ConstraintMap &map);

  // BaseAlgorithm Implementation
  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;

 private:
  std::shared_ptr<llcs::LLCS2_Algorithm> algo_;
};

}// namespace lcs_solver::algorithms::lcs
#endif//LCS_SOLVER_ALGORITHMS_LCS_LCS2_STD_H_