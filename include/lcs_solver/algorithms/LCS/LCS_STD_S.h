#ifndef LCS_SOLVER_ALGORITHMS_LCS_LCS_STD_FL_HPP
#define LCS_SOLVER_ALGORITHMS_LCS_LCS_STD_FL_HPP

#include <stack>
#include <vector>

#include "algorithms/BaseAlgorithm.h"
#include "algorithms/LLCS/LLCS_STD_FL.h"
#include "algorithms/solutions/BaseCollector.h"
#include "algorithms/solutions/BaseIterator.h"
#include "structures/Matrix.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::lcs {
class LCS_STD_S final : public BaseAlgorithm {
 public:
  using AnyIterator = ::lcs_solver::algorithms::solutions::AnyIterator;
  using BaseIterator = ::lcs_solver::algorithms::solutions::BaseIterator;
  using BaseCollector = ::lcs_solver::algorithms::solutions::BaseCollector;
  using LLCS_STD_FL = ::lcs_solver::algorithms::llcs::LLCS_STD_FL;
  using uint = LLCS_STD_FL::uint;
  using Matrix = LLCS_STD_FL::Matrix;
  using StrViewVector = std::vector<util::StringView>;

  static constexpr const char* name = "LCS_STD_S";

  class StackSolution final : public BaseCollector {
   public:
    class StackIterator final : public BaseIterator {
     public:
      using Points = ::lcs_solver::algorithms::solutions::Points;

      explicit StackIterator(
          const Matrix& m, const StrPtrVector& spv, bool end);

      void advance() override;
      [[nodiscard]] std::unique_ptr<BaseIterator> clone() const override;
      [[nodiscard]] std::string DebugString() const override;

     private:
      static bool matchingSymbol(
          const Matrix& mat, const Matrix::uint& pos, const StrViewVector& svv);

      [[nodiscard]] std::string describe(
          const std::string& msg, Matrix::uint pos) const;

      const StrPtrVector& spv;
      const Matrix& dp;
      const StringViewVector& s;

      std::unique_ptr<Points> sol;
      std::stack<Matrix::uint> stack;
    };

    // Constructor
    explicit StackSolution(const StrPtrVector& spv, const Matrix& matrix);

    // BaseCollector Implementations
    [[nodiscard]] AnyIterator begin() const override;
    [[nodiscard]] AnyIterator end() const override;
    [[nodiscard]] BaseSolution* clone() const override;

   private:
    const StrPtrVector& spv;
    const Matrix& m;
  };

  LCS_STD_S(const StrPtrVector& vec, const ConstraintMap& map);

  // BaseAlgorithm Implementation
  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel l) override;

 private:
  std::unique_ptr<LLCS_STD_FL> algo;
};

}  // namespace lcs_solver::algorithms::lcs
#endif  // LCS_SOLVER_ALGORITHMS_LCS_LCS_STD_FL_HPP
