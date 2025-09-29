#ifndef LCS_SOLVER_ALGORITHMS_LCS_LCS2_STD_H_H_
#define LCS_SOLVER_ALGORITHMS_LCS_LCS2_STD_H_H_

// #include <gtest/internal/gtest-port.h>

#include <generator>
#include <memory>

#include "algorithms/BaseAlgorithm.h"            // Required for BaseAlgorithm
#include "algorithms/solutions/BaseCollector.h"  // Required for Iterator

namespace lcs_solver::algorithms::llcs {
class LLCS2_Algorithm;
}

namespace lcs_solver::algorithms::solutions {
class BaseCollector;
class BaseIterator;
class AnyIterator;
}// namespace lcs_solver::algorithms::solutions

namespace lcs_solver::algorithms::lcs {

/*******************************************************************************
 * @brief LCS2 algorithm implementation using the hirschberg trick
 * @details The solution is a `BaseCollector`, which is a wrapper around
 *          `std::generator<std::vector<std::pair<util::uint,util::uint>>>`
 ******************************************************************************/
class LCS2_STD_H final : public BaseAlgorithm {
 public:
  /// @brief Alias for an index pair used in the algorithm.
  using Pair = std::pair<util::uint, util::uint>;
  using PairVec = std::vector<Pair>;

  /// @brief Generator type producing index pairs.
  using Generator = std::generator<std::vector<Pair>>;

  /// @brief Solution class for LCS2_STD_H based on std::generator.
  class GenSolution;

  /// @brief String for identifying the algorithm
  static constexpr auto name = "LCS2_STD_H";

  /// @brief Construct LCS2_STD_H from string pointer vector and constraints.
  LCS2_STD_H(const StrPtrVector &spv, const ConstraintMap &map);

  // Inherited methods
  [[nodiscard]] bool isValid() const override;
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string_view getDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override;
  void doPreprocessing() override;
  void reset(ResetLevel lvl) override;
};

/*******************************************************************************
 * @brief Solution type for LCS2_STD_H based on a generator.
 * @details Provides iteration over LCS2 results using a generator pattern.
 *          Implements the BaseCollector interface.
 ******************************************************************************/
class LCS2_STD_H::GenSolution final : public solutions::BaseCollector {
public:
  using AnyIterator = solutions::AnyIterator;
  using BaseIterator = solutions::BaseIterator;

  /// @brief Iterator over pairs in GenSolution.
  class GenIterator;

  /// @brief Construct a GenSolution with a generator.
  explicit GenSolution(util::StrPtrVector spv);

  /// @brief Begin iterator.
  [[nodiscard]] AnyIterator begin() const override;

  /// @brief End iterator.
  [[nodiscard]] AnyIterator end() const override;

  /// @brief Clone the solution.
  [[nodiscard]] BaseSolution *clone() const override;

private:
  util::StrPtrVector v_;
};

/*******************************************************************************
 * @brief Iterator for GenSolution using std::generator.
 ******************************************************************************/
class LCS2_STD_H::GenSolution::GenIterator final : public BaseIterator {
public:
  /// @brief Construct a GenIterator representing the end
  explicit GenIterator();

  /// @brief Construct a GenIterator from a generator reference.
  explicit GenIterator(const util::StrPtrVector& spv, Generator&& gen);

  /// @brief Move to the next element in the generator.
  void advance() override;

  /// @brief Clone the iterator.
  [[nodiscard]] std::unique_ptr<BaseIterator> clone() const override;

private:
  Generator generator_;
  std::unique_ptr<solutions::Points> current_;
  decltype(generator_.begin()) iter_;
  const util::StrPtrVector& spv_;
  static util::StrPtrVector empty_spv_;
  util::uint n_ = 0;
};

}// namespace lcs_solver::algorithms::lcs
#endif//LCS_SOLVER_ALGORITHMS_LCS_LCS2_STD_H_H_