#ifndef LCS_SOLVER_ALGORITHMS_BASEALGORITHM_H_
#define LCS_SOLVER_ALGORITHMS_BASEALGORITHM_H_

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/AlgoCategory.h"
#include "algorithms/AlgoType.h"
#include "algorithms/BaseSolution.h"
#include "constraints/BaseConstraint.h"
#include "constraints/ConstraintMap.h"
#include "constraints/ConstraintType.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms {

class BaseAlgorithm {
 public:
  using AlgorithmPtr = std::unique_ptr<BaseAlgorithm>;
  using BaseConstraint = constraints::BaseConstraint;
  using ConstraintType = constraints::ConstraintType;
  using ConstraintMap = constraints::ConstraintMap;

  using StringViewVector = std::vector<::lcs_solver::util::StringView>;
  using StrPtr = std::shared_ptr<const ::lcs_solver::util::String>;
  using StrPtrVector = std::vector<StrPtr>;
  using Symbol = ::lcs_solver::util::Symbol;
  using uint = util::uint;

  enum class State : unsigned char {
    Invalid = 0x00,
    Constructed = 0x01,
    Initialized = 0x02,
    Preprocessed = 0x03,
    Processed = 0x04,
    Updated = 0x05
  };
  enum class ResetLevel : unsigned char {
    Full = 0x01,
    Query = 0x02
  };

  virtual ~BaseAlgorithm() = default;

  [[nodiscard]] static std::string toString(const StringViewVector& svv);
  [[nodiscard]] static std::string toString(const StrPtrVector& spv);

  [[nodiscard]] virtual bool isValid() const = 0;
  [[nodiscard]] virtual std::string_view getName() const = 0;
  [[nodiscard]] virtual std::string_view getDescription() const = 0;
  [[nodiscard]] virtual std::string DebugString() const = 0;
  [[nodiscard]] virtual std::unique_ptr<BaseSolution> query() = 0;
  virtual void doPreprocessing() = 0;
  virtual void reset(ResetLevel resetLevel) = 0;

  [[nodiscard]] AlgoType getType() const;// Used in factory
  [[nodiscard]] AlgoCategory getCategory() const;
  [[nodiscard]] const ConstraintMap &getConstraintMap() const;
  [[nodiscard]] bool isFinished() const;
  [[nodiscard]] State getState() const;
  [[nodiscard]] StringViewVector getStringViewVec() const;
  [[nodiscard]] const StrPtrVector &getStrPtrVec() const;
  [[nodiscard]] StrPtrVector genSortedStrPtrVec(bool increasing = true) const;
  void reset();

  //=== Input and Output of Problems ===========================================
  friend std::istream &operator>>(std::istream &in, AlgorithmPtr &ptr);
  friend std::ostream &operator<<(std::ostream &os, const AlgorithmPtr &p);

 protected:
  BaseAlgorithm(AlgoCategory category, AlgoType type, const StrPtrVector &spv, const ConstraintMap &map);

  [[nodiscard]] bool isEachConstraintIndividuallyValid() const;
  [[nodiscard]] bool isFilled() const;///< No StringView has length 0
  [[nodiscard]] bool isSorted() const;///< ascending by StringView length
  [[nodiscard]] static size_t getLongestStrLen(const StrPtrVector &vec);
  void setState(State s);
  void SortStringViews(bool increasing = true);

  const AlgoType type;               ///< unique identifier of the algorithm used while defining problems
  const AlgoCategory category;       ///< describes the kind of algorithm
  const ConstraintMap &constraintMap;///< container used to search for a constraint
  const StrPtrVector &strPtrVec;     ///< std::vector of pointers to the Strings in a problem
  StringViewVector s;          ///< std::vector of StringViews for the Strings given to the algorithm

 private:
  State state;
};

}// namespace lcs_solver::algorithms
#endif// LCS_SOLVER_ALGORITHMS_BASEALGORITHM_H_