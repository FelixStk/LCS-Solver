#ifndef LCS_SOLVER_CONSTRAINTS_GLOBAL_CONSTRAINT_LLCS_KNOWN_H_
#define LCS_SOLVER_CONSTRAINTS_GLOBAL_CONSTRAINT_LLCS_KNOWN_H_

#include "constraints/BaseConstraint.h"
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints::global {

class Constraint_LLCS_Known : public BaseConstraint {
 public:
  static constexpr ConstraintType kType = ConstraintType::LLCS_KNOWN;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::GLOBAL;
  static constexpr const char *kName = "LLCS_Known";
  static constexpr const char *kDescription =
      "The Length of the LCS is a constant `llcs`";
  static constexpr ConstraintType StaticType(){return kType;}

  explicit Constraint_LLCS_Known(Unsigned lengthOfLongestSubsequence);
  Constraint_LLCS_Known(): Constraint_LLCS_Known(static_cast<size_t>(-1)){}
  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &strPtrVec) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  [[nodiscard]] size_t GetLlcs() const;
  void SetLlcs(Unsigned llcs);
 private:
  Unsigned llcs_; ///< LLCS <= k. There are at most k-1 gaps

 public:
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_LLCS_Known, llcs_)
};


}// end of namespace
#endif //LCS_SOLVER_CONSTRAINTS_GLOBAL_CONSTRAINT_LLCS_KNOWN_H_
