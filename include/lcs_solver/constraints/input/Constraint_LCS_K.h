#ifndef LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_LCS_K_H_
#define LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_LCS_K_H_

#include "constraints/BaseConstraint.h"
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints::input {

class Constraint_LCS_K : public BaseConstraint {
 public:
  static constexpr ConstraintType kType = ConstraintType::STRINGS_K;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::INPUT;
  static constexpr const char *kName = "LCS_K";
  static constexpr const char *kDescription =
      "The LCS Problem has K input strings";
  static constexpr ConstraintType StaticType(){return kType;}

  explicit Constraint_LCS_K(Unsigned k);
  Constraint_LCS_K() : Constraint_LCS_K(static_cast<size_t>(-1)){}
  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &strPtrVec) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  [[nodiscard]] Unsigned GetK() const;
  void SetK(Unsigned k);

 protected:
  explicit Constraint_LCS_K(ConstraintType t, Unsigned k);
 private:
  Unsigned k_{};  ///< the length of input sequences
public:
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_LCS_K, k_)
};

} // end of namespace
#endif //LCS_SOLVER_INCLUDE_CONSTRAINTS_INPUT_CONSTRAINT_LCS_K_H_
