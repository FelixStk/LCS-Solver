#ifndef LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_CONST_SIGMA_H_
#define LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_CONST_SIGMA_H_

#include "constraints/BaseConstraint.h"
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints::input {

class Constraint_Const_Sigma : public BaseConstraint {
 public:
  static constexpr ConstraintType kType = ConstraintType::CONST_SIG;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::INPUT;
  static constexpr const char *kName = "LCS_Const_Sigma";
  static constexpr const char *kDescription =
      "The alphabet of the strings in the problem has size sigma";
  static constexpr ConstraintType StaticType() { return kType; }

  Constraint_Const_Sigma() : Constraint_Const_Sigma(static_cast<size_t>(-1)) {}
  explicit Constraint_Const_Sigma(const StrPtrVector &spv);
  explicit Constraint_Const_Sigma(size_t sigma);
  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  [[nodiscard]] size_t GetSigmaSize() const;
  void SetSigmaSize(size_t sigma);

 private:
  static size_t CountDifferentSymbolsInVector(const StrPtrVector &spv);
  size_t sigma_;///< Number of different symbols in the strings
 public:
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_Const_Sigma, sigma_)
};

}// namespace lcs_solver::constraints::input
#endif//LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_CONST_SIGMA_H_
