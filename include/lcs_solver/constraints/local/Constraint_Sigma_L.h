#ifndef LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_SIGMA_L_H_
#define LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_SIGMA_L_H_

#include "constraints/local/Constraint_Sigma.h"

namespace lcs_solver::constraints::local {
class Constraint_Sigma_L : public Constraint_Sigma {
 public:
  static constexpr ConstraintType kType = ConstraintType::SIGMA_L;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_Sigma_L";
  static constexpr const char *kDescription =
      "The constraint is defined by: Let gap_e(w,j) = w[e[j]+1 :e[j+1]-1], where e is an embedding.\n"
      "For every gap between two embedded symbols. If e[i]=sym and left[sym] = (l,u) then for all w : l <= |gap_e(w,j)| <= u.\n";
  static constexpr ConstraintType StaticType(){return kType;}

  explicit Constraint_Sigma_L(const SigmaTupleMap &m);
  Constraint_Sigma_L(const SymbolVector &alph, const GapVector &left);
  Constraint_Sigma_L() = default;

  static Constraint_Sigma_L CreateRelaxed(const StrPtrVector &spv);

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_Sigma_L, left_)
};

} // end of namespace
#endif //LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_SIGMA_L_H_