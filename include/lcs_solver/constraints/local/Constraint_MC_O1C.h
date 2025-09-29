#ifndef LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_O1C_H_
#define LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_O1C_H_

#include "constraints/local/Constraint_MC.h"

namespace lcs_solver::constraints::local {
class Constraint_MC_O1C : public Constraint_MC {
 public:
  static constexpr ConstraintType kType = ConstraintType::MC_O1C;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_MC_O1C";
  static constexpr const char *kDescription =
      "The constraint is defined by:\n"
      "1.) If gap[i] = (li,ui) then for ever string w and every gap_e(w,j) = w[e[j]+1 :e[j] -1] : li <= |gap_e(w,j)| <= ui.\n"
      "2.) The t-tuple for the gap constraints is decomposable into O(1) equivalence classes via gap[i]==gap[j] (elementwise)\n";
  static constexpr ConstraintType StaticType(){return kType;}

  explicit Constraint_MC_O1C(const GapVector &gc);
  Constraint_MC_O1C() = default;

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_MC_O1C, Constraint_MC::gap_, t_, gap_lower_bound_, gap_upper_bound_, gap_max_length_, gap_num_uniques_)

 protected:
  Constraint_MC_O1C(const GapVector &gc, ConstraintType type);
};

} // end of namespace
#endif //LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_O1C_H_