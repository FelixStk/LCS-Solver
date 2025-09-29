#ifndef LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_O1C_SYNC_H_
#define LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_O1C_SYNC_H_

#include "constraints/local/Constraint_MC_O1C.h"

namespace lcs_solver::constraints::local {
class Constraint_MC_O1C_SYNC final : public Constraint_MC_O1C {
 public:
  static constexpr ConstraintType kType = ConstraintType::MC_O1C_SYNC;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_MC_O1C_SYNC";
  static constexpr const char *kDescription =
      "The constraint is defined by:\n"
      "1.) If gap[i] = (li,ui) then for ever string w and every gap_e(w,j) = w[e[j]+1 :e[j] -1] : li <= |gap_e(w,j)| <= ui.\n"
      "2.) The t-tuple for the gap constraints is decomposable into O(1) equivalence classes via gap[i]==gap[j] (elementwise)\n"
      "3.) the t-tuple is synchronized: For i,j in [0..t], if gap[i]==gap[j] and i <= j then gap[i+e] is a subset of gap[j+e]\n"
      " for all e >= 0 such that: i+e <= j+e <= t";
  static constexpr ConstraintType StaticType(){return kType;}

  explicit Constraint_MC_O1C_SYNC(const GapVector &gc);
  Constraint_MC_O1C_SYNC() = default;

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_MC_O1C_SYNC, Constraint_MC::gap_, t_, gap_lower_bound_, gap_upper_bound_, gap_max_length_, gap_num_uniques_)

};

} // end of namespace
#endif //LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_O1C_SYNC_H_