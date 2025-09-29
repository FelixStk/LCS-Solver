#ifndef LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_MC_INC_H_
#define LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_MC_INC_H_

#include "constraints/local/Constraint_MC.h"

namespace lcs_solver::constraints::local {
class Constraint_MC_INC final : public Constraint_MC {
 public:
  static constexpr ConstraintType kType = ConstraintType::MC_INC;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_MC_INC";
  static constexpr const char *kDescription =
      "The constraint is defined by: \n"
      "1.) If gap[i] = (li,ui) then for ever string w and every gap_e(w,j) = w[e[j]+1 :e[j] -1] : li <= |gap_e(w,j)| <= ui.\n"
      "2.) The (t-1)-tuple is increasing: the interval given by gc[i] is contained in that of gc[i+1] for all i=0..k-2\n";
  static constexpr ConstraintType StaticType(){return kType;}

  explicit Constraint_MC_INC(const GapVector &gc);
  Constraint_MC_INC() = default;

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_MC_INC, gap_, t_, gap_lower_bound_, gap_upper_bound_, gap_max_length_, gap_num_uniques_)
};

} // end of namespace
#endif //LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_MC_INC_H_
