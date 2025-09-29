#ifndef LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_1C_H_
#define LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_1C_H_

#include "constraints/local/Constraint_MC.h"

namespace lcs_solver::constraints::local {
class Constraint_MC_1C : public Constraint_MC {
 public:
  static constexpr ConstraintType kType = ConstraintType::MC_1C;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_1C";
  static constexpr const char *kDescription =
      "The constraint is defined by: \n"
      "1.) If gap[i] = (li,ui) then for ever string w and every gap_e(w,j) = w[e[j]+1 :e[j+1] -1] : li <= |gap_e(w,j)| <= ui. \"\n"
      "2.) For all i,j in 0..t gap[i]==gap[j] (element-wise)\"";
  static constexpr ConstraintType StaticType() { return kType; }

  explicit Constraint_MC_1C(const GapVector &vec);
  Constraint_MC_1C(const StrPtrVector &, uint lower, uint upper);
  Constraint_MC_1C(uint gap_vector_length, uint lower, uint upper);
  Constraint_MC_1C() = default;

  static GapVector CreateGapVector1C(uint size, uint l, uint u);

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  [[nodiscard]] uint GetLower() const;
  [[nodiscard]] uint GetUpper() const;

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_MC_1C, Constraint_MC::gap_, t_, gap_lower_bound_, gap_upper_bound_, gap_max_length_, gap_num_uniques_)

 private:
  static uint GetMinStrLen(const StrPtrVector &spv);
};

}// namespace lcs_solver::constraints::local
#endif//LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_1C_H_
