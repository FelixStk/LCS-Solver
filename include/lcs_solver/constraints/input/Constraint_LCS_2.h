#ifndef LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_LCS_2_H_
#define LCS_SOLVER_CONSTRAINTS_INPUT_CONSTRAINT_LCS_2_H_

#include "constraints/input/Constraint_LCS_K.h"

namespace lcs_solver::constraints::input {

class Constraint_LCS_2 : public Constraint_LCS_K {
 public:
  static constexpr ConstraintType kType = ConstraintType::STRINGS_2;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::INPUT;
  inline static const std::string kName = "LCS_2"; ///< Name of the Constraint
  static constexpr const char *kDescription =
      "The LCS Problem has 2 input strings";
  static constexpr ConstraintType StaticType(){return kType;}

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  Constraint_LCS_2();
};

} // end of namespace
#endif //LCS_SOLVER_INCLUDE_CONSTRAINTS_INPUT_CONSTRAINT_LCS_2_H_
