#ifndef LCS_SOLVER_INCLUDE_CONSTRAINTS_GLOBAL_CONSTRAINT_BR_H_
#define LCS_SOLVER_INCLUDE_CONSTRAINTS_GLOBAL_CONSTRAINT_BR_H_

#include "constraints/BaseConstraint.h"
#include "util/CommonTypes.h"
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints::global {

class Constraint_BR : public BaseConstraint {
 public:
  static constexpr ConstraintType kType = ConstraintType::BR;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::GLOBAL;
  static constexpr const char *kName = "Constraint BR";
  static constexpr const char *kDescription =
      "The LCS is contained in a factor of length b";
  static constexpr ConstraintType StaticType(){return kType;}

  Constraint_BR():Constraint_BR(static_cast<size_t>(-1)){}
  explicit Constraint_BR(Unsigned b);
  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &strPtrVec) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  [[nodiscard]] Unsigned GetB() const;
  void SetB(Unsigned b);

 private:
  Unsigned b_;
 public:
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_BR, b_)
};


}// end of namespace
#endif //LCS_SOLVER_INCLUDE_CONSTRAINTS_GLOBAL_CONSTRAINT_BR_H_
