/******************************************************************************
 * @file LCS2_MC.cc
 * @author Steinkopp:Felix
 * @brief Implementation of LCS2_MC
 *****************************************************************************/
#include "problems/AdamsonEtAl/LCS2_MC.h"
#include <ranges>
#include "algorithms/AlgoFactory.h"
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "constraints/ConstraintType.h"

namespace lcs_solver::problems::adamson {

LCS2_MC::LCS2_MC(
    std::string_view name,
    std::string_view description,
    const StrPtrVec &&spv,
    const ConstraintMap &&map)
    : BaseProblem(ProblemType::LCS_MC, name, description, spv, map) {
  if(!map_.contains(ConstraintType::STRINGS_2)){
    using ::lcs_solver::constraints::input::Constraint_LCS_2;
    ConstraintPtr p = std::make_shared<Constraint_LCS_2>();
    AddConstraint(p);
  }
  if(!map_.contains(ConstraintType::CONST_SIG)) {
    using ::lcs_solver::constraints::input::Constraint_Const_Sigma;
    ConstraintPtr p = std::make_shared<Constraint_Const_Sigma>(spv);
    AddConstraint(p);
  }
}

constexpr ProblemType LCS2_MC::getType() {
  return ProblemType::LCS_MC;
}

BaseProblem::ConstraintSpan LCS2_MC::GetReqConstraints() const {
  return kRequiredConstraints;
}

BaseProblem::AlgorithmSpan LCS2_MC::GetAlgoSpan() const {
  return kAvailableAlgorithms;
}

}// end of namespace