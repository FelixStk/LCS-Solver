/******************************************************************************
 * @file LCS2_SIGMA.cc
 * @author Steinkopp:Felix
 * @brief Implementation of LCS2_SIGMA
 *****************************************************************************/

#include "problems/AdamsonEtAl/LCS2_SIGMA.h"
#include <ranges>
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "algorithms/AlgoFactory.h"

namespace lcs_solver::problems::adamson {

LCS2_Sigma::LCS2_Sigma(
    const std::string_view name,
    const std::string_view description,
    const StrPtrVec &&spv,
    const ConstraintMap &&map
) : BaseProblem(ProblemType::LCS_Sigma, name, description, spv, map) {
  if(!map.contains(ConstraintType::STRINGS_2)){
    using ::lcs_solver::constraints::input::Constraint_LCS_2;
    ConstraintPtr p = std::make_shared<Constraint_LCS_2>();
    AddConstraint(p);
  }
  if(!map.contains(ConstraintType::CONST_SIG)) {
    using ::lcs_solver::constraints::input::Constraint_Const_Sigma;
    ConstraintPtr p = std::make_shared<Constraint_Const_Sigma>(spv);
    AddConstraint(p);
  }
}

constexpr ProblemType LCS2_Sigma::GetType() {
  return ProblemType::LCS_Sigma;
}

BaseProblem::ConstraintSpan LCS2_Sigma::GetReqConstraints() const {
  return kRequiredConstraints;
}

BaseProblem::AlgorithmSpan LCS2_Sigma::GetAlgoSpan() const {
  return kAvailableAlgorithms;
}

}// end of namespace