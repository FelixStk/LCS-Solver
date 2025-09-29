/******************************************************************************
 * @file LCS2_SIGMA_L.cc
 * @author Steinkopp:Felix
 * @brief Implementation of LCS2_SIGMA_L
 *****************************************************************************/

#include "problems/AdamsonEtAl/LCS2_SIGMA_L.h"
#include <ranges>
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "algorithms/AlgoFactory.h"

namespace lcs_solver::problems::adamson {

LCS2_Sigma_L::LCS2_Sigma_L(
    const std::string_view name,
    const std::string_view description,
    const StrPtrVec &&spv,
    const ConstraintMap &&map
) : BaseProblem(ProblemType::LCS_Sigma_L, name, description, spv, map) {
  if(!map.contains(ConstraintType::STRINGS_2)){
    using constraints::input::Constraint_LCS_2;
    ConstraintPtr p = std::make_shared<Constraint_LCS_2>();
    AddConstraint(p);
  }
  if(!map.contains(ConstraintType::CONST_SIG)) {
    using constraints::input::Constraint_Const_Sigma;
    ConstraintPtr p = std::make_shared<Constraint_Const_Sigma>(spv);
    AddConstraint(p);
  }
}

constexpr ProblemType LCS2_Sigma_L::GetType() {
  return ProblemType::LCS_Sigma_L;
}

BaseProblem::ConstraintSpan LCS2_Sigma_L::GetReqConstraints() const {
  return kRequiredConstraints;
}

BaseProblem::AlgorithmSpan LCS2_Sigma_L::GetAlgoSpan() const {
  return kAvailableAlgorithms;
}

}// end of namespace