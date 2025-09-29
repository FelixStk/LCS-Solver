/******************************************************************************
 * @file LCS2_MC_O1C.cc
 * @author Steinkopp:Felix
 * @brief Implementation of LCS2_MC_O1C
 *****************************************************************************/

#include "problems/AdamsonEtAl/LCS2_MC_O1C.h"
#include <ranges>
#include "constraints/ConstraintType.h"
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "algorithms/AlgoFactory.h"

namespace lcs_solver::problems::adamson {

LCS2_MC_O1C::LCS2_MC_O1C(
    const std::string_view name,
    const std::string_view description,
    const StrPtrVec &&spv,
    const ConstraintMap &&map)
    : BaseProblem(name, ProblemType::LCS_MC, description, spv, map) {
  if(!map_.contains(ConstraintType::STRINGS_2)){
    using ::lcs_solver::constraints::input::Constraint_LCS_2;
    ConstraintPtr p = std::make_shared<Constraint_LCS_2>();
    addConstraint(p);
  }
  if(!map_.contains(ConstraintType::CONST_SIG)) {
    using ::lcs_solver::constraints::input::Constraint_Const_Sigma;
    ConstraintPtr p = std::make_shared<Constraint_Const_Sigma>(spv);
    addConstraint(p);
  }
}

constexpr ProblemType LCS2_MC_O1C::getType() {
  return ProblemType::LCS_MC_O1C;
}

const BaseProblem::ConstraintSpan LCS2_MC_O1C::getReqConstraints() const {
  return kRequiredConstraints;
}

const BaseProblem::AlgorithmSpan LCS2_MC_O1C::getAlgoSpan() const {
  return kAvailableAlgorithms;
}

}// end of namespace