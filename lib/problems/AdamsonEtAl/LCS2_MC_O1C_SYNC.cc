/******************************************************************************
 * @file LCS2_MC_O1C_SYNC.cc
 * @author Steinkopp:Felix
 * @brief Implementation of LCS2_MC_O1C_SYNC
 *****************************************************************************/

#include "problems/AdamsonEtAl/LCS2_MC_O1C_SYNC.h"
#include <ranges>
#include "constraints/ConstraintType.h"
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "algorithms/AlgoFactory.h"

namespace lcs_solver::problems::adamson {

LCS2_MC_O1C_SYNC::LCS2_MC_O1C_SYNC(
    std::string_view name,
    std::string_view description,
    const StrPtrVec &&spv,
    const ConstraintMap &&map
) : BaseProblem(ProblemType::LCS_MC_O1C_SYNC, name, description, spv, map) {
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

constexpr ProblemType LCS2_MC_O1C_SYNC::GetType() {
  return ProblemType::LCS_MC_O1C_SYNC;
}

BaseProblem::ConstraintSpan LCS2_MC_O1C_SYNC::GetReqConstraints() const {
  return kRequiredConstraints;
}

BaseProblem::AlgorithmSpan LCS2_MC_O1C_SYNC::GetAlgoSpan() const {
  return kAvailableAlgorithms;
}

}// end of namespace