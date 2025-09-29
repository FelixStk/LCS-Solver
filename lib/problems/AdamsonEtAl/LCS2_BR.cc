/******************************************************************************
 * @file LCS2_BR.cc
 * @author Steinkopp:Felix
 * @brief Implementation of LCS2_BR
 *****************************************************************************/

#include "problems/AdamsonEtAl/LCS2_BR.h"
#include "constraints/input/Constraint_LCS_2.h"
#include <utility>
#include <ranges>

namespace lcs_solver::problems::adamson {

LCS2_BR::LCS2_BR(
    std::string_view name,
    std::string_view description,
    const StrPtrVec &&spv,
    const ConstraintMap &&map,
    AlgorithmType algorithm_type)
    : BaseProblem(name, ProblemType::BR, description, spv, map, algorithm_type) {
  if(!map.contains(ConstraintType::STRINGS_2)){
    using ::lcs_solver::constraints::input::Constraint_LCS_2;
    ConstraintPtr p = std::make_shared<Constraint_LCS_2>();
    addConstraint(p);
  }
}

const BaseProblem::ConstraintSpan &LCS2_BR::getReqConstraints() const {
  return requiredConstraints;
}

const BaseProblem::AlgorithmSpan &LCS2_BR::getAlgoSpan() const {
  return availableAlgorithms;
}

}// end of namespace