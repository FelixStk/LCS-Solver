/******************************************************************************
 * @file LCS_Classic.cpp
 * @author Steinkopp:Felix
 * @brief Def of the classic lcs problem with two strings (short LCS2_Classic)
 *****************************************************************************/

#include "problems/LCS_Classic.h"

#include "algorithms/AlgoFactory.h"
#include "constraints/ConstraintType.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_LCS_K.h"
#include "algorithms/LLCS/LLCS_STD_FL.h"
// #include "algorithms/LCS/LCS_STD_S.h"
#include "algorithms/LCS/LCS2_STD_S.h"

namespace lcs_solver::problems {

LCS_Classic::LCS_Classic(
    const std::span<const util::String> s,
    const algorithms::AlgoType at,
    const bool calc_llcs)
    : BaseProblem(ProblemType::LCS_Classic, kName, kWhat, {}, {}) {
  spv_.reserve(s.size());
  for (const auto &str : s) {
    spv_.push_back(std::make_shared<util::String>(str));
  }
  if (spv_.size() == 2) {
    using constraints::input::Constraint_LCS_2;
    const ConstraintPtr p = std::make_shared<Constraint_LCS_2>();
    AddConstraint(p);
  }
  else {
    using constraints::input::Constraint_LCS_K;
    const ConstraintPtr p = std::make_shared<Constraint_LCS_K>(spv_.size());
    AddConstraint(p);
  }
  using constraints::input::Constraint_Const_Sigma;
  const ConstraintPtr p = std::make_shared<Constraint_Const_Sigma>(spv_);
  AddConstraint(p);

  using algorithms::AlgoFactory;
  AlgorithmPtr algorithm_ptr;
  if (calc_llcs) {
    algorithm_ptr = AlgoFactory::Create(at, spv_, map_);
  }
  else {
    if (spv_.size() != 2)
      throw std::runtime_error("LCS_Classic: only 2 strings allowed");
    using algorithms::lcs::LCS2_STD_S;
    algorithm_ptr = std::make_unique<LCS2_STD_S>(spv_, map_);
  }
  this->SetAlgorithm(std::move(algorithm_ptr));
}

LCS_Classic::LCS_Classic(
    const std::string_view name,
    const std::string_view description,
    StrPtrVec &&spv,
    ConstraintMap &&map)
    : BaseProblem(ProblemType::LCS_Classic, name, description, spv, map) {

  if(!map_.contains(ConstraintType::STRINGS_K)) {
    using constraints::input::Constraint_LCS_K;
    const ConstraintPtr p = std::make_shared<Constraint_LCS_K>(spv.size());
    AddConstraint(p);
  }
  if(!map_.contains(ConstraintType::CONST_SIG)) {
    using constraints::input::Constraint_Const_Sigma;
    const ConstraintPtr p = std::make_shared<Constraint_Const_Sigma>(spv);
    AddConstraint(p);
  }

}

constexpr ProblemType LCS_Classic::GetType() {
  return ProblemType::LCS_Classic;
}

BaseProblem::ConstraintSpan LCS_Classic::GetReqConstraints() const {
  return kRequiredConstraints;
}

BaseProblem::AlgorithmSpan LCS_Classic::GetAlgoSpan() const {
  return kAvailableAlgorithms;
}

}// end of namespace