/******************************************************************************
 * @file LLCS_Algorithm.cc
 * @author Steinkopp:Felix
 * @version 1.2
 * @brief Specialization of BaseAlgorithm for LLCS Problems
 *****************************************************************************/
#include "algorithms/LLCS/LLCS_Algorithm.h"

#include "algorithms/AlgoCategory.h"
#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"

namespace lcs_solver::algorithms::llcs {

/*******************************************************************************
 * LLCS_Algorithm - Constructor
 * @param algo AlgoType identifies the algorithm uniquely
 * @param vec StrPtrVector of shared points to the constant strings of a problem
 * @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
 ******************************************************************************/
LLCS_Algorithm::LLCS_Algorithm(
    AlgoType algo,
    const BaseAlgorithm::StrPtrVector &vec,
    const ConstraintMap &map
) : BaseAlgorithm(AlgoCategory::LLCS, algo, vec, map) {
  SortStringViews();
}

} // end of namespace