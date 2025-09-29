/*******************************************************************************
 * @file EmptySolution.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Impl. of a EmbeddingsSolution for LCS Solutions.
 ******************************************************************************/

#include "algorithms/solutions/EmptySolution.h"

namespace lcs_solver::algorithms::solutions {

EmptySolution::EmptySolution() = default;

BaseSolution::SolutionType EmptySolution::getType() const {
  return BaseSolution::SolutionType::Empty;
}

bool EmptySolution::isEqual(const BaseSolution &rhs) const {
  return rhs.getType() == SolutionType::Empty;
}

bool EmptySolution::isLessThan(const BaseSolution &rhs) const {
  return rhs.getType() < SolutionType::Empty;
}

bool EmptySolution::isLessEqualThan(const BaseSolution &rhs) const {
  return rhs.getType() <= SolutionType::Empty;
}

std::string EmptySolution::DebugString() const {
  return "Empty Solution";
}
BaseSolution *EmptySolution::clone() const {
  return new EmptySolution();
}
bool EmptySolution::empty() const {
  return true;
}

}// end of namespace