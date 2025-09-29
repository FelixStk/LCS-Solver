/*******************************************************************************
 * @file BaseIterator.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief BaseClass for Iterators used in LCS Algorithms
 ******************************************************************************/

#include "algorithms/solutions/BaseCollector.h"
#include <sstream>
#include <utility>

namespace lcs_solver::algorithms::solutions {

BaseIterator::reference BaseIterator::operator*() const {
  return *current;
}

BaseIterator::pointer BaseIterator::operator->() const {
  return current;
}

BaseIterator &BaseIterator::operator++() {
  advance();
  if (current->empty()) {
    is_end = true;
  }
  return *this;
}

bool BaseIterator::operator==(const BaseIterator &other) const {
  bool endFlagEq = (is_end == other.is_end);
  if (endFlagEq && is_end) {
    return true;
  }
  if (endFlagEq && !is_end) {
    /* Compare solutions::Points, because they are unique and the Points are
     * sometimes modified in place. So that a simple address comparison is not
     * possible. This could potentially lead to problems if a container has the
     * same Points twice in it. This is possible by abusing (for example)
     * Vector3DSolution.*/
    return *current == *(other.current);
  }
  return false;
}

bool BaseIterator::operator!=(const BaseIterator &other) const {
  return !(*this == other);
}

std::string BaseIterator::DebugString() const {
  if (is_end) {
    return "BaseIteratorEnd";
  }
  std::ostringstream oss;
  oss << "current:" << current->DebugString() << std::endl;
  oss << "is_end:" << is_end << std::endl;
  return oss.str();
}

}// namespace lcs_solver::algorithms::solutions
