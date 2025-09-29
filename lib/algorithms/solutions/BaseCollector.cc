/*******************************************************************************
 * @file BaseCollector.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief BaseClass for collecting the longest common subsequences
 ******************************************************************************/

#include "algorithms/solutions/BaseCollector.h"

// #include <iostream>
#include <sstream>

namespace lcs_solver::algorithms::solutions {

//// Public ////////////////////////////////////////////////////////////////////
BaseSolution::SolutionType BaseCollector::getType() const {
  return BaseSolution::SolutionType::Collector;
}

std::string BaseCollector::DebugString() const {
  if (empty()) {
    return "BaseCollector is Empty";
  }
  const auto set = genSet(this);
  std::ostringstream oss;
  int i = 0;
  bool new_line = false;
  for (const auto &val : set) {
    if (new_line)
      oss << "\n";
    oss << "value[" << i++ << "]: " << val << " ";
    new_line = true;
  }
  return oss.str();
}

bool BaseCollector::empty() const {
  return begin() == end();
}

//// Private ///////////////////////////////////////////////////////////////////

std::set<BaseIterator::value_type> BaseCollector::genSet(const BaseCollector *p) {
  std::set<BaseIterator::value_type> set;
  for (const auto &val : *p) {
    set.insert(val);
  }
  // auto iter1 = set.begin();
  // auto iter2 = set.begin();
  // ++iter2;
  // std::cout << "New Set" << std::endl;
  // while (iter2 != set.end()) {
  //   std::cout << "iter1: " << iter1->DebugString() << std::endl;
  //   std::cout << "iter2: " << iter2->DebugString() << std::endl;
  //   std::cout << "iter1<iter2: " << std::boolalpha << iter1->isLessThan(*iter2);
  //   std::cout << "\n" << std::endl;
  //   ++iter1;
  //   ++iter2;
  // }
  return set;
}

bool BaseCollector::isEqual(const BaseSolution &rhs) const {
  if (const auto rhs_ptr = dynamic_cast<const BaseCollector*>(&rhs)) {
    const auto left_set = genSet(this);
    auto right_set = genSet(rhs_ptr);
    if (left_set.size() != right_set.size()) {
      return false;
    }
    auto it = right_set.begin();
    for (const auto &val : left_set) {
      if (!val.isEqual(*it)) {
        return false;
      }
      ++it;
    }
    return true;
  }
  return false;
}

bool BaseCollector::isLessThan(const BaseSolution &rhs) const {
    if (const auto rhs_ptr = dynamic_cast<const BaseCollector *>(&rhs)) {
      const auto left_set = genSet(this);
      const auto right_set = genSet(rhs_ptr);
      return left_set.size() < right_set.size();
    }
  return getType() < rhs.getType();
}

bool BaseCollector::isLessEqualThan(const BaseSolution &rhs) const {
  if (const auto rhs_ptr = dynamic_cast<const BaseCollector *>(&rhs)) {
    const auto left_set = genSet(this);
    auto right_set = genSet(rhs_ptr);
    if (left_set.size() < right_set.size()) {
      return true;
    }
    if (left_set.size() > right_set.size()) {
      return false;
    }
    auto it = right_set.begin();
    for (const auto &val : left_set) {
      if (val.isEqual(*it)) {
        return false;
      }
      ++it;
    }
    return true;
  }
  return getType() < rhs.getType();
}

}// namespace lcs_solver::algorithms::solutions
