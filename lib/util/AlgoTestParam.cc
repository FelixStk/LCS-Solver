/*******************************************************************************
 * @file AlgoTestParam.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief defines a structure to run a lcs algorithm
 ******************************************************************************/

#include "util/AlgoTestParam.h"
#include <iostream>
#include <sstream>
#include <utility>

namespace lcs_solver::util {

AlgoParam::AlgoParam(std::string s, const StringVec &v, Map m, BaseSolution *p)
    : name(std::move(s)), pointers(convertToPtr(v)), map(std::move(m)), sol(p) {
//  pointers.reserve(v.size());
//  for (auto &str : v) {
//    pointers.push_back(std::make_shared<const String>(str));
//  }
}

AlgoParam::AlgoParam(std::string s, StrPtrVector vec,
                     Map m, BaseSolution *p)
  : name(std::move(s)), pointers(std::move(vec)), map(std::move(m)), sol(p) {
}

AlgoParam::AlgoParam(const AlgoParam &rhs)
    : name(rhs.name), pointers(rhs.pointers), map(rhs.map), sol(nullptr) {
  if (rhs.sol != nullptr) {
    sol =
        rhs.sol->clone();
  }
}

AlgoParam &AlgoParam::operator=(const AlgoParam &rhs) {
  if (this != &rhs) {
    name = rhs.name;
    pointers = rhs.pointers;
    map = rhs.map;

    BaseSolution *newSol = nullptr;
    if (rhs.sol != nullptr) {
      newSol = rhs.sol->clone();
    }
    delete sol; // Delete current sol
    sol = newSol; // Assign the new clone
  }
  return *this;
}

AlgoParam &AlgoParam::operator=(AlgoParam &&rhs) noexcept {
  if (this != &rhs) {
    name = std::move(rhs.name);
    pointers = std::move(rhs.pointers);
    map = std::move(rhs.map);

    delete sol; // Delete the current sol
    sol = rhs.sol; // Transfer ownership
    rhs.sol = nullptr; // Nullify the rhs pointer to prevent deletion
  }
  return *this;
}

AlgoParam::~AlgoParam() {
  if (sol != nullptr) {
    delete sol;
    sol = nullptr;
  }
}
std::string AlgoParam::DebugString() const {
  std::ostringstream oss;
  oss << "name: " << name << std::endl;
  for (size_t i = 0; i < pointers.size(); ++i) {
    oss << "s[" << i << "]: " << util::to_string(*pointers[i]) << std::endl;
  }
  for (const auto &[key, constraint] : map) {
    oss << constraint->DebugString();
  }
  oss << sol->DebugString();
  return oss.str();
}

AlgoParam::StrPtrVector AlgoParam::convertToPtr(const AlgoParam::StringVec &strings) {
  lcs_solver::util::AlgoParam::StrPtrVector vec;
  vec.reserve(strings.size());
  for (auto &str : strings) {
    vec.push_back(std::make_shared<const String>(str));
  }
  return vec;
}

} // end of namespace