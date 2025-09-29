/******************************************************************************
 * @file BaseAlgorithm.cc
 * @author Steinkopp:Felix
 * @version 1.3
 * @brief Pure virtual classes for algorithms
 *****************************************************************************/
#include "algorithms/BaseAlgorithm.h"

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <ranges>
#include <string>
#include <string_view>
//#include <utility>
#include <vector>

#include "algorithms/AlgoCategory.h"
#include "algorithms/AlgoFactory.h"
#include "algorithms/AlgoType.h"
#include "util/CommonTypes.h"
#include "util/Logger.hpp"
#include "util/StrPtrVector.h"

#include <sstream>

namespace lcs_solver::algorithms {
//=== Public Functions ==========================================================
/*******************************************************************************
 * getType
 * @return type of the BaseAlgorithm instance
 ******************************************************************************/
AlgoType BaseAlgorithm::getType() const {
  return this->type;
}

/*******************************************************************************
 * getCategory
 * @return type of the BaseAlgorithm instance
 ******************************************************************************/
AlgoCategory BaseAlgorithm::getCategory() const {
  return category;
}

/*******************************************************************************
 * getConstraintMap
 * @return Ref. to the ConstraintMap of the problem that the algorithm solves
 ******************************************************************************/
const BaseAlgorithm::ConstraintMap &BaseAlgorithm::getConstraintMap() const {
  return constraintMap;
}

/*******************************************************************************
 * isFinished
 * @return whether state == State::Processed
 ******************************************************************************/
bool BaseAlgorithm::isFinished() const {
  return state == State::Processed;
}

/*******************************************************************************
 * getState
 * @return state of algorithm
 ******************************************************************************/
BaseAlgorithm::State BaseAlgorithm::getState() const {
  return state;
}

/*******************************************************************************
 * getStringViewVec
 * @return std::vector with StringViews of the problems Strings
 ******************************************************************************/
BaseAlgorithm::StringViewVector BaseAlgorithm::getStringViewVec() const {
  return s;
}

/*******************************************************************************
 * getStrPtrVec
 * @return std::vector with std::share_ptr to the problems Strings
 ******************************************************************************/
const BaseAlgorithm::StrPtrVector& BaseAlgorithm::getStrPtrVec() const {
  return strPtrVec;
}

/*******************************************************************************
 * Generates a data-owning and sorted vector of the problem's strings
 * @param increasing Flag whether the result shall be sorted increasingly
 * @return StrPtrVector (aka `std::vector<std::shared_ptr<const util::String>>`)
 ******************************************************************************/
BaseAlgorithm::StrPtrVector BaseAlgorithm::genSortedStrPtrVec(const bool increasing) const {
  StrPtrVector spv;
  spv.reserve(strPtrVec.size());
  for (const auto& ptr : strPtrVec) {
    spv.push_back(ptr);
  }
  if (increasing) {
    std::ranges::sort(spv, [](const auto &a, const auto &b) {
      return a->size() < b->size();
    });
  }
  else {
    std::ranges::sort(spv, [](const auto &a, const auto &b) {
     return a->size() > b->size();
    });
  }
  return spv;
}

/*******************************************************************************
 * @brief Helper to convert a Vector of StringViews into a std::string
 * @param svv StringViewVector (aka std::vector<util::StringView>)
 * @return std::string
 ******************************************************************************/
std::string BaseAlgorithm::toString(const StringViewVector &svv) {
  std::ostringstream oss;
  uint i = 0;
  for (const auto &string_view : svv) {
    oss << "s[" << i++ << "]: " << util::to_string(string_view);
    if (string_view != svv.back()) {
      oss << "\n";
    }
  }
  return oss.str();
}

/*******************************************************************************
 * @brief Helper to convert a Vector of StrPtr into a std::string
 * @param spv StrPtrVector (aka std::vector<shared_ptr<const String>)
 * @return std::string
 ******************************************************************************/
std::string BaseAlgorithm::toString(const StrPtrVector &spv) {
  std::ostringstream oss;
  uint i = 0;
  auto toString = [](const std::shared_ptr<const util::String> &p) {
    return lcs_solver::util::to_string(p.get());
  };
  for (const auto &s : spv | std::views::transform(toString)) {
    oss << "s[" << i++ << "]: " << s;
    if (i != spv.size()) {
      oss << "\n";
    }
  }
  // for (const auto &string_view : spv) {
  //   oss << "s[" << i++ << "]: " << *string_view;
  //   if (i != spv.size()) {
  //     oss << "\n";
  //   }
  // }
  return oss.str();
}

/*******************************************************************************
 * getName
 * @details is to overwritten statically to give a name to an algorithm
 * @return std::string_view, the name of the algorithm
 ******************************************************************************/
std::string_view BaseAlgorithm::getName() const {
  return {"Empty Algorithm"};
}

/*******************************************************************************
 * getDescription
 * @return std::string_view, the description of the algorithm
 ******************************************************************************/
std::string_view BaseAlgorithm::getDescription() const {
  return {"Empty Description"};
}

/*******************************************************************************
 * reset - does a full reset of the algorithm
 ******************************************************************************/
void BaseAlgorithm::reset() {
  reset(ResetLevel::Full);
}

/*******************************************************************************
 * operator>> reads a AlgoType and changes the algorithm
 * @param in std::istream used for reading a AlgoType of p
 * @param ptr
 * @return std::istream written to
 ******************************************************************************/
std::istream &operator>>(std::istream &in, BaseAlgorithm::AlgorithmPtr &ptr) {
  std::string input;
  in >> input;
  AlgoType t;
  if (in.rdbuf() == std::cin.rdbuf()) {
    t = AlgoFactory::ReadAlgoType();
  }
  else {
    std::ostringstream oss;
    t = AlgoFactory::ReadAlgoType(in,oss);
  }

  const auto spv = ptr->getStrPtrVec();
  const auto map = ptr->getConstraintMap();
  ptr = AlgoFactory::Create(t, spv, map);
  return in;
}

/*******************************************************************************
 *
 * @param os std::ostream used for writing the AlgoType of p
 * @param p AlgorithmPtr to the algorithm
 * @return std::ostream read from
 ******************************************************************************/
std::ostream &operator<<(std::ostream &os, const BaseAlgorithm::AlgorithmPtr &p) {
  os << AlgoFactory::GetName(p->getType());
  return os;
}

//=== Protected Functions ======================================================
/*******************************************************************************
 * BaseAlgorithm - Constructor
 * @param category AlgoCategory of the algorithm
 * @param type AlgoType of the algorithm
 * @param spv std::vector of std::shared_ptr to the problem Strings
 * @param map std::map from ConstraintType to std::shared_ptr of the constraint
 ******************************************************************************/
BaseAlgorithm::BaseAlgorithm(
    const AlgoCategory category,
    const AlgoType type,
    const StrPtrVector &spv,
    const ConstraintMap &map)
    : type(type),
      category(category),
      constraintMap(map),
      strPtrVec(spv),
      s(util::StringViewVecFrom(strPtrVec)),
      state(State::Invalid) {
}

/*******************************************************************************
 * isEachConstraintIndividuallyValid
 * @return whether all constraints in map return true for isConstraintValid(svv)
 ******************************************************************************/
bool BaseAlgorithm::isEachConstraintIndividuallyValid() const {
  //  ::lcs_solver::util::Logger::Error() << *(getStrPtrVec()[0]) << "\n";
  //  ::lcs_solver::util::Logger::Error() << *(getStrPtrVec()[1]) << "\n";
  //  for (const auto &[key, constraint_ptr] : constraintMap) {
  //    const bool check = constraint_ptr->isConstraintValid(strPtrVec);
  //    if (!check) {
  //      Logger::Warning() << "Constraint " << static_cast<unsigned char>(key)
  //                        << "is not valid!\n";
  //      return false;
  //    }
  //  }
  //  return true;
  const bool all_valid = std::ranges::all_of(constraintMap, [this](const auto &pair) {
    const auto &[key, constraint_ptr] = pair;
    const bool check = constraint_ptr->IsConstraintValid(strPtrVec);
    if (!check) {
      util::Logger::Warning() << "Constraint " << static_cast<unsigned char>(key)
                              << " is not valid!\n";
    }
    return check;
  });
  return all_valid;
}

/*******************************************************************************
 * isFilled
 * @return whether each and every strings given contains at least one symbol
 ******************************************************************************/
bool BaseAlgorithm::isFilled() const {
  for (size_t i = 0; i < s.size(); i++) {
    if (s[i].empty()) {
      util::Logger::Warning() << "String" << i << "should not be empty! (in "
                              << getName() << ")" << '\n';
      return false;
    }
  }
  return true;
}

/*******************************************************************************
 * @brief Checks if the StringViewVector s is sorted
 * @return ture iff the svv is sorted by increasing length otherwise false
 ******************************************************************************/
bool BaseAlgorithm::isSorted() const {
  // Check sorting |svv[0]| <= svv[1] <= ... <= svv[k]
  if (s.size() >= 2) {
    for (size_t i = 0; i < s.size() - 1; i++) {
      if (s[i].size() > s[i + 1].size()) {
        util::Logger::Warning() << "Input is not sorted! (s["
                                << i << "] > s[" << i + 1 << "] in " << getName()
                                << ")" << '\n';
        return false;
      }
    }
  }
  return true;
}

/*******************************************************************************
 * Gets the length of the longest String in a vector of pointers
 * @param vec StrPtrVec a vector of shared pointer to strings
 * @return size_t the length of the longest string in vec
 ******************************************************************************/
size_t BaseAlgorithm::getLongestStrLen(const StrPtrVector &vec) {
  size_t maxLength = 0;
  for (const auto &p : vec) {
    maxLength = std::max(p->size(), maxLength);
  }
  return maxLength;
}

/*******************************************************************************
 * setState
 * @param update new state of the algorithm
 ******************************************************************************/
void BaseAlgorithm::setState(BaseAlgorithm::State update) {
  state = update;
}

void BaseAlgorithm::SortStringViews(const bool increasing) {
  if (increasing) {
    // std::ranges::sort(s);
    std::ranges::sort(s, [](const auto &a, const auto &b) {
      return a.size() < b.size();
    });
  } else {
    // std::ranges::sort(s, std::ranges::greater());
    std::ranges::sort(s, [](const auto &a, const auto &b) {
      return a.size() > b.size();
    });
  }

}

}  // namespace lcs_solver::algorithms