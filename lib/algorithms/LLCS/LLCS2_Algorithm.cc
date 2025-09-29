/******************************************************************************
 * @file LLCS2_Algorithm.cc
 * @author Steinkopp:Felix
 * @version 1
 * @brief Specialization of BaseAlgorithm for LLCS2 Problems
 * @details LLCS2_Algorithm sets the type of indices used in 2D matrices. In
 *  addition, it implements some useful functions for working with them.
 *  It also provides an interface for requesting read-only access to the summary
 *  matrix and for tracking keyPoints.
 *****************************************************************************/
#include "algorithms/LLCS/LLCS2_Algorithm.h"

#include "algorithms/AlgoCategory.h"
#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include <cassert>
#include <sstream>
#include <unordered_map>

namespace lcs_solver::algorithms::llcs {

/**
 * LLCS2_Algorithm - Constructor
 * @param algo AlgoType identifies the algorithm uniquely
 * @param vec StrPtrVector of shared points to the constant strings of a problem
 * @param map std::map<std::string, share_ptr<BaseConstraint> for constraints
 * @param doTracking true, iff the LLCS2_Algorithm algorithm is supposed to do
 *  bookkeeping of points that are extendable: If Pair x is in keyPoints[l],
 *  then there exists a shared subsequence of length l that ends at x (meaning
 *  it is a subsequence of s1[0:x.first-1] and s2[0:x.second-1]
*  @param oneBased Whether the indices used in the matrices of the algorithm are
 *  oneBased (in [1:s1.size()]x[1:s2.size()] or zeroBased
 */
LLCS2_Algorithm::LLCS2_Algorithm(
    const AlgoType algo,
    const BaseAlgorithm::StrPtrVector &vec,
    const ConstraintMap &map,
    const bool doTracking,
    const bool oneBased)
    : BaseAlgorithm(AlgoCategory::LLCS, algo, vec, map),
      areIndicesOneBased(oneBased),
      trackKeyPairs(doTracking) {
  if (s.size() == 2) {
    keyPairs.resize(std::min(s[0].size(), s[1].size()) + 1);
  }
  SortStringViews();
}

/**
 * isMatched
 * @param pair Index (x,y) where to check the strings for a match
 * @param oneBased Whether indices are oneBased or zeroBased
 * @return s1[x] == s1[y]
 * @note s is always zero Based. But some algorithms use oneBased indices
 */
bool LLCS2_Algorithm::isMatched(const Pair &pair, const bool oneBased) const {
  if (s.size() != 2)
    return false;
  if (oneBased) {
    assert(pair.first != 0 && pair.first <= s[0].size() && "Valid Index");
    assert(pair.second != 0 && pair.second <= s[1].size() && "Valid Index");
    return s[0].at(pair.first - 1) == s[1].at(pair.second - 1);
  }
  assert(pair.first < s[0].size() && "Valid Index");
  assert(pair.second < s[1].size() && "Valid Index");
  return s[pair.first] == s[pair.second];
}

/**
 * isMatched
 * @param pair Index (x,y) where to check the strings for a match
 * @return s1[x] == s1[y]
 * @see LLCS2_Algorithm::isMatched(const Pair &pair, const bool oneBased)
 */
bool LLCS2_Algorithm::isMatched(const Pair &pair) const {
  return isMatched(pair, areIndicesOneBased);
}

/**
 * genAllKeyPairsNaive
 * @param oneBased Whether indices are oneBased or zeroBased
 * @return std::vector<LLCS2_Algorithm::Pair> of Pairs p such that
 * @note time complexity is always O(l1 x l2)
 */
std::vector<LLCS2_Algorithm::Pair> LLCS2_Algorithm::genAllKeyPairsNaive(const bool oneBased) const {
  if (s.size() != 2) {
    return {};
  }
  uint l1 = s.size() < 2 ? 0 : s[0].size();
  uint l2 = s.size() < 2 ? 0 : s[1].size();
  std::vector<Pair> keyPoints;
  if (l1 == 0 || l2 == 0) return keyPoints;
  // Normal approach: The complexity is always O(l1 x l2)
  Pair p = {0, 0};
  if (oneBased) {
    p = {1, 1};
    l1++;
    l2++;
  }
  for (uint &i = p.first; i < l1; i++) {
    for (uint j = p.second; j < l2; j++) {
      if (isMatched(p, oneBased)) {
        keyPoints.emplace_back(i, j);
      }
    }
  }
  return keyPoints;
}
/**
 * lenKPtM
 * @param kpm KeyPointMatrix p in kpm[i] <=> isMatched(, true) and llcs(p) == i
 * @return uint the llcs based on the kpm object
 */
BaseAlgorithm::uint LLCS2_Algorithm::lenKPtM(const KeyPointMatrix &kpm) {
  if (kpm.empty()) return 0;
  uint length = 0;
  for (uint i = kpm.size() - 1; i > 0; i--) {
    if (kpm[i].empty()) {
      length = i - 1;
    } else {
      length = std::max(length, i);
    }
  }
  return length;
}

BaseAlgorithm::uint LLCS2_Algorithm::max(const Matrix &matrix) {
  uint max = 0;
  for (const auto &row : matrix) {
    for (const auto &elem : row) {
      max = std::max(max, elem);
    }
  }
  return max;
}

/**
 * genAllKeyPairs
 * @param oneBased Whether indices are oneBased or zeroBased
 * @return std::vector<LLCS2_Algorithm::Pair> of Pairs p such that
 * @note This function optimizes for sparsity. The average case complexity is:
 *  O(l1* l2* k1 * Sum_s(|indices1_s| * |indices2_s|)) where li are the lengths
 *  of the strings, k1 is the number of different symbols in the first string
 *  and indices_s are the positions of a symbol s in either the first or second
 *  string.
 */
std::vector<LLCS2_Algorithm::Pair> LLCS2_Algorithm::genAllKeyPairs(const bool oneBased) const {
  if (s.size() != 2) {
    return {};
  }
  std::unordered_map<util::Symbol, std::vector<uint>, util::SymbolPerfectHash, util::SymbolEqual> map1, map2;
  const uint l1 = s.size() < 2 ? 0 : s[0].size();
  const uint l2 = s.size() < 2 ? 0 : s[1].size();
  std::vector<Pair> keyPoints;
  if (l1 == 0 || l2 == 0) return keyPoints;

  for (uint i = 0; i < l1; ++i)
    map1[s[0][i]].push_back(i);
  for (uint j = 0; j < l2; ++j)
    map2[s[1][j]].push_back(j);
  for (const auto &[symbol, indices1] : map1) {
    if (auto it = map2.find(symbol); it != map2.end()) {
      const auto &indices2 = it->second;
      for (uint i : indices1) {
        for (uint j : indices2) {
          if (oneBased) {
            keyPoints.emplace_back(i + 1, j + 1);
          } else {
            keyPoints.emplace_back(i, j);
          }
        }//j
      }//i
    }//if
  }//each symbol
  return keyPoints;
}

/**
 * setTracking is a setter for trackKeyPairs in LLCS2_Algorithm
 * @param tracking new value of trackKeyPairs
 */
void LLCS2_Algorithm::setTracking(const bool tracking) {
  trackKeyPairs = tracking;
}

/**
 * doesTracking is a getter for trackKeyPairs in LLCS2_Algorithm
 * @return bool trackKeyPairs
 */
bool LLCS2_Algorithm::getTrackingFlag() const {
  return trackKeyPairs;
}

/**
 * track is a helper to do bookkeeping of important points
 * @details An important point p = (i,j) is defined by two properties:
 *  - isMatched(p, oneBased) must be true
 *  - there is a common subsequence of s[0] and s[1] with length llcs
 * @param pair Index Pair to be added
 * @param llcs uint to specify the associated length of the Pair
 * @note stored Pairs in trackKeyPairs shall be oneBased - always
 */
void LLCS2_Algorithm::track(const Pair &pair, const uint llcs) {
  assert(llcs < keyPairs.size() && "Valid keyPairs Index ");
  keyPairs[llcs].emplace_back(pair);
}

/**
 * track is a helper to do bookkeeping of important points
 * @param pair Index Pair to be added
 * @param llcs uint to specify the associated length of the Pair
 * @see LLCS2_Algorithm::track(const Pair& a, const uint llcs)
 */
void LLCS2_Algorithm::track(const Pair &&pair, const uint llcs) {
  assert(llcs < keyPairs.size() && "Valid keyPairs Index ");
  keyPairs[llcs].emplace_back(pair);
}

/**
 * toString
 * @param m Matrix as std::vector<std::vector<uint>
 * @param s StringViewVector
 * @param trim Whether to skip elements m[i][j] where either i or j is zero
 * @param name String to give the matrix a name (optional)
 * @param printStrings Whether to add symbols to the rows and columns
 * @return std::string to describe m
 */
std::string LLCS2_Algorithm::toString(
    const Matrix &m,
    StringViewVector s,
    const bool trim,
    const std::string &name,
    const bool printStrings) {
  std::basic_ostringstream<char> oss;
  if (m.empty()) {
    return name + ": empty";
  }
  else {
    oss << name + ":\n";
  }

  const bool addSymbols = (s.size() == 2 && printStrings);
  if (addSymbols) {
    if (!trim) {
      oss << "  ";// empty space for first matrix column
    }
    oss << "  ";// empty space for symbols
    for (const auto &i : s[1]) {
      oss << util::to_string(i) << " ";
    }
  }
  oss << "\n";
  for (uint i = trim ? 1 : 0; i < m.size(); i++) {
    if (addSymbols) {
      if (i == 0) {
        oss << "  ";
      } else {
        if (i <= s[0].size()) {
          oss << util::to_string(s[0].at(i - 1)) << " ";
        } else {
          oss << " " << " ";
        }
      }
    }
    for (uint j = trim ? 1 : 0; j < m[i].size(); j++) {
      oss << m[i][j] << " ";
    }
    if (const bool lastLine = (i == m.size() - 1); m.size() > 1 && !lastLine)
      oss << "\n";
  }
  return oss.str();
}

/**
 * @brief Helper: toString for KeyPointMatrix
 * @param keyPairs KeyPointMatrix (filled by track in LLCS2_Algorithm)
 * @param trackKeyPairs Weather tracking is enabled or not
 * @return std::string to describe the tracking done with keyPairs and trackKeyPairs
 */
std::string LLCS2_Algorithm::toString(const KeyPointMatrix &keyPairs, const bool trackKeyPairs) {
  std::ostringstream oss;
  oss << (trackKeyPairs ? "\ntrackKeyPairs: true" : "trackKeyPairs: false");
  if (trackKeyPairs) {
    bool firstLineOutput = true;
    bool lineOutput = false;
    for (uint i = 0; const auto &pairs : keyPairs) {
      if (!pairs.empty()) {
        if (firstLineOutput || lineOutput) {
          firstLineOutput = false;
          oss << "\n";
        }
        oss << "i== " << ++i << ": ";
        for (const auto &pair : pairs) {
          oss << "(" << pair.first << " " << pair.second << ") ";
        }
        lineOutput = true;
      }
    }
  }
  return oss.str();
}

/**
 * getKeyPairs - Getter for keyPairs
 * @return keyPairs
 */
const LLCS2_Algorithm::KeyPointMatrix &LLCS2_Algorithm::getKeyPairs() const {
  return keyPairs;
}

void LLCS2_Algorithm::reset(ResetLevel lvl) {
  for (auto &elem : keyPairs) {
    elem.clear();
  }
  setState(State::Constructed);
}

}// namespace lcs_solver::algorithms::llcs