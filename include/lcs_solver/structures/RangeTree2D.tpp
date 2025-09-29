#ifndef LCS_SOLVER_STRUCTURES_RANGETREE2D_TPP_
#define LCS_SOLVER_STRUCTURES_RANGETREE2D_TPP_

#include <cassert>
//#include <iostream>
#include <string>
#include <sstream>
namespace lcs_solver::structures {
///// Result Stuff /////////////////////////////////////////////////////////////
//=== Constructor ==============================================================
template<typename T1, typename T2>
RangeTree2D<T1, T2>::Result::Result(const NodeVector &nodes, const Range2 &r2)
    : nodes(nodes), r2(r2) {}

// //=== Copy Constructor =========================================================
template<typename T1, typename T2>
RangeTree2D<T1, T2>::Result::Result(const Result &other) : nodes(other.nodes), r2(other.r2) {}
//
// //=== Move Constructor =========================================================
template<typename T1, typename T2>
RangeTree2D<T1, T2>::Result::Result(Result &&other) noexcept : nodes(std::move(other.nodes)), r2(std::move(other.r2)) {}

// === Copy assignment ==========================================================
 template<typename T1, typename T2>
 typename RangeTree2D<T1, T2>::Result &RangeTree2D<T1, T2>::Result::operator=(const Result &other) {
   if (this != &other) {
     nodes = other.nodes;
     r2 = other.r2;
   }
   return *this;
 }

//=== Move assignment ==========================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result &RangeTree2D<T1, T2>::Result::operator=(Result &&other) noexcept {
  if (this != &other) {
    nodes = std::move(other.nodes);
    r2 = std::move(other.r2);
  }
  return *this;
}

//=== cbegin ==== ==============================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result::Iterator RangeTree2D<T1, T2>::Result::cbegin() const {
  if (nodes.empty()) {
    return Iterator(nodes.cend(), nodes.cend(), r2);
  }

  return Iterator(nodes.cbegin(), nodes.cend(), r2);
}

//=== cend =====================================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result::Iterator RangeTree2D<T1, T2>::Result::cend() const {
  return Iterator(nodes.cend(), nodes.cend(), r2);
}

///// Iterator Stuff ///////////////////////////////////////////////////////////
//=== Constructor ==============================================================
template<typename T1, typename T2>
RangeTree2D<T1, T2>::Result::Iterator::Iterator(
    NodeIter start, NodeIter end, const Range2 &r2)
    : currentPairIter(start),
      endPairIter(end),
      pCurrentElem(nullptr),
      pLastElem(nullptr),
      r2(r2) {
  if (currentPairIter == endPairIter) {
    return;// no Pairs in vector
  }
  pCurrentElem = currentPairIter->second;
  pLastElem = &currentPairIter->first->arr.back();
  if (!pCurrentElem || !validPosition(pCurrentElem->y, r2)) {
    setNextValidPosition();
  }
  const bool isNull = pCurrentElem == nullptr;
  assert(isNull || (r2.first <= pCurrentElem->y && pCurrentElem->y <= r2.second));
}

template<typename T1, typename T2>
RangeTree2D<T1, T2>::Result::Iterator::Iterator(const Iterator &other)
    : currentPairIter(other.currentPairIter),
      endPairIter(other.endPairIter),
      pCurrentElem(other.pCurrentElem),
      pLastElem(other.pLastElem),
      r2(other.r2) {}

template<typename T1, typename T2>
RangeTree2D<T1, T2>::Result::Iterator::Iterator(Iterator &&other) noexcept
: currentPairIter(other.currentPairIter),
     endPairIter(other.endPairIter),
     pCurrentElem(other.pCurrentElem),
     pLastElem(other.pLastElem),
     r2(other.r2){ }

template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result::Iterator &RangeTree2D<T1, T2>::Result::Iterator::operator=(const Iterator &other) {
  if (this != &other) {
    currentPairIter = other.currentPairIter;
    endPairIter = other.endPairIter;
    pCurrentElem = other.pCurrentElem;
    pLastElem = other.pLastElem;
    r2 = other.r2;
  }
  return *this;
}
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result::Iterator &RangeTree2D<T1, T2>::Result::Iterator::operator=(Iterator &&other) noexcept {
  if (this != &other) {
    currentPairIter = std::move(other.currentPairIter);
    endPairIter = std::move(other.endPairIter);
    pCurrentElem = std::move(other.pCurrentElem);
    pLastElem = std::move(other.pLastElem);
    // r2 = other.r2;
    r2 = std::move(other.r2);
  }
  return *this;
}

template<typename T1, typename T2>
bool RangeTree2D<T1, T2>::Result::Iterator::validPosition(T2 y, const Range2 &r2) {
  const auto &[lower, upper] = r2;
  return lower <= y && y <= upper;
}

//=== Pre Increment Operator ===================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result::Iterator &RangeTree2D<T1, T2>::Result::Iterator::operator++() {
  advance();
  const bool isNull = pCurrentElem == nullptr;
  assert(isNull || (r2.first <= pCurrentElem->y && pCurrentElem->y <= r2.second));
  return *this;
}

//=== operator!= ===============================================================
template<typename T1, typename T2>
bool RangeTree2D<T1, T2>::Result::Iterator::operator!=(const Iterator &other) const {
  return currentPairIter != other.currentPairIter || pCurrentElem != other.pCurrentElem;
}

//=== operator== ===============================================================
template<typename T1, typename T2>
bool RangeTree2D<T1, T2>::Result::Iterator::operator==(const Iterator &other) const {
  return !(*this != other);
}

//=== operator* ================================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::PointType RangeTree2D<T1, T2>::Result::Iterator::operator*() const {
  if (pCurrentElem == nullptr)
    return {T1(), T2()};
  return {pCurrentElem->x, pCurrentElem->y};
}

//=== advance ==================================================================
template<typename T1, typename T2>
void RangeTree2D<T1, T2>::Result::Iterator::advance() {
  if (currentPairIter == endPairIter) {// Nothing to advance to
    pCurrentElem = nullptr;
    pLastElem = nullptr;
    return;
  }

  // Try to advance by moving to the next Element (in the same node)
  if (pCurrentElem != pLastElem && (++pCurrentElem)->y <= r2.second) {
    return;
  }

  // Need to look for valid Pair with currentPairIter
  if (!setNextValidPosition()) {
    pCurrentElem = nullptr;
    pLastElem = nullptr;
  }
}

//=== findValidPosition ========================================================
template<typename T1, typename T2>
bool RangeTree2D<T1, T2>::Result::Iterator::setNextValidPosition() {
  //std::cout << "setNextValidPosition" << std::endl;
  const Element *pOldElem = pCurrentElem;
  while (currentPairIter != endPairIter) {
    if (currentPairIter->second == nullptr) {// Skip irrelevant pair
      ++currentPairIter;
      if (currentPairIter != endPairIter) {
        pCurrentElem = currentPairIter->second;
        pLastElem = &currentPairIter->first->arr.back();
      } else {
        pCurrentElem = nullptr;
        pLastElem = nullptr;
      }
      continue;
    }

    // Go through Element Vector in a node
    while (pCurrentElem != pLastElem) {
      if (r2.first <= pCurrentElem->y && pCurrentElem->y <= r2.second) {
        return true;
      }
      ++pCurrentElem;
    }
    // Last element: pElement == endElem
    if (pOldElem != pCurrentElem && validPosition(pCurrentElem->y, r2)) {
      return true;
    }

    // Move to next Pair (in the vector of pairs)
    ++currentPairIter;
    if (currentPairIter != endPairIter) {
      pCurrentElem = currentPairIter->second;
      pLastElem = &currentPairIter->first->arr.back();
    } else {
      pCurrentElem = nullptr;
      pLastElem = nullptr;
    }
  }// end while
  return false;
}

//=== Iterator =================================================================
template<typename T1, typename T2>
std::string RangeTree2D<T1, T2>::Result::Iterator::DebugString() const {
  std::stringstream ss;
  if (currentPairIter == endPairIter) {
    ss << "currentPairIter == endPairIter \n";
  } else {
    Node *node = currentPairIter->first;
    ss << "currentPairIter != endPairIter\n";
    ss << "\t node = ";
    std::stringstream input(node->DebugString());
    std::string line;
    while (std::getline(input, line)) {
      ss << "\n \t"
         << line;
    }
    ss << "\t pElement = " << currentPairIter->second;
  }
  ss << "pCurrentElem = " << pCurrentElem << " (" << pCurrentElem->x << ", " << pCurrentElem->y << ") \n";
  ss << "pLastElem = " << pLastElem << " (" << pLastElem->x << ", " << pLastElem->y << ") \n";
  ss << "r2 = (" << r2.first << ", " << r2.second << ")";
  return ss.str();
}

///// RangeTree2D //////////////////////////////////////////////////////////////
//=== Constructor ==============================================================
template<typename T1, typename T2>
RangeTree2D<T1, T2>::RangeTree2D(std::span<PointType const> points) {
  std::vector<PointType> sortedPoints(points.begin(), points.end());
  std::sort(sortedPoints.begin(), sortedPoints.end());
  const size_t nPoints = sortedPoints.size();
  size_t nx = 1;
  switch (nPoints) {
    case 0: nx = 0; break;
    case 1: nx = 1; break;
    default:
      for (size_t i = 0; i < nPoints - 1; ++i) {
        if (sortedPoints[i].first != sortedPoints[i + 1].first) {
          ++nx;
        }
      }
  }
  switch (nPoints) {
    case 0: root = nullptr; break;
    case 1: {
      T1 x = points[0].first;
      T2 y = points[0].second;
      Element element = Element(x, y, nullptr, nullptr);                                // no fractional cascading pointers
      root = std::make_unique<Node>(x, std::vector<Element>{element}, nullptr, nullptr);// no children in primitive tree
      break;
    }
    default: {
      std::vector<T1> xs = {sortedPoints[0].first};
      xs.reserve(nx);
      std::vector<std::vector<T2>> ys = {{sortedPoints[0].second}};
      for (size_t i = 0; i < nPoints - 1; ++i) {
        if (sortedPoints[i].first != sortedPoints[i + 1].first) {
          xs.push_back(sortedPoints[i + 1].first);
          auto tempVec = std::vector<T2>{sortedPoints[i + 1].second};
          ys.push_back(tempVec);
        } else {
          ys.back().push_back({sortedPoints[i + 1].second});
        }
      }
      assert(xs.size() == ys.size());
      assert(xs.size() == nx);
      root = buildTree(xs, ys, 0, xs.size() - 1);
    }
  }
}

//=== buildTree ================================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::NodePtr RangeTree2D<T1, T2>::buildTree(
    const std::vector<T1> &xs,
    const std::vector<std::vector<T2>> &ys,
    size_t start,
    size_t end) {
  if (start > end) return nullptr;
  if (start == end) {
    const T1 x = xs[start];
    const std::vector<T2> &yVec = ys[start];
    std::vector<Element> elements;
    elements.reserve(yVec.size());
    for (const auto &y : yVec) {
      elements.emplace_back(x, y, nullptr, nullptr);// no fractional cascading pointers
    }
    return std::make_unique<Node>(x, elements, nullptr, nullptr);// no children
  }

  size_t middle = (start + end) / 2;
  auto node = std::make_unique<Node>(xs[middle]);

  // Build left and right children
  node->lchild = buildTree(xs, ys, start, middle);
  node->rchild = buildTree(xs, ys, middle + 1, end);

  // Merge sorted arrays for fractional cascading and y arrays
  mergeElements(node.get());

  return node;
}

//=== buildTree ================================================================
// template<typename T1, typename T2>
// typename RangeTree2D<T1, T2>::NodePtr RangeTree2D<T1, T2>::buildTree(
//     const std::span<PointType> &points,
//     size_t start,
//     size_t end) {
//   if (start > end) return nullptr;
//   if (start == end) {
//     T1 x = points[start].first;
//     T2 y = points[start].second;
//     Element element = Element(x, y, nullptr, nullptr);// no fractional cascading pointers
//     return std::make_unique<Node>(x, std::vector<Element>{element}, nullptr, nullptr);
//   }
//
//   size_t middle = (start + end) / 2;
//   auto node = std::make_unique<Node>(points[middle].first);
//
//   // Build left and right children
//   node->lchild = buildTree(points, start, middle);
//   node->rchild = buildTree(points, middle + 1, end);
//
//   // Merge sorted arrays for fractional cascading and y arrays
//   mergeElements(node.get());
//
//   return node;
// }

//=== mergeElements ============================================================
template<typename T1, typename T2>
void RangeTree2D<T1, T2>::mergeElements(Node *node) {
  if (!node->lchild && !node->rchild) { return; }
  if (!node->rchild) {
    node->arr = node->lchild->arr;
    return;
  }
  if (!node->lchild) {
    node->arr = node->rchild->arr;
    return;
  }

  // Normal Merge (node has two children)
  const auto &leftVec = node->lchild->arr;
  const auto &rightVec = node->rchild->arr;
  auto leftIter = leftVec.cbegin();
  auto rightIter = rightVec.cbegin();

  size_t i = 0, j = 0;
  while (i < leftVec.size() && j < rightVec.size()) {
    if (leftVec[i].y <= rightVec[j].y) {
      const T1 &x = leftVec[i].x;
      const T2 &y = leftVec[i].y;
      while (rightIter != rightVec.cend() && rightIter->y < y) {
        ++rightIter;
      }
      node->arr.emplace_back(x, y, &leftVec[i], &(*rightIter));
      ++i;
    } else {
      const T1 &x = rightVec[j].x;
      const T2 &y = rightVec[j].y;
      while (leftIter != leftVec.cend() && leftIter->y < y) { ++leftIter; }
      node->arr.emplace_back(x, y, &(*leftIter), &rightVec[j]);
      ++j;
    }
  }

  while (i < leftVec.size()) {
    const T1 &x = leftVec[i].x;
    const T2 &y = leftVec[i].y;
    node->arr.emplace_back(x, y, &leftVec[i], nullptr);
    ++i;
  }

  while (j < rightVec.size()) {
    const T1 &x = rightVec[j].x;
    const T2 &y = rightVec[j].y;
    node->arr.emplace_back(x, y, nullptr, &rightVec[j]);
    ++j;
  }
}

//=== findSplitNode ============================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Node *RangeTree2D<T1, T2>::findSplitNode(const Range1 r1) const {
  const auto [start, end] = r1;
  Node *v = root.get();
  if (v == nullptr) return nullptr;
  while (!isLeaf(v) && (end <= v->x || start > v->x)) {
    if (v->lchild != nullptr && end <= v->x) {
      v = v->lchild.get();// end is in the left subtree of node
    } else {
      v = v->rchild.get();// start is in the right subtree of node
    }
  }
  return v;
}

//=== findFirstElem ============================================================
template<typename T1, typename T2>
const typename RangeTree2D<T1, T2>::Element *RangeTree2D<T1, T2>::findFirstElem(const Node *node, const Range2 &r2) const {
  if (node == nullptr || node->arr.empty()) return nullptr;
  auto it = std::lower_bound(
      node->arr.begin(),
      node->arr.end(), r2.first,
      [](const Element &elem, const T2 &value) {
        return elem.y < value;
      });
  return (it != node->arr.end()) ? &(*it) : nullptr;
}

//=== Query ====================================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result RangeTree2D<T1, T2>::query(Range1 d1, Range2 d2) const {
  // std::cout << "New Query: (" << d1.first << ", " << d1.second << ") x ("
  //           << d2.first << ", " << d2.second << ") \nTree: \n";
  // std::string temp = DebugString();
  // std::cout << temp << "\n"
  //           << std::endl;
  auto nodes = collectNodesInXRange(d1, d2);
  return Result(std::move(nodes), d2);
}

//=== collectNodesInXRange =====================================================
template<typename T1, typename T2>
typename RangeTree2D<T1, T2>::Result::NodeVector
RangeTree2D<T1, T2>::collectNodesInXRange(const Range1 r1, const Range2 r2) const {

  // Report Split Note if necessary
  const Node *vSplit = findSplitNode(r1);//node where paths to r1.first and r1.second split
  if (vSplit == nullptr) return {};
  typename Result::NodeVector nodes;                      // nodes which subtrees together are the solution
  const Element *fstSplitElem = findFirstElem(vSplit, r2);// finds first valid point in vSplit->arr
  if (isLeaf(vSplit)) {
    if (isValueIn(vSplit->x, r1)) {
      nodes.emplace_back(vSplit, fstSplitElem);
    }
    return nodes;
  }

  // Follow path to range.first and report nodes in subtree right of the path
  const Node *v = vSplit->lchild.get();
  const Element *p = fstSplitElem ? fstSplitElem->lc : nullptr;
  while (!isLeaf(v)) {
    if (r1.first <= v->x) {
      if (v->rchild.get() != nullptr) {
        nodes.emplace_back(v->rchild.get(), p ? p->rc : nullptr);
      }
      v = v->lchild.get();
      p = p ? p->lc : nullptr;
    } else {
      v = v->rchild.get();
      p = p ? p->rc : nullptr;
    }
  }
  if (isValueIn(v->x, r1)) {
    nodes.emplace_back(v, &v->arr.front());// Add leaf of search path if necessary
  }

  // Follow path to range.second and report nodes in subtree left of the path
  v = vSplit->rchild.get();
  p = fstSplitElem ? fstSplitElem->rc : nullptr;
  while (!isLeaf(v)) {
    if (v->x <= r1.second) {
      if (v->lchild.get() != nullptr) {
        nodes.emplace_back(v->lchild.get(), p ? p->lc : nullptr);
      }
      v = v->rchild.get();
      p = p ? p->rc : nullptr;
    } else {
      v = v->lchild.get();
      p = p ? p->lc : nullptr;
    }
  }
  if (isValueIn(v->x, r1)) {
    nodes.emplace_back(v, &v->arr.front());// Add leaf of search path if necessary
  }
  return nodes;
}

//=== isLeaf ===================================================================
template<typename T1, typename T2>
bool RangeTree2D<T1, T2>::isLeaf(const Node *p) {
  // if (p == nullptr)
  //   return true;
  return !p->lchild && !p->rchild;
}

//=== isValueIn ================================================================
template<typename T1, typename T2>
bool RangeTree2D<T1, T2>::isValueIn(T1 x, Range1 r) {
  return r.first <= x && x <= r.second;
}

///// DebugStrings /////////////////////////////////////////////////////////////
template<typename T1, typename T2>
std::string RangeTree2D<T1, T2>::DebugString() const {
  std::stringstream ss;
  if (root == nullptr) {
    ss << "empty";
    return ss.str();
  }

  std::queue<std::pair<const Node *, int>> q;
  q.push(std::make_pair(root.get(), 0));
  int bfsIndex = 0;
  while (!q.empty()) {
    const auto [node, level] = q.front();
    // Node * node = q.front().first.get();
    // int level = q.front().second;
    ss << "Node: " << bfsIndex << " ";
    ss << "(l=" << level << ", addr=" << node << ") ";
    ss << "Data: ";
    std::stringstream input(node->DebugString());
    std::string line;
    while (std::getline(input, line)) {
      ss << "\n"
         << line;
    }
    ss << std::endl;
    ++bfsIndex;
    q.pop();
    if (node->lchild != nullptr) {
      q.push({node->lchild.get(), level + 1});
    }
    if (node->rchild != nullptr) {
      q.push({node->rchild.get(), level + 1});
    }
  }
  return ss.str();
}

//=== Element ==================================================================
template<typename T1, typename T2>
std::string RangeTree2D<T1, T2>::Element::DebugString() const {
  std::stringstream ss;
  ss << "x=" << x << " ";
  ss << "y=" << y << " ";
  ss << "addr=" << this << " ";
  ss << "lc=" << lc << " ";
  ss << "rc=" << rc;
  return ss.str();
}
//=== Node =====================================================================
template<typename T1, typename T2>
std::string RangeTree2D<T1, T2>::Node::DebugString() const {
  std::stringstream ss;
  ss << "x=" << x << " ";
  ss << "addr=" << this << " ";
  ss << "lchild=" << lchild << " ";
  ss << "rchild=" << rchild << " ";
  ss << "arr= [\n";
  for (const auto &element : arr) {
    ss << "\t" << element.DebugString();
    if (&element != &arr.back())
      ss << "\n";
  }
  ss << "\n]";
  return ss.str();
}

}// namespace lcs_solver::structures

#endif// LCS_SOLVER_STRUCTURES_RANGETREE2D_TPP_