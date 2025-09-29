/******************************************************************************
 * @file MaxSegTreeND.cc
 * @author Steinkopp:Felix
 * @version 1.2
 * @brief Impl. of SegmentTree with range max updated and point queries
 * @note The code has not yet been tested for N>2 and when different range
 * operations are mixed. The implementation uses a stack to avoid recursive
 * functions and a dense ND-Matrix.
 ******************************************************************************/

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stack>
#include <tuple>
#include <vector>

#include "structures/MaxSegTreeND.h"
#include "util/Logger.hpp"

using ::lcs_solver::util::Logger;

namespace lcs_solver::structures {

/*******************************************************************************
 * Constructor: MaxSegTreeND
 * @param d the dimensionality of the data
 * @param list container of (span, value) pairs
 ******************************************************************************/
MaxSegTreeND::MaxSegTreeND(std::span<const size_t> d, FuncSet list)
    : dim(d.begin(), d.end()),
      tree(calcTreeSize(d), neutralElement),
      lazy(calcTreeSize(d), neutralElement),
      treeChange(tree.size(), false) {
  for (const auto &[idx, val] : list) {
    update(idx, val);
  }
}

/*******************************************************************************
 * Constructor: MaxSegTreeND
 * @param mat ND-Matrix of underlying data
 ******************************************************************************/
MaxSegTreeND::MaxSegTreeND(const Matrix &mat)
    : dim(mat.dim),
      tree(calcTreeSize(dim)),
      lazy(calcTreeSize(dim)),
      treeChange(tree.size(), false) {
  if (mat.empty()) return;

  // Push initial node onto the stack
  std::stack<std::tuple<uint, uint, std::vector<Pair>, bool>> stack;
  uint node = 0;
  uint d = 0;
  std::vector<Pair> r;
  r.reserve(dim.size());
  for (const auto elements : dim)
    r.emplace_back(0, elements - 1);
  bool isSecondVisit = false;
  stack.emplace(node, d, r, isSecondVisit);

  while (!stack.empty()) {
    std::tie(node, d, r, isSecondVisit) = stack.top();
    stack.pop();

    if (r[d].first == r[d].second) { // Base Case for splitting
      if (isLastDim(d)) {
        // Handle leaf node that refers to a single value
        std::vector<size_t> idx;
        idx.reserve(dim.size());
        for (const auto &range : r)
          idx.emplace_back(range.first);
        tree[node] = mat[idx];
      } else {
        // Handle leaf in a higher dimension
        stack.emplace(node, d + 1, r, false);
      }
    } else if (!isSecondVisit) {
      // Handle first visit of truly nd-ranged nodes (add children to stack)
      uint pos = tree.vectorizeIndex(node, d);
      auto [mid, nLeft, nRight] = getIdx(pos, r[d].first, r[d].second);
      uint leftChild = tree.getIdxAbs(node, d, nLeft);
      uint rightChild = tree.getIdxAbs(node, d, nRight);

      stack.emplace(node, d, r, true);        // push parent for the 2nd time
      auto old = r[d];
      r[d] = {old.first, mid};
      stack.emplace(leftChild, d, r, false);  // push leftChild for the 1st time
      r[d] = {mid + 1, old.second};
      stack.emplace(rightChild, d, r, false); // push rightChild for the 1st time

    } // Handle second visit of truly nd-ranged nodes (calculate max now)
    else {
      if (isLastDim(d)) {
        uint pos = tree.vectorizeIndex(node, d);
        auto [mid, nLeft, nRight] = getIdx(pos, r[d].first, r[d].second);
        uint left = tree.getIdxAbs(node, d, nLeft);
        uint right = tree.getIdxAbs(node, d, nRight);
        tree[node] = std::max(tree[left], tree[right]); // Build range nodes
      } else {
        uint posInD = tree.vectorizeIndex(node, d);
        auto [mid, nLeft, nRight] = getIdx(posInD, r[d].first, r[d].second);
        uint leftInD = tree.getIdxAbs(node, d, nLeft);
        uint rightInD = tree.getIdxAbs(node, d, nRight);
        for (size_t i = 0; i < tree.dim[d + 1]; ++i) {
          uint left = tree.getIdxAbs(leftInD, d + 1, i);
          uint right = tree.getIdxAbs(rightInD, d + 1, i);
          uint rangeNodeIdx = tree.getIdxAbs(node, d + 1, i);
          tree[rangeNodeIdx] = std::max(tree[left], tree[right]);
        }
      }
    }

  } // end of while(!stack.empty())
}

/*******************************************************************************
 * MaxSegTreeND::query (point query)
 * @param index Index of a data element stored in tree leaf. Row-major order.
 * @return Result of the point query
 ******************************************************************************/
MaxSegTreeND::uint MaxSegTreeND::query(std::span<size_t> index) {
  assert(index.size() == dim.size());
  std::vector<Pair> ranges;
  for (auto val : index)
    ranges.emplace_back(val, val);
  return query(ranges);
}

/*******************************************************************************
 * MaxSegTreeND::query (range query)
 * @param r Ranges along the different dimensions for a range query
 * @return Result of the range query (maximum value in ranges)
 ******************************************************************************/
MaxSegTreeND::uint MaxSegTreeND::query(std::span<const Pair> r) {
  assert(r.size() == dim.size());
  if (tree.empty()) return neutralElement;
  if (r.size() != tree.dim.size()) return neutralElement;

  // Push initial node onto the stack
  std::stack<std::tuple<uint, uint, uint, uint, uint, uint, bool>> stack;
  uint node = 0, d = 0, start = 0, end = dim[0] - 1;
  bool w = true; // true iff node has permission to write its value as result
  stack.emplace(node, d, start, end, r[0].first, r[0].second, w);

  uint result = neutralElement, ql, qr;
  auto localChild = std::vector<Pair>(dim.size());
  auto ranged = std::vector<bool>(dim.size());
  while (!stack.empty()) {
    std::tie(node, d, start, end, ql, qr, w) = stack.top();
    stack.pop();
    uint pos = tree.vectorizeIndex(node, d);
    auto [mid, nLeft, nRight] = getIdx(pos, start, end);
    uint leftChild = tree.getIdxAbs(node, d, nLeft);
    uint rightChild = tree.getIdxAbs(node, d, nRight);
    localChild[d] = {nLeft, nRight};
    ranged[d] = start != end;

    if (isLastDim(d))
      propagateLazy(node, d, start, end); // no effects if d
    if (isLastDim(d) && treeChange[node])
      translate(node, d, localChild, ranged);

    if (start == ql && end == qr) {
      if (isLastDim(d) && w) // Last dimension and node is in range r
        result = std::max(result, tree[node]);
    } else if (isSubsetEq({ql, qr}, {start, mid})) {
      // Completely fits within the left child
      stack.emplace(leftChild, d, start, mid, ql, std::min(qr, mid), w);

    } else if (isSubsetEq({ql, qr}, {mid + 1, end})) {
      // Completely fits within the right child
      stack.emplace(rightChild, d, mid + 1, end, std::max<uint>(mid + 1, ql), qr, w);

    } else {
      stack.emplace(leftChild, d, start, mid, ql, std::min(qr, mid), w);
      stack.emplace(rightChild, d, mid + 1, end, std::max<uint>(mid + 1, ql), qr, w);
    }

    if (!isLastDim(d)) {
      if (isSubsetEq({start, end}, {ql, qr})) {
        stack.emplace(node, d + 1, 0, dim[d + 1] - 1, r[d + 1].first, r[d + 1].second, w);
      } else {
        stack.emplace(node, d + 1, 0, dim[d + 1] - 1, r[d + 1].first, r[d + 1].second, false);
      }
    }
  } // end while(!stack.empty())

  return result;
}

/*******************************************************************************
 * MaxSegTreeND::update (point update)
 * @param index Index of a data element stored in tree leaf. Row-major order.
 * @param value used in the update tree[index] = max(Sets tree[index], value)
 ******************************************************************************/
void MaxSegTreeND::update(std::span<const size_t> index, uint value) {
  assert(index.size() == dim.size());
  std::vector<Pair> ranges;
  for (auto val : index)
    ranges.emplace_back(val, val);
  update(ranges, value);
}

/*******************************************************************************
 * MaxSegTreeND::update (range update)
 * @param r span of closed intervals along dimensions in row-major order
 * @param value used in the update tree[index] = max(Sets tree[index], value)
 ******************************************************************************/
void MaxSegTreeND::update(std::span<const Pair> r, uint value) {
  assert(r.size() == dim.size());
  if (tree.empty()) return;
  if (r.size() != tree.dim.size()) return;

  // Push initial node onto the stack
  std::stack<std::tuple<uint, uint, uint, uint, uint, uint>> stack;
  uint node = 0, d = 0;
  uint start = 0, end = dim[0] - 1, ul = r[0].first, ur = r[0].second;
  stack.emplace(node, d, start, end, ul, ur);

  while (!stack.empty()) {
    std::tie(node, d, start, end, ul, ur) = stack.top();
    stack.pop();
    uint pos = tree.vectorizeIndex(node, d);
    auto [mid, nLeft, nRight] = getIdx(pos, start, end);
    uint leftChild = tree.getIdxAbs(node, d, nLeft);
    uint rightChild = tree.getIdxAbs(node, d, nRight);

    if (isLastDim(d))
      propagateLazy(node, d, start, end);
    updateTree(node, value); // Always update the current node in the tree

    if (ul == start && end == ur) {
      // Range are identical: Write lazy[node] or reinterpret d via stack push
      if (isLastDim(d)) {
        if (start != end) lazy[node] = std::max(lazy[node], value);
      } else {
        stack.emplace(node, d + 1, 0, dim[d + 1] - 1, r[d + 1].first, r[d + 1].second);
      }
    } else if (isSubsetEq({ul, ur}, {start, mid})) {
      // Only left child is relevant
      stack.emplace(leftChild, d, start, mid, ul, ur);
    } else if (isSubsetEq({ul, ur}, {mid + 1, end})) {
      // Only right child is relevant
      stack.emplace(rightChild, d, mid + 1, end, ul, ur);
    } else {
      // Not a complete overlap and not a leaf node, push children
      stack.emplace(leftChild, d, start, mid, ul, std::min(ur, mid));
      stack.emplace(rightChild, d, mid + 1, end, std::max<uint>(mid + 1, ul), ur);
    }
  } //  end while(!stack.empty())
}

/*******************************************************************************
 * MaxSegTreeND::propagateLazy
 * @param node uint - the index of the current node
 * @param d uint the dimension in which the children are looked up
 * @param start uint begin of the node's range in dimension d
 * @param end  uint end of the node's range in dimension d
 * @note If node is not a leaf the value of lazy[node] is used for updating
 *  the children of node in dim d. The tree is updated with updateTree(node,
 *  lazy[node]). Finally lazy[node] will be reset to the neutralElement.
 ******************************************************************************/
void MaxSegTreeND::propagateLazy(const uint node, const uint d, const uint start, const uint end) {
  if (lazy[node] == neutralElement) return;

  uint pos = tree.vectorizeIndex(node, d);
  auto [mid, nLeft, nRight] = getIdx(pos, start, end);
  if (start == end) {
    updateTree(node, lazy[node]);
    lazy[node] = neutralElement;
  }
  if (start != end) {
    // Calculate Indices
    uint leftChild = tree.getIdxAbs(node, d, nLeft);
    uint rightChild = tree.getIdxAbs(node, d, nRight);
    // Push values down to children
    updateTree(node, lazy[node]);
    lazy[rightChild] = std::max(lazy[rightChild], lazy[node]);
    lazy[leftChild] = std::max(lazy[leftChild], lazy[node]);
    lazy[node] = neutralElement;
  }
}

/*******************************************************************************
 * MaxSegTreeND::operator[]
 * @param indices std::span<const uint> for accessing data via row major order
 * @return query(index)
 ******************************************************************************/
MaxSegTreeND::uint MaxSegTreeND::operator[](std::span<size_t> indices) {
  return query(indices);
}

/*******************************************************************************
 * MaxSegTreeND::calcTreeSize
 * @param n uint the number of data elements managed by in a 1D-Segment tree
 * @return uint Number of nodes in a segment tree with n elements
 ******************************************************************************/
size_t MaxSegTreeND::calcTreeSize(size_t n) {
  return 2 * n - 1;
}

/*******************************************************************************
 * MaxSegTreeND::getIdx
 * @details Getter for middle of two uint, leftChild(n) and rightChild(n)
 * @param node uint index of a node in a 1D Segment tree
 * @param start uint the lowest index of the node's range
 * @param end  uint the largest index of the node's range
 * @return Triple ( (start + end) / 2, leftChild(n), rightChild(n) )
 ******************************************************************************/
MaxSegTreeND::Triple MaxSegTreeND::getIdx(const uint node, const uint start, const uint end) {
  uint mid = (start + end) / 2;
  uint leftChild = node + 1; // left tree for [start, mid] has 2*(mid-start+1)-1 nodes
  uint rightChild = node + 2 * (mid - start + 1);
  return std::make_tuple(mid, leftChild, rightChild);
}

/*******************************************************************************
 * MaxSegTreeND::calcTreeSize
 * @details Getter for a vector with number of nodes for each dimension
 * @param d std::span<uint> with count of leaves in the different dimensions
 * @return std::vector<uint> with count of nodes for each dimension
 ******************************************************************************/
std::vector<size_t> MaxSegTreeND::calcTreeSize(std::span<const size_t> d) {
  if (d.empty())
    return {};
  std::vector<size_t> res;
  res.reserve(d.size());
  std::ranges::transform(d, std::back_inserter(res),
      [](const uint val) {
        return MaxSegTreeND::calcTreeSize(val);
      });
  return res;
}

/*******************************************************************************
 * MaxSegTreeND::size()
 * @return number of nodes in tree
 ******************************************************************************/
MaxSegTreeND::uint MaxSegTreeND::size() const {
  if (dim.empty())
    return 0;
  const auto &vec = tree.dim;
  return std::accumulate(vec.begin(), vec.end(), 1, std::multiplies<>());
}

/*******************************************************************************
 * MaxSegTreeND::translate
 * @details During a query, we ascend down along dimensions 0, ..., d-1 until we
 * reach a node interpreted in dimension d. If one of the nodes along the path to
 * the destination was ranged, we need to fix the value propagation. For example,
 * lazy[{0:1}{3:4}] needs to be propagated to lazy[{0}{3:4}] and lazy[{1}{3:4}].
 * @param dest uint number of a node (must be a leaf in the last dimension)
 * @param d uint dim in which dest is interpreted (should be dim.size()-1)
 * @param child vector<Pair> with numbers of children along different axes on
 * the path to dest. The number of the children is taken with respect to the
 * tree in their dimension.
 * @param ranged vector<bool> ranged[i] is true if and only if the node for the
 * ith dimension on the path to dest represents a range in the ith dimension.
* @note Time Complexity O(2^(d-1))
******************************************************************************/
void MaxSegTreeND::translate(
    uint dest,
    uint d,
    const std::vector<Pair> &child,
    const std::vector<bool> &ranged) {
  // Setup idx = {child[0].first, child[1].first, ..., child[d-1].first, dest}
  if (stopTranslateEarly(dest, d, ranged)) return;
  auto idx = firstPath(dest, d, child, ranged);

  bool validPathToChild = true;
  while (validPathToChild) {
    uint node = lazy.linearIndex(idx); // calc relevant node idx
    if (ranged[d])
      lazy[node] = std::max(lazy[node], lazy[dest]); // dest is to a range in d
    else
      updateTree(node, tree[dest]);

    validPathToChild = incTranslationPath(idx, dest, d, child, ranged);
  }
  treeChange[dest] = false; // treeChange has been pushed up into
}

/*******************************************************************************
 * MaxSegTreeND::isLastDim
 * @param d uint of a dimension
 * @return d == largest dimension of the matrix tree, lazy
 ******************************************************************************/
bool MaxSegTreeND::isLastDim(uint d) const {
  return d + 1u == dim.size(); // dim = tree.dim = lazy.dim
}

/*******************************************************************************
 * MaxSegTreeND::isSubsetEq
 * @param a Pair representing a closed interval
 * @param b Pair representing a closed interval
 * @return true if \f$a \subseteq b\f$, otherwise false
 ******************************************************************************/
bool MaxSegTreeND::isSubsetEq(const MaxSegTreeND::Pair &a, const MaxSegTreeND::Pair &b) {
  return b.first <= a.first && a.second <= b.second;
}

/*******************************************************************************
 * MaxSegTreeND::updateTree
 * @param node uint representing the index of a node in the tree matrix
 * @param value uint used in the update function
 * @note Updates tree[node] to max(tree[node], value)
 * @note Sets treeChange[node] = true
 ******************************************************************************/
void MaxSegTreeND::updateTree(uint node, uint value) {
  if (tree[node] < value) {
    tree[node] = value;
    treeChange[node] = true;
  }
}

/*******************************************************************************
 * MaxSegTreeND::firstPath
 * @param dest uint representing the number of a node in the tree matrix
 * @param d uint dim in that the node dest is interpreted
 * @param child vector of size d+1 with local children number on path to dest
 * @param ranged vector<bool> with ranged[i] true iff the node i (for dim i) on
 *  the path to dest is a range
 * @return std::vector<uint> of size d to compute the first child number of
 *  the nodes that need to be updated when dest was changed. If no such child
 *  exists, the vector will be a vectorization of the node number dest.
 ******************************************************************************/
std::vector<size_t> MaxSegTreeND::firstPath(
    uint dest,
    uint d,
    const std::vector<Pair> &child,
    const std::vector<bool> &ranged) const {
  uint pos = tree.vectorizeIndex(dest, d);
  std::vector<size_t> vec;
  vec.reserve(d + 1);

  for (uint i = 0; i < d; ++i) {
    if (ranged[i])
      vec.emplace_back(child[i].first);
    else
      vec.emplace_back(tree.vectorizeIndex(dest, i));
  }
  vec.emplace_back(pos);
  return vec;
}

/*******************************************************************************
 * MaxSegTreeND::stopTranslateEarly
 * @param dest uint representing a node in tree matrix
 * @param d uint dimension in that the node dest is interpreted (zero-based)
 * @param ranged vector<bool> with ranged[i] true iff the node i (for dim i) on
 *  the path to dest is a range
 * @return false if a translation is necessary, otherwise return true
 ******************************************************************************/
bool MaxSegTreeND::stopTranslateEarly(
    uint dest,
    uint d,
    const std::vector<bool> &ranged) const {
  if (dim.size() == 1) return true;
  if (tree[dest] == neutralElement) return true; // for quicker debugging

  for (uint i = 0; i < d; ++i)
    if (ranged[i])
      return false;
  return true; // There are no children that need to receive the value from dest
}

/*******************************************************************************
 * MaxSegTreeND::incTranslationPath
 * @param path std::vector<uint> vectorized version of the node number of a
 *  valid translation's destination candidate (a child on the path to source)
 * @param source uint index of a node in the tree matrix
 * @param d uint dimension in which the node dest is interpreted (zero-based)
 * @param child vector of size d+1 with local children numbers on path to dest
 * @param ranged vector<bool> with ranged[i] true iff the node i (for dim i) on
 *  the path to dest is a range
 * @note path is set to the next child. If there's no next candidate, the
 *  path switches back and forth between the last candidate and its predecessor
 * @return false if and only if the path is set to the last child (indicates
 *  that all children have been processed)
 ******************************************************************************/
bool MaxSegTreeND::incTranslationPath(std::vector<size_t> &path,
                                      const uint source,
                                      const uint d,
                                      const std::vector<Pair> &child,
                                      const std::vector<bool> &ranged) const {
  uint i = d; // path[d] is const = node in one of the last dims seg trees
  while (i--) {
    if (!ranged[i]) {
      path[i] = tree.vectorizeIndex(source, i); // no children in dim i
      if (i == 0)
        return false;
    } else {
      if (path[i] == child[i].first) {
        path[i] = child[i].second;
        break; // successful incrementation
      } else {
        path[i] = child[i].first; // and need to change path[x] with x<i
      }
      if (i == 0) return false;
    }
  } // end of path change loop
  return true;
}

}  // end of namespace lcs_solver::structures