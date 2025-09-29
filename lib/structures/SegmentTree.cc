/******************************************************************************
 * @file SegmentTree.cc
 * @author Steinkopp:Felix
 * @version 1.2
 * @brief Impl. of SegmentTree
 * @note SegmentTree uses function pointers: q:TxT->T is the query function
 * and u:TxT->T in U are update functions. Updates are done elementwise
 * node[i] = u(data[i], x). (T,q) should be a commutative semigroup this means
 * q(a,b) = q(b,a). Additionally, its assumed that, functions in U are be closed
 * under composition and distribute over q: u(q(a,b)) = q(u(a), u(b))
 * This notion is based on: https://doi.org/10.4230/LIPIcs.ITCS.2024.35
 ******************************************************************************/

#include <algorithm>
#include <bit>
#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

#include "structures/SegmentTree.h"
#include "util/CommonTypes.h"
#include "util/Logger.hpp"

namespace lcs_solver::structures {

/*******************************************************************************
 * Constructor (without vector)
 *
 * @tparam T typename of data in the segment tree
 * @param query q:TxT->T the query function for the segment tree (e.g. max)
 * @param range the number of leaves of the tree
 * @param order MemoryOrder used in the std::vector<T> for the data in the tree
 * @param init T value used to initialize the leaves.
 * @param nElem neutral element of the semi group (T, q).
 ******************************************************************************/
template<typename T>
SegmentTree<T>::SegmentTree(FunctionPtr query, size_t range,
                            MemoryOrder order, T init, T nElem)
    : q(query),
      nElements(range),
      treeSize(calcTreeSize(range, order)),
      traversal(order),
      neutralElement(nElem) {
  tree.resize(treeSize, init);
  lazy.resize(treeSize, false);
}

/*******************************************************************************
 * Constructor (with vector)
 *
 * @tparam T typename of data in the segment tree
 * @param query q:TxT->T the query function for the segment tree (e.g. max)
 * @param vec used to initialize the leaves of the tree
 * @param order MemoryOrder used in the std::vector<T> for the data in the tree
 * @param neutralElement neutral element of the semi group (T, q)
 ******************************************************************************/
template<typename T>
SegmentTree<T>::SegmentTree(FunctionPtr query, const std::vector<T> &vec,
                            MemoryOrder order, T neutralElement)
    : q(query),
      nElements(vec.size()),
      treeSize(calcTreeSize(nElements, order)),
      traversal(order),
      neutralElement(neutralElement) {
  tree.resize(treeSize);
  lazy.resize(treeSize);
  if (nElements > 0)
    build(vec, 0, 0, nElements - 1);
}

/*******************************************************************************
 * Query (Simple)
 *
 * @tparam T typename of data in the segment tree
 * @param index size_t in [0:size()-1] of the number of the leaf
 * @return the value in the index-th leaf of the tree
 * @note Time complexity: O(log(n))
 ******************************************************************************/
template<typename T>
T SegmentTree<T>::query(size_t index) {
  if (nElements == 0)
    return neutralElement;
  assert(index < nElements && "SegmentTree<T>::query invalid index");
  return query(0, 0, nElements - 1, index, index);
}

/*******************************************************************************
 * Query (Range)
 *
 * @tparam T typename of data in the segment tree
 * @param left index of the first element in the range (must be in 0..size()-1)
 * @param right index of the last element in the range (must be in 0..size()-1)
 * @return q evaluated over the disjunct ranges that add up to [left,right]
 ******************************************************************************/
template<typename T>
T SegmentTree<T>::query(size_t left, size_t right) {
  if (nElements == 0)
    return neutralElement;
  assert(left < nElements && "SegmentTree<T>::query invalid left index");
  assert(right < nElements && "SegmentTree<T>::query invalid right index");
  return query(0, 0, nElements - 1, left, right);
}

/*******************************************************************************
 * Getter: size()
 *
 * @tparam T typename of data in the segment tree
 * @return size_t number of leaves in the tree
 ******************************************************************************/
template<typename T>
size_t SegmentTree<T>::size() const {
  return nElements;
}


/*******************************************************************************
 * Update (Simple)
 *
 * @tparam T typename of data in the segment tree
 * @param u update function TxT->T. Will be used as data[i] = u(data[i],val)
 * @param index number of the leaf to update. Must be in [0: size()-1]
 * @param val the value used in the update function: data[i] = u(data[i],val)
 ******************************************************************************/
template<typename T>
void SegmentTree<T>::update(SegmentTree::FunctionPtr u, size_t index, T val) {
  if (nElements == 0)
    return;
  assert(index < nElements && "SegmentTree<T>::update invalid index");
  updateRange(u, 0, 0, nElements - 1, index, index, val);
}


/*******************************************************************************
 * Update (Range): Sets data[i] with i in [left, right]
 *
 * @tparam T typename of data in the segment tree
 * @param u update function TxT->T. Will be used as data[i] = u(data[i],val)
 * @param left start of range
 * @param right end of the range
 * @param val value used in the update: data[i] = u(data[i],val)
 ******************************************************************************/
template<typename T>
void SegmentTree<T>::update(SegmentTree::FunctionPtr u, size_t left,
                            size_t right, T val) {
  if (nElements == 0)
    return;
  assert(left < nElements && "SegmentTree<T>::update invalid left index");
  assert(right < nElements && "SegmentTree<T>::update invalid right index");
  updateRange(u, 0, 0, nElements - 1, left, right, val);
}


/*******************************************************************************
 * calcTreeSize
 *
 * @tparam T typename of data in the segment tree
 * @param s number leaves in the segment tree (equal to size() )
 * @param o MemoryOrder that is used in the tree
 * @return upper bound for number of nodes in the segment tree
 ******************************************************************************/
template<typename T>
size_t SegmentTree<T>::calcTreeSize(size_t s, MemoryOrder o) {
  if (s == 0)
    return 0;
  switch (o) {
    case MemoryOrder::Array: {
//      size_t treeHeight = ceil(log2(s));
//      size_t treeSize = pow(2, treeHeight + 1) - 1;   // is < than 4 * s
      size_t treeHeight = std::bit_width(s);
      size_t treeSize = (1 << (treeHeight + 1)) - 1;
      return treeSize;
    }
    case MemoryOrder::Euler: {
      return 2 * s - 1;
    }
    default: {
      util::Logger::Error() << "SegmentTree<T>::calcTreeSize: switch default";
      throw std::invalid_argument("Unimplemented item");
    }
  }
//  return 0;
}


/*******************************************************************************
 * Helper: Getter for important
 *
 * @tparam T typename of data in the segment tree
 * @param node index of a node in the segment tree
 * @param start the first index of a range
 * @param end the last index of a range
 * @return [floor( (start+end)/2 ), idx(leftchild(node), idx(rightChild(node)]
 ******************************************************************************/
template<typename T>
SegmentTree<T>::Triple SegmentTree<T>::getIdx(size_t node,
                                              size_t start, size_t end) {
  size_t mid = 0, leftChild = 0, rightChild = 0;
  mid = (start + end) / 2;
  if (traversal == MemoryOrder::Array) {
    leftChild = 2 * node + 1;
    rightChild = 2 * node + 2;
  }
  if (traversal == MemoryOrder::Euler) {
    // the tree rooted at index leftChild is responsible for range [start, mid]
    // it has 2 * ( mid - start + 1 ) - 1 vertices
    leftChild = node + 1;
    rightChild = node + 2 * (mid - start + 1);
  }
  return std::make_tuple(mid, leftChild, rightChild);
}


/*******************************************************************************
 * Helper: build - initializes the tree with data recursively
 *
 * @tparam T typename of data in the segment tree
 * @param vec vector used for initialization: leaf[i] := vec[i]
 * @param node current node to visite
 * @param start lower index of a nodes range
 * @param end upper index of a nodes range
 ******************************************************************************/
template<typename T>
void SegmentTree<T>::build(const std::vector<T> &vec, size_t node,
                           size_t start, size_t end) {
  if (start == end) {  // Base Case: Set a leaf to a value in vec
    tree[node] = vec[start];
  } else {  // Post-Order DFS Recursion
    auto [mid, leftChild, rightChild] = getIdx(node, start, end);
    build(vec, leftChild, start, mid);
    build(vec, rightChild, mid + 1, end);
    tree[node] = q(tree[leftChild], tree[rightChild]);
  }

}

/*******************************************************************************
 * Helper: updateRage - updates a range in the segment tree with a function u
 *
 * @tparam T T typename of data in the segment tree
 * @param u update function TxT->T. Will be used as data[i] = u(data[i],val)
 * @param node index of node to visite
 * @param start first value of the range of a node
 * @param end second value of the range of a node (>=start)
 * @param ul lower index of range in the update
 * @param ur upper index of range in the update
 * @param val value to used in the update: data[i] = u(data[i],val)
 ******************************************************************************/
template<typename T>
void SegmentTree<T>::updateRange(SegmentTree::FunctionPtr u, size_t node,
                                 size_t start, size_t end,
                                 size_t ul, size_t ur, T val) {
  propagateLazy(node, start, end);

  auto [mid, leftChild, rightChild] = getIdx(node, start, end);
  if (start > ur || end < ul) {
    // Range is Empty
    return;
  } else if (ul <= start && end <= ur) {
    // Current node represents the range
    tree[node] = u(tree[node], val);
    if (start != end) {
      lazy[node] = true;
    }

  } else {
    updateRange(u, leftChild, start, mid, ul, ur, val);
    updateRange(u, rightChild, mid + 1, end, ul, ur, val);
    tree[node] = q(tree[leftChild], tree[rightChild]);
  }
}


/*******************************************************************************
 * Helper: query
 *
 * @tparam T typename segment tree's data
 * @param node number of of node at is visited (in [0:nElements - 1])
 * @param nl the lower index of the range of the node
 * @param nr the upper index of the range of the node
 * @param ql lower argument of the range query
 * @param qr upper argument of the range query
 * @return T query value of the a node's range: q([nl, nr])
 ******************************************************************************/
template<typename T>
T SegmentTree<T>::query(size_t node, size_t s, size_t e,
                              size_t ql, size_t qr) {
  propagateLazy(node, s, e);
  auto [mid, leftChild, rightChild] = getIdx(node, s, e);
  if(s == ql && e == qr){
    return tree[node];
  } else if( s <= ql && qr <= mid) {
    return query(leftChild, s, mid, ql, std::min(qr,mid));
  } else if(mid + 1 <= ql && qr <= e ){
    return query(rightChild, mid+1, e, std::max(mid+1,ql), qr);
  } else {
    T queryLeft = neutralElement, queryRight = neutralElement;
    queryLeft = query(leftChild, s, mid, ql, std::min(qr,mid));
    queryRight = query(rightChild, mid + 1, e, std::max(mid+1,ql), qr);
    return q(queryLeft, queryRight);
  }
//  Logger::Warning() << "SegmentTree<T>::Bad query - "
//                      << "q(" << qLeft << ":" << qRight << ") asked about"
//                      << "[" << start << "," << end<< "] \n" ;
}


/*******************************************************************************
 * Helper: propagateLazy - pushes the query value of a node to its children
 *
 * @tparam T typename segment tree's data
 * @param node index of the node for which the value is used for pushing
 * @param start smallest leaf number node's range (in [0: size()-1]
 * @param end largest leaf number of node's range  (in 0:size()-1]
 ******************************************************************************/
template<typename T>
void SegmentTree<T>::propagateLazy(size_t node, size_t start, size_t end) {
  auto [mid, leftChild, rightChild] = getIdx(node, start, end);
  if (lazy[node] && start!= end) {
    // Push values down to children
    tree[rightChild] = tree[node];
    lazy[rightChild] = true;
    tree[leftChild] = tree[node];
    lazy[leftChild] = true;
    lazy[node] = false;
  }
}

//== Explict Initialization ====================================================
template class SegmentTree<int>;
template class SegmentTree<util::uint>;

} // namespace lcs_solver::structures
