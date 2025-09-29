/*******************************************************************************
 * @file RMQ.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Impl. of Range Minimum/Maximum Queries (see: The LCA Problem Revisited
 * https://link.springer.com/chapter/10.1007/10719839_9)
 ******************************************************************************/
#include "util/RMQ.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <bitset>
#include <stack>

//#include "util/Logger.hpp"

namespace lcs_solver::util {

//== RMQ_ONN <O(n^2),O(1)> =====================================================

/*******************************************************************************
 * Constructor for a range minimum (maximum) query data structure instance.
 * @tparam T typename of data
 * @param type range query type (RMQ_TYPE::MIN or RMQ_TYPE::MIN)
 * @param arr vector with data for which the queries are asked
 * @note does the preprocessing in O(n*n)
 ******************************************************************************/
template<typename T>
RMQ_ONN<T>::RMQ_ONN(RMQ_TYPE type, const std::vector<T> &arr)
    : mType(type), h(arr), n(arr.size()), M(Matrix(n, std::vector<T>(n, 0))) {
  preprocess();
}

/*******************************************************************************
 * Does the Preprocessing: Precompute answers with dynamic Programming
 * @tparam T typename of data
 * @note Time complexity: O(N^2)
 ******************************************************************************/
template<typename T>
void RMQ_ONN<T>::preprocess() {
  //Base case
  for (size_t i = 0; i < n; ++i) {
    M[i][i] = i;
  }
  //Fill with Dynamic Programming O(N^2)
  if (mType == RMQ_TYPE::MIN) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (h[j] < h[M[i][j - 1]])
          M[i][j] = j;
        else
          M[i][j] = M[i][j - 1];
      }
    }
  } else {
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        if (h[j] > h[M[i][j - 1]])
          M[i][j] = j;
        else
          M[i][j] = M[i][j - 1];
      }
    }
  }
}

/*******************************************************************************
 * query the RMQ data structures. Time complexity: O(1)
 * @tparam T typename of data
 * @param left Position of the first element in range
 * @param right The length of the first string v = m.
 * @return Index of the min or max of { h[left], h[left + 1], ..., h[right]}
 ******************************************************************************/
template<typename T>
size_t RMQ_ONN<T>::query(size_t left, size_t right) {
  // Do some checks
  // assert(left <n && "Bad RMQ_ONN Query");
  // assert(right < n && "Bad RMQ_ONN Query");
  // assert(left <= right && "Bad RMQ_ONN Query");
  if (left > right)
    std::swap(left, right);
  if (mType == RMQ_TYPE::MIN && (left > n || right > n))
    return std::numeric_limits<size_t>::max();
  if (mType == RMQ_TYPE::MAX && (left > n || right > n))
    return std::numeric_limits<size_t>::max();
  return M[left][right]; // Return precomputed answer
}

//== RMQ with O(n log n) preprocessing =========================================

/*******************************************************************************
 * Constructor for a range minimum (maximum) query data structure instance
 * @tparam T typename of data
 * @param type Type of queries
 * @param arr The array for which the queries are asked
 ******************************************************************************/
template<typename T>
RMQ_nlogn<T>::RMQ_nlogn(RMQ_TYPE type, const std::vector<T> &arr)
    : mType(type), h(arr), n(arr.size()),
      LOG2(n + 1),
      POW2(LOG2[n + 1] + 2),
      M(arr.empty() ? Matrix()
                    : Matrix(n,
                             std::vector<size_t>(
                                 std::max<size_t>(LOG2[n] + 1, 2),
                                 0
                             ))
      ) {
  if (n != 0)
    preprocess();
}

/*******************************************************************************
 * preprocess. Time complexity: O(n log n)
 * @tparam T typename of data
 * @note Effect: M[i][j] = index of min { h[i], h[i+1], ..., h[i+pow(2,j)-1] }
 ******************************************************************************/
template<typename T>
void RMQ_nlogn<T>::preprocess() {
  /* Note: If you use `<=` in the comparisons the first occurrence will be
   * selected. If you change `<=` to `<` the last occurrence will be selected
   * If you change it, don't forget to change it in query too! */
  if (n == 0) return;
  if (mType == RMQ_TYPE::MIN) {
    // k = 1
    for (size_t i = 0; i < n; ++i) {
      M[i][0] = i;
      size_t sndPos = (i + 1) < n ? i + 1 : n - 1;
      if (h.at(i) <= h.at(sndPos)) {
        M[i][1] = i;
      } else {
        M[i][1] = sndPos;
      }
    }
    //k > 1
    for (size_t k = 2; k <= LOG2[n]; ++k) { // if you use floating log2 write <
      for (size_t i = 0; i < n; ++i) {
        // size_t sndPos = i + pow(2, k - 1) < n ? i + pow(2, k - 1) : n - 1;
        size_t sndPos = i + POW2[k - 1] < n ? i + POW2[k - 1] : n - 1;
        if (h.at(M[i][k - 1]) <= h.at(M[sndPos][k - 1]))
          M[i][k] = M[i][k - 1];
        else
          M[i][k] = M[sndPos][k - 1];
      }
    }
  } else { // mType == RMQ_TYPE::MAX
    // k = 1
    for (size_t i = 0; i < n; ++i) {
      M[i][0] = i;
      if (h.at(i) >= h.at((i + 1) < n ? i + 1 : n - 1)) {
        M[i][1] = i;
      } else {
        M[i][1] = (i + 1) < n ? i + 1 : n - 1;
      }
    }
    //k > 1
    for (size_t k = 2; k <= LOG2[n]; ++k) { // if you use floating log2 write <
      for (size_t i = 0; i < n; ++i) {
        // size_t sndPos = i + pow(2, k - 1) < n ? i + pow(2, k - 1) : n - 1;
        size_t sndPos = i + POW2[k - 1] < n ? i + POW2[k - 1] : n - 1;
        if (h.at(M[i][k - 1]) >= h.at(M[sndPos][k - 1]))
          M[i][k] = M[i][k - 1];
        else
          M[i][k] = M[sndPos][k - 1];
      }
    }
  }
}

/*******************************************************************************
 * query the RMQ data structures. Time complexity: O(1)
 * @tparam T typename of data
 * @param left Position of the first element in range
 * @param right The length of the first string v = m.
 * @return Index of the min or max of { h[left], h[left + 1], ..., h[right]}
 ******************************************************************************/
template<typename T>
size_t RMQ_nlogn<T>::query(size_t left, size_t right) {
  // assert(left <n && "Bad RMQ_ONN Query");
  // assert(right < n && "Bad RMQ_ONN Query");
  // assert(left <= right && "Bad RMQ_ONN Query");
  if (left > right)
    std::swap(left, right);
  if (n == 0)
    return std::numeric_limits<size_t>::max();
  if (mType == RMQ_TYPE::MIN && (left > n || right > n))
    return std::numeric_limits<size_t>::max();
  if (mType == RMQ_TYPE::MAX && (left > n || right > n))
    return std::numeric_limits<size_t>::max();

  // Answer in O(1) with (overlapping) windows of size log2(right-left)
  size_t &i = left;
  size_t &j = right;
  if (i == j) // need because log2(0) = -inf
    return i;
  size_t s = LOG2[j - i];
  if (mType == RMQ_TYPE::MIN)
    return h.at(M[i][s]) <= h.at(M[j - POW2[s] + 1][s])
           ? M[i][s]
           : M[j - POW2[s] + 1][s];
  // mType == RMQ_TYPE::MAX
  return h.at(M[i][s]) >= h.at(M[j - POW2[s] + 1][s])
         ? M[i][s]
         : M[j - POW2[s] + 1][s];
}

//== RMQ for plus minus one arrays =============================================

/*******************************************************************************
 * Constructor for a range minimum (maximum) query data structure instance
 * @tparam T typename of data
 * @param type type of queries (RMQ_TYPE::MIN or RMQ_TYPE::MAX)
 * @param arr data for which the queries are asked
 ******************************************************************************/
template<typename T>
RMQ_pm1<T>::RMQ_pm1(RMQ_TYPE type, const std::vector<T> &arr)
    : mType(type),
      h(arr),
      n(arr.size()),
      s(n > 3 ? log2(n) / 2 : 1),
      hp(std::vector<T>(n / s + (n % s != 0), 0)), // vector of size ceil(n/s)
      pos(std::vector<size_t>(n / s + (n % s != 0), 0)) {
  if (n < 4)return;
  preprocess();
  rmq = std::make_unique<RMQ_nlogn<T>>(mType, hp);
  precompute();
}

/*******************************************************************************
 * Helper used to fill hp and pos (minimus values of blocks and their positions)
 * @tparam T typename of data
 * @note Time complexity: O(n), where n is the length of the input data vector
 ******************************************************************************/
template<typename T>
void RMQ_pm1<T>::preprocess() {
  // Fill hp and pos
  for (size_t ip = 0; ip * s < n; ++ip) {
    size_t t = ip * s; // Reset temporary pos at the head of each block
    for (size_t i = ip * s; i < ip * s + s && i < n; ++i) {
      if ((mType == RMQ_TYPE::MIN && h[i] < h[t]) ||
          (mType == RMQ_TYPE::MAX && h[i] > h[t])) {
        t = i; // Found new minimum (maximum) in a block
      }
    }
    pos[ip] = t;
    hp[ip] = h[t];
  }
}

/*******************************************************************************
 * precompute lookup tables to answer RMQ queries for normalized arrays
 * @tparam T typename of data
 * @note There are O(2**s)= O(sqrt(n)) arrays of length s. Using RMQ_ONN this
 * adds up to a time complexity of: \f$\mathcal{O}(\sqrt(n)*s^2) = \mathcal{O}( \sqrt{n}(\log n)^2) )
 * \subseteq \mathcal{O}(n)\f$. Example of Normalisation: \f$ [3, 4, 5, 4, 3, 4]
 * \mapsto (2; [1,1,1,-1,-1,1])\f$  * with getNormBlockNum([1,1,1,-1,-1,1], 0, 5) = 0b111001
 ******************************************************************************/
template<typename T>
void RMQ_pm1<T>::precompute() {
  //BlockM.reserve(pow(2,s));
  BlockM.resize(pow(2, s));
  size_t total = 1 << s; // total number of vectors in {-1,1}^s
  for (int i = 0; i < static_cast<int>(total); i++) {
    // Convert i to binary and store in a vector
    std::bitset<MAXLOGS> bs(i); // assuming MAXLOGS is maximum s you will have
    std::vector<T> binaryVector(s);
    for (size_t j = 0; j < s; j++) {
      binaryVector[j] = bs[s - j - 1];
    }
    std::vector<size_t> nVec(s);
    nVec[0] = s; // initialize with s to avoid underflow when creating nVec
    for (size_t j = 1; j < s; j++) {
      if (binaryVector[j] == 1) {
        nVec[j] = nVec[j - 1] + 1;
      } else { // binaryVector[j]==0
        nVec[j] = nVec[j - 1] - 1;
      }
    }
    size_t nx = getNormBlockNum(nVec, 0, s - 1);
    BlockM[nx] = std::make_unique<RMQ_ONN<size_t >>(mType, nVec);
    //BlockM.push_back(std::make_unique<RMQ_ONN>(mType,nVec));
  }
  lookUpBlockNum.resize(h.size() / s + 1);
  for (size_t ip = 0; ip * s < h.size(); ++ip) {
    size_t is = ip * s;
    size_t ie = ip * s + s - 1;
    if (ie < h.size()) {
      lookUpBlockNum[ip] = getNormBlockNum(h, is, ie);
    } else {
      std::vector<T> vec;
      vec.reserve(s);
      for (size_t i = is; i < h.size(); ++i)
        vec.push_back(h[i]);
      while (vec.size() < s)
        vec.push_back(vec.back() + 1);
      lookUpBlockNum[ip] = getNormBlockNum(vec, 0, s - 1);
    }
  }
}

/*******************************************************************************
 * Helper: getNormBlockNum (convert to arr[left:right] to a number)
 * @tparam T typename of data
 * @param arr data that contains blocks
 * @param left start of the block
 * @param right end of the block
 * @return size_t r = 0b1 concatenate_{i=left}^{right-1} (arr[i+1]-arr[i+1]==1)
 * @note time complexity O(right-left)
 ******************************************************************************/
template<typename T>
template<typename U>
size_t RMQ_pm1<T>::getNormBlockNum(const std::vector<U> &arr,
                                   size_t left, size_t right) {
  // Check that a block beginning at left and ending at right has the right size
  assert(right - left + 1 == s && "Block Size matches [left:right]");
  assert(left <= right && "invalid arguments");

  size_t binary = 1;
  for (size_t i = left + 1; i <= right; i++) {
    binary <<= 1;
    // if the normalized value of a block at index i is 1, set the least significant bit of binary to 1
    if (arr[i] - arr[i - 1] == 1) { // normalizedArray[i] == 1
      binary |= 1;
    }
  }
  return binary; // Example: Binary( {1,0,1,1} ) = 2^3 + 2^1+ 2^0
}

/*******************************************************************************
 * query the RMQ data structures. Time complexity: O(1)
 * @tparam T typename of data
 * @param left first position in the range to search in
 * @param right last position in the range to search in
 * @return left most position of the minimum or maximum value
 ******************************************************************************/
template<typename T>
size_t RMQ_pm1<T>::query(size_t left, size_t right) {
  // assert(left <n && "Bad RMQ_ONN Query");
  // assert(right < n && "Bad RMQ_ONN Query");
  // assert(left <= right && "Bad RMQ_ONN Query");
  if (left > right)
    std::swap(left, right);
  if (mType == RMQ_TYPE::MIN && (left > n || right > n))
    return std::numeric_limits<size_t>::max();
  if (mType == RMQ_TYPE::MAX && (left > n || right > n))
    return std::numeric_limits<size_t>::max();

  if (n < 4) { // if division by s is not well-defined (s=log2(n)/2)
    if (n == 0) return std::numeric_limits<size_t>::max();
    if (n == 1) return 0;
    if (n == 2) {
      if (left == right) return left;
      else if (mType == RMQ_TYPE::MIN) return h[0] <= h[1] ? 0 : 1;
      else return h[0] >= h[1] ? 0 : 1; // return max pos
    }
    if (n == 3) {
      if (left == right) return left;
      else if (right - left == 1) {
        if (mType == RMQ_TYPE::MIN)
          return h[left] <= h[right] ? left : right; // min pos of two elements
        else
          return h[left] >= h[right] ? left : right; // max pos of two elements
      } else { // left=0 and right = 2
        if (mType == RMQ_TYPE::MIN) { // calc min
          if (h[0] <= h[1] && h[0] <= h[2]) return 0;
          else if (h[0] > h[1] && h[1] <= h[2]) return 1;
          else return 2;
        } else { // calc max
          if (h[0] >= h[1] && h[0] >= h[2]) return 0;
          else if (h[0] < h[1] && h[1] >= h[2]) return 1;
          else return 2;
        }
      }
    }
  }

  const size_t i = left, j = right;     // want to find RMQ[i:j]_h
  const size_t ip = i / s, jp = j / s;  // number blocks that contain pos i, j
  const size_t is = ip * s;             // first position_h in the ip-th block
  const size_t js = jp * s;             // first position_h in the jp-th block
  const size_t ie = is + s - 1;     // last position_h in block ip
  // const size_t je = js + s - 1;     // last position_h in block jp
  const size_t lBlock = lookUpBlockNum[ip]; //getNormBlockNum(h, is, ie);
  const size_t rBlock = lookUpBlockNum[jp]; //getNormBlockNum(h, js, je);
  size_t l, m, r, result;
  /*
   * l = RMQ[i:min(j,ie)]_h (first block)
   * m = RMQ[ie+ 1: js-1]_h = RMQ[ip:jp]_hp (corresponds to multiple blocks)
   * r = RMQ[max(js,i):j]_h (last block)
   */
  l = BlockM[lBlock]->query(i % s, j <= ie ? j % s : s - 1) + ip * s;
  m = (jp - ip > 1) ? pos[rmq->query(ip + 1, jp - 1)] : l;
  r = BlockM[rBlock]->query(i < js ? 0 : i % s, right % s) + jp * s;
  T val;
  if (mType == RMQ_TYPE::MIN)
    val = std::min({h[l], h[m], h[r]});
  else
    val = std::max({h[l], h[m], h[r]});
  result = (val == h[l]) ? l : ((val == h[m]) ? m : r);
  return result;
}


//== CartesianTree & Lowest Common Ancestor ====================================


template<typename T>
CartesianTree<T>::CartesianTree(RMQ_TYPE type, const std::vector<T> &arr)
    : mType(type), h(arr), n(arr.size()), tree(std::vector<Node<T>>(0)) {
  if (arr.empty()) return;
  buildTree();
  fillEulerArrays();
  rmq = new RMQ_pm1<size_t>(RMQ_TYPE::MIN, level); // Always search lowest level
  // rmq = std::make_unique<RMQ_pm1<size_t>>(RMQ_TYPE::MIN, level);
}
template<typename T>
Node<T> *CartesianTree<T>::getRoot() {
  return root;
}
template<typename T>
void CartesianTree<T>::buildTree() {
  tree.clear();
  tree.resize(n);

  if (h.empty()) {
    return;
  }
  std::stack<Node<T> *> nodeStack;
  for (size_t i = 0; i < n; ++i) {

    Node<T> *node = &tree.at(i);
    node->label = i;
    node->value = h[i];

    Node<T> *lastPopped = nullptr;
    if (mType == RMQ_TYPE::MIN) {
      while (!nodeStack.empty() && nodeStack.top()->value > h[i]) {
        lastPopped = nodeStack.top();
        nodeStack.pop();
      }
    } else { //mType==RMQ_TYPE::MAX
      while (!nodeStack.empty() && nodeStack.top()->value < h[i]) {
        lastPopped = nodeStack.top();
        nodeStack.pop();
      }
    }

    if (!nodeStack.empty()) {
      nodeStack.top()->right = node;
      node->parent = nodeStack.top();
    }

    if (lastPopped != nullptr) {
      node->left = lastPopped;
      lastPopped->parent = node;
    }

    nodeStack.push(node);
  }

  while (!nodeStack.empty()) {
    root = nodeStack.top();
    nodeStack.pop();
  }
}

template<typename T>
void CartesianTree<T>::dfs(Node<T> *node, size_t depth) {
  if (node == nullptr) {
    return;
  }

  // visit the node: record the node value and level
  euler.push_back(node->label);
  level.push_back(depth);
  // if(node->label < nodeRepresentative[node->value])   // firstInEuler[i] = argmin_{1<=j<2n} (euler[j]==i)
  //     firstInEuler[node->value] = node->label;

  // visit the left child
  if (node->left != nullptr) {
    dfs(node->left, depth + 1);
    euler.push_back(node->label);  // visit node again after returning from left subtree
    level.push_back(depth);
  }

  // visit the right child
  if (node->right != nullptr) {
    dfs(node->right, depth + 1);
    euler.push_back(node->label);  // visit node again after returning from right subtree
    level.push_back(depth);
  }
}

template<typename T>
void CartesianTree<T>::fillEulerArrays() {
  euler.clear();
  level.clear();

  dfs(root, 0); // fills euler & levels
  nodeRepresentative.clear();
  nodeRepresentative.resize(n, 2 * n);
  for (size_t j = 0; j < euler.size(); j++) {
    if (j < nodeRepresentative[euler[j]])
      nodeRepresentative[euler[j]] = j;
  }
}

template<typename T>
size_t CartesianTree<T>::getLCA(Node<T> *u, Node<T> *v) {
  return rmq->query(nodeRepresentative[u->value], nodeRepresentative[v->value]);
}

template<typename T>
size_t CartesianTree<T>::query(size_t i, size_t j) {
  // assert(left <n && "Bad RMQ_ONN Query");
  // assert(right < n && "Bad RMQ_ONN Query");
  // assert(left <= right && "Bad RMQ_ONN Query");
  if (i > j)
    std::swap(i, j);
  if (n == 0)
    return std::numeric_limits<size_t>::max();
  if (mType == RMQ_TYPE::MIN && (i > n || j > n))
    return std::numeric_limits<size_t>::max();
  if (mType == RMQ_TYPE::MAX && (i > n || j > n))
    return std::numeric_limits<size_t>::max();

  size_t ansInLevels;
  if (nodeRepresentative[i] > nodeRepresentative[j])
    ansInLevels = rmq->query(nodeRepresentative[j], nodeRepresentative[i]);
  else
    ansInLevels = rmq->query(nodeRepresentative[i], nodeRepresentative[j]);
  return euler[ansInLevels];
}


//== RMQ with O(n) preprocessing ===============================================

/*******************************************************************************
 * Constructor for a range minimum (maximum) query data structure instance
 * @tparam T typename of data
 * @param type type of queries (RMQ_TYPE::MIN or RMQ_TYPE::MAX)
 * @param arr data for which the queries are asked
 ******************************************************************************/
template<typename T>
RMQ_ON<T>::RMQ_ON(RMQ_TYPE type, const std::vector<T> &arr)
    :tree(type, arr) {

}

/*******************************************************************************
 * query the RMQ data structures. Time complexity: O(1)
 * @tparam T typename of data
 * @param i first index in the range to search in
 * @param j last index in the range to search in
 * @return  left most index of the minimum or maximum value in the data array
 ******************************************************************************/
template<typename T>
size_t RMQ_ON<T>::query(size_t i, size_t j) {
  return tree.query(i, j);
}

//=== Explicit instantiation ===================================================

template class RMQ_ONN<int>;
template class RMQ_nlogn<int>;
template class RMQ_pm1<int>;
template class Node<int>;
template class CartesianTree<int>;
template class RMQ_ON<int>;

template class RMQ_ONN<util::uint>;
template class RMQ_nlogn<util::uint>;
template class RMQ_pm1<util::uint>;
template class CartesianTree<util::uint>;
template class RMQ_ON<util::uint>;

} /*end of namespace RMQ*/

