#ifndef LCS_SOLVER_UTIL_RMQ_H_
#define LCS_SOLVER_UTIL_RMQ_H_

#include "Logger.hpp"
#include<vector>
#include <memory>
#include <tuple>
#include <unordered_map>
#include "util/MetaTableCalc.h"

// The LCA Problem Revisited: https://link.springer.com/chapter/10.1007/10719839_9

namespace lcs_solver::util {

/*******************************************************************************
 * Types of Ranged Queries
 ******************************************************************************/
enum class RMQ_TYPE {
  MIN = 0,
  MAX
};

/*******************************************************************************
 * Range Maximum or Minimum. Complexity: <O(n^2),O(1)>
 * @tparam T typename of data
 ******************************************************************************/
template<typename T>
class RMQ_ONN {
 public:
  RMQ_ONN(RMQ_TYPE type, const std::vector<T> &arr);
  size_t query(size_t left, size_t right);
 private:
  using Matrix = std::vector<std::vector<T>>;
  RMQ_TYPE mType = RMQ_TYPE::MIN;
  const std::vector<T> &h;
  const size_t n;
  Matrix M;
  void preprocess();
};

/*******************************************************************************
 * Range Maximum or Minimum Query with complexity: <O(n*log(n)) , O(1)>
 * @tparam T typename of data
 ******************************************************************************/
template<typename T>
class RMQ_nlogn {
  using Matrix = std::vector<std::vector<size_t>>;

  RMQ_TYPE mType;                  ///< type of queries
  const std::vector<T> &h;         ///< reference to vector with data
  const size_t n;                  ///< size of data vector
  Log2_table<size_t> LOG2;         ///< LOG2[x] = floor(log2(x)) for x >0
  Pow2_table<size_t> POW2;         ///< POW[x] = 2**x for x=>0
  Matrix M;                        ///< dp table

//  static constexpr size_t MAX_INDEX = 256;
//  static constexpr size_t MAX_POWER = 10;
//  static constexpr auto LOG2 = Meta_Log2_table<MAX_INDEX>(); ///< [log2(1), log2(2), ..., log2(MAX_INDEX)] rounded down to integer value
//  static constexpr auto POW2 = Meta_Pow2_table<MAX_POWER>(); ///< [2**0, 2**1, ..., 2**MAX_POWER]

 public:
  RMQ_nlogn(RMQ_TYPE type, const std::vector<T> &arr);
  size_t query(size_t left, size_t right);

 private:
  void preprocess();
};

/*******************************************************************************
 * RMQ for arrays with pm1 property and complexity: <O(n), O(1)>
 * @note pm1 property: neighboring data elements in the array differ by +-1
 * @tparam T typename of data (needs to be number like)
 ******************************************************************************/
template<typename T>
class RMQ_pm1 {
  using Matrix = std::vector<std::vector<size_t>>;

  static constexpr int MAXLOGS = 8 * sizeof(int); // length of int bitmask >= s

  const RMQ_TYPE mType;     ///< Type of Query (MIN or MAX)
  const std::vector<T> &h;  ///< Reference to data vector
  const size_t n;           ///< Number of elements in data vector
  const size_t s;           ///< Length of a block s := log2(n)/2 <=> n<4^s
  std::vector<T> hp;        ///< hp[i] is the minimum element in the i-th block defined by a division of h by s
  std::vector<size_t> pos;  ///< pos[i] is the minimum element in the i-th block, relative to the number in h, where hp[i] is found
  Matrix M;                 ///< dp table with Indices
  std::unique_ptr<RMQ_nlogn<T>> rmq;
  std::vector<std::unique_ptr<RMQ_ONN<size_t>>> BlockM;
  std::vector<size_t> lookUpBlockNum;
 public:
  RMQ_pm1(RMQ_TYPE type, const std::vector<T> &arr);
  size_t query(size_t left, size_t right);

 private:
  void preprocess();
  void precompute();
  template<typename U>
  size_t getNormBlockNum(const std::vector<U> &arr, size_t left, size_t right);
};

/*******************************************************************************
 * Node in CartesianTree
 * @tparam T datatype of the node value
 ******************************************************************************/
template<typename T>
class Node {
 public:
  T value;
  size_t label;
  Node *left;
  Node *right;
  Node *parent;
  Node() : value(0), label(0), left(nullptr), right(nullptr), parent(nullptr) {}
  Node(size_t label, size_t val)
      : value(val),
        label(label),
        left(nullptr),
        right(nullptr),
        parent(nullptr) {}
  Node(size_t label, T val, Node *l, Node *r, Node *p)
      : value(val), label(label), left(l), right(r), parent(p) {}
};

/*******************************************************************************
 * RMQ for arrays via solving the lowest common ancestor problem.
 * @tparam T typename of data
 * @note Time Complexity: <O(n), O(1)>
 ******************************************************************************/
template<typename T>
class CartesianTree {
  const RMQ_TYPE mType;
  const std::vector<T> &h;
  const size_t n;

  Node<T> *root;
  std::vector<Node<T>> tree;  // n nodes for h
  std::vector<size_t> val;    // n nodes for h
  std::vector<size_t> label;  // n nodes for h
  std::vector<size_t> level;  ///< Levels of the nodes in the euler array. Root has level 0. Its children have level 1.
  std::vector<size_t> euler;  ///< Values of nodes visited in an euler tour, length of euler tour = 2 * edges + 1 = 2*(n-1)+1.
  //std::unordered_map<Number,Index> firstInEuler; // firstInEuler[i] = argmin{E[j]=i} for all 0<=j<=2n-1
  std::vector<size_t> nodeRepresentative;
  //  std::unique_ptr<RMQ_pm1<size_t>> rmq;
  RMQ_pm1<size_t> *rmq = nullptr;

 public:
  CartesianTree(RMQ_TYPE type, const std::vector<T> &arr);
  ~CartesianTree() {
    delete rmq;
  }
  Node<T> *getRoot();
  size_t getLCA(Node<T> *u, Node<T> *v);
  size_t query(size_t i, size_t j);

 private:
  void buildTree();
  void fillEulerArrays();
  void dfs(Node<T> *node, size_t depth);
};

/*******************************************************************************
 * RMQ_ON depends on all classes in this file!
 * @tparam T typename of the data
 ******************************************************************************/
template<typename T>
class RMQ_ON {
  CartesianTree<T> tree;
 public:
  RMQ_ON(RMQ_TYPE type, const std::vector<T> &arr);
  size_t query(size_t i, size_t j);
};

}
#endif  // LCS_SOLVER_UTIL_RMQ_H_