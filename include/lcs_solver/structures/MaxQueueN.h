/*********************************************************************
 * @file MaxQueueN.h
 * @author Steinkopp:Felix
 * @brief Impl. of MaxQueue (for problems with symbol dependent gap lengths)
 * @note Based on https://www.mimuw.edu.pl/~kubica/publications/2007-cpm/vlcs.pdf
 ********************************************************************/
#ifndef LCS_SOLVER_STRUCTURES_MAX_QUEUE_H_
#define LCS_SOLVER_STRUCTURES_MAX_QUEUE_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <deque>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "util/RMQ.h"

namespace lcs_solver::structures {

/*******************************************************************************
 * MaxQueue is a kind of priority queue that provides the maximum of the last l
 * elements put into the queue (for a fixed l)
 * @tparam T typename of data in the queue
 ******************************************************************************/
template<typename T>
struct MaxQueue {
 public:
  using uint = std::size_t;
  using Pair = std::pair<uint, T>;

 private:
  std::deque<Pair> q; ///< Two-Linked Queue for Insertions Pair{index, value}
  const uint l;       ///< Context Length
  uint c{0};          ///< Counter: indexes consecutive insertions
  const T minValue;   ///< Value used to return from max for an empty Queue

 public:
  /*****************************************************************************
   * Constructor of MaxQueue objects
   * @param contextLen uint the capacity of the MaxQueue
   * @param min T used for returning the maximum of an empty queue
   ****************************************************************************/
  explicit MaxQueue(
      uint contextLen,
      T min = std::numeric_limits<T>::min()
  ) : q(), l(contextLen), minValue(min) {}

  /*****************************************************************************
   * @brief Checks if a MaxQueue object is empty
   * @return true iff no insertion happened in the past of the object
   ****************************************************************************/
  bool empty() { return q.empty() && c == 0; }

  /*****************************************************************************
   * @brief Getter for the Capacity of a MaxQueue object
   * @return uint `l`
   ****************************************************************************/
  [[nodiscard]] uint size() const { return l; }

  /*****************************************************************************
   * @brief Inserts a value into the object
   * @param val T the value to be inserted
   ****************************************************************************/
  void insert(T val) {
    /* Remove such pairs (i, value) from q with value <= val */
    while (!q.empty() && q.back().second <= val) q.pop_back();
    ++c;
    q.emplace_back(Pair(c, val));

    /* Remove such pairs (i, val), that i<=c-l */
    while (!q.empty() && q.front().first + l < c + 1) q.pop_front();
    assert(q.size() <= l);
    assert(!q.empty() || l == 0);
  }

  /*****************************************************************************
   * @brief Retrieves insertion with the largest Value (=head of the queue)
   * @return pair<uint, T> with maximum T value out of the last `l` insertions
   ****************************************************************************/
  Pair max() {
    if (q.empty()) return {std::numeric_limits<uint>::max(), minValue};
    return q.front();
  }

};

template<typename T>
struct MaxQueue2D {
  /*****************************************************************************
   * Using Definitions
   * Matrix M in this class has dimensions m x n
   * Tuple t = [l1,d1,l2,d2] ,where:
   * l1,l2 = length of window along axis 1 and along axis 2
   * d1,d2 = offset of the window when going over all indices in a Matrix M such
   *         that M[i-d1][M-d2] is the top left (first) element in the window
   ****************************************************************************/
  using uint = std::size_t;
  using Matrix = std::vector<std::vector<T>>;
  using RMQ = ::lcs_solver::util::RMQ_ON<T>;
  using Tuple = std::tuple<uint, uint, uint, uint>;

 private:
  Matrix &M;        ///< Ref to 2D data matrix
  const uint m;   ///< m = NumberOfRows(M)
  const uint n;   ///< n = NumberOfColumns(M) = NumberOfQueues
  const T minValue; ///< neutral element of `<` in the domain of T
  uint yBuffer;
  std::vector<MaxQueue<T>> Q;///< maximum of the last l2 elements put into cols
  std::vector<T> max;        ///< stores the max[f] of window[ _ : f]
  std::unique_ptr<RMQ> rmq;  ///< RMQ data structure build on max after each row

 public:
  /*****************************************************************************
   * @brief Constructor (without setting the capacity)
   * @param matrix over that the window is slided
   * @param min neutral element of `<`
   ****************************************************************************/
  explicit MaxQueue2D(Matrix &matrix, T min = std::numeric_limits<T>::min())
      : M(matrix),
        m(matrix.size()),
        n(matrix.empty() ? 0 : matrix[0].size()),
        minValue(min),
        yBuffer(1),
        Q(std::vector<MaxQueue<T>>()),
        max(std::vector<T>(n, min)) {

  }

  /******************************************************************************
   * @brief Constructor (with setting the capacity)
   * @param matrix over that the window is slided
   * @param yBufferLength buffer length of 1d MaxQueue
   * @param min neutral element of `<`
   *****************************************************************************/
  explicit MaxQueue2D(Matrix &matrix, const uint yBufferLength, T min)
      : MaxQueue2D(matrix, min) {
    setYBuffer(yBufferLength);
  }

//  /*****************************************************************************
//    * @brief Slides a Window of `t` over `M` and executes a llcs gap algorithm
//    * @details slide a window of size l1 x l2 over a matrix of size m x n
//    * relative to M[i,j] with i in[0,m-1] and j in [0:n-1]. After each movement
//    * the maximum of M[i-d1:i-d1+l1][j-d2:j-d2+l2] is found via a rmq structure
//    * If Phi(i,j) is one M[i][j] is updated to the found maximum plus one.
//    * @param Phi functional predicate that take a index {i,j} in M
//    * @param t [s0,e0,s1,e1] defines the window to M[i-s0:i-e0][j-s1:j-e1]
//    ***************************************************************************/
//  void slide(const std::function<bool(uint, uint)> &Phi, Tuple t) {
//    auto &[s0, e0, s1, e1] = t;
//    uint const n0 = s0 - e0 + 1; // number of elements in [i-s0:i-e0]
//    uint const n1 = s1 - e1 + 1; // number of elements in [i-s0:i-e0]
//    setYBuffer(std::max<uint>(1, n1 - 1));// set min capacity to 1 so that Q[f] can store a cell in the current row
//    for (uint i = 0; i < m; ++i) {
//      doLineSetup(i, n0, s0); // Generate RMQ for Max[a:b] over updated window
//      for (uint j = 0; j < n; ++j) {
//        doComputingSetup(i, j, t);  // insert M[i-e0][j-e1]
//        T currMax = minValue;
//        if (!(j < s1 && j < e1)) {
//          uint start = j - s1;
//          uint end = j - e1;
//          if (j < s1) start = 0; // j - s1 would be negative (avoid overflow)
//          if (j < e1) end = 0;  // j - e1 would be negative (avoid overflow)
//          uint maxPos = rmq->query(start, end);
//          currMax = Q.empty() ? minValue : Q[maxPos].max().second;
//        }
//
//        // Update Matrix if Phi(i,j)==1
//        if (static_cast<int>(Phi(i, j)) == 1)
//          M[i][j] = currMax + 1;
//        else if (i > 0 && j > 0)
//          M[i][j] = std::max(M[i - 1][j], M[i][j - 1]);
//
//      } // j (row loop end)
//    } // i (colum loop end)
//  }

  /*****************************************************************************
   * doLineSetup: Creates a rmq data structure for a new line
   * @note Sideeffect: When a new row M[i-e0][...] changes the window, the
   * vector max and the rmq data structure based on it are updated such that:
   * rmq(a,b) = argmax_{x in [0:n-1]} M[i-s0:i-e0][x] with a,b in [0:m-1]
   * @note the elements of the new row are not added to the MaxQueues
   * @param i the number of the reference row
   * @param s0 offset to specify the first in the window (M[i-s0][...])
   * @param e0 offset to specify the last in the window (M[i-e0][...])
   ****************************************************************************/
  void doLineSetup(const uint i, const uint s0, const uint e0) {
    const uint nElements = s0 - e0 + 1; // = (i-e0) - (i-s0) + 1
    for (uint f = 1; f < n; ++f) {
      if (nElements == 1) {
        if (inRange(i - e0, {0, m - 1}))
          max[f] = M[i - e0][f];
        else
          max[f] = 0;
      } else {
        const uint oldValue = Q[f].max().second; // max M[i-s0:i-(e0+1)][f]
        uint candidateValue = minValue;
        if (inRange(i - e0, {0, m - 1})) // avoids i - e0 < 0
          candidateValue = M[i - e0][f];
        max[f] = std::max(oldValue, candidateValue);
      }
    }
    rmq = std::make_unique<RMQ>(lcs_solver::util::RMQ_TYPE::MAX, max);
  }

  /*****************************************************************************
   * @brief Checks if val is in closed interval
   * @param val uint to be checked
   * @param interval std::pair<uint,uint> the interval to be used  in the check
   * @return true iff `val` is in the closed interval `interval`
   ****************************************************************************/
  static bool inRange(const uint val, const std::pair<uint, uint> &interval) {
    return interval.first <= val && val <= interval.second;
  }

  /*****************************************************************************
   * doComputingSetup inserts M[i - e0][ j - e1] into Q[j - e1]
   * @param i index along the first dimension of M
   * @param j index along the second dimension of M
   * @param t =[s0, e0, s1, e1] tuple to specify the window offset. The window
   * ist set relative to a point (i,j) to M[i-s0:i-s1][j-s1:j-e1]
   ****************************************************************************/
  void doComputingSetup(const uint i, const uint j, Tuple t) {
    if (auto &[s0, e0, s1, e1] = t; i >= e0 && j >= e1) {
      const uint f = j - e1;
      const T &value = M[i - e0][f];
      Q[f].insert(value);
    }
  }

  /*****************************************************************************
   * doComputingSetup inserts M[i - e0][ j - e1] into Q[j - e1]
   * @param i Index along the first dimension of M
   * @param j Index along the second dimension of M
   * @param e0 Offset to specify the last row in the window (M[i-e0][...])
   * @param  e1 Offset to specify the last colum in the window (M[...][j-e1])
   * @note t =[s0, e0, s1, e1] tuple to specify the window offset. The window
   * ist set relative to a point (i,j) to M[i-s0:i-s1][j-s1:j-e1]
   ****************************************************************************/
  void doComputingSetup(const uint i, const uint j, const uint e0, const uint e1) {
    if (i >= e0 && j >= e1) {
      const uint f = j - e1;
      const T &value = M[i - e0][f];
      Q[f].insert(value);
    }
  }

  /*****************************************************************************
   * setYBuffer initializes the 1D MaxQueue by setting their maximum capacity
   * and sets the yBuffer to the new bufferSize.
   * @param bufferSize maximum number of values stored in a MaxQueue (field l).
   * @note if bufferSize is equal to zero, the MaxQueues will always be empty
   ****************************************************************************/
  void setYBuffer(uint bufferSize) {
    yBuffer = bufferSize;
    Q = std::vector<MaxQueue<T>>(n, MaxQueue<T>(yBuffer));
  }

  /*****************************************************************************
   * @brief query finds the maximum in M[i-s0:i-s1][j-s1:j-e1] based on Q.max
   * @param i index along the first dimension of M
   * @param j index along the second dimension of M
   * @param t =[s0, e0, s1, e1] tuple to specify the window offset. The window
   * ist set relative to a point (i,j) to M[i-s0:i-s1][j-s1:j-e1]
   * @return Maximum Values in M[i-s0:i-s1][j-s1:j-e1] if the window is not
   * well-defined return minValue
   ****************************************************************************/
  T query(const uint i, const uint j, Tuple t) {
    const auto &[s0, e0, s1, e1] = t;
    T res = minValue;
    if (!(j < s1 && j < e1)) { // llcs = RMQ window(i,j)
      uint start = j - s1;
      uint end = j - e1;
      if (j < s1) start = 0; // j - s1 might be negative (avoid overflow)
      if (j < e1) end = 0;  // j - e1 might be negative (avoid overflow)
      uint maxPos = rmq->query(start, end);
      //res = Q.empty() ? minValue : Q[maxPos].max().second;
      res = max[maxPos];
    }
    return res;
  }
};

} // end of namespace
#endif // LCS_SOLVER_STRUCTURES_MAX_QUEUE_H_