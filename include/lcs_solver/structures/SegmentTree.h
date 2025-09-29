#ifndef LCS_SOLVER_STRUCTURES_SEGMENT_TREE_H_
#define LCS_SOLVER_STRUCTURES_SEGMENT_TREE_H_

#include <cstddef>
#include <tuple>
#include <vector>

namespace lcs_solver::structures {

/*******************************************************************************
 * MemoryOrder defines the memory layout of a tree. Let's say we have this tree:
 *         A
 *        / \
 *       B   C
 *      / \   \
 *     D  E    F
 * Then
 * Nodes:  A   B   C   D   E   F
 * Index:  0   1   2   3   4   6 with MemoryOrder::Array
 * Index:  0   1   2   3   4   5 with MemoryOrder::Euler
 *
 * The getter for the nodes' children indices is defined in SegmentTree::
 * getIdx(node, start, end) where [start:end] is the range that the node manages
 ******************************************************************************/
enum class MemoryOrder { Array, Euler };

/*******************************************************************************
 * SegmentTree class that uses lazy propagation
 * @tparam T type of the data in the Segment Tree
 * @details Works when (T, q) is a semi group & update functions are closed
 * under composition and distributes over q
 * @todo LazyPropagation is not implemented correctly.
 * - It needs to also push the update operation to children
 * - updates need to take into account the update operation!!!
 ******************************************************************************/
template<typename T>
class SegmentTree {
 public:
  using uint = std::size_t;
  using FunctionPtr = T (*)(T, T);

  SegmentTree(
      FunctionPtr q,
      uint range,
      MemoryOrder order = MemoryOrder::Euler,
      T init = T(),
      T neutralElement = T()
  );

  SegmentTree(
      FunctionPtr q,
      const std::vector<T> &values,
      MemoryOrder order = MemoryOrder::Euler,
      T neutralElement = T()
  );

  T query(uint index);
  T query(uint front, uint back);
  [[nodiscard]] uint size() const;
  void update(FunctionPtr u, uint index, T val);
  void update(FunctionPtr u, uint front, uint back, T val);

 private:
  using Triple = std::tuple<uint, uint, uint>;

  static uint calcTreeSize(uint s, MemoryOrder o);
  inline Triple getIdx(uint node, uint start, uint end);

  void build(const std::vector<T> &vec, uint node, uint start, uint end);
  void updateRange(FunctionPtr u, uint node, uint start, uint end, uint ul, uint ur, T val);
  T query(uint node, uint start, uint end, uint ql, uint qr);
  void propagateLazy(uint node, uint start, uint end);

  const FunctionPtr q;
  const uint nElements;
  const uint treeSize;
  const MemoryOrder traversal;
  std::vector<T> tree;
  std::vector<bool> lazy;  ///< lazy[i] <=> node(i) has to push its value to children
  const T neutralElement;
};

} // namespace lcs_solver::structures

#endif // LCS_SOLVER_STRUCTURES_SEGMENT_TREE_H_