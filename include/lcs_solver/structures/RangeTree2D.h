#ifndef LCS_SOLVER_STRUCTURES_RANGETREE2D_H_
#define LCS_SOLVER_STRUCTURES_RANGETREE2D_H_

#include <algorithm>
#include <functional>
#include <memory>
#include <queue>
#include <span>
#include <vector>

namespace lcs_solver::structures {

template<typename T1, typename T2>
struct RangeTree2D {
 public:
  using Range1 = std::pair<T1, T1>;
  using Range2 = std::pair<T2, T2>;
  using Point = std::pair<T1, T2>;
  using PointType = std::pair<T1, T2>;

  struct Element {
    T1 x;
    T2 y;
    const Element *lc = nullptr;// Left child fractional cascading pointer
    const Element *rc = nullptr;// Right child fractional cascading pointer
    // explicit Element(T2 y) : y(y) {}
    Element(T1 x, T2 y) : x(x), y(y), lc(nullptr), rc(nullptr) {}
    Element(const T1 x, const T2 y, const Element *lc, const Element *rc)
        : x(x), y(y), lc(lc), rc(rc) {}
    [[nodiscard]] std::string DebugString() const;
  };

  struct Node {
    T1 x;
    std::vector<Element> arr;// Fractional cascading array
    std::unique_ptr<Node> lchild = nullptr;
    std::unique_ptr<Node> rchild = nullptr;
    explicit Node(T1 x) : x(x) {}
    Node(T1 x, std::vector<Element> arr,
         std::unique_ptr<Node> lchild,
         std::unique_ptr<Node> rchild)
        : x(x), arr(std::move(arr)), lchild(std::move(lchild)), rchild(std::move(rchild)) {}

    [[nodiscard]] std::string DebugString() const;
  };
  using NodePtr = std::unique_ptr<Node>;

  class Result {
   public:
    using ElementIter = typename std::vector<Element>::const_iterator;
    using NodeVector = std::vector<std::pair<const Node *, const Element *>>;
    using NodeIter = typename NodeVector::const_iterator;
    using Predicate = std::function<bool(const PointType &)>;

    Result() = default;
    explicit Result(const NodeVector &nodes, const Range2 &r2);
    Result(const Result& other);
    Result(Result&& other) noexcept;

    Result &operator=(const Result &other);
    Result &operator=(Result &&other) noexcept;

    class Iterator {
     public:
      Iterator() = default;
      Iterator(NodeIter start, NodeIter end, const Range2 &r2);
      Iterator(const Iterator &other);
      Iterator(Iterator &&other) noexcept;
      Iterator &operator=(const Iterator &other);
      Iterator &operator=(Iterator &&other) noexcept;

      Iterator &operator++();
      bool operator!=(const Iterator &other) const;
      bool operator==(const Iterator &other) const;
      PointType operator*() const;
      [[nodiscard]] std::string DebugString() const;

     private:
      void advance();
      bool setNextValidPosition();
      static bool validPosition(T2 y, const Range2 &r2);

      NodeIter currentPairIter, endPairIter;
      const Element *pCurrentElem, *pLastElem;
      Range2 r2;
    };// class Iterator

    Iterator cbegin() const;
    Iterator cend() const;
    // These allow range-based for loops
    auto begin() const { return cbegin(); }
    auto end() const { return cend(); }

   private:
    NodeVector nodes;
    Range2 r2;
  };// class Result

  explicit RangeTree2D(std::span<PointType const> points);
  Result query(Range1 d1, Range2 d2) const;
  [[nodiscard]] std::string DebugString() const;

 private:
  // NodePtr buildTree(const std::span<PointType> &points, size_t start, size_t end);
  NodePtr buildTree(const std::vector<T1> &xs, const std::vector<std::vector<T2>> &ys, size_t start, size_t end);
  void mergeElements(Node *node);
  Node *findSplitNode(Range1 r1) const;
  const Element *findFirstElem(const Node *node, const Range2 &r2) const;
  typename Result::NodeVector collectNodesInXRange(Range1 r1, Range2 r2) const;
  static bool isLeaf(const Node *p);
  static bool isValueIn(T1 x, Range1 r);

  NodePtr root;
};
}// namespace lcs_solver::structures

#include "structures/RangeTree2D.tpp"

#endif// LCS_SOLVER_STRUCTURES_RANGETREE2D_H_