#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_EMBEDDINGSVector_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_EMBEDDINGSVector_H_

#include <span>

#include "algorithms/BaseSolution.h"

namespace lcs_solver::algorithms::solutions {

/**
 * @details Storage Class to represent an embedding
 * @details Stores:
 *          - vector for shared_ptr to strings associated with an embedding
 *          - A vector of marked positions (as std::vector<std::vector<uint>>)
 *          - Flag reversed_ to indicate that positions are stored in reverse
 */
class Points final : public BaseSolution {
 public:
  using uint = std::size_t;
  using Point = std::vector<uint>;
  using Matrix2D = std::vector<Point>;
  using StrPtr = std::shared_ptr<const util::String>;
  using StrPtrVec = std::vector<StrPtr>;

  Points();
  explicit Points(StrPtrVec spv, bool reversed = false);
  explicit Points(StrPtrVec spv, const std::span<Embedding> &embs, bool reversed = false);
  explicit Points(StrPtrVec spv, Matrix2D matrix, bool one_based = false, bool reversed = false);
  [[nodiscard]] SolutionType getType() const override;
  [[nodiscard]] bool isEqual(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessThan(const BaseSolution &rhs) const override;
  [[nodiscard]] bool isLessEqualThan(const BaseSolution &rhs) const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] BaseSolution *clone() const override;
  [[nodiscard]] bool empty() const override;
  [[nodiscard]] const Point& access(uint pos) const;
  [[nodiscard]] uint GetLLCS() const;
  [[nodiscard]] uint GetNumOfStrings() const;
  [[nodiscard]] const StrPtrVec &GetStrPtrVec() const;
  [[nodiscard]] Embedding GenEmbedding(uint str_idx, bool one_based) const;
  [[nodiscard]] const std::vector<uint> &operator[](uint i, bool reverse) const;
  [[nodiscard]] const std::vector<uint> &operator[](uint i) const;

  template<typename... Args>
  void emplace_back(Args &&...args);
  void push_back(const std::vector<uint> &point);
  void push_back(const std::vector<uint> &&point);
  void push_back(const std::pair<uint, uint> &point);
  void push_back(const std::pair<uint, uint> &&point);
  void mod(uint pos, const std::pair<uint, uint> &point, bool one_based);
  std::vector<uint>& at(uint pos, bool one_based=true);
  Matrix2D& data();
  void resize(uint size);
  void clear();
  void reverse();
  void set_embedding(const std::vector<std::pair<uint, uint>> &embedding);
  [[nodiscard]] size_t size() const;
  [[nodiscard]] const Point &back() const;
  [[nodiscard]] std::pair<uint, uint> BackPair() const;
  void pop_back();
  void modify(uint lcs_idx, uint str_idx, uint pos, bool reverse);
  friend std::ostream &operator<<(std::ostream &os, const Points &points);

 private:
  static inline bool smaller(const std::vector<uint> &a, const std::vector<uint> &b);
  static inline bool SmallerOrEqual(const std::vector<uint> &a, const std::vector<uint> &b);
  const StrPtrVec spv_;
  Matrix2D mat_;
  const bool reversed_ = false;
};

}// namespace lcs_solver::algorithms::solutions
#endif//LCS_SOLVER_ALGORITHMS_SOLUTIONS_EMBEDDINGSVector_H_
