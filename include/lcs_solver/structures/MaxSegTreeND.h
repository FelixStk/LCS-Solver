#ifndef LCS_SOLVER_STRUCTURES_MAXSEGTREEND_H_
#define LCS_SOLVER_STRUCTURES_MAXSEGTREEND_H_

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <span>
#include <utility>
#include <vector>

#include "structures/Matrix.h"
#include "util/CommonTypes.h"

namespace lcs_solver::structures {

class MaxSegTreeND {
 public:
  using uint = util::uint;
  using Pair = std::pair<size_t, size_t>;
  using Matrix = lcs_solver::structures::Matrix<uint>;
  using FuncSet = std::span<std::pair<std::span<const size_t>, size_t>>;

  explicit MaxSegTreeND(std::span<const size_t> d, FuncSet list = {});
  explicit MaxSegTreeND(const Matrix &mat);
  uint operator[](std::span<size_t> indices);
  uint query(std::span<size_t> index);
  uint query(std::span<const Pair> ranges);
  void update(std::span<const size_t> index, uint value);
  void update(std::span<const Pair> index, uint value);
  [[nodiscard]] uint size() const;

  const std::vector<size_t> dim;

 private:
  using Function = const uint& (*)(const uint&, const uint&);
  using Triple = std::tuple<uint, uint, uint>;
  static constexpr Function q = std::max<uint>;
  static constexpr Function u = std::max<uint>;
  static constexpr uint neutralElement = 0;

  static inline size_t calcTreeSize(size_t d);
  static std::vector<size_t> calcTreeSize(std::span<const size_t> d);
  static inline Triple getIdx(uint node, uint start, uint end);
  void propagateLazy(uint node, uint d, uint start, uint end);


  [[nodiscard]] bool isLastDim(uint d) const;
  [[nodiscard]] static bool isSubsetEq(const Pair &a, const Pair &b) ;
  void updateTree(uint node, uint value);

  void translate(uint dest, uint d, const std::vector<Pair>& child, const std::vector<bool> &ranged);
  [[nodiscard]] bool stopTranslateEarly(uint dest, uint d, const std::vector<bool> &ranged) const;
  [[nodiscard]] std::vector<size_t> firstPath(uint dest, uint d, const std::vector<Pair> &child, const std::vector<bool> &ranged) const;
  bool incTranslationPath(std::vector<size_t>& path, uint source, uint d, const std::vector<Pair> &child, const std::vector<bool> &ranged) const;

  Matrix tree;
  Matrix lazy;
  std::vector<bool> treeChange;

};

}  // namespace lcs_solver::structures
#endif  // LCS_SOLVER_STRUCTURES_MAXSEGTREEND_H_