#ifndef LCS_SOLVER_ALGORITHMS_SOLUTIONS_VECTOR3DSOLUTION_H_
#define LCS_SOLVER_ALGORITHMS_SOLUTIONS_VECTOR3DSOLUTION_H_

#include <vector>
#include "algorithms/BaseSolution.h"
#include "algorithms/solutions/Points.h"
#include "algorithms/solutions/BaseIterator.h"
#include "algorithms/solutions/BaseCollector.h"

namespace lcs_solver::algorithms::solutions {

class Vector3DSolution final : public BaseCollector {
 public:
  using uint = size_t;
  using Matrix3D = std::vector<std::vector<std::vector<uint>>>;
  using Points = ::lcs_solver::algorithms::solutions::Points;
  using StrPtr = std::shared_ptr<const util::String>;
  using StrPtrVec = std::vector<StrPtr>;

  class VecIterator final : public BaseIterator {
   public:
    using VecIter = std::vector<Points>::const_iterator;
    explicit VecIterator(VecIter iter, VecIter end);
    void advance() override;
    [[nodiscard]] std::unique_ptr<BaseIterator> clone() const override;
    [[nodiscard]] std::string DebugString() const override;

   private:
    std::vector<Points>::const_iterator mCurrent, mEnd;
  };// VecIterator

  Vector3DSolution();
  explicit Vector3DSolution(const StrPtrVec &spv, const Matrix3D &matrix);
  explicit Vector3DSolution(const std::vector<Points> &vec);

  [[nodiscard]] AnyIterator begin() const override;
  [[nodiscard]] AnyIterator end() const override;
  [[nodiscard]] BaseSolution *clone() const override;

  template<typename... Args>
  void emplace_back(Args &&...args);
  void push_back(const Points &element);
  void push_back(const Points &&element);
  void modify(uint solIdx, uint lcsIdx, uint strIdx, uint pos, bool reverse);

 private:
  std::vector<Points> mSolutions;
};

}// namespace lcs_solver::algorithms::solutions
#endif //LCS_SOLVER_ALGORITHMS_SOLUTIONS_VECTOR3DSOLUTION_H_
