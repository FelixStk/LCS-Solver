#ifndef LCS_SOLVER_ALGORITHMS_BASESOLUTION_H_
#define LCS_SOLVER_ALGORITHMS_BASESOLUTION_H_

#include "structures/Embedding.h"
#include "util/CommonTypes.h"
#include <ostream>

//=== BaseSolution =============================================================
namespace lcs_solver::algorithms {
using ::lcs_solver::util::uint;
using ::lcs_solver::structures::Embedding;
class BaseSolution {
 public:
  virtual ~BaseSolution() = default;

  // used to define an order of SolutionTypes
  enum class SolutionType : unsigned char {
    Empty = 0x01,
    Unsigned = 0x02,
    Points = 0x03,
    Collector = 0x04,
    UFT8String = 0x05,
    DiffOutput = 0x06
  };

  [[nodiscard]] virtual bool empty() const = 0;
  bool operator==(const BaseSolution &rhs) const;
  bool operator!=(const BaseSolution &rhs) const;
  bool operator<(const BaseSolution &rhs) const;
  bool operator>(const BaseSolution &rhs) const;
  bool operator<=(const BaseSolution &rhs) const;
  bool operator>=(const BaseSolution &rhs) const;

  [[nodiscard]] virtual SolutionType getType() const = 0;
  [[nodiscard]] virtual std::string DebugString() const =0;
  void print() const;
  void print(std::ostream &os) const;
  [[nodiscard]] virtual BaseSolution * clone() const = 0;

  friend std::ostream &operator<<(std::ostream &os, const BaseSolution &sol){
    sol.print(os);
    return os;
  }

 private:
  [[nodiscard]] virtual bool isEqual(const BaseSolution &rhs) const = 0;
  [[nodiscard]] virtual bool isLessThan(const BaseSolution &rhs) const = 0;
  [[nodiscard]] virtual bool isLessEqualThan(const BaseSolution &rhs) const = 0;

};

}
#endif //LCS_SOLVER_ALGORITHMS_BASESOLUTION_H_
