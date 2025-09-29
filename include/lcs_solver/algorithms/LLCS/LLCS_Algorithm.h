#ifndef LCS_SOLVER_INCLUDE_ALGORITHMS_LLCS_LLCS_ALGORITHM_H_
#define LCS_SOLVER_INCLUDE_ALGORITHMS_LLCS_LLCS_ALGORITHM_H_

#include "algorithms/BaseAlgorithm.h"
#include "algorithms/AlgoType.h"

namespace lcs_solver::algorithms::llcs {

class LLCS_Algorithm : public BaseAlgorithm {
 public:
  // using uint = util::uint;
  // using Matrix = structures::Matrix<util::uint>;
  // using Point = std::vector<uint>;
  // using Pair = std::pair<uint, uint>;
  // using Window = std::vector<Pair>;

  ~LLCS_Algorithm() override = default;
  LLCS_Algorithm(AlgoType algo,
                 const StrPtrVector& vec,
                 const ConstraintMap& map);

 //  [[nodiscard]] virtual const Matrix& getMatrix() const = 0;
 //  [[nodiscard]] virtual bool IsExtensible(const Point& a,
 //                                          const Point& b,
 //                                          uint llcs_of_a) const = 0;
 //  [[nodiscard]] virtual Window GetPrevRange(const Pair& pair,
 //                                            uint llcs) const = 0;
 //
 //  void reset(ResetLevel lvl) override;
 //
 //  [[nodiscard]] bool IsMatched(const Point& p) const;
 //  [[nodiscard]] bool IsMatched(const Point& p, bool one_based) const;
 //  [[nodiscard]] std::vector<Pair> GenAllKeyPairs(bool one_based) const;
 //  [[nodiscard]] std::vector<Pair> GenAllKeyPairsNaive(bool one_based) const;
 //
 //  [[nodiscard]] static uint LenKPtM(const KeyPointMatrix& kpm);
 //  [[nodiscard]] static uint max(const Matrix& matrix);
 //
 //  void track(const Pair& pair, uint llcs);
 //  void track(const Pair&& pair, uint llcs);
 //  void setTracking(bool tracking);  ///< Vector of zero-based matching Indices
 //  [[nodiscard]] bool getTrackingFlag() const;
 //  [[nodiscard]] const KeyPointMatrix& getKeyPairs() const;
 //
 // protected:
 //  [[nodiscard]] static std::string toString(const Matrix& m,
 //                                            StringViewVector s,
 //                                            bool trim = true,
 //                                            const std::string& name = "",
 //                                            bool printStrings = true);
 //
 //  [[nodiscard]] static std::string toString(const KeyPointMatrix& kpm,
 //                                            bool trackKeyPairs);
 //
 //  const bool areIndicesOneBased;
 //  bool trackKeyPairs;
 //  KeyPointMatrix keyPairs;
};

}  // namespace lcs_solver::algorithms::llcs

#endif  // LCS_SOLVER_INCLUDE_ALGORITHMS_LLCS_LLCS_ALGORITHM_H_
