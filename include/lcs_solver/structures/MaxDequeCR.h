#ifndef LCS_SOLVER_STRUCTURES_MAXDEQUECR_H_
#define LCS_SOLVER_STRUCTURES_MAXDEQUECR_H_

#include <cstddef>
#include <deque>
#include <functional>
#include <string>
#include <util/ParamGenerator.h>
#include <utility>
#include <vector>

namespace lcs_solver::structures {
/*******************************************************************************
 * @brief MaxDequeCR allow the calculate of the max Phi(i,j) where Phi(i,j)
 * := {M[ip,jp] | ip + Tl(M[ip,jp]) + 1 <= i <= ip + Tu(M[ip,jp] + 1 and
 *                jp + Tl(M[ip,jp]) + 1 <= j <= ip + Tu(M[ip,jp] + 1 and
 *                with ip in [i] and jp in [j] }
 * @details Tl and Tu are predefined getter for the gap bounds and M stores the
 * llcs of an mc inc problem.
 * Operations of MaxDequeCR maintain the following invariant properties
 * C[f] = [(t_{1,f}, S[i_{1,f},f]), ... , (t_{e,f},S[i_{e,f},f]) ]
 * then,
 * - i_{1,f} < i_{2,f} < ... < i_{e,f}
 * - S[i_{1,f},f] > ... > S[i_{e,f},f]
 * - S[i_{g,f},f] > M[h,f] for all i_{g-1,f} < h < i_{g,f} for g in [2:e]
 * - t_{1,f} >= i
 * similarly let:
 * R[r] = [ (t_{1,r}, S[j_{1,r},f] ), ... , (t_{e,r}, S[j_{e,f},f] )]
 * then,
 * - j_{1,r} < j_{2,r} < ... < j_{e,r}
 * - S[i_{1,r},j_{1,r}] > ... > S[i_{e,r},j{e,r}]
 * - S[r,j_{g,r}] > S[r,h] for all i_{g-1,r} < h < j_{g,r} for g in [2:e]
 * - t_{1,r} >= j
 ******************************************************************************/
class MaxDequeCR {
 public:
  using uint = util::uint;
  using Pair = std::pair<uint, uint>;
  using Matrix = std::vector<std::vector<uint>>;

  MaxDequeCR(
      const Matrix& mat,
      std::function<uint(uint, uint)> Phi,
      std::function<uint(uint)> Tl,
      std::function<uint(uint)> Tu
  );

  void update(uint i, uint j);
  void insert(uint i, uint j);
  uint extractMax(uint i, uint j);
  std::string DebugString();

 private:
  using DequeItem = std::pair<uint, uint *>;
  using Deques = std::vector<std::deque<DequeItem>>;

  const Matrix & M;
  const uint m; //Number of rows in S
  const uint n; //Number of cols in S
  const std::function<uint(uint, uint)> Phi;
  const std::function<uint(uint)> Tl;
  const std::function<uint(uint)> Tu;

  Matrix S;
  Deques C;
  Deques R;
};
}  // namespace lcs_solver::structures

#endif //LCS_SOLVER_STRUCTURES_MAXDEQUECR_H_
