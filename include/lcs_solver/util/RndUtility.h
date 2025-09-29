#ifndef LCS_SOLVER_RND_UTILITY_HPP
#define LCS_SOLVER_RND_UTILITY_HPP

#include "CommonTypes.h"
#include <vector>
#include <random>
namespace lcs_solver::util {
class RndUtility {
 public:
  using uint = std::size_t;
  using Pair = std::pair<uint, uint>;
  static void setSeed();
  static void setSeed(unsigned int seed);

  static uint binomialCoeff(uint n, uint k);
  static std::vector<uint> ballsAndBins(uint nBalls, uint nBins);
  static std::vector<uint> combNoRep(uint nObjects, uint nSet);
  static std::vector<uint> combNoRepInRange(uint nObjects,
                                                Pair interval);
  static uint rndNumInRange(const Pair &interval);
  static Pair nrdPairInRange(const Pair &interval, const Pair &lengthBound);
  static Pair nrdPairInRange(const Pair &interval);
  static uint numbersInRange(const Pair &interval);

 private:
  using Matrix2D = std::vector<std::vector<uint>>;
  static Matrix2D nCr;
  static std::random_device rd;  ///< Uniformly-distributed non-deterministic integer random number generator
  static std::mt19937 rng;       ///< Efficient pseudo-random number generator (32bit Mersenne twister algorithm)
};

}
#endif //LCS_SOLVER_RND_UTILITY_HPP
