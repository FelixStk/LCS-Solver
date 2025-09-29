/*******************************************************************************
 * @file RndUtility.cc
 * @author Steinkopp:Felix
 * @version 1.2
 * @brief Impl. of various common function used in combinatorics problems
 ******************************************************************************/

#include "util/RndUtility.h"

#include <algorithm>
#include <cassert>
#include <set>
#include <stdexcept>

namespace lcs_solver::util {

// Initialize static variables
RndUtility::Matrix2D RndUtility::nCr = {{1}, {1, 1}, {1, 2, 1}}; ///< Initialize nCr dp-array with the base case
std::random_device RndUtility::rd;
std::mt19937 RndUtility::rng(RndUtility::rd());

/**
 * @brief Sets the seed of the static random number generator rng via std::random_device
 */
void RndUtility::setSeed() {
  rng.seed(rd());
}

/**
 * @brief Sets the seed of the static random number generator rng via a seed
 * @param seed const unsigned int that is used to set the seed
 */
void RndUtility::setSeed(const unsigned int seed) {
  rng.seed(seed);
}

/**
 * @brief Computes the binomial coefficient C(n, k), also known as "n choose k".
 *
 * The binomial coefficient C(n, k) represents the number of ways to choose k elements
 * from a set of n elements without regard to the order of selection. This function
 * uses dynamic programming and memoization to efficiently compute the result.
 *
 * @param n The total number of elements. Must be a non-negative integer.
 * @param k The number of elements to choose. Must be a non-negative integer and <= n.
 * @return The binomial coefficient C(n, k) as an unsigned integer.
 *         Returns 0 if k > n, and 1 if k == 0 or k == n.
 *
 * @note This function uses the symmetry property C(n, k) = C(n, n-k) to optimize
 *       computation by reducing the problem size when k > n/2.
 * @warning The function does not handle integer overflow. For large values of n and k,
 *          the result may exceed the maximum value representable by the return type.
 *
 * @see https://en.wikipedia.org/wiki/Binomial_coefficient for more details on binomial coefficients.
 */
RndUtility::uint RndUtility::binomialCoeff(uint n, uint k) {
  // Check trivial cases
  if (k > n) {
    return 0;
  }
  if (k == 0 || n == k) {
    return 1;
  }

  if (2 * k > n)//Uses symmetry property nCr(n,k) = nCr(n, n-k) to reduce the problem size
    k = n - k;

  // Resize if needed
  if (nCr.size() <= n) {
    const uint old_n = nCr.size() - 1;
    nCr.resize(n + 1);
    for (auto i = old_n + 1; i <= n; i++) {
      nCr[i].resize(i + 1, 0);
      for (uint j = 0; j <= i; j++) {
        if (j == 0 || j == i) {
          nCr[i][j] = 1;
        } else {
          nCr[i][j] = nCr[i - 1][j - 1] + nCr[i - 1][j];
        }
      }
    }
  }
  return nCr[n][k];
}

/**
 * @brief Randomly distributes indistinguishable balls into distinguishable bins.
 *
 * This function solves the "balls and bins" problem by generating a random distribution
 * of `nBalls` indistinguishable balls into `nBins` distinguishable bins. The solution
 * is based on the "stars and bars" theorem, which uses combinatorial methods to
 * randomly select dividers between balls.
 *
 * @param nBalls The number of indistinguishable balls to distribute. Must be >= 0.
 * @param nBins The number of distinguishable bins to distribute the balls into. Must be >= 1.
 * @return A vector of size `nBins`, where each element represents the number of balls
 *         in the corresponding bin. The sum of all elements in the vector equals `nBalls`.
 *
 * @note The function uses a uniform random distribution to select dividers between balls.
 *       The algorithm ensures that the distribution is uniformly random across all possible
 *       distributions.
 * @warning If `nBins` is 0, the behavior is undefined. If `nBalls` is 0, the function
 *          returns a vector of zeros with size `nBins`.
 *
 * @see https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics) for more details on the
 *      "stars and bars" theorem.
 */
std::vector<RndUtility::uint> RndUtility::ballsAndBins(const uint nBalls, const uint nBins) {
  // Handle edge cases
  if (nBins == 0) {
    throw std::invalid_argument("Number of bins must be greater than 0.");
  }
  if (nBalls == 0) {
    return std::vector<uint>(nBins, 0);
  }

  // Generate nBins - 1 unique dividers in the range [0, nBalls + nBins - 1]
  std::uniform_int_distribution<uint> dist(0, nBalls + nBins - 1);
  std::set<uint> dividers;
  while (dividers.size() + 1< nBins) {
    uint random_value = dist(rng);
    dividers.insert(random_value);
  }
  // std::vector<uint> dividers;
  // dividers.reserve(nBins - 1);
  // while (dividers.size() < nBins - 1) {
  //   if (uint random_value = dist(rng); std::ranges::find(dividers, random_value) == dividers.end()) {
  //     dividers.push_back(random_value);
  //   }
  // }
  // std::ranges::sort(dividers);

  // Calculate the number of balls in each bin
  std::vector<uint> k_values;
  k_values.reserve(nBins);
  uint prevDividerEnd = 0; // Points to the start of the first bin
  for (const auto &divider : dividers) {
    k_values.push_back(divider - prevDividerEnd);
    prevDividerEnd = divider + 1; // Move past the divider to the start of the next bin
  }
  k_values.push_back(nBalls + nBins - 1 - prevDividerEnd); // Add the last bin

  return k_values;
}

/**
 * @brief Generates a random combination without repetition using Fisher-Yates shuffle
 * @param nObjects Number of elements to select (must be less or equal than nSet)
 * @param nSet Cardinality of the set {0, 1, ..., nSet-1}
 * @throws std::invalid_argument if nObjects > nSet
 * @return Vector of nObjects unique random elements in random order
 * @note time and space complexity are of order \f$/mathcal{O}(nObjects)\f$. It
 *  works best if nObjects and nSet are of the same order of magnitude. So if
 *  nSet > 1e7 or nObjects << nSet, consider using an alternative approach with
 *  std::unordered_set
 *  @see https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 */
std::vector<RndUtility::uint> RndUtility::combNoRep(const uint nObjects, const uint nSet) {
  if (nObjects > nSet) {
    throw std::invalid_argument("nObjects cannot exceed nSet");
  }
  std::vector<uint> elements(nSet);
  std::iota(elements.begin(), elements.end(), 0);

  for (uint i = 0; i < nObjects; ++i) {
    std::uniform_int_distribution<uint> dist(i, nSet - 1);
    uint j = dist(rng);
    std::swap(elements[i], elements[j]);
  }
  elements.resize(nObjects);

  // Alternative reservoir sampling approach
  // std::unordered_set<uint> selected;
  // while (selected.size() < nObjects) {
  //   selected.insert(uniform_int_distribution(0, nSet - 1)(rng));
  // }
  // return {selected.begin(), selected.end()};

  return elements;
}

/**
 * @brief Generates unique random numbers within a closed interval [lower_bound, upper_bound]
 * @param nObjects Number of distinct values to select (must satisfy 1 <= nObjects <= upper_bound - lower_bound + 1)
 * @param interval Closed interval as pair [lower_bound, upper_bound] (requires upper_bound ≥ lower_bound)
 * @return Vector of unique random numbers in random order, uniformly distributed within [lower_bound, upper_bound]
 * @throws std::invalid_argument if input constraints are violated
 *
 * @note Delegates to combNoRep() for core algorithm, then adjusts values to target interval
 * @see combNoRep() for implementation details of the base random selection algorithm
 */
std::vector<RndUtility::uint> RndUtility::combNoRepInRange( const uint nObjects, std::pair<uint,uint> interval){
  const auto [lower, upper] = interval;
  if (upper < lower) {
    throw std::invalid_argument("Upper bound must be >= lower bound");
  }

  const auto nSet = upper - lower + 1;
  auto comb = combNoRep(nObjects, nSet); // Throws if nObjects > nSet
  for (auto& el : comb) {
    el += lower; // Offset to target interval
  }

  return comb;
}

/**
 * @brief Utility function for choosing randomly and uniformly a number from an interval
 * @param interval Closed interval as pair [lower_bound, upper_bound]
 * @return uint within [lower_bound, upper_bound]
 */
RndUtility::uint RndUtility::rndNumInRange(const std::pair<uint, uint> &interval) {
  return std::uniform_int_distribution(interval.first, interval.second)(rng);
}

/**
 * @brief Generates a random sub-interval within specified bounds using uniform distribution
 * @param interval The containing interval as [big_lower, big_upper] (must be valid: big_upper >= big_lower)
 * @param lengthBound Length constraints [min_len, max_len] where 1 <= min_len <= max_len <= (big_upper - big_lower + 1)
 * @return Random sub-interval [a,b] satisfying:
 *         - big_lower ≤ a ≤ b ≤ big_upper
 *         - min_len ≤ (b - a + 1) ≤ max_len
 * @throws std::invalid_argument for invalid input constraints
 * @note Time Complexity \f$/mathcal{O}(1)\f$
 */
RndUtility::Pair RndUtility::nrdPairInRange(const Pair &interval, const Pair &lengthBound) {
  const auto [big_lower, big_upper] = interval;
  const auto [min_len, max_len] = lengthBound;
  const auto super_length = big_upper - big_lower + 1;

  if(big_upper < big_lower) throw std::invalid_argument("Invalid super interval");
  if(min_len < 1 || max_len < min_len) throw std::invalid_argument("Invalid length bounds");
  if(max_len > super_length) throw std::invalid_argument("Max length exceeds super interval");

  std::uniform_int_distribution<uint> len_dist(min_len, max_len);
  const uint len = len_dist(rng);

  std::uniform_int_distribution<uint> start_dist(big_lower, big_upper - len + 1);
  const uint a = start_dist(rng);

  return {a, a + len - 1};
}

/**
 * @brief Calculates the count of numbers in a valid closed interval
 * @param interval Closed interval as pair [lower, upper] (requires lower <= upper)
 * @return The count of integers: upper - lower + 1
 * @throws std::invalid_argument if interval is invalid (upper < lower)
 */
RndUtility::uint RndUtility::numbersInRange(const Pair &interval) {
  const auto [lower, upper] = interval;
  if(upper < lower) throw std::invalid_argument("Invalid interval bounds");
  return  upper - lower + 1;
}

/**
 * @brief Generates a uniformly random sub-interval within a super interval
 * @details Uniform distribution across all possible sub-intervals of size ≥1
 * @param interval Valid parent interval [pl, pr] (requires pl <= pr)
 * @return Random sub-interval [a,b] where pl <= a <= b <= pr
 * @throws std::invalid_argument if parent interval is invalid
 * @note Time Complexity \f$/mathcal{O}(1)\f$
 */
RndUtility::Pair RndUtility::nrdPairInRange(const Pair &interval) {
  const auto [pl, pr] = interval;
  if(pr < pl) throw std::invalid_argument("Invalid parent interval");

  const uint n = pr - pl + 1;
  const uint totalSubintervals = n * (n + 1) / 2; // Triangular number

  std::uniform_int_distribution<uint> dist(0, totalSubintervals - 1);
  const uint k = dist(rng); // Represents the index of the subinterval when all are listed in lexicographic order.

  // Map k to (length, start) using combinatorial inversion: m*(m+1)/2 <= k <= (m+1)*(m+2)/2
  const uint m = static_cast<uint>((sqrt(static_cast<double>(8 * k + 1)) - 1) / 2); // Solution of m^2+m-2k = 0
  const uint remaining = k - m*(m + 1)/2;
  const uint length = m + 1;
  const uint start = pl + remaining;

  return {start, start + length - 1};
}

}// namespace lcs_solver::util