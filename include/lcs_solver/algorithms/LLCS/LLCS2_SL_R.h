/******************************************************************************
 * @file LLCS2_SL_R.h
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Wrapper: Uses a LCS2_S-Algorithm to solve a LCS2_SL Problem
 * @details Time Complexity: max( template Algo, n+m)
 *****************************************************************************/
#ifndef LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SL_R_H_
#define LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SL_R_H_

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/BaseSolution.h"
#include "algorithms/LLCS/LLCS2_SIG_Algorithm.h"
#include "constraints/BaseConstraint.h"
#include "constraints/ConstraintType.h"
#include "constraints/local/Constraint_Sigma_R.h"
#include "util/CommonTypes.h"

namespace lcs_solver::algorithms::llcs {

template<typename T>
concept IsDerivedFromLLCS2_SIG_Algorithm  = std::is_base_of_v<LLCS2_SIG_Algorithm, T>;

/*******************************************************************************
 * concat
 * Helper to generate meaningful names via templates
 * @tparam N Length of first string
 * @tparam M Length of second string
 * @param a Pointer to first cstyle string
 * @param b Pointer to second cstyle string
 * @return std::array<char, M+N> a+b
 ******************************************************************************/
template<std::size_t N, std::size_t M>
static constexpr auto concat(const char(&a)[N], const char(&b)[M]) {
  std::array<char, N + M - 1> result{};
  for (std::size_t i = 0; i < N - 1; ++i)
    result[i] = a[i];
  for (std::size_t i = 0; i < M; ++i)
    result[N - 1 + i] = b[i];
  return result;
}

/*******************************************************************************
 * LLCS2_SL_R
 * @details Wrapper solves LLCS_SL by using LLCS2_SR. It reverses the problems
 *  strings and passes them together with right(symb) := left(symb) to Algo.
 * @tparam Algo LLCS2_SIG_Algorithm that is able to solve LLCS2_SR
 ******************************************************************************/
template<IsDerivedFromLLCS2_SIG_Algorithm Algo>
struct LLCS2_SL_R : public LLCS2_SIG_Algorithm {
public:
  static constexpr const char prefix[] = "LLCS2_SL_R_";
  static constexpr auto name_array = concat(prefix, Algo::name);
  static constexpr const char *name = name_array.data();

  /*****************************************************************************
   * LLCS2_SL_R - Constructor
   * @tparam Algo LLCS_SR algorithm
   * @param vec std::vector of shared points to the constant strings of a problem
   * @param map map<std::string, share_ptr<BaseConstraint> with the constraints
   ****************************************************************************/
  LLCS2_SL_R(const StrPtrVector &vec,const ConstraintMap &map)
      : LLCS2_SIG_Algorithm(AlgoType::LLCS2_SR_MQ, vec, map,
                            LLCS2_SIG_Algorithm::getSigLMap(
                                map,
                                constraints::ConstraintType::SIGMA_L),
                            LLCS2_SIG_Algorithm::getSigRMap(
                                map,
                                constraints::ConstraintType::SIGMA_L)),
        rs(reverseEachIn(vec)),
        rMap(genRightMap(map)),
        algo(rs, rMap) {

  };

  [[nodiscard]] bool isValid() const override { return algo.isValid(); }
  [[nodiscard]] std::string_view getName() const override { return name; }
  [[nodiscard]] std::string DebugString() const override { return algo.DebugString(); }
  [[nodiscard]] std::string_view getDescription() const override {
    return {"Uses " + std::string(algo.getName()) + " algorithm to solve a LLCS2_SL Problem"};
  };
  [[nodiscard]] std::unique_ptr<BaseSolution> query() override { return algo.query(); }
  void doPreprocessing() override {
    algo.setTracking(trackKeyPairs);
    algo.doPreprocessing();
    keyPairs = algo.getKeyPairs();
    setState(State::Preprocessed);
  }
  void reset(ResetLevel l) override { algo.reset(l); }

  // LLCS2 Algorithm Interface
  [[nodiscard]] const Matrix &getMatrix() const override {
    return algo.getMatrix();
  };
  [[nodiscard]] bool isExtensible(Pair a, Pair b, uint llcsOfA) const override {
    return algo.isExtensible(a, b, llcsOfA);
  }
  [[nodiscard]] Window getPrevRange(const Pair &pair, uint llcs) const override {
    return algo.getPrevRange(pair, llcs);
  }

 private:
  const StrPtrVector rs;
  ConstraintMap rMap;
  Algo algo;

  /*****************************************************************************
   * reverseEachIn
   * @param vec vector of pointers to input string
   * @return vector of pointers to reversed input strings
   ****************************************************************************/
  StrPtrVector reverseEachIn(const StrPtrVector &vec) {
    using ::lcs_solver::util::String;
    StrPtrVector result;
    for (const auto &p : vec) {
      String s = String(p->rbegin(), p->rend());
      result.push_back(std::make_shared<String>(s));
    }
    return result;
  }

  /*****************************************************************************
   * genRightMap
   * @return ConstraintMap m with m[Sigma_R] = left
   ****************************************************************************/
  ConstraintMap genRightMap(const ConstraintMap & /*map*/) {
    ConstraintMap m;
    using ::lcs_solver::constraints::ConstraintType;
    using ::lcs_solver::constraints::local::Constraint_Sigma_R;
    m[ConstraintType::SIGMA_R] = std::make_shared<Constraint_Sigma_R>(left);
    return m;
  }
};

} // end of namespace

#endif // LCS_SOLVER_ALGORITHMS_LLCS_LLCS2_SL_R_H_