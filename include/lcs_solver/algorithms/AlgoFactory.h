#ifndef LCS_SOLVER_ALGORITHMS_ALGOFACTORY_H_
#define LCS_SOLVER_ALGORITHMS_ALGOFACTORY_H_

#include <cstddef>// size_t
#include <map>
#include <memory>
#include <span>
#include <string_view>
#include <iostream>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "util/StrPtrVector.h"

namespace lcs_solver::algorithms {

class AlgoFactory {
public:
  using uint = size_t;
  using AlgorithmPtr = std::unique_ptr<BaseAlgorithm>;
  using AlgorithmSpan = std::span<const AlgoType>;
  using StrPtrVector = util::StrPtrVector;
  using ConstraintMap = constraints::ConstraintMap;

  virtual ~AlgoFactory()= default;

  static std::unique_ptr<BaseAlgorithm> Create(
      AlgoType t,
      const StrPtrVector &spv,
      const ConstraintMap &map,
      AlgoType c = AlgoType::Unknown
  );

  static std::unique_ptr<BaseAlgorithm> Create(
      std::pair<AlgoType, AlgoType> conf_pair,
      const StrPtrVector &spv,
      const ConstraintMap &map
  );

  static std::unique_ptr<BaseAlgorithm> ReadAlgorithm(
      const StrPtrVector &spv,
      const ConstraintMap &,
      AlgorithmSpan options = GetAvailable(),
      std::istream &in=std::cin,
      std::ostream &os=std::cout
  );

  [[nodiscard]] static AlgoType ReadAlgoType(std::istream &in=std::cin, std::ostream &os=std::cout, std::string_view msg = "");
  [[nodiscard]] static AlgoType ReadAlgoType(AlgorithmSpan algorithm_span, AlgoType default_value = AlgoType::Unknown, std::istream &in=std::cin,std::ostream &os=std::cout, std::string_view msg = "");
  [[nodiscard]] static std::string_view GetName(AlgoType t);
  [[nodiscard]] static std::string_view GetDescription(AlgoType t);
  [[nodiscard]] static AlgorithmSpan GetAvailable();
  [[nodiscard]] static AlgorithmSpan GetLLCSAvailable();
  [[nodiscard]] static const std::map<std::string_view, AlgoType> &GetMapNameToAlgoType();
  [[nodiscard]] static const std::map<AlgoType, std::string_view> &GetMapAlgoTypeToName();
  static void WriteAlgorithm(const BaseAlgorithm *p, std::ostream &os);

  // private:
  //   static constexpr const auto available_ = std::to_array<AlgoType>(
  //       {AlgoType::LLCS_STD_FL,
  //        AlgoType::LLCS2_STD_FL,
  //        AlgoType::LLCS2_MC,
  //        AlgoType::LLCS2_MC_INC,
  //        AlgoType::LLCS2_MC_1C,
  //        AlgoType::LLCS2_MC_O1_SYNC,
  //        AlgoType::LLCS2_SR_MQ,
  //        AlgoType::LLCS2_SR_RMQ,
  //        AlgoType::LLCS2_SL_MQ,
  //        AlgoType::LLCS2_SL_RMQ,
  //        AlgoType::LLCS2_SA_MQ,
  //        AlgoType::LLCS2_SA_RMQ,
  //        AlgoType::LCS2_Stack,
  //        AlgoType::LCS2_RT});
  // };
};
}// namespace lcs_solver::algorithms
#endif//LCS_SOLVER_ALGORITHMS_ALGOFACTORY_H_
