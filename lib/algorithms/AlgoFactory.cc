/******************************************************************************
 * @file AlgoFactory.cc
 * @author Steinkopp:Felix
 * @version 2.0
 * @brief Provides a factory for generating algorithms. The main function is
 * `Create` it takes a given AlgoType t, a ConstraintMap map and vector
 *  of shared pointers. A child Algorithm can be specified with c. For example,
 *  LCS2_RT can be created with all LLCS2 Algorithms. Additionally, it provides
 *  static functions to get the descriptions and names of the implemented
 *  algorithms by looking up the AlgoType key.
 *****************************************************************************/
#include "algorithms/AlgoFactory.h"

#include <array>
#include <map>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/LCS/LCS2_RT.h"
#include "algorithms/LCS/LCS2_STD_S.h"
#include "algorithms/LLCS/LLCS2_MC.h"
#include "algorithms/LLCS/LLCS2_MC_1C.h"
#include "algorithms/LLCS/LLCS2_MC_INC.h"
#include "algorithms/LLCS/LLCS2_MC_INC_E.h"
#include "algorithms/LLCS/LLCS2_MC_O1_SYNC.h"
#include "algorithms/LLCS/LLCS2_SA_MQ.h"
#include "algorithms/LLCS/LLCS2_SA_RMQ.h"
#include "algorithms/LLCS/LLCS2_SL_R.h"
#include "algorithms/LLCS/LLCS2_SR_MQ.h"
#include "algorithms/LLCS/LLCS2_SR_RMQ.h"
#include "algorithms/LLCS/LLCS2_STD_FL.h"
#include "algorithms/LLCS/LLCS_STD_FL.h"
#include "util/InOutHelper.h"

namespace {

using lcs_solver::algorithms::llcs::LLCS2_MC;
using lcs_solver::algorithms::llcs::LLCS2_MC_1C;
using lcs_solver::algorithms::llcs::LLCS2_MC_INC;
using lcs_solver::algorithms::llcs::LLCS2_MC_INC_E;
using lcs_solver::algorithms::llcs::LLCS2_MC_O1_SYNC;
using lcs_solver::algorithms::llcs::LLCS2_SR_MQ;
using lcs_solver::algorithms::llcs::LLCS2_SR_RMQ;
using lcs_solver::algorithms::llcs::LLCS2_STD_FL;
using lcs_solver::algorithms::llcs::LLCS_STD_FL;
using LLCS2_SL_MQ = lcs_solver::algorithms::llcs::LLCS2_SL_R<LLCS2_SR_MQ>;
using LLCS2_SL_RMQ = lcs_solver::algorithms::llcs::LLCS2_SL_R<LLCS2_SR_RMQ>;
using lcs_solver::algorithms::lcs::LCS2_RT;
using lcs_solver::algorithms::lcs::LCS2_STD_S;
using lcs_solver::algorithms::llcs::LLCS2_SA_MQ;
using lcs_solver::algorithms::llcs::LLCS2_SA_RMQ;
}  // namespace

namespace lcs_solver::algorithms {
/*******************************************************************************
 * Create
 * @param t AlgoType identifier of the algorithm
 * @param spv std::vector of std::shared_ptr to the problem Strings
 * @param map std::map from ConstraintType to std::shared_ptr of the constraint
 * @param c AlgoType of a sub-algorithm. E.g., to generate a dp-table or a lcs
 * @return std::unique_ptr<BaseAlgorithm> to a new algorithm instantiation
 ******************************************************************************/
std::unique_ptr<BaseAlgorithm> AlgoFactory::Create(
    const AlgoType t,
    const BaseAlgorithm::StrPtrVector& spv,
    const BaseAlgorithm::ConstraintMap& map,
    AlgoType c) {
  switch (t) {
    case AlgoType::LLCS_STD_FL: return std::make_unique<LLCS_STD_FL>(spv, map);
    case AlgoType::LLCS2_STD_FL: return std::make_unique<LLCS2_STD_FL>(spv, map);
    case AlgoType::LLCS2_MC: return std::make_unique<LLCS2_MC>(spv, map);
    case AlgoType::LLCS2_MC_INC: return std::make_unique<LLCS2_MC_INC>(spv, map);
    case AlgoType::LLCS2_MC_INC_E: return std::make_unique<LLCS2_MC_INC_E>(spv, map);
    case AlgoType::LLCS2_MC_1C: return std::make_unique<LLCS2_MC_1C>(spv, map);
    case AlgoType::LLCS2_MC_O1_SYNC: return std::make_unique<LLCS2_MC_O1_SYNC>(spv, map);
    case AlgoType::LLCS2_SR_MQ: return std::make_unique<LLCS2_SR_MQ>(spv, map);
    case AlgoType::LLCS2_SR_RMQ: return std::make_unique<LLCS2_SR_RMQ>(spv, map);
    case AlgoType::LLCS2_SL_MQ: return std::make_unique<LLCS2_SL_MQ>(spv, map);
    case AlgoType::LLCS2_SL_RMQ: return std::make_unique<LLCS2_SL_RMQ>(spv, map);
    case AlgoType::LLCS2_SA_MQ: return std::make_unique<LLCS2_SA_MQ>(spv, map);
    case AlgoType::LLCS2_SA_RMQ: return std::make_unique<LLCS2_SA_RMQ>(spv, map);
    case AlgoType::LCS2_STD_S: return std::make_unique<LCS2_STD_S>(spv, map);
    case AlgoType::LCS2_RT: return std::make_unique<LCS2_RT>(spv, map, c);
    case AlgoType::Unknown:
    default: throw std::runtime_error("Algorithm Creation was not defined for a type");
  }
  // return nullptr;
}

/*******************************************************************************
 * Create
 * @details Expands the recommended algorithm from `GetLcsAlgoTypes` (in
 * ConstraintMap.h) and calls Create on it.
 * @param conf_pair Pair of AlgoType informs construction
 * @param spv std::vector of std::shared_ptr to the problem Strings
 * @param map std::map from ConstraintType to std::shared_ptr of the constraint
 * @return std::unique_ptr<BaseAlgorithm> to a new algorithm instantiation
 ******************************************************************************/
std::unique_ptr<BaseAlgorithm> AlgoFactory::Create(
    std::pair<AlgoType, AlgoType> conf_pair, const StrPtrVector &spv,
    const ConstraintMap &map) {
  const auto [main, child] = conf_pair;
  return Create(main, spv, map, child);
}

/*******************************************************************************
 * readAlgorithm
 * @param in std::istream for reading the problem in
 * @param os std::ostream for prompting inputs
 * @param options std::span<const AlgoType> for the possible AlgoTypes
 * @param spv std::vector<std::share_ptr<String>> of the problem Strings
 * @param map std::map from ConstraintType to std::shared_ptr of the constraint
 * @return std::unique_ptr<BaseAlgorithm> to the created BaseAlgorithm
 ******************************************************************************/
std::unique_ptr<BaseAlgorithm> AlgoFactory::ReadAlgorithm(
    const BaseAlgorithm::StrPtrVector &spv,
    const BaseAlgorithm::ConstraintMap &map,
    const AlgorithmSpan options,
    std::istream &in,
    std::ostream &os) {
  if (options.empty())
    return nullptr;
  // const auto type = ReadAlgoType(in, os);
  const auto type = ReadAlgoType(options, options.front(), in, os);
  auto base = AlgoType::Unknown;
  switch (type) {
    case AlgoType::LCS2_RT:{
      // Need information about which algorithm to use for reconstruction
      std::vector<AlgoType> next;
      for (const auto &llcs_available: GetLLCSAvailable()) {
        for (const auto &option: options) {
          if (option == llcs_available) {
            next.push_back(option);
          }
        }
      }
      if (next.empty()) {
        std::ostringstream oss;
        oss << "No base algorithm available for LCS2_RT of type ";
        oss << GetName(type);
        throw std::runtime_error(oss.str());
      }
      const std::string prompt = "Enter base algorithm for LCS2_RT: ";
      base = ReadAlgoType(next, next.front(), in, os, prompt);
      break;
    }
    default:
      break;
  }
  auto algorithm = Create(type, spv, map, base);
  return algorithm;
}

/*******************************************************************************
 * GetName
 * @param t AlgoType unique identifier of an implemented BaseAlgorithm
 * @return the name associated with the BaseAlgorithm of type t
 ******************************************************************************/
std::string_view AlgoFactory::GetName(const AlgoType t) {
  const auto it = GetMapAlgoTypeToName().find(t);
  if (it == GetMapAlgoTypeToName().end()) {
    throw std::runtime_error("Unknown algorithm type");
  }
  return it->second;
}

/*******************************************************************************
 * GetDescription
 * @param t AlgoType unique identifier of an implemented BaseAlgorithm
 * @return Description of Algorithm with type t
 * @todo complete switch
 ******************************************************************************/
std::string_view AlgoFactory::GetDescription(AlgoType t) {
  switch (t) {
    case AlgoType::Unknown:
      return "";
    case AlgoType::LLCS2_STD_FL:
      return LLCS2_STD_FL::description;
    // case AlgoType::LLCS2_MC:
    //   break;
    // case AlgoType::LLCS2_MC_INC:
    //   break;
    // case AlgoType::LLCS2_MC_INC_E:
    //   break;
    // case AlgoType::LLCS2_MC_1C:
    //   break;
    // case AlgoType::LLCS2_MC_O1_SYNC:
    //   break;
    // case AlgoType::LLCS2_SR_MQ:
    //   break;
    // case AlgoType::LLCS2_SR_RMQ:
    //   break;
    // case AlgoType::LLCS2_SL_MQ:
    //   break;
    // case AlgoType::LLCS2_SL_RMQ:
    //   break;
    // case AlgoType::LLCS2_SA_MQ:
    //   break;
    // case AlgoType::LLCS2_SA_RMQ:
    //   break;
    // case AlgoType::LLCS_STD_FL:
    //   break;
    // case AlgoType::LCS2_Stack:
    //   break;
    // case AlgoType::LCS2_RT:
    //   break;
    default:
      throw std::runtime_error(
          "AlgoFactory::GetDescription type not implemented t = " +
          std::string(GetName(t)));
  }
}

/*******************************************************************************
 * Helper that prompts via streams for an AlgoType
 * @param algorithm_span Span of a Container from which to choose an `AlgoType`
 * @param default_value Value in the algorithm_span. If `AlgoType::Unknown`
 *        (which is also the default argument), the head of algorithm_span is
 *        used.
 * @param in std::istream from which to read
 * @param os std::ostream reference for writing prompts and inform the user
 * @param msg std::string_view to change the prompt message
 * @return AlgoType (read from in)
 * @throws std::runtime_error if the algorithm_span is empy or an invalid input
 *  was read
 ******************************************************************************/
AlgoType AlgoFactory::ReadAlgoType(
    const AlgorithmSpan algorithm_span,
    AlgoType default_value,
    std::istream &in,
    std::ostream &os,
    const std::string_view msg) {
  if (algorithm_span.empty())
    throw std::runtime_error("No Algorithm available");
  if (default_value == AlgoType::Unknown) {
    default_value = algorithm_span.front();
  }

  // Print options
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << (!msg.empty() ? msg : "\nThe following algorithms are available:");
    os << "\n";
    uint i = 0;
    for (const AlgoType& t : algorithm_span) {
      i++;
      os << i << ") " << GetName(t) << "\n";
    }
    os  << "Enter the number or the name of the algorithm to use "
        << "(default=" << GetName(default_value) << "): ";
  }

  // Read algorithm identifier form istream and save it to name
  std::string input;
  std::getline(in, input);

  if (input.empty()) {
    return default_value;
  }

  // Convert identifier string to AlgoType
  auto type = AlgoType::Unknown;

  if (util::IsNumber(input)) {
    // Interpret identifier as a number. Cast: string => size_t => AlgoType
    size_t const num = std::stoull(input) - 1;
    if (num > algorithm_span.size()) {
      throw std::runtime_error(
          "AlgoFactory (ReadAlgoType): Input number was to large");
    }
    type = algorithm_span[num];
  } else {
    // Interpret input as a string
    const auto& look_up_name_map = GetMapNameToAlgoType();
    if (const auto it = look_up_name_map.find(input);
        it != look_up_name_map.end()) {
      type = it->second;
    } else {
      throw std::runtime_error(
          "AlgoFactory (ReadAlgoType): Name is not in defined.");
    }
  }
  return type;
}

/*******************************************************************************
 * Calls ReadAlgoType with the GetAvailable() AlgoSpan
 * @param in std::istream from which to read
 * @param os std::ostream pointer for writing prompts and inform the user
 * @param msg std::string_view to change the prompt message
 * @return AlgoType
 ******************************************************************************/
AlgoType AlgoFactory::ReadAlgoType(std::istream& in, std::ostream& os,
                                   const std::string_view msg) {
  return ReadAlgoType(GetAvailable(), AlgoType::LLCS_STD_FL, in, os, msg);
}

AlgoFactory::AlgorithmSpan AlgoFactory::GetAvailable() {
  static constexpr std::array available =
      std::to_array<AlgoType>({AlgoType::LLCS_STD_FL,
                               AlgoType::LLCS2_STD_FL,
                               AlgoType::LLCS2_MC,
                               AlgoType::LLCS2_MC_INC,
                               AlgoType::LLCS2_MC_INC_E,
                               AlgoType::LLCS2_MC_1C,
                               AlgoType::LLCS2_MC_O1_SYNC,
                               AlgoType::LLCS2_SR_MQ,
                               AlgoType::LLCS2_SR_RMQ,
                               AlgoType::LLCS2_SL_MQ,
                               AlgoType::LLCS2_SL_RMQ,
                               AlgoType::LLCS2_SA_MQ,
                               AlgoType::LLCS2_SA_RMQ,
                               AlgoType::LCS2_STD_S,
                               AlgoType::LCS2_RT});
  return available;
}
AlgoFactory::AlgorithmSpan AlgoFactory::GetLLCSAvailable() {
  static constexpr std::array available = {
      AlgoType::LLCS2_STD_FL,
      AlgoType::LLCS2_MC,
      AlgoType::LLCS2_MC_INC,
      AlgoType::LLCS2_MC_INC_E,
      AlgoType::LLCS2_MC_1C,
      AlgoType::LLCS2_MC_O1_SYNC,
      AlgoType::LLCS2_SR_MQ,
      AlgoType::LLCS2_SR_RMQ,
      AlgoType::LLCS2_SL_MQ,
      AlgoType::LLCS2_SL_RMQ,
      AlgoType::LLCS2_SA_MQ,
      AlgoType::LLCS2_SA_RMQ,
  };
  return {available};
}

/*******************************************************************************
 * Getter: GetMapNameToAlgoType
 * @return reference to a locally and statically created lookup map
 ******************************************************************************/
const std::map<std::string_view, AlgoType>&
AlgoFactory::GetMapNameToAlgoType() {
  static std::map<std::string_view, AlgoType> map;
  if (map.empty()) {
    for (const auto& [key, val] : GetMapAlgoTypeToName()) {
      map.insert({val, key});
    }
  }
  return map;
}

/*******************************************************************************
 * Getter: GetMapNameToAlgoType
 * @return Reference to the hard-coded static member `name_map_` in AlgoFactory
 ******************************************************************************/
const std::map<AlgoType, std::string_view>&
AlgoFactory::GetMapAlgoTypeToName() {
  static std::map<AlgoType, std::string_view> name_map = {
      {AlgoType::Unknown, "NULL"},
      {AlgoType::LLCS_STD_FL, LLCS_STD_FL::name},
      {AlgoType::LLCS2_STD_FL, LLCS2_STD_FL::name},
      {AlgoType::LLCS2_MC, LLCS2_MC::name},
      {AlgoType::LLCS2_MC_INC, LLCS2_MC_INC::name},
      {AlgoType::LLCS2_MC_INC_E, LLCS2_MC_INC_E::name},
      {AlgoType::LLCS2_MC_1C, LLCS2_MC_1C::name},
      {AlgoType::LLCS2_MC_O1_SYNC, LLCS2_MC_O1_SYNC::name},
      {AlgoType::LLCS2_SR_MQ, LLCS2_SR_MQ::name},
      {AlgoType::LLCS2_SR_RMQ, LLCS2_SR_RMQ::name},
      {AlgoType::LLCS2_SL_MQ, LLCS2_SL_MQ::name},
      {AlgoType::LLCS2_SL_RMQ, LLCS2_SL_RMQ::name},
      {AlgoType::LLCS2_SA_MQ, LLCS2_SA_MQ::name},
      {AlgoType::LLCS2_SA_RMQ, LLCS2_SA_RMQ::name},
      {AlgoType::LCS2_STD_S, LCS2_STD_S::name},
      {AlgoType::LCS2_RT, LCS2_RT::name},
  };
  return name_map;
}

/*******************************************************************************
 * Writes identifying of a BaseAlgorithm to an out stream
 * @param p Pointer to a BaseAlgorithm
 * @param os std::ostream pointer for writing prompts and inform the user
 ******************************************************************************/
void AlgoFactory::WriteAlgorithm(const BaseAlgorithm* p, std::ostream& os) {
  os << p->getName();
  if (const AlgoType t = p->getType(); t == AlgoType::LCS2_RT) {
    const auto algo_ptr = dynamic_cast<const LCS2_RT*>(p);
    os << algo_ptr->getBase()->getName();
  }
}

}  // namespace lcs_solver::algorithms