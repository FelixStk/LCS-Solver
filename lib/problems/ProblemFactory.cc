/******************************************************************************
 * @file ProblemFactory.cc
 * @author Steinkopp:Felix
 * @version 2.0
 * @brief Implementation of the Problem Factory
 *****************************************************************************/

#include "problems/ProblemFactory.h"

#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <stdexcept>

#include "algorithms/AlgoFactory.h"
#include "algorithms/LCS/LCS2_RT.h"
#include "constraints/ConstraintFactory.h"
#include "constraints/local/Constraint_MC.h"
#include "constraints/local/Constraint_Sigma.h"
#include "problems/AdamsonEtAl/LCS2_MC.h"
#include "problems/AdamsonEtAl/LCS2_MC_1C.h"
#include "problems/AdamsonEtAl/LCS2_MC_INC.h"
#include "problems/AdamsonEtAl/LCS2_MC_O1C.h"
#include "problems/AdamsonEtAl/LCS2_MC_O1C_SYNC.h"
#include "problems/AdamsonEtAl/LCS2_SIGMA.h"
#include "problems/AdamsonEtAl/LCS2_SIGMA_L.h"
#include "problems/AdamsonEtAl/LCS2_SIGMA_R.h"
#include "problems/BaseProblem.h"
#include "problems/LCS_Classic.h"
#include "problems/ProblemType.h"
#include "util/InOutHelper.h"

namespace lcs_solver::problems {

using adamson::LCS2_MC;
using adamson::LCS2_MC_1C;
using adamson::LCS2_MC_INC;
using adamson::LCS2_MC_O1C;
using adamson::LCS2_MC_O1C_SYNC;
using adamson::LCS2_Sigma;
using adamson::LCS2_Sigma_L;
using adamson::LCS2_Sigma_R;
using algorithms::AlgoFactory;
using algorithms::AlgoType;
using constraints::ConstraintFactory;
using constraints::ConstraintType;

// ----------------------------------------------------------------------------
// Factory Method Implementations
// ----------------------------------------------------------------------------

/**
 * @copydoc ProblemFactory::CreateFromDialog()
 *
 * @implementation
 * Implements the interactive creation flow through these steps:
 * 1. Problem type selection via ReadProblemType()
 * 2. Metadata collection using ReadDescription()
 * 3. String input handling with ReadStrPtrVec()
 * 4. Constraint configuration via ConstraintFactory
 * 5. Algorithm selection through AlgoFactory
 *
 * @note Uses standard cin/cout streams by default. For file-based input,
 *       use FromFile() instead.
 */
ProblemFactory::ProblemPtr ProblemFactory::CreateFromDialog() {
  // Define interaction streams
  auto& in = std::cin;
  auto& out = std::cout;

  // Aks user for name, description, strings and constraints of the problem
  const ProblemType t = ReadProblemType(in, out);
  auto [name, description] = ReadDescription(t);
  StrPtrVec spv;
  ConstraintMap map;
  if (t == ProblemType::LCS_Base) {
    spv = ReadStrPtrVec(in, out);
    ConstraintFactory::ConstraintSpan all = ConstraintFactory::GetAvailable();
    map = ConstraintFactory::ReadConstraintMap(all, spv, nullptr, in, out);
  }
  else {
    spv = util::ReadStrPtrVec(2,"", nullptr, in, out);
  }

  ProblemPtr ptr = Create(t,
      std::move(name),
      std::move(description),
      std::move(spv),
      std::move(map));

  // Add required Constraints
  const ConstraintMap& map_ref = ptr->GetConstraints();
  const StrPtrVec& spv_ref = ptr->GetStrPtrVector();
  for (const auto& constraint_type : ptr->GetReqConstraints()) {
    if (!map_ref.contains(constraint_type)) {
      auto constraint_ptr = ConstraintFactory::ReadConstraint(spv_ref, constraint_type, nullptr, in, out);
      ptr->AddConstraint(constraint_ptr);
    }
  }

  // Aks user which algorithm to use for solving
  const AlgoFactory::AlgorithmSpan opts = ptr->GetAlgoSpan();
  auto algo_ptr = AlgoFactory::ReadAlgorithm(spv_ref, map_ref, opts, in, out);
  ptr->SetAlgorithm(std::move(algo_ptr));
  return ptr;
}

ProblemFactory::ProblemPtr ProblemFactory::CreateFromJson(const Path& p) {
  std::unique_ptr<BaseProblem> ptr = std::make_unique<BaseProblem>();
  nlohmann::json j;
  if (std::ifstream in(p); !in)
    throw std::runtime_error("Could not open file " + p.string() + " as input");
  else {
    j = nlohmann::json::parse(in);
  }
  from_json(j, ptr.get());
  return ptr;
}

/**
 * @copydoc ProblemFactory::Create(ProblemType,span<const String>,AlgoType,bool)
 *
 * @implementation
 * Handles constraint auto-generation based on the problem type:
 * - For MC variants: Generates relaxed gap constraints using GenRelaxedGapVector()
 * - For Sigma variants: Calculates alphabet and generates sigma constraints
 */
ProblemFactory::ProblemPtr ProblemFactory::Create(
    const ProblemType t,
    std::string&& name,
    std::string&& description,
    StrPtrVec&& spv,
    ConstraintMap&& map) {
  ProblemPtr ptr;
  switch (t) {
    case ProblemType::LCS_Base:
      ptr = std::make_unique<BaseProblem>(name, description, std::move(spv), std::move(map));
      // throw std::runtime_error(
      //     "Construction of LCS_Base Object is not allowed.");
      break;
    case ProblemType::LCS_Classic:
      ptr = std::make_unique<LCS_Classic>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_MC:
      ptr = std::make_unique<LCS2_MC>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_MC_INC:
      ptr = std::make_unique<LCS2_MC_INC>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_MC_1C:
      ptr = std::make_unique<LCS2_MC_1C>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_MC_O1C_SYNC:
      ptr = std::make_unique<LCS2_MC_O1C_SYNC>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_Sigma_L:
      ptr = std::make_unique<LCS2_Sigma_L>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_Sigma_R:
      ptr = std::make_unique<LCS2_Sigma_R>(
          name, description, std::move(spv), std::move(map));
      break;
    case ProblemType::LCS_Sigma:
      ptr = std::make_unique<LCS2_Sigma>(
          name, description, std::move(spv), std::move(map));
      break;
  }
  return ptr;
}

ProblemFactory::ProblemPtr ProblemFactory::Create(const ProblemType t,
    const std::span<const String> s,
    const AlgoType at,
    bool calc_llcs) {
  GapVector gc;
  SymbolVector alph;
  switch (t) {
    case ProblemType::LCS_Classic:
      return std::make_unique<LCS_Classic>(s, at, calc_llcs);
    case ProblemType::LCS_MC:
    case ProblemType::LCS_MC_INC:
    case ProblemType::LCS_MC_1C:
    case ProblemType::LCS_MC_O1C_SYNC:
      gc = constraints::local::Constraint_MC::GenRelaxedGapVector(s);
      return Create(t, s, gc, at, calc_llcs);
    case ProblemType::LCS_Sigma_L:
    case ProblemType::LCS_Sigma_R:
      alph = util::CalcAlphabet(s);
      gc = constraints::local::Constraint_Sigma::GenRelaxedGapVector(alph, s);
      return Create(t, s, alph, gc, at, calc_llcs);
    case ProblemType::LCS_Sigma:
      alph = util::CalcAlphabet(s);
      gc = constraints::local::Constraint_Sigma::GenRelaxedGapVector(alph, s);
      return Create(t, s, alph, gc, gc, at, calc_llcs);
    default:
      return nullptr;
  }
}

ProblemFactory::ProblemPtr ProblemFactory::Create(const ProblemType t,
    const std::span<const String> s,
    const GapSpan gap_tuples,
    const AlgoType at,
    const bool calc_llcs) {
  // Create a vector of string Pointers
  StrPtrVec spv;
  spv.reserve(s.size());
  for (const auto& str : s) {
    spv.push_back(std::make_shared<util::String>(str));
  }

  ProblemPtr ptr;

  const GapVector gc(gap_tuples.begin(),gap_tuples.end());
  const std::string_view name = GetName(t);
  const std::string_view description = GetDescription(t);
  std::pair<AlgoType, AlgoType> conf_pair;
  switch (t) {
    case ProblemType::LCS_Classic: {
      ConstraintMap map;
      ptr = std::move(std::make_unique<LCS_Classic>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_MC: {
      constexpr auto ct = ConstraintType::MC;
      ConstraintMap map;
      map[ct] = ConstraintFactory::Create(ct, gc);
      ptr = std::move(std::make_unique<LCS2_MC>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_MC_INC: {
      constexpr auto ct = ConstraintType::MC_INC;
      ConstraintMap map;
      map[ct] = ConstraintFactory::Create(ct, gc);
      ptr = std::move(std::make_unique<LCS2_MC_INC>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_MC_1C: {
      constexpr auto ct = ConstraintType::MC_1C;
      ConstraintMap map;
      map[ct] = ConstraintFactory::Create(ct, gc);
      ptr = std::move(std::make_unique<LCS2_MC_1C>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_MC_O1C_SYNC: {
      constexpr auto ct = ConstraintType::MC_O1C_SYNC;
      ConstraintMap map;
      map[ct] = ConstraintFactory::Create(ct, gc);
      ptr = std::move(std::make_unique<LCS2_MC_O1C_SYNC>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    default: {
      return nullptr;
    }
  }

  // Add algorithm
  if (!calc_llcs) {
    conf_pair = {AlgoType::LCS2_RT, at};
  } else {
    conf_pair = {at, AlgoType::Unknown};
  }
  AlgorithmPtr algorithm_ptr = AlgoFactory::Create(
      conf_pair, ptr->GetStrPtrVector(), ptr->GetConstraints());
  ptr->SetAlgorithm(std::move(algorithm_ptr));
  return ptr;
}

ProblemFactory::ProblemPtr ProblemFactory::Create(const ProblemType t,
    const std::span<const String> s,
    const SymbolSpan alph,
    const GapSpan gc,
    const AlgoType at,
    const bool calc_llcs) {
  // Create A String container
  StrPtrVec spv;
  spv.reserve(s.size());
  for (const auto& str : s) {
    spv.push_back(std::make_shared<util::String>(str));
  }

  // Create Problem
  ProblemPtr ptr;
  const std::string_view name = GetName(t);
  const std::string_view description = GetDescription(t);
  std::pair<AlgoType, AlgoType> conf_pair;
  switch (t) {
    case ProblemType::LCS_Classic: {
      break;
    }
    case ProblemType::LCS_Sigma_L: {
      constexpr auto ct = ConstraintType::SIGMA_L;
      ConstraintMap map;
      using constraints::local::Constraint_Sigma;
      const auto gc_map = Constraint_Sigma::GenSigmaTupleMap(alph, gc);
      map[ct] = ConstraintFactory::Create(ct, gc_map);
      ptr = std::move(std::make_unique<LCS2_Sigma_L>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_Sigma_R: {
      constexpr auto ct = ConstraintType::SIGMA_R;
      ConstraintMap map;
      using constraints::local::Constraint_Sigma;
      const auto gc_map = Constraint_Sigma::GenSigmaTupleMap(alph, gc);
      map[ct] = ConstraintFactory::Create(ct, gc_map);
      ptr = std::move(std::make_unique<LCS2_Sigma_R>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    default:
      return nullptr;
  }

  // Add algorithm
  if (!calc_llcs) {
    conf_pair = {AlgoType::LCS2_RT, at};
  } else {
    conf_pair = {at, AlgoType::Unknown};
  }
  AlgorithmPtr algorithm_ptr = AlgoFactory::Create(
      conf_pair, ptr->GetStrPtrVector(), ptr->GetConstraints());
  ptr->SetAlgorithm(std::move(algorithm_ptr));
  return ptr;
}

ProblemFactory::ProblemPtr ProblemFactory::Create(
    const ProblemType t,
    const std::span<const String> s,
    const SymbolSpan alph,
    const GapSpan l,
    const GapSpan r,
    AlgoType at,
    bool calc_llcs) {
  // Create A String container
  StrPtrVec spv;
  spv.reserve(s.size());
  for (const auto& str : s) {
    spv.push_back(std::make_shared<util::String>(str));
  }

  // Create Problem
  ProblemPtr ptr;
  const std::string_view name = GetName(t);
  const std::string_view description = GetDescription(t);
  std::pair<AlgoType, AlgoType> conf_pair;
  switch (t) {
    case ProblemType::LCS_Classic: {
      break;
    }
    case ProblemType::LCS_Sigma_L: {
      constexpr auto ct = ConstraintType::SIGMA_L;
      ConstraintMap map;
      using constraints::local::Constraint_Sigma;
      const auto gc_map = Constraint_Sigma::GenSigmaTupleMap(alph, l);
      map[ct] = ConstraintFactory::Create(ct, gc_map);
      ptr = std::move(std::make_unique<LCS2_Sigma_L>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_Sigma_R: {
      constexpr auto ct = ConstraintType::SIGMA_R;
      ConstraintMap map;
      using constraints::local::Constraint_Sigma;
      const auto gc_map = Constraint_Sigma::GenSigmaTupleMap(alph, r);
      map[ct] = ConstraintFactory::Create(ct, gc_map);
      ptr = std::move(std::make_unique<LCS2_Sigma_R>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    case ProblemType::LCS_Sigma: {
      constexpr auto ct = ConstraintType::SIGMA;
      ConstraintMap map;
      using constraints::local::Constraint_Sigma;
      const auto gc_l_map = Constraint_Sigma::GenSigmaTupleMap(alph, l);
      const auto gc_r_map = Constraint_Sigma::GenSigmaTupleMap(alph, r);
      map[ct] = ConstraintFactory::Create(ct, gc_l_map, gc_r_map);
      ptr = std::move(std::make_unique<LCS2_Sigma>(
          name, description, std::move(spv), std::move(map)));
      break;
    }
    default:
      return nullptr;
  }

  // Add algorithm
  if (!calc_llcs) {
    conf_pair = {AlgoType::LCS2_RT, at};
  } else {
    conf_pair = {at, AlgoType::Unknown};
  }
  AlgorithmPtr algorithm_ptr = AlgoFactory::Create(
      conf_pair, ptr->GetStrPtrVector(), ptr->GetConstraints());
  ptr->SetAlgorithm(std::move(algorithm_ptr));
  return ptr;
}
ProblemFactory::ProblemPtr ProblemFactory::Create(
    const ProblemType t,
    const std::span<const String> s,
    const GapSpan gc,
    const SymbolSpan alph,
    const GapSpan l,
    const GapSpan r,
    const AlgoType at,
    const bool calc_llcs) {
  ProblemPtr ptr;
  switch (t) {
    case ProblemType::LCS_Base: {
      return nullptr;
    }
    case ProblemType::LCS_Classic:{
      return Create(t, s, at, calc_llcs);
    }
    case ProblemType::LCS_MC:
    case ProblemType::LCS_MC_INC:
    case ProblemType::LCS_MC_1C:
    case ProblemType::LCS_MC_O1C_SYNC:{
      if (gc.empty()) {
        Create(t,s,at,calc_llcs);
      }
      return Create(t, s, gc, at, calc_llcs);
    }
    case ProblemType::LCS_Sigma_L: {
      if (l.empty()) {
        return Create(t, s, at, calc_llcs);
      }
      return Create(t, s, alph, l, at, calc_llcs);
    }
    case ProblemType::LCS_Sigma_R:{
      if (r.empty()) {
        return Create(t, s, at, calc_llcs);
      }
      return Create(t, s, alph, r, at, calc_llcs);
    }
    case ProblemType::LCS_Sigma: {
      if (r.empty() && l.empty()) {
        return Create(t, s, at, calc_llcs);
      }
      return Create(t, s, alph, l, r, at, calc_llcs);
    }
    default: {
      std::cerr << "ProblemFactory::Create: Unknown Problem Type" << GetName(t)
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

// -----------------------------------------------------------------------------
// Helper Implementations and Meta Data Getter
// -----------------------------------------------------------------------------

ProblemType ProblemFactory::ReadProblemType(
    std::istream& in, std::ostream& os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    uint n = 0;
    os << "The available problem types are:\n";
    for (const auto& type : kAvailable) {
      n++;
      os << n << ") " << GetName(type) << "\n";
    }
    os << "Enter a number or name of the problem type (default=1):";
  }
  std::string input;
  std::getline(in, input);
  if (input.empty()) return kAvailable[0];
  if (util::IsNumber(input)) {
    const uint idx = std::stoul(input) - 1;
    if (idx >= kAvailable.size())
      throw std::runtime_error(
          "ConstraintFactory: readConstraintType input number was out of "
          "range");
    return kAvailable[idx];
  }
  if (const auto it = GetNameTypeMap().find(input);
      it != GetNameTypeMap().end()) {
    return it->second;
  }
  throw std::runtime_error(
      "ConstraintFactory: readConstraintType type is not available");
}

std::pair<std::string, std::string> ProblemFactory::ReadDescription(
    const ProblemType t, std::istream& in, std::ostream& os) {
  std::string name, description;
  if (t == ProblemType::LCS_Base || os.rdbuf() == std::cout.rdbuf()) {
    // always read in name and description from a file
    name = util::ReadStdString("Enter a name for the problem", "", in, os);
    description = util::ReadStdString("Enter a description", "", in, os);
    return std::make_pair(name, description);
  }
  name = std::string(GetName(t).data());
  description = std::string(GetDescription(t).data());
  return std::make_pair(name, description);
}

ProblemFactory::StrPtrVec ProblemFactory::ReadStrPtrVec(
    std::istream& in, std::ostream& os) {
  return util::ReadStrPtrVec("", nullptr, in, os);
}

std::string_view ProblemFactory::GetName(const ProblemType type) {
  switch (type) {
    case ProblemType::LCS_Base:
      return "General Problem";
    case ProblemType::LCS_Classic:
      return LCS_Classic::kName;
    case ProblemType::LCS_MC:
      return LCS2_MC::kName;
    case ProblemType::LCS_MC_INC:
      return LCS2_MC_INC::kName;
    case ProblemType::LCS_MC_1C:
      return LCS2_MC_1C::kName;
    case ProblemType::LCS_MC_O1C_SYNC:
      return LCS2_MC_O1C_SYNC::kName;
    case ProblemType::LCS_Sigma_L:
      return LCS2_Sigma_L::kName;
    case ProblemType::LCS_Sigma_R:
      return LCS2_Sigma_R::kName;
    case ProblemType::LCS_Sigma:
      return LCS2_Sigma::kName;
  }
  return "";
}

std::string_view ProblemFactory::GetDescription(const ProblemType type) {
  switch (type) {
    case ProblemType::LCS_Base:
      return "Problem without suitability checks";
    case ProblemType::LCS_Classic:
      return LCS_Classic::kWhat;
    case ProblemType::LCS_MC:
      return LCS2_MC::kWhat;
    case ProblemType::LCS_MC_INC:
      return LCS2_MC_INC::kWhat;
    case ProblemType::LCS_MC_1C:
      return LCS2_MC_1C::kWhat;
    case ProblemType::LCS_MC_O1C_SYNC:
      return LCS2_MC_O1C_SYNC::kWhat;
    case ProblemType::LCS_Sigma_L:
      return LCS2_Sigma_L::kWhat;
    case ProblemType::LCS_Sigma_R:
      return LCS2_Sigma_R::kWhat;
    case ProblemType::LCS_Sigma:
      return LCS2_Sigma::kWhat;
  }
  return "";
}

const std::map<std::string_view, ProblemType>&
ProblemFactory::GetNameTypeMap() {
  static std::map<std::string_view, ProblemType> map = {};
  if (map.empty()) {
    for (const auto type : kAvailable) {
      map.insert({GetName(type), type});
    }
  }
  return map;
}

void ProblemFactory::SaveToJson(const Path& file, const BaseProblem* ptr) {
  nlohmann::json j;
  to_json(j, ptr);
  std::ofstream os(file);
  if (!os) {
    throw std::runtime_error("Could not open file " + file.string());
  }
  os << j.dump(4);
}

ProblemFactory::ConstraintSpan ProblemFactory::GetReqConstrains(
    const ProblemType type) {
  switch (type) {
    case ProblemType::LCS_Base:
      return {};
    case ProblemType::LCS_Classic:
      return LCS_Classic::kRequiredConstraints;
    case ProblemType::LCS_MC:
      return LCS2_MC::kRequiredConstraints;
    case ProblemType::LCS_MC_INC:
      return LCS2_MC_INC::kRequiredConstraints;
    case ProblemType::LCS_MC_1C:
      return LCS2_MC_1C::kRequiredConstraints;
    case ProblemType::LCS_MC_O1C_SYNC:
      return LCS2_MC_O1C_SYNC::kRequiredConstraints;
    case ProblemType::LCS_Sigma_L:
      return LCS2_Sigma_L::kRequiredConstraints;
    case ProblemType::LCS_Sigma_R:
      return LCS2_Sigma_R::kRequiredConstraints;
    case ProblemType::LCS_Sigma:
      return LCS2_Sigma::kRequiredConstraints;
  }
  return {};
}

ProblemFactory::ConstraintSpan ProblemFactory::GetOptConstrains(
    const ProblemType type) {
  switch (type) {
    case ProblemType::LCS_Base:
      return ConstraintFactory::GetAvailable();
    // case ProblemType::LCS_Classic:;
    // case ProblemType::LCS_MC:;
    // case ProblemType::LCS_MC_INC:;
    // case ProblemType::LCS_MC_1C:;
    // case ProblemType::LCS_MC_O1C_SYNC:;
    // case ProblemType::LCS_Sigma_L:;
    // case ProblemType::LCS_Sigma_R:;
    // case ProblemType::LCS_Sigma:;
    default:
      break;
  }
  return {};
}

ProblemFactory::AlgorithmSpan ProblemFactory::GetAlgorithms(
    const ProblemType type) {
  switch (type) {
    case ProblemType::LCS_Base:
      return AlgoFactory::GetAvailable();
    case ProblemType::LCS_Classic:
      return LCS_Classic::kAvailableAlgorithms;
    case ProblemType::LCS_MC:
      return LCS2_MC::kAvailableAlgorithms;
    case ProblemType::LCS_MC_INC:
      return LCS2_MC_INC::kAvailableAlgorithms;
    case ProblemType::LCS_MC_1C:
      return LCS2_MC_1C::kAvailableAlgorithms;
    case ProblemType::LCS_MC_O1C_SYNC:
      return LCS2_MC_O1C_SYNC::kAvailableAlgorithms;
    case ProblemType::LCS_Sigma_L:
      return LCS2_Sigma_L::kAvailableAlgorithms;
    case ProblemType::LCS_Sigma_R:
      return LCS2_Sigma_R::kAvailableAlgorithms;
    case ProblemType::LCS_Sigma:
      return LCS2_Sigma::kAvailableAlgorithms;
  }
  return AlgoFactory::GetAvailable();
}

// ----------------------------------------------------------------------------
// Deprecated Serialization Implementations
// ----------------------------------------------------------------------------

/**
 * @copydoc ProblemFactory::FromFile(const Path&)
 *
 * @implementation
 * Legacy file format structure:
 * 1. Problem type name
 * 2. Metadata (name/description)
 * 3. String count followed by strings
 * 4. Constraint definitions
 * 5. Algorithm configuration
 *
 * @deprecated Prefer JSON-based serialization using CreateFromJson()
 */
ProblemFactory::ProblemPtr ProblemFactory::FromFile(const Path& file) {
  std::ifstream in(file);
  if (!in)
    throw std::runtime_error(
        "Could not open file " + file.string() + " as input");

  std::ostringstream os_dummy;
  const ProblemType t = ReadProblemType(in, os_dummy);
  auto [name, description] = ReadDescription(t, in, os_dummy);
  StrPtrVec spv = util::ReadStrPtrVec("", nullptr, in, os_dummy);
  auto map = ConstraintFactory::ReadConstraintMap(spv, nullptr, in, os_dummy);
  ProblemPtr ptr = Create(t,
      std::move(name),
      std::move(description),
      std::move(spv),
      std::move(map));
  AlgoType main = AlgoFactory::ReadAlgoType(in, os_dummy);
  auto child = AlgoType::Unknown;
  if (main == AlgoType::LCS2_RT)
    child = AlgoFactory::ReadAlgoType(in, os_dummy);
  AlgorithmPtr algorithm_ptr = AlgoFactory::Create(
      {
          main,
          child,
      },
      ptr->GetStrPtrVector(),
      ptr->GetConstraints());
  ptr->SetAlgorithm(std::move(algorithm_ptr));
  return ptr;
}


void ProblemFactory::ToFile(const Path& file, const BaseProblem* p) {
  if (!p->IsProblemValid()) {
    throw std::runtime_error("Attempted to write a invalid problem.");
  }
  std::ofstream fs;
  fs.open(file, std::ofstream::out | std::ofstream::trunc);
  if (!fs) {
    throw std::runtime_error("Could not open file " + file.string());
  }
  WriteProblem(fs, p);
  fs.close();
}

void ProblemFactory::WriteProblem(std::ostream& os, const BaseProblem* p) {
  if (!p) return;
  const auto type = p->GetProblemType();
  const std::string_view type_name = GetName(type);
  const std::string_view& name = p->GetName();
  const std::string_view& description = p->GetDescription();
  const uint num_strings = p->GetStrPtrVector().size();
  const StrPtrVec& sting_ptr_vec = p->GetStrPtrVector();
  os << type_name << std::endl;
  os << name << std::endl;
  os << description << std::endl;
  os << num_strings << std::endl;
  for (const auto& sp : sting_ptr_vec) {
    os << util::to_string(*sp) << std::endl;
  }
  ConstraintFactory::WriteConstraintMap(p->GetConstraints(), os);
  // Save AlgoType of the main algorithm that calculates the solution
  AlgoFactory::WriteAlgorithm(p->GetAlgorithm(), os);

  // Save AlgoType of the child algorithm that is used by the main algorithm
  if (p->GetAlgorithm()->getType() == AlgoType::LCS2_RT) {
    auto* main = dynamic_cast<const algorithms::lcs::LCS2_RT*>(p);
    AlgoFactory::WriteAlgorithm(main->getBase(), os);
  }
  os << std::endl;
}

}  // namespace lcs_solver::problems
