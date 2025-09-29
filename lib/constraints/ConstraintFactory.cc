/******************************************************************************
 * @file ConstraintFactory.cc
 * @author Steinkopp:Felix
 * @version 2.0
 * @brief Implementation of a factor to create Constraints
 *****************************************************************************/

#include "constraints/ConstraintFactory.h"

#include <map>
#include <ranges>
#include <sstream>

#include "constraints/global/Constraint_BR.h"
#include "constraints/global/Constraint_LLCS_Known.h"
#include "constraints/input/Constraint_Const_Sigma.h"
#include "constraints/input/Constraint_LCS_2.h"
#include "constraints/input/Constraint_LCS_K.h"
#include "constraints/local/Constraint_MC.h"
#include "constraints/local/Constraint_MC_1C.h"
#include "constraints/local/Constraint_MC_INC.h"
#include "constraints/local/Constraint_MC_O1C.h"
#include "constraints/local/Constraint_MC_O1C_SYNC.h"
#include "constraints/local/Constraint_Sigma.h"
#include "constraints/local/Constraint_Sigma_L.h"
#include "constraints/local/Constraint_Sigma_R.h"
#include "util/InOutHelper.h"
#include "util/StrPtrVector.h"

namespace lcs_solver::constraints {

/*******************************************************************************
 * Getter for globally enabled Constraints
 * @return std::span<const ConstraintType>
 ******************************************************************************/
std::span<const ConstraintType> ConstraintFactory::GetAvailable() {
  static constexpr std::array available_types = {
      ConstraintType::MC,          ConstraintType::MC_INC,
      ConstraintType::MC_1C,       ConstraintType::MC_O1C,
      ConstraintType::MC_O1C_SYNC, ConstraintType::SIGMA,
      ConstraintType::SIGMA_R,     ConstraintType::SIGMA_L};
  return available_types;
}

/*******************************************************************************
 * @brief Creates a constraint object based on the specified constraint type and
 * gap vector.
 * @details This overload of Create constructs and returns a shared pointer to a
 * constraint object that operates on a gap vector. The type of constraint to be
 * created is determined by the provided ConstraintType enum value.
 * @param t The type of constraint to create.
 * @param gc A constant reference to a std::vector<std::pair<size_t, size_t>>
 * object used to initialize the constraint.
 * @return A shared pointer to the created constraint object, or nullptr if the
 * type is Empty.
 * @throws std::runtime_error If the specified constraint type is unknown or
 * construction with a GapVector is unsupported.
 ******************************************************************************/
ConstraintFactory::ConstraintPtr ConstraintFactory::Create(
    const ConstraintType t, const GapVector &gc) {
  switch (t) {
    case ConstraintType::Empty:
      return nullptr;
    case ConstraintType::MC:
      return std::make_shared<local::Constraint_MC>(gc);
    case ConstraintType::MC_INC:
      return std::make_shared<local::Constraint_MC_INC>(gc);
    case ConstraintType::MC_1C:
      return std::make_shared<local::Constraint_MC_1C>(gc);
    case ConstraintType::MC_O1C:
      return std::make_shared<local::Constraint_MC_O1C>(gc);
    case ConstraintType::MC_O1C_SYNC:
      return std::make_shared<local::Constraint_MC_O1C_SYNC>(gc);
    default:
      throw std::runtime_error("ConstraintFactory::Create: unknown constraint type");
  }
}

/*******************************************************************************
 * @brief Creates a constraint object based on the specified constraint type and
 * a SigmaTupleMap.
 * @details This overload of Create constructs and returns a shared pointer to a
 * constraint object that uses a SigmaTupleMap for initialization. Typically
 * used for sigma-based constraints that operate on a single `SigmaTupleMap`.
 * @param t The type of constraint to create.
 * @param m A constant reference to a SigmaTupleMap object used to initialize
 * the constraint.
 * @return A shared pointer to the created constraint object, or nullptr if the
 * type is Empty.
 * @throws std::runtime_error If the specified constraint type is unknown or
 *  construction with a one SigmaTupleMap is unsupported.
 ******************************************************************************/
ConstraintFactory::ConstraintPtr ConstraintFactory::Create(
    const ConstraintType t,
    const SigmaTupleMap &m) {
  switch (t) {
    case ConstraintType::Empty:
      return nullptr;
    case ConstraintType::SIGMA_R:
      return std::make_shared<local::Constraint_Sigma_R>(m);
    case ConstraintType::SIGMA_L:
      return std::make_shared<local::Constraint_Sigma_L>(m);
    default:
      throw std::runtime_error("ConstraintFactory::Create: unknown constraint type");
  }
}

/*******************************************************************************
 * @brief Creates a constraint object based on the specified constraint type and
 * two SigmaTupleMaps.
 * @details This overload of Create constructs and returns a shared pointer to a
 * constraint object that uses two SigmaTupleMaps (aka std::unordered_map<
 * util::Symbol, Pair, util::SymbolPerfectHash, util::SymbolEqual>) for
 * initialization. It is typically used for constraints that compare or
 * synchronize data between two sigma maps.
 * @param t The type of constraint to create.
 * @param l A constant reference to the left SigmaTupleMap. Such that left[c] =
 * (l,u) implies l ≤ (length of the gap before Symbol c) ≤ u
 * @param r A constant reference to the right SigmaTupleMap. Such that right[c]
 * = (l,u) implies l ≤ (length of the gap after Symbol c) ≤ u
 * @return A shared pointer to the created constraint object, or nullptr if the
 * type is Empty.
 * @throws std::runtime_error If the specified constraint type is unknown or
 * construction with a two SigmaTupleMap is unsupported.
 ******************************************************************************/

ConstraintFactory::ConstraintPtr ConstraintFactory::Create(
    const ConstraintType t,
    const SigmaTupleMap &l,
    const SigmaTupleMap &r) {
  switch (t) {
    case ConstraintType::Empty:
      return nullptr;
    case ConstraintType::SIGMA:
      return std::make_shared<local::Constraint_Sigma>(l, r);
    case ConstraintType::SIGMA_L:
      return std::make_shared<local::Constraint_Sigma_L>(l);
    case ConstraintType::SIGMA_R:
      return std::make_shared<local::Constraint_Sigma_R>(r);
    default:
      throw std::runtime_error("ConstraintFactory::Create: unknown constraint type");
  }
}

/*******************************************************************************
 * Getter for names associated with ConstraintTypes
 * @param t ConstraintType
 * @return std::string_view containing the name associated with t
 ******************************************************************************/
std::string_view ConstraintFactory::GetName(const ConstraintType t) {
  return GetMapTypeToName().at(t);
}

/*******************************************************************************
 * Getter for Descriptions associated with ConstraintTypes
 * @param t ConstraintType
 * @return string_view containing the description associated with t
 ******************************************************************************/
std::string_view ConstraintFactory::GetDescription(const ConstraintType t) {
  return GetMapTypeToDesc().at(t);
}

/*******************************************************************************
 * IO-Function: ReadGapVector
 * @param spv vector with shared_ptr of Strings
 * @param name std::string_view that constraints the name of the GapVector to be
 *  read (is displayed after calling ReadGapVector if os != nullptr)
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 * @return std::vector<std::pair<size_t, size_t>> represents the gc tuple
 ******************************************************************************/
ConstraintFactory::GapVector ConstraintFactory::ReadGapVector(
    const util::StrPtrVector &spv,
    const std::string_view name,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << "Enter " << name << (name.empty() ? "" : " ") << "a GapVector\n";
  }
  GapVector gc;
  const uint min_length = util::CalcMinStrLen(spv);
  const uint max_length = util::CalcMaxStrLen(spv);
  uint number_of_gaps = min_length - 1;
  if (os.rdbuf() != std::cout.rdbuf())
    number_of_gaps = util::ReadUnsigned("",in,os); // For file format reasons
  gc.reserve(number_of_gaps);
  for (uint i = 0; i < number_of_gaps; ++i) {
    std::string gap_name = "gap[" + std::to_string(i) + "]";
    auto pair = util::ReadSubinterval(gap_name, {0, max_length}, {0, max_length}, in ,os);
    gc.push_back(pair);
  }
  return gc;
}

/*******************************************************************************
 * IO-Function: ReadSigmaTupleMap
 * @param spv vector with shared_ptr of Strings
 * @param in std::istream from which to read
 * @param r_map std::unordered_map<util::Symbol, std::string> to display symbols
 * @param os std::ostream for prompting and outputting information
 * @param name std::string_view that constraints the name of the GapVector to be
 *  read (is displayed after calling ReadGapVector if os != nullptr)
 * @return  std::unordered_map<util::Symbol, Pair>
 ******************************************************************************/
ConstraintFactory::SigmaTupleMap ConstraintFactory::ReadSigmaTupleMap(
    const util::StrPtrVector &spv,
    const std::string_view name,
    const ReverseMap* r_map,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    os << "Enter " << name << (name.empty() ? "" : " ") << "a SigmaTupleMap\n";
  }
  SigmaTupleMap map;
  const uint max_length = util::CalcMaxStrLen(spv);
  if (os.rdbuf() != std::cout.rdbuf()) {  // Read from a file (e.g., a text file)
    const uint n = util::ReadUnsigned("",in, os);
    for (uint i = 0; i < n; ++i) {
      std::string temp;
      std::getline(in, temp);
      util::String input = util::to_String<util::Symbol>(std::move(temp));
      util::Symbol c;
      uint l, u;
      if (std::basic_istringstream<util::Symbol> iss(input); !(iss >> c >> l >> u))
        throw std::runtime_error("ConstraintFactory: pair input had invalid format");
      map.insert({c, {l, u}});
    }
  } else { // Read from the terminal interactively
    for (const auto &c : util::CalcAlphabet(spv)) {
      std::string gap_name;
      if (r_map == nullptr) {
        gap_name = std::string(name) + "[" + util::to_string(c) + "]";
      }else {
        gap_name = std::string(name) + " for \n````\n" + r_map->at(c) + "\n```\nThis is";
      }
      auto pair = util::ReadSubinterval(gap_name, {0, max_length}, {0, max_length},in, os);
      map.insert({c, pair});
    }
  }
  return map;
}

/*******************************************************************************
 * IO-Function: ReadConstraintType
 * @param v ConstraintSpan from which to choose a ConstraintType
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 * @return ConstraintType chosen
 ******************************************************************************/
ConstraintType ConstraintFactory::ReadConstraintType(
    const ConstraintSpan v,
    std::istream &in,
    std::ostream &os) {
  if (os.rdbuf() == std::cout.rdbuf()) {
    uint n = 0;
    os << "The available ConstraintTypes are:\n";
    for (const auto &type : v) {
      n++;
      os << n << ") " << GetName(type) << std::endl;
    }
    os << "Enter a number or name of a constraint:\n";
  }
  std::string input;

  std::getline(in, input);
  if (util::IsNumber(input)) {
    const uint idx = std::stoul(input) - 1;
    if (idx >= v.size())
      throw std::runtime_error("ConstraintFactory: readConstraintType input number was out of range");
    return v[idx];
  }
  ConstraintType type;
  const auto & look_up_name_map = GetMapNameToType();
  if (const auto it = look_up_name_map.find(input); it == look_up_name_map.end()) {
    throw std::runtime_error("ConstraintFactory: readConstraintType type is not available");
  } else {
    type = it->second;
  }
  return type;
}

/*******************************************************************************
 * IO-Function: ReadConstraint
 * @param spv vector with shared_ptr of Strings
 * @param t ConstraintType of the Constraint to be generated
 * @param r_map std::unordered_map<util::Symbol, std::string> to display symbols
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 * @return std::shared_ptr<BaseConstraint>
 ******************************************************************************/
ConstraintFactory::ConstraintPtr ConstraintFactory::ReadConstraint(
    const util::StrPtrVector &spv,
    const ConstraintType t,
    const ReverseMap* r_map,
    std::istream &in,
    std::ostream &os) {
  std::shared_ptr<BaseConstraint> ptr;
  switch (t) {
    case ConstraintType::Empty: {
      return nullptr;
    }
    case ConstraintType::MC: {
      auto gc_vector = ReadGapVector(spv, "" , in, os);
      return std::make_shared<local::Constraint_MC>(gc_vector);
    }
    case ConstraintType::MC_INC: {
      auto gc_vector = ReadGapVector(spv, "", in, os);
      return std::make_shared<local::Constraint_MC_INC>(gc_vector);
    }
    case ConstraintType::MC_1C: {
      const Pair def = {0, util::CalcMinStrLen(spv)};
      auto [l, u] = util::ReadSubinterval("1C", def, def, in, os);
      return std::make_shared<local::Constraint_MC_1C>(spv, l, u);
    }
    case ConstraintType::MC_O1C: {
      auto gc_vector = ReadGapVector(spv,"", in, os);
      return std::make_shared<local::Constraint_MC_O1C>(gc_vector);
    }
    case ConstraintType::MC_O1C_SYNC: {
      auto gc_vector = ReadGapVector(spv, "", in, os);
      return std::make_shared<local::Constraint_MC_O1C_SYNC>(gc_vector);
    }
    case ConstraintType::SIGMA: {
      SigmaTupleMap left = ReadSigmaTupleMap(spv, "left", r_map, in, os);
      SigmaTupleMap right = ReadSigmaTupleMap(spv, "right", r_map, in, os);
      return std::make_shared<local::Constraint_Sigma>(left, right);
    }
    case ConstraintType::SIGMA_R: {
      SigmaTupleMap right = ReadSigmaTupleMap(spv, "right", r_map, in, os);
      return std::make_shared<local::Constraint_Sigma_R>(right);
    }
    case ConstraintType::SIGMA_L: {
      SigmaTupleMap left = ReadSigmaTupleMap(spv, "left", r_map, in, os);
      return std::make_shared<local::Constraint_Sigma_L>(left);
    }
    case ConstraintType::BR: {
      const std::pair<size_t,size_t> range = {0, util::CalcMaxStrLen(spv)};
      size_t b = util::ReadUnsigned( "b", range.second, range, in, os);
      return std::make_shared<global::Constraint_BR>(b);
    }
    case ConstraintType::LLCS_KNOWN: {
      const Pair range = {0, util::CalcMinStrLen(spv)};
      size_t llcs = util::ReadUnsigned("llcs", range.second, range, in, os);
      return std::make_shared<global::Constraint_LLCS_Known>(llcs);
    }
    case ConstraintType::STRINGS_K: {
      // Pair range = {0, spv.size()};
      // size_t llcs = util::ReadUnsigned(&in, os, "k", range.second, range);
      return std::make_shared<input::Constraint_LCS_K>(spv.size());
    }
    case ConstraintType::STRINGS_2: {
      return std::make_shared<input::Constraint_LCS_2>();
    }
    case ConstraintType::CONST_SIG: {
      const std::vector<util::Symbol> alphabet = util::CalcAlphabet(spv);
      // Pair range = {0, alphabet.size()};
      // size_t sig = util::ReadUnsigned(&in, os, "sigma", range.second, range);
      return std::make_shared<input::Constraint_Const_Sigma>(alphabet.size());
    }
    default: return nullptr;
  }
}

/*******************************************************************************
 * IO-Function: ReadConstraintMap
 * @param v ConstraintSpan from which to choose ConstraintTypes
 * @param spv vector with shared_ptr of Strings
 * @param r_map std::unordered_map<util::Symbol, std::string> to display symbols
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 * @return std::map<ConstraintType, std::shared_ptr<BaseConstraint>>
 * @throws std::runtime_error If the input could not be read successfully
 ******************************************************************************/
ConstraintMap ConstraintFactory::ReadConstraintMap(
    const ConstraintSpan v,
    const util::StrPtrVector &spv,
    const ReverseMap* r_map,
    std::istream &in,
    std::ostream &os) {
  ConstraintMap map;
  if (v.empty()) {
    return map;
  }
  constexpr uint default_value = 1;
  const auto prompt =
      "How many Constraints do you want to add? Enter a number";
  const size_t n = util::ReadUnsigned(prompt, default_value, in, os);
  for (uint i = 0; i < n; ++i) {
    ConstraintType t = ReadConstraintType(v,in, os);
    if (map.contains(t)) {
      throw std::runtime_error("ConstraintFactory::ReadConstraintMap: duplicate constraint");
    }
    auto constraint_ptr = ReadConstraint(spv, t, r_map, in, os);
    map.insert({t, std::move(constraint_ptr)});
  }
  return map;
}

/*******************************************************************************
 * IO-Function: ReadConstraintMap
 * @param req ConstraintSpan for Constraints that must be added
 * @param opt ConstraintSpan for Constraints that are optional
 * @param spv vector with shared_ptr of Strings
 * @param r_map std::unordered_map<util::Symbol, std::string> to display symbols
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 * @return std::map<ConstraintType, std::shared_ptr<BaseConstraint>>
 * @throws std::runtime_error If the input could not be read successfully
 ******************************************************************************/
ConstraintMap ConstraintFactory::ReadConstraintMap(
    const ConstraintSpan req,
    const ConstraintSpan opt,
    const util::StrPtrVector &spv,
    const ReverseMap* r_map,
    std::istream &in,
    std::ostream &os
) {
  ConstraintMap map;
  for (auto t : req) {
    auto constraint_ptr = ReadConstraint(spv, t, r_map, in ,os);
    map.insert({t, std::move(constraint_ptr)});
  }
  for (auto &[constraintType, constraint_ptr]: ReadConstraintMap(opt, spv, r_map, in, os) ) {
    map.insert({constraintType, std::move(constraint_ptr)});
  }
  return map;
}

/*******************************************************************************
 * IO-Function: WriteConstraint
 * @param ptr BaseConstraint pointer
 * @param os std::ostream for writing information about the constraint
 ******************************************************************************/
void ConstraintFactory::WriteConstraint(
    const BaseConstraint *ptr,
    std::ostream &os) {
  os << ptr->DebugString() << std::endl;
}

/*******************************************************************************
 * IO-Function: WriteConstraintMap
 * @param map ConstraintMap reference
 * @param os std::ostream for writing information about the constraint
 ******************************************************************************/
void ConstraintFactory::WriteConstraintMap(
    const ConstraintMap &map,
    std::ostream &os) {
  os << map.size() << std::endl;
  for (const auto & c: map | std::views::values) {
    WriteConstraint(c.get(), os);
  }
}

/*******************************************************************************
 * IO-Function: ReadConstraintMap
 * @param spv vector with shared_ptr of Strings
 * @param r_map std::unordered_map<util::Symbol, std::string> to display symbols
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 * @return std::map<ConstraintType, std::shared_ptr<BaseConstraint>>
 * @throws std::runtime_error if the input cannot be parsed
 ******************************************************************************/
ConstraintMap ConstraintFactory::ReadConstraintMap(
    const util::StrPtrVector& spv,
    const ReverseMap* r_map,
    std::istream& in,
    std::ostream& os) {
  const auto& v = GetAvailable();
  return ReadConstraintMap(v, spv, r_map, in, os);
}

/*******************************************************************************
 * @brief Adds a ConstraintPtr to a Constraint Map
 * @param[in,out] map ConstraintMap in which to add a constraint of type `t`
 * @param[in] t ConstraintType of the Constraint to be generated
 * @param[in] spv vector with shared_ptr of Strings
 * @param[in] r_map unordered_map for translating Symbols to string (for `os`)
 * @param[in] in std::istream from which to read
 * @param[in] os std::ostream for prompting and outputting information
 ******************************************************************************/
void ConstraintFactory::AddConstraint(ConstraintMap& map,
                                      const ConstraintType t,
                                      const util::StrPtrVector& spv,
                                      const ReverseMap* r_map,
                                      std::istream& in,
                                      std::ostream& os) {
  map.at(t) = ReadConstraint(spv, t, r_map, in, os);
}

/*******************************************************************************
 * Adds a Constraint Pointers to a Constraint Map
 * @param map ConstraintMap in which to add a constraint of type `t`
 * @param v Span of ConstraintTypes to consider
 * @param force_add Whether to force a new insertion for every element of `v`
 * @param spv vector with shared_ptr of Strings
 * @param r_map unordered_map for translating Symbols to string (for `os`)
 * @param in std::istream from which to read
 * @param os std::ostream for prompting and outputting information
 ******************************************************************************/
void ConstraintFactory::AddConstraint(ConstraintMap& map,
                                      ConstraintSpan v,
                                      bool force_add,
                                      const util::StrPtrVector& spv,
                                      const ReverseMap* r_map,
                                      std::istream& in,
                                      std::ostream& os) {
  if (force_add) {
    for (const auto& type : v) {
      if (os.rdbuf() == std::cout.rdbuf()) {
        os << "Adding Constraint of type" << GetName(type) ;
      }
      AddConstraint(map, type, spv, r_map, in, os);
    }
  }else {
    // force_add == false
    for (ConstraintMap temp = ReadConstraintMap(v, spv, r_map, in, os);
         const auto& [key, val] : temp) {
      map[key] = val; // override the constraint if it exists
    }
  }
}

}  // namespace lcs_solver::constraints