#ifndef LCS_SOLVER_CONSTRAINTS_CONSTRAINTFACTORY_H_
#define LCS_SOLVER_CONSTRAINTS_CONSTRAINTFACTORY_H_

#include "constraints/ConstraintMap.h"
#include "constraints/ConstraintType.h"
#include "util/CommonTypes.h"
#include "util/StrPtrVector.h"
#include <cstddef>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <iostream>

namespace lcs_solver::constraints {

class BaseConstraint;

class ConstraintFactory {
 public:
  using uint = size_t;
  using ConstraintPtr = std::shared_ptr<BaseConstraint>;
  using ConstraintSpan = std::span<const ConstraintType>;
  using Pair = std::pair<uint, uint>;
  using GapVector = std::vector<Pair>;
  using SigmaTupleMap = std::unordered_map<util::Symbol, Pair, util::SymbolPerfectHash, util::SymbolEqual>;
  using ReverseMap = std::unordered_map<util::Symbol, std::string>;

  virtual ~ConstraintFactory()= default;

  static ConstraintPtr Create(ConstraintType t, const GapVector &gc);
  static ConstraintPtr Create(ConstraintType t, const SigmaTupleMap &m);
  static ConstraintPtr Create(ConstraintType t, const SigmaTupleMap &l, const SigmaTupleMap &r);

  static std::string_view GetName(ConstraintType t);
  static std::string_view GetDescription(ConstraintType t);
  static ConstraintSpan GetAvailable();

  static GapVector ReadGapVector(const util::StrPtrVector &spv, std::string_view name = "", std::istream &in=std::cin, std::ostream &os = std::cout);
  static SigmaTupleMap ReadSigmaTupleMap(const util::StrPtrVector &spv, std::string_view name = "", const ReverseMap* r_map = nullptr, std::istream &in=std::cin, std::ostream &os = std::cout);
  static ConstraintType ReadConstraintType(ConstraintSpan v = GetAvailable(), std::istream &in = std::cin, std::ostream& os = std::cout);
  static ConstraintPtr ReadConstraint(const util::StrPtrVector &spv, ConstraintType t, const ReverseMap* r_map, std::istream &in = std::cin, std::ostream &os = std::cout);
  static ConstraintMap ReadConstraintMap(ConstraintSpan v, const util::StrPtrVector &spv, const ReverseMap* r_map, std::istream &in, std::ostream &os=std::cout);
  static ConstraintMap ReadConstraintMap(ConstraintSpan req, ConstraintSpan opt, const util::StrPtrVector &spv, const ReverseMap* r_map, std::istream &in = std::cin, std::ostream& os=std::cout);
  static ConstraintMap ReadConstraintMap(const util::StrPtrVector &spv, const ReverseMap* r_map, std::istream &in=std::cin, std::ostream &os=std::cout);

  static void AddConstraint(ConstraintMap & map, ConstraintType t, const util::StrPtrVector &spv, const ReverseMap* r_map, std::istream &in=std::cin, std::ostream &os=std::cout);
  static void AddConstraint(ConstraintMap & map, ConstraintSpan v, bool force_add, const util::StrPtrVector &spv, const ReverseMap* r_map, std::istream &in=std::cin, std::ostream &os=std::cout);

  static void WriteConstraint(BaseConstraint const *ptr, std::ostream &os = std::cout);
  static void WriteConstraintMap(const ConstraintMap &map, std::ostream &os = std::cout);

};

}// namespace lcs_solver::constraints

#endif//LCS_SOLVER_CONSTRAINTS_CONSTRAINTFACTORY_H_
