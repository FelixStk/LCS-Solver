#ifndef LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_SIGMA_H_
#define LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_SIGMA_H_

#include <unordered_map>
#include "constraints/BaseConstraint.h"
#include "util/CommonTypes.h"
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints::local {
class Constraint_Sigma : public BaseConstraint {
 public:
  static constexpr ConstraintType kType = ConstraintType::SIGMA;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_Sigma";
  static constexpr const char *kDescription =
      "The constraint is defined by: Let gap_e(w,j) = w[e[j]+1 :e[j+1]-1], where e is an embedding.\n"
      "1.) For every gap between two embedded symbols. If e[i]=sym and left[sym] = (l,u) then for all w : l <= |gap_e(w,j)| <= u.\n"
      "2.) For every gap between two embedded symbols. If e[i+1]=sym and right[sym] = (l,u) then for all w : l <= |gap_e(w,j)| <= u.\n";
  static constexpr ConstraintType StaticType(){return kType;}

  using uint = size_t;
  using Pair = std::pair<uint, uint>;
  using GapVector = std::vector<Pair>;
  using GapSpan = std::span<const Pair>;
  using Symbol = util::Symbol;
  using SymbolEqual = util::SymbolEqual;
  using SymbolHash = util::SymbolPerfectHash;
  using SigmaTupleMap = std::unordered_map<Symbol, Pair, SymbolHash, SymbolEqual>;
  using SymbolVector = std::vector<Symbol>;
  using SymbolSpan = std::span<const Symbol>;

  Constraint_Sigma(const SymbolVector &symbol_vec, GapVector left_gaps, GapVector right_gaps);
  Constraint_Sigma(SigmaTupleMap l, SigmaTupleMap r);
  Constraint_Sigma() : Constraint_Sigma(SigmaTupleMap(), SigmaTupleMap()) {}

  static SigmaTupleMap GenSigmaTupleMap(SymbolSpan alph, GapSpan gc);
  static SigmaTupleMap GenSigmaTupleMap(const SymbolVector& alph, const GapVector& gc);
  static GapVector GenRelaxedGapVector(const SymbolVector& alph, const StrPtrVector &spv);
  static GapVector GenRelaxedGapVector(const SymbolVector& alph, const std::span<const util::String> &s);
  static Constraint_Sigma CreateRelaxed(const StrPtrVector &spv);

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  const SigmaTupleMap &GetLeft() const;
  const SigmaTupleMap &GetRight() const;
  void SetLeft(const SigmaTupleMap &l);
  void SetRight(const SigmaTupleMap &r);

 protected:
  Constraint_Sigma(SigmaTupleMap l, SigmaTupleMap r, ConstraintType t);
  Constraint_Sigma(const SymbolVector &symbol_vec, GapVector left_gaps, GapVector right_gaps,  ConstraintType t);
  static bool IsSigmaTupleMapValid(const SigmaTupleMap &map, const StrPtrVector &spv);
  static bool DoLeftCheck(const SigmaTupleMap &l,const Embedding &e);
  static bool DoRightCheck(const SigmaTupleMap &r, const Embedding &e);

  SigmaTupleMap left_;  ///< left[c]  = (l,u) => [l ≤ length(gap before Symbol c) ≤ u]
  SigmaTupleMap right_; ///< right[c] = (l,u) => [l ≤ length(gap after Symbol c) ≤ u]

public:
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_Sigma, left_, right_)
};

} // end of namespace
#endif //LCS_SOLVER_INCLUDE_CONSTRAINTS_LOCAL_CONSTRAINT_SIGMA_H_