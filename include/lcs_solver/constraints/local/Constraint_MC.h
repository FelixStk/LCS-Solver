#ifndef LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_MC_H_
#define LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_MC_H_

#include "constraints/BaseConstraint.h"
#include <nlohmann/json.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lcs_solver::constraints::local {

/*******************************************************************************
 * @brief Base Class for LCS-1C, LCS-MC, LCS-MC-INC, LCS-O(1), LCS-O(1)-SYNC
 *  Constraints.
 *
 * @details This struct inherits from the BaseConstraint class and constraints
 * what a valid embedding is. The gap[i] represents in a pair the minimum and
 * the maximum number of uncommon symbols that can come after the i-th symbol in
 * the longest common subsequence.
 ******************************************************************************/
class Constraint_MC : public BaseConstraint {
 public:
  using uint = std::size_t;
  using Pair = std::pair<uint, uint>;
  using GapVector = std::vector<Pair>;

  static constexpr ConstraintType kType = ConstraintType::MC;
  static constexpr ConstraintCategory kCategory = ConstraintCategory::LOCAL;
  static constexpr const char *kName = "Constraint_MC";
  static constexpr const char *kDescription =
      "If gap[i] = (li,ui) then li <= |gap_e(w,j)| <= ui for every string w "
      "and every gap_e(w,j) = w[e[j]+1 : e[j+1]-1] (where e is an embedding)";
  static constexpr ConstraintType StaticType(){return kType;}

  Constraint_MC() : Constraint_MC(GapVector{}) {};
  explicit Constraint_MC(const GapVector &gc);

  [[nodiscard]] std::string_view GetName() const override;
  [[nodiscard]] std::string_view GetDescription() const override;
  [[nodiscard]] std::string DebugString() const override;
  [[nodiscard]] bool IsConstraintValid(const StrPtrVector &spv) const override;
  [[nodiscard]] bool IsEmbeddingValid(const Embedding &e) const override;
  [[nodiscard]] uint GetGapsLowerBound() const;
  [[nodiscard]] uint GetGapsUpperBound() const;
  [[nodiscard]] uint GetGapsMaxLength() const;
  [[nodiscard]] auto GetGapsNumUniques() const -> uint;
  [[nodiscard]] const GapVector &GetGapVector() const;

  void SetGapVector(const GapVector &gc);

  static GapVector GenRelaxedGapVector(const std::span<const util::String> &spv);
  static GapVector GenRelaxedGapVector(const StrPtrVector &spv);

  // bool operator==(const BaseConstraint &other) const override;

 protected:
  Constraint_MC(const GapVector &gc, ConstraintType type);
  [[nodiscard]] bool IsGapVectorValid(const StrPtrVector &sv) const;

  static uint CalcLowerBound(const GapVector &gc);
  static uint CalcUpperBound(const GapVector &gc);
  static uint CalcMaxLength(const GapVector &gc);
  static uint CalcNumUniques(const GapVector &gc);

  std::vector<Pair> gap_;///< Gap Constraint vector: gap[i] = (l,u) => [l ≤ length(gap_i) ≤ u]
  uint t_;               ///< upper bound on the llcs, k-1 is the maximum possible number of gaps in the lcs
  uint gap_lower_bound_; ///< min of gap[i].first for i in 0..k-2 (there are k-1 gaps)
  uint gap_upper_bound_; ///< max of gap[i].first for i in 0..k-2 (there are k-1 gaps)
  uint gap_max_length_;  ///< equal to gapUpperBound - gapLowerBound
  uint gap_num_uniques_; ///< number of different gaps in the vector gap (value wise)

 public:
  // NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_MC, kName)
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(Constraint_MC, gap_, t_, gap_lower_bound_, gap_upper_bound_,gap_max_length_,gap_num_uniques_)
};

}// namespace lcs_solver::constraints::local
#endif//LCS_SOLVER_CONSTRAINTS_LOCAL_CONSTRAINT_MC_H_
