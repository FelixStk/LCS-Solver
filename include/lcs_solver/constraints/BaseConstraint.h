#ifndef LCS_SOLVER_CONSTRAINTS_BASECONSTRAINT_H_
#define LCS_SOLVER_CONSTRAINTS_BASECONSTRAINT_H_

#include <vector>
#include <string>
#include <string_view>
#include "constraints/ConstraintType.h"
#include "constraints/ConstraintCategory.h"
#include "structures/Embedding.h"
#include "util/CommonTypes.h"
#include "util/StrPtrVector.h"
#include <nlohmann/json.hpp>

namespace lcs_solver::constraints {

class BaseConstraint {
 public:
  using Unsigned = util::uint;
  using StrPtrVector = util::StrPtrVector;
  using Embedding = structures::Embedding;

  BaseConstraint(ConstraintCategory t, ConstraintType identifier);
  virtual ~BaseConstraint() = default;

  static std::shared_ptr<BaseConstraint> FromJson(const std::filesystem::path &path);
  static void ToJson(const std::filesystem::path &path, const BaseConstraint *constraints);

  [[nodiscard]] virtual bool IsConstraintValid(const StrPtrVector &strPtrVec) const = 0;
  [[nodiscard]] virtual bool IsEmbeddingValid(const Embedding &e) const = 0;
  [[nodiscard]] virtual std::string_view GetName() const = 0;
  [[nodiscard]] virtual std::string_view GetDescription() const = 0;
  [[nodiscard]] virtual std::string DebugString() const = 0;
  [[nodiscard]] ConstraintCategory GetCategory() const;
  [[nodiscard]] ConstraintType GetType() const;

  void ToJson(const std::filesystem::path &path) const;
  // [[nodiscard]] virtual bool operator==(const BaseConstraint &other) const = 0;

 private:
  ConstraintCategory category_;
  ConstraintType type_;
};

} // end of namespace
#endif /* LCS_SOLVER_CONSTRAINTS_BASECONSTRAINT_H_ */