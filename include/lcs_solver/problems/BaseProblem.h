#ifndef LCS_SOLVER_PROBLEMS_BASE_PROBLEM_H_
#define LCS_SOLVER_PROBLEMS_BASE_PROBLEM_H_

#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "algorithms/AlgoType.h"
#include "algorithms/BaseAlgorithm.h"
#include "algorithms/BaseSolution.h"
#include "constraints/BaseConstraint.h"
#include "constraints/ConstraintMap.h"
#include "constraints/ConstraintType.h"
#include "problems/ProblemType.h"
#include "util/CommonTypes.h"

namespace lcs_solver::problems {
class BaseProblem;

/*******************************************************************************
 * @brief Interface for LCS problems.
 ******************************************************************************/
class BaseProblem {
 public:
  using AlgorithmType = algorithms::AlgoType;
  using AlgorithmSpan = std::span<const AlgorithmType>;
  using AlgorithmPtr = std::unique_ptr<algorithms::BaseAlgorithm>;
  using BaseAlgorithm = algorithms::BaseAlgorithm;
  using BaseSolution = algorithms::BaseSolution;
  using BaseConstraint = constraints::BaseConstraint;

  using ConstraintPtr = std::shared_ptr<BaseConstraint>;
  using ConstraintType = constraints::ConstraintType;
  using ConstraintSpan = std::span<const ConstraintType>;
  using ConstraintMap = constraints::ConstraintMap;
  using ProblemPtr = std::unique_ptr<BaseProblem>;
  using StrPtrVec = util::StrPtrVector;

  BaseProblem() = default;
  BaseProblem(std::string_view name,
    std::string_view description,
    StrPtrVec &&spv,
    ConstraintMap &&map);
  virtual ~BaseProblem() = default;
  BaseProblem(const BaseProblem &);
  BaseProblem &operator=(const BaseProblem &);
  BaseProblem(BaseProblem &&) noexcept;
  BaseProblem &operator=(BaseProblem &&) noexcept;

  static std::unique_ptr<BaseProblem> FromJson(const std::filesystem::path &path);

  [[nodiscard]] virtual ConstraintSpan GetReqConstraints() const;
  [[nodiscard]] virtual AlgorithmSpan GetAlgoSpan() const;

  [[nodiscard]] ProblemType GetProblemType() const;               ///< Gets type_
  [[nodiscard]] std::string_view GetName() const;                 ///< Gets name_
  [[nodiscard]] std::string_view GetDescription() const;          ///< Gets description_
  [[nodiscard]] const StrPtrVec &GetStrPtrVector() const;         ///< Gets spv_
  [[nodiscard]] const ConstraintMap &GetConstraints() const;      ///< Gets map_
  [[nodiscard]] const BaseAlgorithm *GetAlgorithm() const;        ///< Gets algo_

  void SetProblemType(ProblemType type);               ///< Sets type_
  void SetName(std::string_view name);                 ///< Sets name_
  void SetDescription(std::string_view description);   ///< Sets description_
  void SetStrPtrVector(const StrPtrVec &spv);          ///< Sets spv_
  void SetConstraintMap(const ConstraintMap &map);     ///< Sets map_
  void SetAlgorithm(AlgorithmPtr p);                   ///< Sets algo_

  [[nodiscard]] bool IsProblemValid() const;                      ///< Checks if the problem is well-defined
  [[nodiscard]] std::unique_ptr<BaseSolution> ExecuteAlgo() const;///< Executes Algorithm and returns its solution
  [[nodiscard]] std::string DebugString() const;                  ///< String about name, sp, map, and algorithmPtr

  void AddConstraint(const ConstraintPtr &p);  ///< Adds a constraint owned by the problem
  void AddString(const util::String &&s);       ///< Adds a shared ptr to spv_
  void ToJson(const std::filesystem::path &path) const;

  /// Retrieves a pointer to a constraint of type `T` (calls `constraints::Get`)
  template<constraints::ConstraintConcept T>
  T *Get() {
    return constraints::Get<T>(map_);
  }

  /// Retrieves a pointer to an algorithm to of type `T`
  template<typename T>
  requires std::is_base_of_v<algorithms::BaseAlgorithm,T>
  T *Get() const{
    return dynamic_cast<T*>(algo_.get());
  }

  template<typename T>
  requires std::is_base_of_v<algorithms::BaseSolution, T>
  std::unique_ptr<T> ExecuteAlgo(std::type_identity_t<T>* = nullptr) const {
    BaseSolution *p = ExecuteAlgo().release();
    auto result = std::unique_ptr<T>(dynamic_cast<T*>(p));
    if (!result) {
      delete p;
      throw std::runtime_error("ExecuteAlgo<T>: dynamic_cast failed");
    }
    return result;
  }

 protected:
  BaseProblem(
      ProblemType type,
      std::string_view name,
      std::string_view description
  );

  BaseProblem(
      ProblemType type,
      std::string_view name,
      std::string_view description,
      StrPtrVec spv,
      ConstraintMap map
  );

  void Clear();                   ///< Removes resources allocate for the problem

  ProblemType type_ = ProblemType::LCS_Base; ///< Identifier of problem kind
  std::string name_;              ///< Name given to a Problem Instance
  std::string description_;       ///< The Problem's Description
  StrPtrVec spv_;                 ///< Input Strings of the Problem (original)
  ConstraintMap map_;             ///< Map of LCS Constraint Instances
  AlgorithmPtr algo_;             ///< Algorithm to Solve the Problem

};// end of BaseProblem

void to_json(nlohmann::json& j, const BaseProblem* p);
void from_json(const nlohmann::json& j, BaseProblem* p);

}// namespace lcs_solver::problems
#endif /* LCS_SOLVER_PROBLEMS_BASE_PROBLEM_H_ */