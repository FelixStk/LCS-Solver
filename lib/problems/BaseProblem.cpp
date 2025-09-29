/******************************************************************************
 * @file BaseProblem.cpp
 * @author Steinkopp:Felix
 * @version 3
 * @brief Implementation of fundamental methods for handling problems
 *****************************************************************************/

#include "problems/BaseProblem.h"

#include <algorithm>
#include <fstream>
#include <ranges>
#include <sstream>
#include <utility>

#include "algorithms/AlgoFactory.h"
#include "algorithms/LCS/LCS2_RT.h"
#include "constraints/ConstraintFactory.h"
#include "constraints/ConstraintMap.h"
#include "problems/ProblemFactory.h"

namespace lcs_solver::problems {

/*******************************************************************************
 * BaseProblem - Minimal Constructor
 * @param type ProblemType of the LCS Problem
 * @param name string_view to define the name of the problem
 * @param description string_view to define a description
 ******************************************************************************/
BaseProblem::BaseProblem(const ProblemType type,
                         const std::string_view name,
                         const std::string_view description)
    : type_(type),
      name_(name),
      description_(description){}

/*******************************************************************************
 * BaseProblem - Constructor with util::StrPtrVector and map
 * @param type ProblemType of the LCS Problem
 * @param name string_view to define the name of the problem
 * @param description string_view to define a description
 * @param spv vector with shared pointers to the strings of the problem
 * @param map ConstraintMap: ConstraintType -> ConstraintPtr
 ******************************************************************************/
BaseProblem::BaseProblem(const ProblemType type,
                         const std::string_view name,
                         const std::string_view description,
                         StrPtrVec spv,
                         ConstraintMap map)
    : type_(type),
      name_(name),
      description_(description),
      spv_(std::move(spv)),
      map_(std::move(map)),
      algo_(nullptr) {
}

/*******************************************************************************
 * BaseProblem - Constructor with util::StrPtrVector and map
 * @param name string_view to define the name of the problem
 * @param description string_view to define a description
 * @param spv vector with shared pointers to the strings of the problem
 * @param map ConstraintMap: ConstraintType -> ConstraintPtr
 ******************************************************************************/
BaseProblem::BaseProblem(const std::string_view name,
                         const std::string_view description,
                         StrPtrVec&& spv,
                         ConstraintMap&& map)
    : name_(name),
     description_(description),
     spv_(std::move(spv)),
     map_(std::move(map)),
     algo_(nullptr){
}
/*******************************************************************************
 * Copy Constructor
 * @param other BaseProblem to be copied
 ******************************************************************************/
BaseProblem::BaseProblem(const BaseProblem &other)
    : type_(other.type_),
      name_(other.name_),
      description_(other.description_),
      spv_(other.spv_),
      map_(other.map_)
{
    const auto algo_type = other.algo_->getType();
    algo_ = algorithms::AlgoFactory::Create(algo_type, spv_, map_);
}

/*******************************************************************************
 * Copy Assignment Operator
 * @param other BaseProblem to be used in the allocation
 * @return this object
 ******************************************************************************/
BaseProblem &BaseProblem::operator=(const BaseProblem &other) {
  if (this == &other) {
    return *this;
  }
  // Copy Simple Members
  type_ = other.type_;
  name_ = other.name_;
  description_ = other.description_;
  spv_ = other.spv_;
  map_ = other.map_;

  // Clone the algorithmPtr if it is not null, otherwise set to nullptr
  const auto algo_type = other.algo_->getType();
  if (other.algo_) {
    algo_ = algorithms::AlgoFactory::Create(algo_type, spv_, map_);
  } else {
    algo_ = nullptr;
  }

  return *this;
}

/*******************************************************************************
 * Move Constructor
 * @param other BaseProblem to be moved
 ******************************************************************************/
BaseProblem::BaseProblem(BaseProblem &&other) noexcept
    : type_(other.type_),
      name_(std::move(other.name_)),
      description_(std::move(other.description_)),
      spv_(std::move(other.spv_)),
      map_(std::move(other.map_)),
      algo_(std::move(other.algo_))
{
  Clear();
}

/*******************************************************************************
 * Move Assignment Operator
 * @param other BaseProblem to be moved
 * @return this object
 ******************************************************************************/
BaseProblem &BaseProblem::operator=(BaseProblem &&other) noexcept {
  if (this == &other) {
    return *this;  // Protect against self-assignment
  }
  type_ = other.type_;
  name_ = std::move(other.name_);
  description_ = std::move(other.description_);
  spv_ = std::move(other.spv_);
  map_ = std::move(other.map_);
  algo_ = std::move(other.algo_);
  other.Clear();
  return *this;
}
std::unique_ptr<BaseProblem> BaseProblem::FromJson(
    const std::filesystem::path &path) {
  std::ifstream f(path);
  if (!f) {
    throw std::runtime_error("Could not open file: " + path.string());
  }
  const nlohmann::json j = nlohmann::json::parse(f);
  std::unique_ptr<BaseProblem> ptr;
  from_json(j, ptr.get());
  return ptr;
}

BaseProblem::ConstraintSpan BaseProblem::GetReqConstraints() const {
  return {};
  // return constraints::ConstraintFactory::GetAvailable();
}

BaseProblem::AlgorithmSpan BaseProblem::GetAlgoSpan() const {
  return algorithms::AlgoFactory::GetAvailable();
}

/*******************************************************************************
 * Getter for the type of the problem
 * @return ProblemType type_
 ******************************************************************************/
ProblemType BaseProblem::GetProblemType() const {
  return type_;
}

/*******************************************************************************
 * Getter for name_
 * @return std::string_view of the problem name
 ******************************************************************************/
std::string_view BaseProblem::GetName() const {
  return std::string_view(name_);
}

/*******************************************************************************
 * Getter for description_
 * @return std::string_view of the problem description
 ******************************************************************************/
std::string_view BaseProblem::GetDescription() const {
  return std::string_view(description_);
}

/*******************************************************************************
 * Getter for sorted_spv_
 * @return StrPtrVec Reference of the problem strings
 ******************************************************************************/
const BaseProblem::StrPtrVec &BaseProblem::GetStrPtrVector() const{
  return spv_;
}

/*******************************************************************************
 * Getter for map_
 * @return ConstraintMap to the problem constraints
 ******************************************************************************/
const BaseProblem::ConstraintMap &BaseProblem::GetConstraints() const {
  return map_;
}

/*******************************************************************************
 * Getter for algo_
 * @return std::unique_ptr<BaseAlgorithm> of the problem algorithm
 ******************************************************************************/
const algorithms::BaseAlgorithm *BaseProblem::GetAlgorithm() const {
  return algo_.get();
}

/*******************************************************************************
 * SetProblemType is a setter for type_
 * @param type ProblemType new type of the problem
 ******************************************************************************/
void BaseProblem::SetProblemType(const ProblemType type) {
  type_ = type;
}

/*******************************************************************************
 * SetProblemType is a setter for name_
 * @param name std::string_view new name of the problem
 ******************************************************************************/
void BaseProblem::SetName(const std::string_view name) {
  name_ = name;
}

/*******************************************************************************
 * SetProblemType is a setter for type_
 * @param description  std::string_view new description of the problem
 ******************************************************************************/
void BaseProblem::SetDescription(const std::string_view description) {
  description_ = description;
}

/*******************************************************************************
 * SetStrPtrVector is a setter for original_spv_
 * @details If the `do_sorting` is true it triggers SortAndGenReversePerm()
 * @param spv util::StrPtrVector new String for the problem
 ******************************************************************************/
void BaseProblem::SetStrPtrVector(const StrPtrVec &spv) {
  spv_ = spv;
}

/*******************************************************************************
 * SetConstraintMap is a setter for `map_`
 * @param map ConstraintMap new constraints for the problem
 ******************************************************************************/
void BaseProblem::SetConstraintMap(const ConstraintMap &map) {
  map_ = map;
}

/*******************************************************************************
 * setAlgorithm
 * @param p AlgorithmPtr for the algorithm added to this BaseProblem object
 ******************************************************************************/
void BaseProblem::SetAlgorithm(AlgorithmPtr p) {
  algo_ = std::move(p);
}

/*******************************************************************************
 * clear
 * @details Sets name_ and description_ to an empty string. Additionally, it
 *  clears files_, original_spv_, sorted_spv_, perm_, perm_r_, map_. It also
 *  pushes the unique smart pointer for the algorithm out of scope and replaces
 *  it with the null pointer. The ProblemType type_ is set to LCS_Base
 ******************************************************************************/
void BaseProblem::Clear() {
  type_ = ProblemType::LCS_Base;
  name_ = "";
  description_ = "";
  spv_.clear();
  map_.clear();
  algo_ = nullptr;
}

/*******************************************************************************
 * Predicate: IsProblemValid
 * @return false, iff at least one constraint or the algorithm is invalid
 ******************************************************************************/
bool BaseProblem::IsProblemValid() const {
  // Existence of Required Constraints
  for (const auto &constraint_type : GetReqConstraints()) {
    if (!map_.contains(constraint_type)) {
      return false;
    }
  }
  // Constraints
  for (auto &values : std::views::values(map_)) {
    if (!values->IsConstraintValid(spv_))
      return false;
  }
  // Algorithm
  if (!algo_)
    return false;
  return algo_->isValid();
}

/*******************************************************************************
 * ExecuteAlgo
 * @details Checks whether the problem is valid. If so, it accesses
 *  algorithmPtr to reset the algorithm. Then it continues with pre-processing.
 *  Finally, the method calls `query` on the algorithm and returns its solution.
 * @throws std::runtime_error if the problem is invalid or the algorithm is null
 ******************************************************************************/
std::unique_ptr<BaseProblem::BaseSolution> BaseProblem::ExecuteAlgo() const {
  if (!IsProblemValid()) {
    throw std::runtime_error("Attempted to solve a invalid problem.");
  }
  if (!algo_) {
    throw std::runtime_error("Solution algorithm is not defined");
  }
  algo_->reset();
  algo_->doPreprocessing();
  return algo_->query();
}

/*******************************************************************************
 * Getter for DebugString
 * @return
 ******************************************************************************/
std::string BaseProblem::DebugString() const {
  std::ostringstream oss;
  oss << "Problem name: " << name_ << std::endl;
  oss << "Problem desc: " << description_ << std::endl;
  if (algo_ == nullptr) {
    oss << "Algorithm: " << "NULL" << std::endl;
  } else {
    oss << "Algorithm: " << algo_->getName() << std::endl;
  }

  oss << "Constraints: " << map_.size() << std::endl;
  for (auto &value : std::views::values(map_)) {
    oss << value->DebugString() << std::endl;
  }
  if (spv_.empty()) oss << "No Strings in Problem" << std::endl;
  else
    oss << "Strings:" << std::endl;
  for (size_t i = 0; i < spv_.size(); i++)
    oss << "s[" << i << "] = " << util::to_string(*spv_[i])
        << std::endl;
  return oss.str();// Convert the stream contents to a std::string
}

/*******************************************************************************
 * AddConstraint
 * @param p ConstraintPtr for the constraint added to this BaseProblem object
 ******************************************************************************/
void BaseProblem::AddConstraint(const ConstraintPtr &p) {
  if (p == nullptr) {
    return;
  }
  ConstraintType const key = p->GetType();
  map_[key] = p;
}

/*******************************************************************************
 * AddString
 * @param s util::String && to be moved as a shared_ptr into original_spv_
 ******************************************************************************/
void BaseProblem::AddString(const util::String &&s) {
  spv_.push_back(std::make_shared<util::String>(s));
}

void BaseProblem::ToJson(const std::filesystem::path &path) const {
  std::ofstream fs;
  fs.open(path, std::ofstream::out | std::ofstream::trunc);
  if (!fs) {
    throw std::runtime_error("Could not open file " + path.string() + " for writing");
  }
  nlohmann::json j;
  to_json(j, this);
  fs << j.dump(4);
  fs.close();
}

void to_json(nlohmann::json& j, const BaseProblem* p) {
  j["type_"] = p->GetProblemType();
  j["name_"] = p->GetName();
  j["description_"] = p->GetDescription();
  auto strings = util::StrPtrVectorTo<std::vector<std::string>>(p->GetStrPtrVector());
  j["spv_"] = strings;
  j["map_"] = p->GetConstraints();
  // auto temp = p->GetAlgorithm()->getType();
  auto temp2 = p->GetAlgorithm()->getName();
  j["algo_"] = temp2;
  if (p->GetAlgorithm()->getType() == algorithms::AlgoType::LCS2_RT) {
    using algorithms::lcs::LCS2_RT;
    const auto * temp_ptr = p->Get<const LCS2_RT>();
    j["base_"] = temp_ptr->getBase()->getType();
  }
}

void from_json(const nlohmann::json& j, BaseProblem* p) {
  if (j.contains("type_")) {
    const ProblemType type = j["type_"];
    p->SetProblemType(type);
  }
  else {
    p->SetProblemType(ProblemType::LCS_Base);
  }

  if (j.contains("name_")) {
    const std::string name = j["name_"];
    p->SetName(name);
  }
  else {
    p->SetName("Unnamed Problem");
  }
  if (j.contains("name_")) {
    const std::string description = j["description_"];
    p->SetDescription(description);
  }
  else {
    p->SetDescription("No description");
  }
  if (j.contains("spv_")) {
      const std::vector<std::string> vec = j["spv_"];
    const util::StrPtrVector spv = util::StrPtrVectorFrom(vec);
    p->SetStrPtrVector(spv);
  }
  if (j.contains("map_")) {
    constraints::ConstraintMap map;
    nlohmann::json j2 = j["map_"];
    constraints::from_json(j2,map);
    // constraints::ConstraintMap map = j["map_"].get<constraints::ConstraintMap>();
    p->SetConstraintMap(map);
  }
  if (j.contains("algo_")) {
    const algorithms::AlgoType algo_type = j["algo_"];
    auto base_type = algorithms::AlgoType::Unknown;
    if (j.contains("base_")) {
      base_type = j["base_"];
    }
    p->SetAlgorithm(algorithms::AlgoFactory::Create(algo_type, p->GetStrPtrVector(), p->GetConstraints(), base_type));
  }
}

}// namespace lcs_solver::problems
