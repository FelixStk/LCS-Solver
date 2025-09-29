// #ifndef LCS_SOLVER_LCS2_STD_INC_HPP
// #define LCS_SOLVER_LCS2_STD_INC_HPP

// #include "../../include/algorithms/Strategy.hpp"
// #include "../../include/constraints/BaseConstraint.hpp"
// #include <map>
// #include <memory>
//namespace lcs_solver::algorithms::lcs {
// class Alg_Dag : public Strategy {
//     using StrViewVector = std::vector<const StringView>;
//     using Constraints = std::map<String, std::unique_ptr<BaseConstraint>>;
//     using EmbeddingVector = std::vector<Embedding>;
// public:
//     Alg_Dag(const StrViewVector& svv, Constraints c);
//     void doPreprocessing() override;
//     void resetQuerying() override;
//     Integer getLLCS() const override;
//     void printDataStructure(std::string msg) override;
// private:
//     const Embedding genNext() override;
// };
// }
// #endif /* LCS_SOLVER_LCS2_STD_INC_HPP */