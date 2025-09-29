/*******************************************************************************
 * @file test_solver.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests the basic functionality of the solver demonstration program
 ******************************************************************************/

#include "algorithms/LCS/LCS2_RT.h"
#include "demo/solver/solver.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"
#include "problems/ProblemType.h"

namespace {

using lcs_solver::problems::ProblemType;

using testing::PrintToString;
using testing::ValuesIn;

//==== Generate Test Parameters=================================================

using TestParam = std::tuple<std::vector<const char*>,  // argv
                             std::vector<const char*>,  // interactions
                             std::string>;              // name

std::string JoinWithNewlines(const std::vector<const char*>& inputs) {
  std::ostringstream oss;
  for (const char* line : inputs) {
    oss << line << '\n';
  }
  return oss.str();
}

std::vector<TestParam> GetParams(ProblemType t) {
  std::vector<TestParam> params;
  std::vector command = {
      "solver",  // placeholder for the binary name (without meaning here)
      "-i",
  };
  switch (t) {
    // case ProblemType::LCS_Base: {break;}
    case ProblemType::LCS_Classic: {
      std::vector classic_1 = {
          "n",                // Create a new problem
          "LCS_Classic",      // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "LLCS_STD_FL",      // Algorithm to use in the problem
      };
      std::vector classic_2 = {
          "n",                // Create a new problem
          "LCS_Classic",      // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "LLCS2_STD_FL",     // Algorithm to use in the problem
      };
      std::vector classic_3 = {
          "n",                // Create a new problem
          "LCS_Classic",      // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "LCS2_STD_S",       // Algorithm to use in the problem
      };
      std::vector classic_4 = {
          "n",                // Create a new problem
          "LCS_Classic",      // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "LCS2_RT",          // Algorithm to use in the problem
          "LLCS2_STD_FL",     // Algorithm to use in LCS2_RT
      };
      params.emplace_back(command, classic_1, classic_1.back());
      params.emplace_back(command, classic_2, classic_2.back());
      params.emplace_back(command, classic_3, classic_3.back());
      params.emplace_back(command, classic_4, "LCS2_RT");
      break;
    }
    case ProblemType::LCS_MC: {
      std::vector llcs2_mc = {
          "n",                // Create a new problem
          "LCS2_MC",          // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set gap[0]
          "0 4",              // Set gap[1]
          "0 4",              // Set gap[2]
          "LLCS2_MC",         // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_mc, llcs2_mc.back());
      break;
    }
    case ProblemType::LCS_MC_INC: {
      std::vector llcs2_mc_inc = {
          "n",                // Create a new problem
          "LCS2_MC_INC",      // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set gap[0]
          "0 4",              // Set gap[1]
          "0 4",              // Set gap[2]
          "LLCS2_MC_INC",     // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_mc_inc, llcs2_mc_inc.back());
      break;
    }
    case ProblemType::LCS_MC_1C: {
      std::vector llcs2_mc_1c = {
          "n",                // Create a new problem
          "LCS2_MC_1C",       // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set gap[i] for all i
          "LLCS2_MC_1C",      // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_mc_1c, llcs2_mc_1c.back());
      break;
    }
    case ProblemType::LCS_MC_O1C_SYNC: {
      std::vector llcs2_mc_o1c_sync = {
          "n",                 // Create a new problem
          "LCS2_MC_O1C_SYNC",  // Type of the problem
          "my_problem_name",   // Name of the problem
          "my_problem_desc",   // Description of the problem
          "test",              // First string in the problem
          "test",              // Second string in the problem
          "0 4",               // Set gap[0]
          "0 4",               // Set gap[1]
          "0 4",               // Set gap[2]
          "LLCS2_MC_O1_SYNC",  // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_mc_o1c_sync, llcs2_mc_o1c_sync.back());
      break;
    }
    case ProblemType::LCS_Sigma_L: {
      std::vector llcs2_sl_mq = {
          "n",                       // Create a new problem
          "LCS2_Sigma_L",            // Type of the problem
          "my_problem_name",         // Name of the problem
          "my_problem_desc",         // Description of the problem
          "test",                    // First string in the problem
          "test",                    // Second string in the problem
          "0 4",                     // Set left[e]
          "0 4",                     // Set left[s]
          "0 4",                     // Set left[t]
          "LLCS2_SL_R_LLCS2_SR_MQ",  // Algorithm to use in the problem
      };
      std::vector llcs2_sl_rmq = {
          "n",                        // Create a new problem
          "LCS2_Sigma_L",             // Type of the problem
          "my_problem_name",          // Name of the problem
          "my_problem_desc",          // Description of the problem
          "test",                     // First string in the problem
          "test",                     // Second string in the problem
          "0 4",                      // Set left[e]
          "0 4",                      // Set left[s]
          "0 4",                      // Set left[t]
          "LLCS2_SL_R_LLCS2_SR_RMQ",  // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_sl_mq, llcs2_sl_mq.back());
      params.emplace_back(command, llcs2_sl_rmq, llcs2_sl_rmq.back());
      break;
    }
    case ProblemType::LCS_Sigma_R: {
      std::vector llcs2_sr_mq = {
          "n",                // Create a new problem
          "LCS2_Sigma_R",     // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set right[e]
          "0 4",              // Set right[s]
          "0 4",              // Set right[t]
          "LLCS2_SR_MQ",      // Algorithm to use in the problem
      };
      std::vector llcs2_sr_rmq = {
          "n",                // Create a new problem
          "LCS2_Sigma_R",     // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set right[e]
          "0 4",              // Set right[s]
          "0 4",              // Set right[t]
          "LLCS2_SR_RMQ",     // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_sr_mq, llcs2_sr_mq.back());
      params.emplace_back(command, llcs2_sr_rmq, llcs2_sr_rmq.back());
      break;
    }
    case ProblemType::LCS_Sigma: {
      std::vector llcs2_sa_mq = {
          "n",                // Create a new problem
          "LCS2_SIGMA",       // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set left[e]
          "0 4",              // Set left[s]
          "0 4",              // Set left[t]
          "0 4",              // Set right[e]
          "0 4",              // Set right[s]
          "0 4",              // Set right[t]
          "LLCS2_SA_MQ",      // Algorithm to use in the problem
      };
      std::vector llcs2_sa_rmq = {
          "n",                // Create a new problem
          "LCS2_SIGMA",       // Type of the problem
          "my_problem_name",  // Name of the problem
          "my_problem_desc",  // Description of the problem
          "test",             // First string in the problem
          "test",             // Second string in the problem
          "0 4",              // Set left[e]
          "0 4",              // Set left[s]
          "0 4",              // Set left[t]
          "0 4",              // Set right[e]
          "0 4",              // Set right[s]
          "0 4",              // Set right[t]
          "LLCS2_SA_RMQ",     // Algorithm to use in the problem
      };
      params.emplace_back(command, llcs2_sa_mq, llcs2_sa_mq.back());
      params.emplace_back(command, llcs2_sa_rmq, llcs2_sa_rmq.back());
      break;
    }
    default: {
      std::vector quit_action = {"e"};
      params.emplace_back(command, quit_action, "quit action");
    }
  }
  return params;
}

//==== Initialize Test Suites ==================================================
class IntegrationDemoSolver : public testing::TestWithParam<TestParam> {};

INSTANTIATE_TEST_SUITE_P(
    LCS_Classic,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_Classic)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS_MC,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_MC)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS_MC_INC,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_MC_INC)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS_MC_1C,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_MC_1C)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS_MC_O1C_SYNC,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_MC_O1C_SYNC)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS_Sigma_R,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_Sigma_R)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS_Sigma_L,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_Sigma_L)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(
    LCS2_SIGMA,
    IntegrationDemoSolver,
    ValuesIn(GetParams(ProblemType::LCS_Sigma)),
    [](const testing::TestParamInfo<IntegrationDemoSolver::ParamType>& info) {
      return std::get<2>(info.param);
    });

//=== Definition of Value-Parameterized Tests ==================================
TEST_P(IntegrationDemoSolver, CreateAndPrint) {
  const TestParam& param = GetParam();
  const std::vector<const char*>& command = std::get<0>(param);
  const std::vector<const char*>& interactions = std::get<1>(param);
  const std::string& name = std::get<2>(param);
  const auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  std::string in = JoinWithNewlines(interactions) + "p\n" + "e\n";

  // Redirect in and out streams for testing
  std::ostringstream output_stream;
  std::istringstream input_stream(in);
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::streambuf* cin_buf = std::cin.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  std::cin.rdbuf(input_stream.rdbuf());
  EXPECT_NO_THROW(solver::Run(argc, argv)) << "Algo:" << name;
  std::cout.rdbuf(cout_buf);
  std::cin.rdbuf(cin_buf);
  std::string output = output_stream.str();
  for (const auto& action : interactions) {
    if (action == "n")  // Action for starting the creation process
      continue;
    EXPECT_NE(output.find(action), std::string::npos) << action;
  }
}

TEST_P(IntegrationDemoSolver, CreateAndRun) {
  const TestParam& param = GetParam();
  const std::vector<const char*>& command = std::get<0>(param);
  const std::vector<const char*>& interactions = std::get<1>(param);
  const std::string& name = std::get<2>(param);
  const auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  std::string in = JoinWithNewlines(interactions) + "x\n" + "e\n";

  // Redirect in and out streams for testing
  std::ostringstream output_stream;
  std::istringstream input_stream(in);
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::streambuf* cin_buf = std::cin.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  std::cin.rdbuf(input_stream.rdbuf());
  EXPECT_NO_THROW(solver::Run(argc, argv)) << "Algo:" << name;
  std::cout.rdbuf(cout_buf);
  std::cin.rdbuf(cin_buf);
}

TEST_P(IntegrationDemoSolver, Json) {
  const TestParam& param = GetParam();
  const std::vector<const char*>& command = std::get<0>(param);
  const std::vector<const char*>& actions = std::get<1>(param);
  const auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  std::string in = JoinWithNewlines(actions) + "s\n" + "test.json\n" + "o\n" + "test.json\n" + "e\n";

  // Redirect in and out streams for testing
  std::ostringstream output_stream;
  std::istringstream input_stream(in);
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::streambuf* cin_buf = std::cin.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  std::cin.rdbuf(input_stream.rdbuf());
  EXPECT_NO_THROW(solver::Run(argc, argv));
  std::cout.rdbuf(cout_buf);
  std::cin.rdbuf(cin_buf);
}

}  // namespace