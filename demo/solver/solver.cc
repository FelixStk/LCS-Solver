/*******************************************************************************
 * @file solver.cc
 * @author Felix Steinkopp
 * @version 1.0.0
 * @brief Implementation of the solver demonstration for the lcs_solver library
 * @details Allows an interactive creation of problems and saving/loading them
 ******************************************************************************/

#include "solver.h"

#include <ctime>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "lcs_solver/algorithms/BaseSolution.h"
#include "lcs_solver/problems/BaseProblem.h"
#include "lcs_solver/problems/ProblemFactory.h"
#include "lcs_solver/util/InOutHelper.h"

namespace solver {

void Run(const int argc, char* argv[]) {
  using lcs_solver::problems::BaseProblem;
  using lcs_solver::problems::ProblemFactory;
  using lcs_solver::util::InputParser;

  InputParser const parser(argc, argv);  // for checking "-f filename"
  std::unique_ptr<BaseProblem> problem = nullptr;

  HandleHelp(parser);       // `-h` flag (print usage msg)
  HandleFileMode(parser);   // `-f` prob.json flag (solve problem in prob.json)
  HandleInputMode(parser);  // `-i` flag (interactive problem creation)
}

void HandleHelp(const InputParser& parser) {
  // Normal Help
  if (parser.HasOption("-h") || parser.HasOption("--help")) {
    PrintUsageMessage(parser);
    return;
  }

  // User needs infos (because he doesn't know better)
  if (!parser.HasOption("-f") && !parser.HasOption("-i")) {
    PrintUsageMessage(parser);
  }
}

void HandleFileMode(const InputParser& parser) {
  // Return early if not in File Mode
  if (!parser.HasOption("-f")) {
    return;
  }

  // Solve a problem from a json file
  using lcs_solver::algorithms::BaseSolution;
  using lcs_solver::problems::BaseProblem;
  using lcs_solver::problems::ProblemFactory;
  const auto& file = parser.GetOption("-f");
  std::cout << "Generate a new problem from file " + file << std::endl;
  const auto problem_uniq_ptr = ProblemFactory::CreateFromJson(file);
  const std::unique_ptr<BaseSolution> sol_ptr = problem_uniq_ptr->ExecuteAlgo();
  std::cout << "Solution:\n" << *sol_ptr;
}

void HandleInputMode(const InputParser& parser) {
  // Return early if not in File Mode
  if (!parser.HasOption("-i")) {
    return;
  }

  // Namespace aliasing
  using lcs_solver::problems::BaseProblem;
  using lcs_solver::problems::ProblemFactory;
  using lcs_solver::util::InputParser;
  using lcs_solver::util::ReadChar;
  using lcs_solver::util::ReadStdString;

  char c;
  char default_val = 'n';
  std::string file_name = "test.json";
  std::unique_ptr<BaseProblem> problem;
  PrintOptionsMessage();

  do {
    c = ReadChar("Enter operation", default_val);
    switch (c) {
      case 'n': {
        std::cout << "Generate a new problem by using the terminal."
                  << std::endl;
        problem = ProblemFactory::CreateFromDialog();
        default_val = 'p';
        file_name = SuggestFileName();
        break;
      }
      case 'o': {
        file_name = ReadStdString("Open file", file_name);
        problem = ProblemFactory::CreateFromJson(file_name);
        default_val = 'p';
        break;
      }
      case 's': {
        if (!problem) throw std::runtime_error("Problem is not defined");
        file_name = ReadStdString("Save to file", file_name);
        ProblemFactory::SaveToJson(file_name, problem.get());
        default_val = 'o';
        break;
      }
      case 'p': {
        if (!problem) {
          std::cout << "Empty\n";
        } else {
          std::cout << problem->DebugString() << std::endl;
        }
        default_val = 'x';
        break;
      }
      case 'x': {
        auto solution_ptr = problem->ExecuteAlgo();
        std::cout << "Execute Algorithm. Solution:\n"
                  << *solution_ptr << "\n";
        default_val = 's';
        break;
      }
      case 'e': {
        std::cout << "Exit program." << std::endl;
        break;
      }
      default: {
        PrintOptionsMessage();
      }
    } /*switch*/
    // std::cout << std::endl;
  } while (c != 'e');
}

void PrintUsageMessage(const InputParser& parser) {
  std::ostringstream oss;
  oss << "Usage: " << parser.GetProgramName() << " [-h] [-f file]\n";
  oss << "Options:\n";
  oss << "  -f file  : Solve the problem described in the specified file.\n";
  oss << "  -i       : Interactively create a the problem in the terminal.\n";
  oss << "  -h       : Display this usage message.\n";
  std::cout << oss.str();
}

void PrintOptionsMessage() {
  std::ostringstream oss;
  oss << "The options in the main loop are" << "\n";
  oss << "n: create a new problem in the terminal dialog" << "\n";
  oss << "o: open a problem from a file" << "\n";
  oss << "s: save a problem to a file" << "\n";
  oss << "x: execute algorithm (assumes a valid problem)" << '\n';
  oss << "p: prints the current problem" << '\n';
  oss << "h: displays this help" << '\n';
  oss << "e: exit" << '\n';
  std::cout << oss.str();
}

std::string SuggestFileName() {
  const auto now = std::chrono::system_clock::now();
  std::time_t const now_c = std::chrono::system_clock::to_time_t(now);
  std::tm const local_tm = *std::localtime(&now_c);
  std::ostringstream oss;
  oss << "problem_" << std::put_time(&local_tm, "%Y-%m-%d_%H-%M-%S");
  // oss << "test";
  oss << ".json";
  return oss.str();
}

}  // namespace solver