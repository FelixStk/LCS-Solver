#ifndef DEMO_SOLVER_H_
#define DEMO_SOLVER_H_

#include <string>

#include "lcs_solver/util/InputParser.h"

namespace solver {

using InputParser   = lcs_solver::util::InputParser;

/*******************************************************************************
 * @brief Main Program. Solves a user specific LCS problem.
 * @details Solves an LCS problem based on user input. It either gets the problem
 *          from a file (FromFileGenerator) or from the terminal dialog
 * @param[in] argc Number of command line arguments
 * @param[in] argv Array of command line argument strings
 ******************************************************************************/
void Run(int argc, char *argv[]);

/*******************************************************************************
 * @brief Informs the user about the program
 * @param[in] parser Parses terminal arguments as flags or options
 ******************************************************************************/
void HandleHelp(const InputParser& parser);

/*******************************************************************************
 * @brief Handles the `-f` flag for reading a problem from a file
 * @param[in] parser Parses terminal arguments as flags or options
 ******************************************************************************/
void HandleFileMode(const InputParser& parser);

/*******************************************************************************
 * @brief Handles the `-i` flag for interactive problem creation
 * @param[in] parser Parses terminal arguments as flags or options
 ******************************************************************************/
void HandleInputMode(const InputParser& parser);

/*******************************************************************************
 * @brief Prints a usage message to the standard out stream
 * @param[in] parser Parses terminal arguments as flags or options
 ******************************************************************************/
void PrintUsageMessage(const InputParser& parser);

/*******************************************************************************
 * getTermOptionsStr
 * @details Prints text describing the options in the main loop of the program
 ******************************************************************************/
void PrintOptionsMessage();

/*******************************************************************************
 * @brief Suggests a file name based on the current clock time
 * @return String in the format "problem_%Y-%m-%d_%H-%M-%S.txt"
 ******************************************************************************/
std::string SuggestFileName();

}  // namespace solver

#endif  // DEMO_SOLVER_H_
