/*******************************************************************************
 * @file main.cc
 * @author Steinkopp:Felix
 * @version 2.0
 * @brief Implementation of the main entry point of the solver demonstration
 ******************************************************************************/

#include "solver.h"

#include <iostream>

int main(const int argc, char* argv[]) {
  try {
    solver::Run(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
