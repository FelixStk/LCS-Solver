/*******************************************************************************
 * @file main.cc
 * @author Steinkopp:Felix
 * @version 1.0.0
 * @brief Demonstrates lcs_solver in the context of gene tree inference
 * @see dna::Run and dna::MainLoop
 ******************************************************************************/

#include "dna.h"

int main(const int argc, char* argv[]) {
  try {
    dna::Run(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}