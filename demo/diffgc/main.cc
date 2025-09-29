/*******************************************************************************
 * @file main.cc
 * @author Steinkopp:Felix
 * @version 1.0.1
 * @brief Runs the diff algorithm with a gap lcs algorithm
 ******************************************************************************/

#include "diffgc.h"

namespace fs = std::filesystem;

int main(const int argc, char *argv[]) {
  try {
    diffgc::Run(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}