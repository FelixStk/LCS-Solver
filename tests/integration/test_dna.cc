/*******************************************************************************
 * @file test_dna.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests to check the dna demo program
 ******************************************************************************/

#include "demo/dna/dna.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"

namespace {

TEST(IntegrationDemoDna, Version) {
  auto command = std::to_array({
      "dna",
      "-V",
      "-P",
      "../demo/dna/swisstree",
  });
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  // Redirect std::cout to check whether the version message starts with "dna"
  std::ostringstream captured_output;
  std::streambuf* original_buf = std::cout.rdbuf();
  std::cout.rdbuf(captured_output.rdbuf());
  EXPECT_NO_THROW(dna::Run(argc, argv));
  std::cout.rdbuf(original_buf);
  std::string output = captured_output.str();
  EXPECT_NE(output.find("dna"), std::string::npos);
}

TEST(IntegrationDemoDna, JsonGapVector) {
  auto command = std::to_array({
      "dna",
      "-p", "LCS2_MC",  // ProblemType Identifier
      "-t", "LLCS2_MC", // AlgoType Identifier
      "-g", "./integration/res/gap.json",
      "-P",
      "../demo/dna/swisstree",
  });
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());
  const dna::InputParser parser(argc, argv);
  const std::unique_ptr<dna::GapVector> v = dna::GetGapVector(&parser);
  EXPECT_EQ(v->size(), 3);
}

TEST(IntegrationDemoDna, JsonAlphabet) {
  auto command = std::to_array({
      "dna",
      "-p", "LCS2_MC",  // ProblemType Identifier
      "-t", "LLCS2_MC", // AlgoType Identifier
      "-g", "./integration/res/gap.json",
      "-a", "./integration/res/alph.json",
      "-l", "./integration/res/gap.json",
      "-r", "./integration/res/gap.json",
      "-P",
      "../demo/dna/swisstree",
  });
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());
  const dna::InputParser parser(argc, argv);
  const std::shared_ptr<dna::SymbolVector> alph = dna::GetAlphabet(&parser);
  EXPECT_EQ(alph->size(), 3);
  EXPECT_EQ(alph->at(0), U'A');
  EXPECT_EQ(alph->at(1), U'B');
  EXPECT_EQ(alph->at(2), U'C');
}

TEST(IntegrationDemoDna, Output) {
  auto command = std::to_array({
      "dna",
      "-p", "LCS_Classic",  // ProblemType Identifier
      "-t", "LLCS2_STD_FL", // AlgoType Identifier
      "-f", "./integration/res/test.fasta", "./integration/res/test.fasta",
      "-g", "./integration/res/gap.json",
      "-a", "./integration/res/alph.json",
      "-l", "./integration/res/gap.json",
      "-r", "./integration/res/gap.json",
      "-P",
      "../demo/dna/swisstree",
  });
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  // Redirects cout to check whether the output message contains distance 0
  std::ostringstream captured_output;
  std::streambuf* original_buf = std::cout.rdbuf();
  std::cout.rdbuf(captured_output.rdbuf());
  EXPECT_NO_THROW(dna::Run(argc, argv));
  std::cout.rdbuf(original_buf);
  std::string output = captured_output.str();
  EXPECT_NE(output.find("0"), std::string::npos) << captured_output.str();
}

TEST(IntegrationDemoDna, Algorithms) {
  std::array prob_alo_pairs = {
    std::make_pair("LCS_Classic", "LLCS2_STD_FL"),
    std::make_pair("LCS2_MC", "LLCS2_MC"),
    std::make_pair("LCS2_MC_INC", "LLCS2_MC_INC"),
    std::make_pair("LCS2_MC_1C", "LLCS2_MC_1C"),
    std::make_pair("LCS2_MC_O1C_SYNC", "LLCS2_MC_O1_SYNC"),
    std::make_pair("LCS2_Sigma_R", "LLCS2_SR_MQ"),
    std::make_pair("LCS2_Sigma_R", "LLCS2_SR_RMQ"),
    std::make_pair("LCS2_Sigma_L", "LLCS2_SL_R_LLCS2_SR_MQ"),
    std::make_pair("LCS2_Sigma_L", "LLCS2_SL_R_LLCS2_SR_RMQ"),
    std::make_pair("LCS2_SIGMA", "LLCS2_SA_MQ"),
    std::make_pair("LCS2_SIGMA", "LLCS2_SA_RMQ"),
  };
  for (const auto& [prob, algo] : prob_alo_pairs) {
    auto command = std::to_array({
      "dna",
      "-p", prob,  // ProblemType Identifier
      "-t", algo, // AlgoType Identifier
      "-f", "./integration/res/test.fasta", "./integration/res/test.fasta",
      "-g", "./integration/res/gap.json",
      "-a", "./integration/res/alph.json",
      "-l", "./integration/res/gap.json",
      "-r", "./integration/res/gap.json",
      "-P",
      "../demo/dna/swisstree",
    });
    constexpr auto argc = static_cast<int>(command.size());
    const auto argv = const_cast<char**>(command.data());
    std::ostringstream captured_output;
    std::streambuf* original_buf = std::cout.rdbuf();
    std::cout.rdbuf(captured_output.rdbuf());
    EXPECT_NO_THROW(dna::Run(argc, argv)) << prob << " " << algo;
    std::cout.rdbuf(original_buf);
  }
}

}  // namespace
