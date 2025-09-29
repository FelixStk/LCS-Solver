/*******************************************************************************
 * @file test_diffgc.cc
 * @author Felix Steinkopp
 * @version 1.0
 * @brief GTests to the integration of lcs_solver in the diffgc demo
 ******************************************************************************/

#include <regex>

#include "demo/diffgc/diffgc.h"
#include "gmock/gmock-more-matchers.h"
#include "gtest/gtest.h"

namespace {

TEST(IntegrationDemoDiffgc, Version) {
  auto command = std::to_array({
      "diffgc",  // placeholder for the binary name (without meaning here)
      "--version",
  });
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  // Redirect std::cout
  const std::ostringstream output_stream;
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  EXPECT_NO_THROW(diffgc::Run(argc, argv));
  std::cout.rdbuf(cout_buf);

  std::string output = output_stream.str();
  EXPECT_NE(output.find("diffgc"), std::string::npos);
}

TEST(IntegrationDemoDiffgc, Files) {
  auto command = std::to_array({
      "diffgc",
      "./integration/res/file1.cc",
      "./integration/res/file2.cc",
  });
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  const std::ostringstream output_stream;
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  EXPECT_NO_THROW(diffgc::Run(argc, argv));
  std::cout.rdbuf(cout_buf);

  const std::string expected_output = R"(+ | #include <span>
+ | 
= | int calculateTotal(const std::span<int> items) {
= |   int sum = 0;
- |   for (int number : items) {
- |     sum += number;
- |   }
+ |   if (items.empty()) return sum;
+ |   sum = items.front() + calculateTotal(items.subspan(1));
= |   return sum;
= | })";

  const std::string output = output_stream.str();
  EXPECT_EQ(output, expected_output) << "actual:\n"
                                     << output << "\n"
                                     << "expected:\n"
                                     << expected_output;
}

TEST(IntegrationDemoDiffgc, Load_Constraint_MC) {
  auto command = std::to_array({"diffgc",
                                "./integration/res/file1.cc",
                                "./integration/res/file2.cc",
                                "-j",
                                "./integration/res/config_mc.json"});
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  const std::ostringstream output_stream;
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  EXPECT_NO_THROW(diffgc::Run(argc, argv));
  std::cout.rdbuf(cout_buf);
}

TEST(IntegrationDemoDiffgc, Gnu) {
  auto command = std::to_array({"diffgc",
                                "./integration/res/file1.cc",
                                "./integration/res/file2.cc",
                                "--gnu"});
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  const std::ostringstream output_stream;
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  EXPECT_NO_THROW(diffgc::Run(argc, argv));
  std::cout.rdbuf(cout_buf);

  const std::string raw_output = output_stream.str();
  const std::regex timestamp_regex(R"((---|\+\+\+) .*\t.*\n)");
  const std::string normalized_output = std::regex_replace(raw_output, timestamp_regex, "$1 <path>\n");

  const std::string expected_output =
      R"(--- <path>
+++ <path>
@@ -1,0 +1,2 @@
+#include <span>
+
@@ -3,3 +5,2 @@
-  for (int number : items) {
-    sum += number;
-  }
+  if (items.empty()) return sum;
+  sum = items.front() + calculateTotal(items.subspan(1));
)";

  EXPECT_EQ(normalized_output, expected_output) << "actual:\n"
                                                << normalized_output << "\n"
                                                << "expected:\n"
                                                << expected_output;
}

TEST(IntegrationDemoDiffgc, Interaction) {
  auto command = std::to_array({"diffgc",
                                "./integration/res/file1.cc",
                                "./integration/res/file2.cc",
                                "-i"});
  constexpr auto argc = static_cast<int>(command.size());
  const auto argv = const_cast<char**>(command.data());

  std::ostringstream output_stream;
  std::istringstream input_stream("no\n no\n");
  std::streambuf* cout_buf = std::cout.rdbuf();
  std::streambuf* cin_buf = std::cin.rdbuf();
  std::cout.rdbuf(output_stream.rdbuf());
  std::cin.rdbuf(input_stream.rdbuf());
  EXPECT_NO_THROW(diffgc::Run(argc, argv));
  std::cout.rdbuf(cout_buf);
  std::cin.rdbuf(cin_buf);
  std::string output = output_stream.str();
}

}  // namespace