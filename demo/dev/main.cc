/*******************************************************************************
 * @file devStuff.cc
 * @brief main file of the devStuff programm for testing random code peaces
 ******************************************************************************/

#include <format>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>

#include "algorithms/LCS/LCS2_RT.h"
#include "constraints/BaseConstraint.h"
#include "constraints/ConstraintMap.h"
#include "constraints/local/Constraint_MC.h"
#include "constraints/local/Constraint_MC_INC.h"
#include <lcs_solver/util/CommonTypes.h>
#include <lcs_solver/util/StrPtrVector.h>
#include <lcs_solver/util/UnicodeRange.h>
#include "problems/BaseProblem.h"

int main(int argc, char *argv[]) {
  using lcs_solver::constraints::ConstraintType;
  using lcs_solver::constraints::local::Constraint_MC;
  using lcs_solver::constraints::local::Constraint_MC_INC;
  using lcs_solver::util::Symbol;
  using lcs_solver::util::String;
  using lcs_solver::util::uint;
  using lcs_solver::util::UnicodeRange;
  using lcs_solver::util::StrPtrVector;
  using lcs_solver::constraints::ConstraintMap;
  using lcs_solver::constraints::BaseConstraint;

  using GapVector = std::vector<std::pair<size_t, size_t>>;
  GapVector gap_vector;
  gap_vector.emplace_back(0, 0);
  gap_vector.emplace_back(1, 1);

  //=== Constraint Ptr
  // std::shared_ptr<BaseConstraint> ptr = std::make_shared<Constraint_MC>(gap_vector);
  // nlohmann::json j;
  // lcs_solver::constraints::to_json(j,ptr);
  // std::cout << j.dump(4) << "\n \n =============== \n";
  //
  // std::string name = j["type_"];
  // ConstraintType t = lcs_solver::constraints::GetMapNameToType().at(name);
  // std::cout << "t = " << lcs_solver::constraints::GetMapTypeToName().at(t) << "\n";
  //
  // std::shared_ptr<BaseConstraint> ptr2;
  // lcs_solver::constraints::from_json(j,ptr2);
  // std::cout << ptr2->DebugString() << '\n';


  //=== Constraint Map
  ConstraintMap map;
  map.emplace(ConstraintType::MC, std::make_shared<Constraint_MC>(gap_vector));
  map.emplace(ConstraintType::MC_INC, std::make_shared<Constraint_MC_INC>(gap_vector));
  nlohmann::json j;
  j["map"] = map;
  std::cout << j.dump(4) << "\n \n =============== \n";

  ConstraintMap map2;
  lcs_solver::constraints::from_json(j["map"],map2);
  for (const auto &[key, value]: map2) {
    std::cout << "key = " << lcs_solver::constraints::GetMapTypeToName().at(key) << "\n";
    std::cout << "value = " << value->DebugString() << "\n";
    std::cout << "\n";
  }
  // std::cout << j.dump(4) << "\n \n =============== \n";
  // int n = j["map"].size();
  // std::cout << "n = " << n << '\n';
  // ConstraintMap map2;
  // for (auto &[key, value]: j["map"].items()) {
  //
  //   ConstraintType t = value[0].get<ConstraintType>();
  //   std::cout << "t = " << lcs_solver::constraints::GetMapTypeToName().at(t) << "\n";
  //   std::shared_ptr<BaseConstraint> ptr;
  //
  //   std::shared_ptr<BaseConstraint> ptr2;
  //   lcs_solver::constraints::from_json(value[1],ptr2);
  //   std::cout << ptr2->DebugString() << '\n';
  //   std::cout << "\n";
  //   map2.emplace(t, ptr2);
  // }

  // ConstraintMap map2 = j["map_"].get<ConstraintMap>();
  //

  //=== Files
  // std::vector<std::filesystem::path> v;
  // v.emplace_back("/var/log/folder1/folder2/log.log");
  // v.emplace_back("/var/log/folder1/folder2/log2.log");
  // nlohmann::json j;
  // j["Files"] = v;
  // std::cout << j.dump(4) << '\n';
  // std::cout << std::endl;


  //=== StrPtrVector JSON
  // StrPtrVector spv = lcs_solver::util::StrPtrVectorFrom({"asdf1", "asdf2"});
  // nlohmann::json j;
  // auto strings = lcs_solver::util::StrPtrVectorTo<std::vector<std::string>>(spv);
  // for (const std::string &p : strings) {
  //   std::cout << p << '\n';
  // }
  // j["Strings"] = lcs_solver::util::StrPtrVectorTo<std::vector<std::string>>(spv);
  // std::cout << j.dump(4) << '\n';
  // std::cout << std::endl;



  // nlohmann::json j;
  // std::vector<std::pair<size_t, size_t>> gc1 = {{0,0},{1,1,}};
  // Constraint_MC_INC c1(gc1);
  // ;
  // j["constraint"] = c1;
  // std::cout << j.dump(4) << '\n';
  // Constraint_MC_INC c2 = j["constraint"].get<Constraint_MC_INC>();
  // // std::cout << c2.DebugString() << '\n';
  // // std::cout << c2.GetName() << '\n';
  // std::cout << std::endl;

  // === ConstraintType
  // auto ct = ConstraintType::MC_INC;
  // nlohmann::json j;
  // j["ct"] = ct;
  // ConstraintType ct2 = j["ct"];
  // std::cout << j.dump(4) << '\n';
  // std::cout << lcs_solver::constraints::GetMapTypeToName().at(ct2) << '\n';
  // std::cout << std::endl;

  // std::shared_ptr<Constraint_MC_INC> p = std::make_shared<Constraint_MC_INC>(gc1);
  // ConstraintMap map = {std::make_pair(ConstraintType::MC_INC,p)};

  // auto res = lcs_solver::constraints::Get<Constraint_MC_INC>(map);
  // std::cout << res->DebugString() << std::endl;

  // result: {20, 40, 10, 30}

  // std::vector<int> input = {10, 20, 30, 40};
  // std::vector<std::size_t> perm = {2, 0, 3, 1}; // Move input[0] to position 2, input[1] to position 0, etc.
  //
  // auto result = lcs_solver::util::ApplyPermutation(input, perm);
  // for (const auto &i : result) {
  //   std::cout << i << std::endl; // result: 20, 40, 10, 30
  // }


  // ::lcs_solver::util::ConfigureUTF8Console();
  // StrPtrVector spv = lcs_solver::util::StrPtrVectorFrom({"asdf", "asdf"});
  // for (const auto &p : spv) {
  //   // std::cout << lcs_solver::util::to_string(*p) << std::endl;
  //   using lcs_solver::util::operator<<;
  //   std::cout << (*p) << std::endl;
  // }
  // auto v = lcs_solver::util::StrPtrVectorTo<std::vector<std::string>>(spv); // OK
  // for (const auto &s : v) {
  //   std::cout << s << std::endl;
  // }

  // UnicodeRange alph1(UnicodeRange::Alphabet::kLatinLowerUpper);
  // auto max = alph1.size();
  // auto vec = alph1.uintAlphabetVector();
  // std::cout << "alpha1.size(): " <<  max << std::endl;
  // uint32_t c =alph1[max-1];
  // std::cout << "last Symbol: " << lcs_solver::util::to_string(static_cast<char32_t>(c)) << std::endl;
  // std::cout << "vec: " << vec[max-1] << std::endl;
  //
  // std::vector<String> vec2 = {lcs_solver::util::to_String<Symbol>("TestString")};
  // std::cout << lcs_solver::util::to_string(vec2[0]) << std::endl;

  // ::lcs_solver::util::Symbol c = 'X';
  // const char32_t y = alph1[i-1]; // 66
  // const std::u32string y = U"a😊sd😊f";
  // std::string s = lcs_solver::util::to_string(y);
  // auto x = lcs_solver::util::char32_from_utf8(s);
  // std::cout << s.size() << std::endl;
  // int i = 0;
  // for (char c : s) {
  //   std::cout << i++ << " " << static_cast<int>(c) << std::endl;
  // }
  // fmt::print("Smile: {}",s);

  // ::lcs_solver::util::configureUTF8Console();
  // SetConsoleOutputCP(CP_UTF8);
  //     setvbuf(stdout, nullptr, _IOFBF, 1000);


  // std::cout << "s.c_str(): " << s.c_str() << std::endl;
  // std::cout << static_cast<uint>(x)<< std::endl;

  // std::string s = lcs_solver::util::to_utf8(y);
  // std::array<char, 5> smilie = {-16, -97, -104, -118, '\0'};
  // std::string question = "\xE2\x9D\x93";
  // std::string bee = "\U0001F41D"; // OR char* bee = "\U0001F41D";
  // std::cout << "To " + bee + " or not to " + bee + " that is the " + question<< std::endl;
  // std::cout << smilie.data() << std::endl;


  // printf("😊\n");
  // std::print("{}\n", "ф");

  // const char32_t y = U'😊'; // 128522
  // std::string utf8;
  // if (y <= 0x7F) {
  //   // 1-byte UTF-8
  //   utf8.push_back(static_cast<char>(y));
  // } else if (y <= 0x7FF) {
  //   // 2-byte UTF-8
  //   utf8.push_back(static_cast<char>(0xC0 | ((y >> 6) & 0x1F)));
  //   utf8.push_back(static_cast<char>(0x80 | (y & 0x3F)));
  // } else if (y <= 0xFFFF) {
  //   // 3-byte UTF-8
  //   utf8.push_back(static_cast<char>(0xE0 | ((y >> 12) & 0x0F)));
  //   utf8.push_back(static_cast<char>(0x80 | ((y >> 6) & 0x3F)));
  //   utf8.push_back(static_cast<char>(0x80 | (y & 0x3F)));
  // } else if (y <= 0x10FFFF) {
  //   // 4-byte UTF-8 (for emoji and other supplementary characters)
  //   utf8.push_back(static_cast<char>(0xF0 | ((y >> 18) & 0x07)));
  //   utf8.push_back(static_cast<char>(0x80 | ((y >> 12) & 0x3F)));
  //   utf8.push_back(static_cast<char>(0x80 | ((y >> 6) & 0x3F)));
  //   utf8.push_back(static_cast<char>(0x80 | (y & 0x3F)));
  // }
  //
  // std::cout << utf8 << std::endl;

  return 0;
}