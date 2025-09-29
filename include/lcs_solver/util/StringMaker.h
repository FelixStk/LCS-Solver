// #ifndef LCS_SOLVER_STRING_MAKER_HPP
// #define LCS_SOLVER_STRING_MAKER_HPP
//
// #include <memory>
// #include <set>
// #include <vector>
//
// #include "CommonTypes.h"
// #include "algorithms/solutions/Points.h"
// #include "algorithms/solutions/UnsignedSolution.h"
//
// namespace lcs_solver::util {
//
// class StringMaker {
//  public:
//   // Setters for creation environment.
//   void setLengths(const std::vector<size_t>& lengths);
//   void setLLCS(size_t llcs);
//   void setShrSymbols(const std::vector<Symbol>& shared);
//   void setUnqSymbols(const std::vector<Symbol>& uniques);
//   void setGaps(const std::vector<std::pair<uint, uint>>& gap);
//   void setBoundNLCS(std::pair<uint, uint> bounds);
//
//   // Getters for the creation.
//   std::vector<std::shared_ptr<const String>> getStrings();
//   algorithms::solutions::UnsignedSolution getLLCS();
//   algorithms::solutions::UnsignedSolution getNLCS();
//   std::set<algorithms::solutions::Points> getLCS();
//
//   // Methods for adding additional Longest Common Subsequences into the strings.
//   std::vector<uint> selectVertex();
//   bool checkVertex(std::vector<uint> vertex);
//   void processVertex();
//
//   // Interface for running.
//   void run();
//
//  private:
//   // Aliases for convenience.
//   using Points = ::lcs_solver::algorithms::solutions::Points;
//   using UnsignedSolution = ::lcs_solver::algorithms::solutions::UnsignedSolution;
//   using StrPtrVector = std::vector<std::shared_ptr<const String>>;
//   using Pair = std::pair<uint, uint>;
//
//   // Data members.
//   size_t mLLCS = 0;
//   size_t mNLCS = 0;
//   std::vector<size_t> mStrLength;
//   std::vector<Pair> mGap;
//   std::vector<Symbol> mSharedSymbol;
//   std::vector<Symbol> mUniqueSymbol;
//   std::set<Symbol> mUsedSymbols;
//   std::vector<std::vector<Symbol>> mStr;
//   std::vector<std::vector<std::vector<size_t>>> mEmb;
//   Pair mBoundNLCS = Pair({0, 20});
//   StrPtrVector mStringResult;
//   std::set<Points> mPointsSet;
//
//   // Internal helper methods.
//   void initialize();
//   void addFirstLCS();
//   void construct();
// };
//
// } // namespace lcs_solver::util
//
//
// #endif  // LCS_SOLVER_STRING_MAKER_HPP
