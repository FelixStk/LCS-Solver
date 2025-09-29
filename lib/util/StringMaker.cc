// /*******************************************************************************
//  * @file StringMaker.cc
//  * @author Steinkopp:Felix
//  * @version 1.0
//  * @brief Impl. of a random string generator
//  ******************************************************************************/
// #include "util/StringMaker.h"
//
// #include <algorithm>
// #include <cassert>
// #include <utility>
//
// #include "util/RndUtility.h"
//
// namespace lcs_solver::util {
// //== Setters ===========================================================================================================
// void StringMaker::setLengths(const std::vector<size_t> &lengths) {
//   mStrLength = lengths;
// }
//
// void StringMaker::setGaps(const std::vector<StringMaker::Pair> &gapVector) {
//   mGap = gapVector;
// }
//
// void StringMaker::setLLCS(size_t llcs) {
//   mLLCS = llcs;
// }
//
// void StringMaker::setBoundNLCS(const StringMaker::Pair bounds) {
//   mBoundNLCS = bounds;
// }
//
// void StringMaker::setShrSymbols(const std::vector<Symbol> &shared) {
//   mSharedSymbol = shared;
// }
//
// void StringMaker::setUnqSymbols(const std::vector<Symbol> &uniques) {
//   mUniqueSymbol = uniques;
// }
//
// //== Solution Getters ==================================================================================================
//
// StringMaker::StrPtrVector StringMaker::getStrings() {
//   return mStringResult;
// }
//
// StringMaker::UnsignedSolution StringMaker::getLLCS() {
//   if (mBoundNLCS.second == 0)
//     return  UnsignedSolution(0);
//   return UnsignedSolution(mLLCS);
// }
//
// StringMaker::UnsignedSolution StringMaker::getNLCS() {
//   return UnsignedSolution(mNLCS);
// }
//
// std::set<StringMaker::Points> StringMaker::getLCS() {
//   return mPointsSet;
// }
//
// //== String Construction ===============================================================================================
//
// void StringMaker::initialize() {
//   assert(mUniqueSymbol.size()
//              >= mStrLength.size()); // each string should have a unique symbol available
//   std::sort(mStrLength.begin(), mStrLength.end());
//   for (size_t i = 0; i < mStrLength.size(); i++) {
//     mStr[i] = std::vector<Symbol>(mStrLength[i], mUniqueSymbol[i]);
//     mUsedSymbols.insert(mUniqueSymbol[i]);
//   }
//   mNLCS = 0;
//   mEmb.clear();
// }
//
// void StringMaker::addFirstLCS() {
//   std::vector<algorithms::Embedding> sol;
//   sol.reserve(mStrLength.size());
//
//   std::vector<std::vector<size_t>>
//       solution; // solution[i] represents the embeddings of the LCS into mStr[i]
//   solution.reserve(mStr.size());
//
//   // Go through strings and generate a random first LCS-Embedding.
//   for (auto &s : mStr) {
//     if (mGap.empty()) {
//       solution.emplace_back(RndUtility::combNoRepInRange(mLLCS,
//                                                          {0, s.size() - 1}));
//       size_t k = 0;
//       // Modify Strings accordingly
//       for (const auto &i : solution.back()) {
//         s[i] = mSharedSymbol[k];
//         k++;
//       }
//     } else {
//
//     }
//   }
//
//   // Remember what Symbols are used
//   for (size_t k = 0; k < mLLCS; k++) {
//     mUsedSymbols.insert(mSharedSymbol[k]);
//   }
//   mEmb.push_back(solution);
//   mNLCS++;
// }
//
// void StringMaker::construct() {
//   // create strings in shared pointers
//   mStringResult.clear();
//   mStringResult.reserve(mStrLength.size());
//   for (const auto &str : mStr) {
//     auto ptr = std::make_shared<const String>(str.begin(), str.end());
//     mStringResult.emplace_back(ptr);
//   }
//
//   // Process 3D Matrix into a set of EmbeddingsSolutions
//   mPointsSet.clear();
//   for (const auto &solution : mEmb) {
//     using algorithms::Embedding;
//     std::vector<Embedding> tempSol;
//     tempSol.reserve(mStringResult.size());
//     for (size_t i = 0; i < mStringResult.size(); i++) {
//       auto e = Embedding(std::shared_ptr<String>(), solution[i]);
//       tempSol.emplace_back(e);
//     }
//     mPointsSet.insert(Points(mStringResult, tempSol));
//   }
// }
//
// void StringMaker::run() {
//   initialize();
//   if (mBoundNLCS.second > 0)
//     addFirstLCS();
//   construct();
// }
// }