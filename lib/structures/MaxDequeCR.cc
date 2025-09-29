/*********************************************************************
 * @file MaxDequeCR.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Data Structure of LCS_MC_INC_E. Consists of arrays of Queues.
 ********************************************************************/

#include "structures/MaxDequeCR.h"
// #include "util/Logger.hpp"

#include <algorithm>
#include <deque>
#include <functional>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

namespace {
using uint = ::lcs_solver::structures::MaxDequeCR::uint;
// using Logger = ::lcs_solver::util::Logger;
}

namespace lcs_solver::structures {

/*******************************************************************************
 * Constructor for MaxDequeCR objects
 * @param mat reference to dp matrix
 * @param phi Predicate Phi(i,j) == 1, if s1[i]==s2[j] with (i,j) in [0:m-1]x[0:n-1]
 * @param tl function from [m] to [m]: tl(i) := lower bound for the ith gap
 * @param tu function from [m] to [m]: tu(i) := upper bound for the ith gap
 * @note tl and tu should together fulfill the increasing gap property
 *  - tu(i) > tl(i) for all i in [m]
 *  - tu(i-1) - tl(i-1) <= tu(i) - tl(i) for all i in [2:m]
 ******************************************************************************/
MaxDequeCR::MaxDequeCR(
    const Matrix &mat,
    std::function<uint(uint, uint)> phi,
    std::function<uint(uint)> tl,
    std::function<uint(uint)> tu
) : M(mat),
    m(mat.size()),
    n(!mat.empty() ? mat[0].size() : 0),
    Phi(std::move(phi)), Tl(std::move(tl)), Tu(std::move(tu)) {
  // Logger::Error() << "New MaxDequeCR Object! xxx" << Tl(0) << " \n";
  S = Matrix(m, std::vector<uint>(n, 0));
  for (uint i = 0; i < m; ++i) {
    for (uint j = 0; j < n; j++) {
      if (Phi(i, j) == 1) {  //s1[i] == s2[j] with i in [0:m-1] and j in [0:n-1]
        // uint ip = Tl(0) + i + 1;
        // uint jp = Tl(0) + j + 1;
        // if(ip < m && jp < n){
        //   Logger::Error() << "init insert " << i << " " << j << " (" << ip << " " << jp << ") " << 1 << "\n";
        //   S[ip][jp] = 1;
        // }
        //S[i][j] = 1;
      }
    }
  }
  C = Deques(n, std::deque<DequeItem>());
  R = Deques(m, std::deque<DequeItem>());
  // Logger::Error() << DebugString();
}

/*******************************************************************************
 * update
 * @param i current row in dp loop (i in [0:m-1])
 * @param j current column in dp loop (j in [0:n-1])
 ******************************************************************************/
void MaxDequeCR::update(MaxDequeCR::uint i, MaxDequeCR::uint j) {
  // Expire
  while (!C[j].empty() && C[j].front().first < i) {
    C[j].pop_front();
  }
  while (!R[i].empty() && R[i].front().first < j) {
    R[i].pop_front();
  }

  // Update
  // Logger::Error() << "a ijValues " << i << " " << j << "\n";
  // Logger::Error() << DebugString();
  if (S[i][j] > 0) {
    while (!C[j].empty() && *(C[j].back().second) <= S[i][j]) {
      C[j].pop_back();
    }
    if(S[i][j]-2 < std::min<uint>(m,n))
      C[j].emplace_back(Tu(S[i][j]-2) - Tl(S[i][j]-2) + i, &S[i][j]);
    else
      C[j].emplace_back(i, &S[i][j]); // Tu(S[i][j]-1) is not defined
  }
  if(!C[j].empty()){
    auto *p = C[j].front().second;
    while (!R[i].empty() && *(R[i].back().second) <= *p) {
      R[i].pop_back();
    }
    if(*p - 2 < std::min<uint>(m,n))
      R[i].emplace_back(Tu(*p-2) - Tl(*p-2) + j, p);
  }
  // Logger::Error() << "b ijValues " << i << " " << j << "\n";
  // Logger::Error() << DebugString();
}

/*******************************************************************************
 * insert
 * @param i current row in dp loop (i in [0:m-1])
 * @param j current column in dp loop (j in [0:n-1])
 ******************************************************************************/
void MaxDequeCR::insert(MaxDequeCR::uint i, MaxDequeCR::uint j) {
  const uint p = M[i][j];
  if(p < std::min<uint>(m,n)){
    // const uint ip = std::min<uint>(Tl(p) + i + 1, m-1);
    // const uint jp = std::min<uint>(Tl(p) + j + 1, n-1);
    const uint ip = Tl(p-1) + i + 1;
    const uint jp = Tl(p-1) + j + 1;
    // Logger::Error() << "insert " << i << " " << j << " (" << ip << " " << jp << ") " << p+1 << "\n";
    if(ip < m && jp < n)
      S[ip][jp] = std::max<uint>(S[ip][jp],p+1);
  }
}

/*******************************************************************************
 * extractMax
 * @param i current row in dp loop (i in [0:m-1])
 * @param j current column in dp loop (j in [0:n-1])
 ******************************************************************************/
MaxDequeCR::uint MaxDequeCR::extractMax(MaxDequeCR::uint i, MaxDequeCR::uint j) {
  if(R[i].empty())
    return 1;//Phi(i,j);
  
  // Logger::Error() << "extraction set M " << i << " " << j << " to " << *R[i].front().second << "\n";
  return *R[i].front().second;
}

/*******************************************************************************
 * DebugString
 * @return std::string that contains the Matrix S and every deque in the object
 ******************************************************************************/
std::string MaxDequeCR::DebugString() {
  std::ostringstream oss;
  oss << "Matrix S: \n";
  for (const auto &row : S) {
    for (const auto &val : row) {
      oss << val << " ";
    }
    oss << "\n";
  }
  for(uint x = 0; x< C.size(); ++x){
    oss << "C[" << x << "]: ";
    for(const auto& [t,p]: C[x]){
      oss << t << " " << *p << "|";
    }
    oss << "\n";
  }
  for(uint x = 0; x< R.size(); ++x){
    oss << "R[" << x << "]: ";
    for(const auto& [t,p]: R[x]){
      oss << t << " " << *p << "|";
    }
    oss << "\n";
  }
  return oss.str();
}

}  // namespace lcs_solver::structures