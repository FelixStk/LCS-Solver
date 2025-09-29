#ifndef LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_RT_TPP_
#define LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_RT_TPP_



//// Iterator Stuff ////////////////////////////////////////////////////////////

/*******************************************************************************
 * getIterator
 * @param i uint current row position
 * @param j uint current column position
 * @param l the value of the llcs
 * @return RangeTree2D<uint,uint>::Iterator with points (i,j) such that
 *  svv[0][i] == svv[1][j] == s and s is lth symbol in the lcs
 ******************************************************************************/
template<OutputType T>
std::tuple<typename Reconstructor_RT<T>::uint,
           typename Reconstructor_RT<T>::uint,
           typename Reconstructor_RT<T>::uint,
           typename Reconstructor_RT<T>::uint>
Reconstructor_RT<T>::getWindow(uint i, uint j, uint l) const {
  uint l0, l1, u0, u1;
  if (l == kpv[i].size()) {
    l0 = 0;
    l1 = 0;
  }
  else {

  }
  return std::tuple<uint, uint, uint, uint>(l0, l1, u1, u1);
}

//== Constructor ===============================================================
// template<OutputType T>
// Reconstructor_RT<T>::Iterator::Iterator(const Tree &t) {
//
// }
//
// template<OutputType T>
// typename Reconstructor_RT<T>::Iterator &Reconstructor_RT<T>::Iterator::operator++() {
//   return *this;
// }
//
// //=== operator!= ===============================================================
// template<OutputType T>
// bool Reconstructor_RT<T>::Iterator::operator!=(const Iterator &other) const {
//   return false;
// }
//
// //=== operator== ===============================================================
// template<OutputType T>
// bool Reconstructor_RT<T>::Iterator::operator==(const Iterator &other) const {
//   return !(*this == other); ;
// }
//
// template<OutputType T>
// typename Reconstructor_RT<T>::Iterator::value_type Reconstructor_RT<T>::Iterator::operator*() const {
//  return {};
// }
//
// template<OutputType T>
// std::string Reconstructor_RT<T>::Iterator::DebugString() const {
//  return {};
// }
//
// template<OutputType T>
// void Reconstructor_RT<T>::Iterator::advance() {
//
// }



//// Reconstructor_RT Stuff ////////////////////////////////////////////////////
//== Constructor ===============================================================
template<OutputType T>
Reconstructor_RT<T>::Reconstructor_RT(const AlgorithmPtr &algo)
    : svv(algo->getStringViewVec()),
      kpv(algo->getKeyPointVec()),
      tree(),
      algo(algo), startPoints(queryStartPoints()) {
  const uint llcs = kpv.size();
  tree.reserve(llcs);
  for (uint i = 0; i < llcs; i++) {
    std::vector<Pair> keyPoints = kpv[i];
    tree.emplace_back(keyPoints);
  }
}

//== cbegin ====================================================================
template<OutputType T>
typename Reconstructor_RT<T>::Iterator Reconstructor_RT<T>::cbegin() const {
  return Iterator(startPoints.begin(), startPoints.end(), tree);
}

//== cend ====================================================================
template<OutputType T>
typename Reconstructor_RT<T>::Iterator Reconstructor_RT<T>::cend() const {
  return Iterator(startPoints.end(), startPoints.end(), tree);
}

//== empty =====================================================================
template<OutputType T>
bool Reconstructor_RT<T>::empty() const {
  return tree.empty();
}

//== DebugString ===============================================================
template<OutputType T>
std::string Reconstructor_RT<T>::DebugString() const {
  std::ostringstream oss;
  uint i = 0;
  for (const auto &s : svv) {
    oss << "svv["  << i++ << "] = " << s << "\n";

  }
  i = 0;
  for (const auto &indices : kpv) {
    oss << "kpv["  << i++ << "] = ";
    if (indices.empty()) {
      oss << "empty\n";
    }
    else {
      for (const auto & val : indices) {
        oss << "[" << val.first << ", " << val.second << "]" << " ";
      }
      oss << "\n";
    }
  }
  return oss.str();
}

#endif//LCS_SOLVER_ALGORITHMS_LCS_RECONSTRUCTOR_RT_TPP_
