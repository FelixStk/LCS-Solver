/*******************************************************************************
 * @file Points.cc
 * @author Steinkopp:Felix
 * @version 1.1
 * @brief Solution consisting of a sequence that contains indices of a dp matrix
 ******************************************************************************/
#include "algorithms/solutions/Points.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <ranges>
#include <span>
#include <sstream>
#include <utility>

#include "util/Logger.hpp"

namespace lcs_solver::algorithms::solutions {

/**
 * @brief Default constructor
 */
Points::Points() = default;

/**
 * @brief Constructor with StrPtrVec
 * @param spv Vector of const String Share Pointers
 * @param reversed Whether the lcs points are pushed in reverse order
 */
Points::Points(StrPtrVec spv, const bool reversed) : Points(std::move(spv), std::span<Embedding>(), reversed) {}

/**
 * @brief Constructor with StrPtrVec and Embeddings Container
 * @param spv Vector of const String Share Pointers
 * @param embs Vector of Embedding
 * @param reversed Whether the embeddings shall be added back to front
 */
Points::Points(StrPtrVec spv, const std::span<Embedding> &embs, const bool reversed)
    : spv_(std::move(spv)), reversed_(reversed) {
  if (!reversed) {
    for (const auto &emb : embs) {
      mat_.push_back(emb.getPositions());
    }
  } else {
    for (const auto &emb : std::ranges::reverse_view(embs)) {
      mat_.push_back(emb.getPositions());
    }
  }
}

/**
 * @brief Constructor with StrPtrVec and Matrix
 * @param spv Vector of const String Share Pointers
 * @param matrix 2DMatrix. M[i][s] = pos(the i-th lcs symbol) in the s-th string
 * @param one_based If true, subtract one from all embedding in mat_.
 * @param reversed Whether the lcs points are pushed in reverse order (sets reversed)
 */
Points::Points(StrPtrVec spv, Matrix2D matrix, const bool one_based, const bool reversed)
    : spv_(std::move(spv)), mat_(std::move(matrix)), reversed_(reversed) {
  if (one_based) {
    for (auto &point : mat_) {
      for (auto &val : point) {
        --val;
      }
    }
  }
}

/**
 * @brief Getter of the SolutionType (defined in BaseSolution)
 * @return BaseSolution::SolutionType
 */
BaseSolution::SolutionType Points::getType() const {
  return SolutionType::Points;
}

/**
 * @brief Checks whether this Points object is equal to another BaseSolution.
 * @details If the instant and the rhs have different reverse_ flags, one
 *          solution is accessed in reverse to compare Points independent
 *          of the internal memory layout.
 * @param rhs Reference to the right-hand-side BaseSolution to compare against
 * @return true iff the objects are of the same type and mark the same LCS
 */
bool Points::isEqual(const BaseSolution &rhs) const {
  if (getType() != rhs.getType())
    return false;

  const auto rhs_ptr = dynamic_cast<const Points *>(&rhs);
  if (rhs_ptr == nullptr)
    return false;

  if (mat_.size() != rhs_ptr->mat_.size())
    return false;

  const auto& lhs_data = this->mat_;
  const auto& rhs_data = rhs_ptr->mat_;
  return this->reversed_ == rhs_ptr->reversed_
    ? std::ranges::equal(lhs_data, rhs_data)
    : std::ranges::equal(lhs_data, std::views::reverse(rhs_data));
}

/**
 * @brief Checks whether this Points object is lexicographically less than
 *        another BaseSolution.
 * @details If the types differ, ordering is determined by comparing their type
 *          enums. If both are of type Points, the internal point sequences are
 *          compared using lexicographical order (via access() and smaller()).
 *          The use of access() abstracts the reverse order of the points if
 *          reverse_ is true.
 * @param rhs Reference to the right-hand-side BaseSolution to compare against
 * @return true if and only if either:
 *         - lhs.getType() < rhs.getType(), or
 *         - lhs.getType() == rhs.getType() and lhs.size() < lhs.size(), or
 *         - lhs.getType() == rhs.getType() and lhs.size() == lhs.size() and
 *           the lhs point vector is lexicographically smaller than that on the
 *           rhs
 */
bool Points::isLessThan(const BaseSolution &rhs) const {
  if (getType() != rhs.getType())
    return getType() < rhs.getType();

  const auto rhs_ptr = dynamic_cast<const Points *>(&rhs);
  if (!rhs_ptr)
    return false;

  if (this->size() < rhs_ptr->size())
    return true;
  if (this->size() > rhs_ptr->size())
    return false;

  for (uint i = 0; i < this->size(); ++i) {
    const auto& a = this->access(i);
    const auto& b = rhs_ptr->access(i);
    if (smaller(a, b))
      return true;
    if (smaller(b, a))
      return false;
  }

  return false; // equal
}

/**
 * @brief Checks whether this Points object is lexicographically less than or
 *        equal to another BaseSolution.
 * @details If the types differ, ordering is determined by comparing their type
 *          enums. If both are of type Points, the internal point sequences are
 *          compared using lexicographical order (via access() and smaller()).
 *          The use of access() abstracts the reverse order of the points if
 *          reverse_ is true.
 * @param rhs Reference to the right-hand-side BaseSolution to compare against
 * @return true if and only if either:
 *         - lhs.getType() < rhs.getType(), or
 *         - lhs.getType() == rhs.getType() and lhs.size() < lhs.size(), or
 *         - lhs.getType() == rhs.getType() and lhs.size() == lhs.size() and
 *           the lhs point vector is lexicographically smaller than or equal to
 *           that on the rhs
 */
bool Points::isLessEqualThan(const BaseSolution &rhs) const {
  if (getType() == rhs.getType())
    return getType() < rhs.getType();

  const auto rhs_ptr = dynamic_cast<const Points *>(&rhs);
  if (!rhs_ptr)
    return false;

  if (this->size() < rhs_ptr->size())
    return true;
  if (this->size() > rhs_ptr->size())
    return false;

  for (uint i = 0; i < this->size(); ++i) {
    const auto& a = this->access(i);
    const auto& b = rhs_ptr->access(i);

    if (smaller(a, b))
      return true;
    if (smaller(b, a))
      return false;
  }
  return true; // equal
}

/**
 * @brief Generates a human-readable debug representation of the Points object
 * @details Creates a formatted string showing the point count and coordinates
 *          of all points. Points are displayed in their current order (normal
 *          or reversed based on the reversed_ flag).
 * @return A debug string in the format:
 *         - "Empty Points Object" if no points exist
 *         - "Length: N, Points: | x1 y1 z1 | x2 y2 z2 | ..." for populated objects
 * @note The output respects the reversed_ state. So points will be shown in
 *       reverse order if reversed_ is true
 */
std::string Points::DebugString() const {
  if (mat_.empty()) {
    return "Empty Points Object";
  }
  std::ostringstream oss;
  oss << "Length: " << mat_.size() << ", Points: | ";
  if (reversed_) {
    for (const auto &point : mat_ | std::views::reverse) {
      for (const auto &val : point) {
        oss << val << " ";
      }
      oss << "| ";
    }

  }
  else {
    for (const auto &point : mat_) {
      for (const auto &val : point) {
        oss << val << " ";
      }
      oss << "| ";
    }
  }
  // oss << "Strings: \n";
  // for (size_t i = 0; i < M.size(); i++) {
  //   oss << "s[" << i << "]" << *spv[i] << " \n";
  // }
  return oss.str();
}

/**
 * @brief Implementation of Clone from the BaseSolution interface
 * @return BaseSolution Pointer to a deep copy of the object on the heap
 */
BaseSolution *Points::clone() const {
  return new Points(spv_, mat_, false,  reversed_);
}

/**
 * @brief Tests whether the Points instance is empty
 * @return true iff no points are stored
 */
bool Points::empty() const {
  return mat_.empty();
}

/**
 * @brief Getter for nth Point that is stored in the object
 * @details If reversed_==true, the point container will be accessed in reverse
 * @param pos Unsigned integer n
 * @return Reference to the nth Point
 */
const Points::Point& Points::access(const uint pos) const {
  assert(pos < mat_.size());
  if (reversed_) {
    auto iter = mat_.rbegin();
    std::advance(iter, pos);
    return *iter;
  }
  return mat_[pos];
}

/**
 * @brief Alias for size()
 * @return The number of points stored in the Points object
 */
Points::uint Points::GetLLCS() const {
  return mat_.size();
}

/**
 * @brief Getter for the number of strings that are stored in the instance
 * @return Size of the String Pointer Container: spv_.size()
 */
Points::uint Points::GetNumOfStrings() const {
  return spv_.size();
}

/**
 * @brief Getter for the string pointers that are stored in the instance
 * @return Const reference to the member spv_
 */
const Points::StrPtrVec &Points::GetStrPtrVec() const {
  return spv_;
}

/**
 * @brief Generates an Embedding as defined in the paper by Adamson et al.
 * @param str_idx Number of Problems String with respect to a StrPtrVec for
 *        which the Embedding should be generated.
 * @param one_based If true, subtract one from all embedding positions.
 * @return Embedding object
 */
Embedding Points::GenEmbedding(const uint str_idx, const bool one_based) const {
  if (str_idx >= GetNumOfStrings())
    return {};

  std::vector<uint> temp;
  temp.reserve(mat_.size());
  if (reversed_) {
    for (const auto &p : mat_ | std::views::reverse)
      temp.push_back(p.at(str_idx));
  }else {
    for (const auto &p : mat_)
      temp.push_back(p.at(str_idx));
  }
  if (one_based)
    for (auto& val : temp)
      --val;
  return {spv_[str_idx], temp};
}

/**
 * @brief Access a Point (aka std::vector<uint>) stored in the instance.
 * @param i Index to access the i-th point. Zero-Based Index.
 * @param reverse If true, the index of point is counted from right to left.
 * @return The i-th Point (aka std::vector<Points::uint>) stored
 */
const std::vector<Points::uint> &Points::operator[](const uint i, const bool reverse) const {
  if (reverse) {
    const uint pos = mat_.size() - 1 - i;
    assert(pos < mat_.size());
    return mat_[pos];
  }
  return mat_[i];
}

/**
 * @brief Access a Point (aka std::vector<uint>) stored in the instance.
 * @param i Index to access the i-th point. Zero-Based Index.
 * @return The i-th Point (aka std::vector<Points::uint>) stored
 * @note The access operator respects the internal reversed_ flag.
 */
const std::vector<Points::uint> &Points::operator[](const uint i) const {
  return operator[](i, reversed_);
}

/**
 * @brief Construct and insert a Point at the end
 * @tparam Args template args
 * @param args Arguments forwarded to construct the new element.
 */
template<typename... Args>
void Points::emplace_back(Args &&...args) {
  mat_.emplace_back(std::vector<uint>{std::forward<Args>(args)...});
}

/**
 * @brief Adds a Point at the end via a reference
 * @param point Reference to a vector of unsigned numbers that represent a point
 */
void Points::push_back(const std::vector<uint> &point) {
  mat_.push_back(point);
}

/**
 * @brief @brief Adds a Point at the end via an r-value
 * @param point R-Value, a vector of unsigned numbers that represent a point
 */
void Points::push_back(const std::vector<uint> &&point) {
  mat_.push_back(point);
}

/**
 * @brief Adds a Point at the end via a reference
 * @param point Reference to a pair of unsigned numbers that represent a point
 */
void Points::push_back(const std::pair<uint, uint> &point) {
  mat_.push_back({point.first, point.second});
}

/**
 * @brief @brief Adds a Point at the end via an r-value
 * @param point R-Value, a pair of unsigned numbers that represent a point
 */
void Points::push_back(const std::pair<uint, uint> &&point) {
  mat_.push_back({point.first, point.second});
}

/**
 * @brief Function to change points stored in the object
 * @param pos Number of the point to be changed (does not consider reverse_)
 * @param point Reference to the new value of the pos-th point
 * @param one_based If true, subtract one from pos.
 */
void Points::mod(const uint pos, const std::pair<uint, uint> &point, const bool one_based) {
  if (!one_based) {
    mat_[pos] = std::vector{point.first, point.second};
  }
  mat_[pos - 1] = std::vector{point.first, point.second};
}

/**
 * Access to the n-th element (point) stored in the object
 * @param pos Number n where to access the internal points container mat_
 * @param one_based If true, subtract one from pos.
 * @return
 */
std::vector<Points::uint>& Points::at(uint pos, const bool one_based) {
  if (one_based) {
    assert(pos > 0);
    pos = pos - 1;
  }
  // if (reversed) {
  //   pos = M.size() - 1 - pos;
  // }
  return mat_[pos];
}

/**
 * @brief Read and Write Access to the data of a Point collection
 * @return Reference to the 2D Matrix. Its first dim gives access to the
 *         different points, and the second dim is for the coordinates of a
 *         point.
 */
Points::Matrix2D& Points::data() {
  return mat_;
}

/**
 * @brief Changes the size of the point stored in the instance
 * @param size The new number of points in the container
 */
void Points::resize(const uint size) {
  mat_.resize(size);
}

/**
 * @brief Clear the points stored in the object
 * @note Does not clear the Strings owned by the object
 */
void Points::clear() {
  mat_.clear();
}

/**
 * @brief Transforms the data matrix such that it equivalent to an embedding
 *        that is read from right to left.
 * @details Be careful when using this function: In general, spv is in this
 *          library the unmodified string input (e.g., no sorting). So if you
 *          call it on an algorithm that outputs Points that based on the sorted
 *          StringView, you will get logical errors. In LCS2_RT, this is solved
 *          by swapping the coordinates if the ordering of the StringView is not
 *          equal to that of the StringPtr Vector (see `aligned_` flag there).
 * @see LCS2_RT::reversePre() and LCS2_RT::reversePost().
 */
void Points::reverse() {
  // const uint llcs = size();

  for(Point& point : mat_){
    for(uint i = 0; i < point.size(); ++i) {
      const uint val = point[i];
      const uint max = spv_[i]->size();
      point[i] = max - val + 1;
    }
  }
  std::ranges::reverse(mat_);
}

/**
 * @brief Clears the Point container and sets marks one single embedding
 * @param embedding 2d points position to be marked
 */
void Points::set_embedding(const std::vector<std::pair<uint, uint>>& embedding) {
  clear();
  for (const auto& point : embedding) {
    push_back(point);
  }
}

/**
 * @brief Getter for the number of Points stored in the object
 * @return First dim of the storage container (mat_.size)
 */
size_t Points::size() const {
  return mat_.size();
}

/**
 * @brief Getter for the last marked poisoning stored at the end of the container
 * @return Last Point (aka std::vector<uint>) at the end of the container
 *         (given by mat_.back())
 */
const Points::Point &Points::back() const {
  return mat_.back();
}

/**
 * @brief Getter for the last marked position stored at the end of the container
 * @return A new pair of uint equal to the end first two coordinates in the last
 *         point at the end of the container
 */
std::pair<Points::uint, Points::uint> Points::BackPair() const {
  if (mat_.empty()) {
    return std::make_pair(0, 0);
  }
  return {mat_.back().at(0), mat_.back().at(1)};
}

/**
 * @brief Delete last element pop_back
 * @details Removes the last point added to the object
 */
void Points::pop_back() {
  mat_.pop_back();
}

/**
 * @brief Function to change points stored in the object
 * @param lcs_idx uint number of the point
 * @param str_idx uint Index of a string in spv_ of the object. (zero-based)
 * @param pos uint of character to be marked in the strIdx-th string.
 * @param reverse bool if true positions in the string are counted from right to left
 */
void Points::modify(const uint lcs_idx, const uint str_idx, const uint pos, const bool reverse) {
  assert(str_idx < spv_.size());
  // if (strIdx >= spv.size()) {
  //   util::Logger::Error() << "Points cannot be modified (spv size error)"
  //                         << "lcsIdx: " << lcsIdx
  //                         << "strIdx: " << strIdx << " pos: " << pos
  //                         << "reverse: " << reverse;
  //   return;
  // }
  if (reverse == false) {
    mat_[lcs_idx][str_idx] = pos;
  }
  const uint idx = mat_[lcs_idx].size() - 1 - str_idx;
  assert(idx < mat_.size());
  // if (idx >= M.size()) {
  //   util::Logger::Error() << "Points cannot be modified (idx error)"
  //                         << "lcsIdx: " << lcsIdx
  //                         << "strIdx: " << strIdx << " pos: " << pos
  //                         << "reverse: " << reverse
  //                         << "idx: " << idx;
  //   return;
  // }
  mat_[lcs_idx][idx] = pos;
}

/**
 * @brief Lexicographical less_equal operation for std::vector<uint> objects
 * @param a std::vector<uint> lhs
 * @param b std::vector<uint> rhs
 * @return lhs < rhs (lexicographic)
 */
bool Points::SmallerOrEqual(const std::vector<uint> &a, const std::vector<uint> &b) {
  return std::ranges::lexicographical_compare(a, b, std::less_equal());
}

/**
 * Lexicographical less operation for std::vector<uint> objects
 * @param a std::vector<uint> lhs
 * @param b std::vector<uint> rhs
 * @return lhs < rhs (lexicographic)
 */
bool Points::smaller(const std::vector<uint> &a, const std::vector<uint> &b) {
  return std::ranges::lexicographical_compare(a, b);
}

/**
 * @brief Stream insertion operator for Points objects
 * @details Enables Points objects to be directly output to streams (like std::cout,
 *          file streams, etc.) by delegating to the DebugString() method.
 * @param os The output stream to write to
 * @param points The Points object to output
 * @return Reference to the output stream (for chaining operations)
 * @note This operator uses the same format as DebugString():
 *       - "Empty Points Object" for empty objects
 *       - "Length: N, Points: | x1 y1 z1 | x2 y2 z2 | ..." for populated objects
 */
std::ostream &operator<<(std::ostream &os, const Points &points) {
  os << points.DebugString();
  return os;
}
}// namespace lcs_solver::algorithms::solutions