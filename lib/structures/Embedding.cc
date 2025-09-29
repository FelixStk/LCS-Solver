/*********************************************************************
 * @file Embedding.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Impl. of Embedding class (represents a possible lcs solution)
 ********************************************************************/

#include "structures/Embedding.h"

#include <cassert>
#include <compare>
#include <cstddef>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

//#define DEBUG_MODE_EMBEDDING // Enables print debugging for constructors
//#define DEBUG_MODE_INDEX // Enables printing in modPos(pos,value)

#if defined(DEBUG_MODE_EMBEDDING) || defined(DEBUG_MODE_INDEX)
#include "util/Logger.hpp"
#endif

namespace {
using String = ::lcs_solver::structures::Embedding::String;
using StringView = ::lcs_solver::structures::Embedding::StringView;
using Symbol = ::lcs_solver::structures::Embedding::Symbol;
}// namespace

namespace lcs_solver::structures {

/*******************************************************************************
 * @brief Default Constructor
 ******************************************************************************/
Embedding::Embedding()
    : strPtr(nullptr) {
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding constructed: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Constructor (with StringPtr & without marking of Symbols)
 * @param sp StringPtr points to the string that the embedding is based on
 ******************************************************************************/
Embedding::Embedding(const std::shared_ptr<const String> &sp) : strPtr(sp) {
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding constructed: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Constructor (with StringPtr & length of the embedding)
 * @param sp StringPtr points to the string that the embedding is based on
 * @param length uint number of symbols in marked in the string `*sp`
 ******************************************************************************/
Embedding::Embedding(const std::shared_ptr<const String> &sp, uint length)
    : strPtr(sp), position(std::vector<uint>(length)) {
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding constructed: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Constructor (with StringPtr & marking of symbols)
 * @param sp StringPtr points to the string that the embedding is based on
 * @param vec vector of positions to be marked in teh string `*os`
 * @param oneBased
 ******************************************************************************/
Embedding::Embedding(
    const std::shared_ptr<const String> &sp,
    const std::vector<uint> &vec,
    bool oneBased) : strPtr(sp), position(vec) {
  if (oneBased) {
    for (auto &pos : position) {
      --pos;
    }
  }
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding constructed: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Copy Constructor
 * @param other Embedding to be copied into the object
 ******************************************************************************/
Embedding::Embedding(const Embedding &other)
    : strPtr(other.strPtr), position(other.position) {
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding copied: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Move Constructor
 * @param other Embedding to be moved into the object
 ******************************************************************************/
Embedding::Embedding(Embedding &&other) noexcept
    : strPtr(std::move(other.strPtr)), position(std::move(other.position)) {
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding moved: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Destructor
 ******************************************************************************/
Embedding::~Embedding() {
#ifdef DEBUG_MODE_EMBEDDING
  Logger::Info() << "Embedding destroyed: " << this << std::endl;
#endif
}

/*******************************************************************************
 * @brief Copy Assignment operator
 * @param other Embedding to be used to assign this object to
 * @return Embedding reference of this object
 ******************************************************************************/
Embedding &Embedding::operator=(const Embedding &other) {
  if (this != &other) {// self-assignment check
    strPtr = other.strPtr;
    position = other.position;
#ifdef DEBUG_MODE_EMBEDDING
    Logger::Info() << "Embedding copy-assigned: " << this << std::endl;
#endif
  }
  return *this;
}

/*******************************************************************************
 * @brief Move Assignment Operator
 * @param other Embedding to be moved
 * @return Embedding reference of this object
 ******************************************************************************/
Embedding &Embedding::operator=(Embedding &&other) noexcept {
  if (this != &other) {// Protect against self-assignment
    strPtr = std::move(other.strPtr);
    position = std::move(other.position);
#ifdef DEBUG_MODE_EMBEDDING
    Logger::Info() << "Embedding move-assigned: " << this << std::endl;
#endif
  }
  return *this;
}

/*******************************************************************************
 * @brief Bracket Operator
 * @param index uint that represent the index of the marked symbols
 * @return Position of the `i`-th marked symbol in the string `*sp`
 ******************************************************************************/
const Embedding::uint &Embedding::operator[](const uint index) const {
  assert(index < position.size() && "Embedding: Index out of bounds (too large)!");
  return position[index];
}

/*******************************************************************************
 * @brief Spaceship Operator
 * @param other Embedding& to represent the rhs of the spaceship operator
 * @return std::strong_ordering based on comparing the strings first and then as
 * a second criteria the marked positions, if and only if the string are equal.
 ******************************************************************************/
std::strong_ordering Embedding::operator<=>(const Embedding &other) const {
  // First, compare the strings lexicographically
  if (const auto cmp = *strPtr <=> *(other.strPtr); cmp != nullptr) {
    return cmp;
  }

  // Second (only when the strings are equal), compare the position vector (again lexicographically)
  return position <=> other.position;
}

/*******************************************************************************
 * @brief Not Operator
 * @return true if the Embedding object is empty
 ******************************************************************************/
bool Embedding::operator!() const {
  return empty();
}

/*******************************************************************************
 * @brief Modifies the position to change a symbol marking
 * @details Executes position[indexInSubsequence] = indexInString
 * @param indexInSubsequence uint the index in the subsequence
 * @param indexInString uint the index in the string
 ******************************************************************************/
void Embedding::modPosition(uint indexInSubsequence, uint indexInString) {
  assert(indexInSubsequence < position.size() && "Embedding: Access out of bounds (not enough space in position )!");
  assert(indexInSubsequence < strPtr->size() && "Embedding: Access out of bounds (indexInSubsequence > str length)!");
  assert(indexInString < strPtr->size() && "Embedding: Access out of bounds (indexInString too large)!");
#ifdef DEBUG_MODE_INDEX
  if (indexInSubsequence >= position.size())
    util::Logger::Info() << "Embedding out of range! (tried to set e["
                         << indexInSubsequence << "] to "
                         << indexInString << std::endl;
#endif
  position[indexInSubsequence] = indexInString;
}

/*******************************************************************************
 * @brief Changes the number of marked symbols (resizes position vector)
 * @param l uint new length of the embedding
 ******************************************************************************/
void Embedding::resize(const uint l) {
  position.resize(l);
}

/*******************************************************************************
 * @brief Generates the embedded substring
 * @return the subsequence of marked symbols
 ******************************************************************************/
String Embedding::getEmbeddedStr() const {
  String res;
  for (const auto i : position) {
    res += strPtr->at(i);
  }
  return res;
}

/*******************************************************************************
 * @brief Checks if the embedded subsequence has length 0 (is empty)
 * @return bool position.empty()
 ******************************************************************************/
bool Embedding::empty() const {
  return position.empty();
}

/*******************************************************************************
 * @brief Getter for the length of the subsequence
 * @return size_t position.size()
 ******************************************************************************/
size_t Embedding::size() const {
  return position.size();// length of subsequence
}

/*******************************************************************************
 * @brief Gets a Symbol of the Subsequence
 * @param indexInSubsequence uint index of the Symbol in teh Subsequence
 * @return strPtr->at(position[indexInSubsequence])
 ******************************************************************************/
const Symbol &Embedding::symbolAt(size_t indexInSubsequence) const {
  assert(indexInSubsequence < position.size()
         && "Embedding: Access out of bounds (not enough space in position )!");
  assert(indexInSubsequence < strPtr->size()
         && "Embedding: Access out of bounds (indexInSubsequence > str length)!");
  return strPtr->at(position[indexInSubsequence]);
}

/*******************************************************************************
 * Gets the positions vector of the embedded symbols
 * @return std::vector<uint> represents the position
 ******************************************************************************/
const std::vector<Embedding::uint> &Embedding::getPositions() const {
  return position;
}

/*******************************************************************************
 * Checks if the Embedding object is valid
 * @details check whether the indices in `position` are strictly increasing and
 * its length is plausible
 * @return true if the embedding object is valid
 ******************************************************************************/
bool Embedding::isValidEmbedding() const {
  if (position.empty())
    return true;

  auto frontIndex = position.front();
  size_t i = 1;
  while (i < position.size()) {
    if (frontIndex <= position[i])
      return false;
    frontIndex = position[i];
    i++;
  }
  return position.back() < strPtr->size();
}

/*******************************************************************************
 * Gets a symbol of the reference string
 * @param idx std::size_t the position where to look in the string (0 based)
 * @return idx-th symbol of the reference string stored (idx in [0:len(str)-1]
 ******************************************************************************/
Symbol Embedding::getSymbolAt(std::size_t idx) const {
  return strPtr->at(idx);
}

/*******************************************************************************
 * Gets a StringView that contains the idx-th gap
 * @note the embedded symbols surrounding the gap are the contained borders of
 * the StringView returned
 * @param idx uint number of the gap (zero based)
 * @return StringView of the (*strPtr)[position[i]:position[i+1]]
 ******************************************************************************/
StringView Embedding::getGapStrView(std::size_t idx) const {
  return {strPtr->data() + position[idx], strPtr->data() + position[idx + 1]};
}

/*******************************************************************************
 * @brief Gets the start and end of a gap
 * @param idx std::size_t for indexing the idx-th gap (zero based)
 * @return std::pair<uint, uint> of the borders that surround the idx-th gap
 ******************************************************************************/
std::pair<Embedding::uint, Embedding::uint> Embedding::getGapStartEnd(std::size_t idx) const {
  return {position[idx], position[idx + 1]};
}

/*******************************************************************************
 * @brief Gets the StringView of the String on which the Subsequence is based on
 * @return StringView of the reference string pointer (given to the constructor)
 ******************************************************************************/
StringView Embedding::getStr() const {
  return *strPtr;
}

/*******************************************************************************
 * Gets the position of the i-th embedded symbol
 * @param idx std::size_t index in the subsequence
 * @return uint containing position[idx]
 ******************************************************************************/
Embedding::uint Embedding::getPosition(const std::size_t idx) const {
  return position[idx];
}

/*******************************************************************************
 * Creates a String describing an Embedding object
 * @return std::string with the following parts:
 * - length of the subsequence
 * - the std::string representation of the subsequence
 * - the std::string representation of the reference string
 * - a description of marked symbols (string of `positions` vector)
 ******************************************************************************/
std::string Embedding::DebugString() const {
  std::ostringstream oss;
  oss << "Embedding of length " << position.size();
  oss << "(" << util::to_string(getEmbeddedStr()) << " in " << util::to_string(getStr()) << ")" << std::endl;
  for (size_t i = 0; i < position.size(); ++i) {
    oss << "\t p[" << i << "]="
        << position[i] << "(c=" << util::to_string(getSymbolAt(i))
        << ", hex = " << util::to_string(getSymbolAt(i), true);
    oss << std::endl;
  }
  return oss.str();
}

}// namespace lcs_solver::structures
