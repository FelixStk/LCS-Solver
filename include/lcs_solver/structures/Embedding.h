#ifndef LCS_SOLVER_STRUCTURES_EMBEDDING_H_
#define LCS_SOLVER_STRUCTURES_EMBEDDING_H_

#include <compare>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>
#include <utility>

#include "util/CommonTypes.h"

namespace lcs_solver::structures {
class Embedding {
 public:
  using Symbol = lcs_solver::util::Symbol;
  using String = lcs_solver::util::String;
  using StringView = lcs_solver::util::StringView;
  using StringPtr = std::shared_ptr<const String>;
  using uint = std::size_t;

  Embedding();
  explicit Embedding(const std::shared_ptr<const String> &sp);
  Embedding(const StringPtr &sp, uint length);
  Embedding(const StringPtr &sp, const std::vector<uint> &positions, bool oneBased = false);
  Embedding(const Embedding &rhs);
  Embedding(Embedding &&rhs) noexcept;
  ~Embedding();

  Embedding &operator=(const Embedding &rhs);
  Embedding &operator=(Embedding &&rhs) noexcept;
  const uint &operator[](uint index) const;
  std::strong_ordering operator<=>(const Embedding &other) const;
  bool operator!() const;

  void modPosition(uint indexInSubsequence, uint indexInString);
  void resize(uint l);

  [[nodiscard]] StringView getStr() const;
  [[nodiscard]] String getEmbeddedStr() const;
  [[nodiscard]] Symbol getSymbolAt(uint idx) const;
  [[nodiscard]] StringView getGapStrView(uint idx) const;
  [[nodiscard]] std::pair<uint,uint> getGapStartEnd(uint idx) const;
  [[nodiscard]] uint getPosition(uint idx) const;
  [[nodiscard]] const std::vector<uint> &getPositions() const;
  [[nodiscard]] bool empty() const;
  [[nodiscard]] uint size() const;
  [[nodiscard]] const Symbol &symbolAt(uint indexInSubsequence) const;
  [[nodiscard]] bool isValidEmbedding() const;
  [[nodiscard]] std::string DebugString() const;

 private:
  std::shared_ptr<const String> strPtr;
  std::vector<uint> position; //position[a] == b means the ath symbol of the subsequence is strPtr[b]
};

using Embeddings = std::vector<Embedding>;

} // end of namespace
#endif // LCS_SOLVER_STRUCTURES_EMBEDDING_H_