#ifndef LCS_SOLVER_PROBLEMS_PROBLEMFACTORY_H_
#define LCS_SOLVER_PROBLEMS_PROBLEMFACTORY_H_

#include <array>
#include <filesystem>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <string_view>

#include "algorithms/AlgoType.h"
#include "problems/BaseProblem.h"
#include "problems/ProblemType.h"

namespace lcs_solver::problems {

/**
 * @class ProblemFactory
 * @brief Factory class for creating LCS problem instances and related components
 *
 * Provides functionality to:
 * - Create problems through various methods (interactive, JSON, direct parameters)
 * - Serialize/deserialize problem configurations
 * - Retrieve metadata about available problem types and algorithms
 *
 * @note All factory methods are static. This class cannot be instantiated.
 *
 * @see ProblemType, BaseProblem, AlgoType
 */
class ProblemFactory {
 public:
  /// @name Type aliases
  /// @{
  using uint            = std::size_t;                        ///< Unsigned integer type used for indices and sizes
  using GapVector       = std::vector<std::pair<uint, uint>>; ///< Vector of gap constraint tuples (min, max) pairs
  using GapSpan         = std::span<const std::pair<uint, uint>>; ///< GapSpan for convenient access to a GapVector
  using Path            = std::filesystem::path;              ///< Filesystem path type for I/O operations
  using SymbolVector    = std::vector<util::Symbol>;          ///< Vector of symbols representing an alphabet
  using SymbolSpan      = std::span<const util::Symbol>;      ///< Symbol Span for convenient access to a SymbolVector

  // Forward declarations from BaseProblem
  using AlgoType        = algorithms::AlgoType;               ///< Enum Class Type used for identifying algorithms
  using AlgorithmPtr    = BaseProblem::AlgorithmPtr;          ///< Smart Pointer used to manage BaseProblem instances
  using AlgorithmSpan   = BaseProblem::AlgorithmSpan;         ///< AlgoType Span for giving options to the user
  using ConstraintMap   = BaseProblem::ConstraintMap;         ///< Map from ConstraintType to BaseConstraint SharedPtr
  using ConstraintSpan  = BaseProblem::ConstraintSpan;        ///< ConstraintType Span for giving options to the user
  using ProblemPtr      = BaseProblem::ProblemPtr;            ///< Unique Pointer to manage a problem instance
  using StrPtrVec       = BaseProblem::StrPtrVec;             ///< Vector of shared pointers util::String
  using String          = util::String;                       ///< String type for input sequences
  /// @}

  /// @name Factory Methods
  /// @{

  /**
   * @brief Interactive problem creation through the console dialog
   * @return ProblemPtr Created problem instance ready for solving
   * @throws std::runtime_error If input validation fails or required components
   *         missing
   *
   * Guides user through:
   * 1. Problem type selection
   * 2. Metadata entry (name/description)
   * 3. String input
   * 4. Constraint configuration
   * 5. Algorithm selection
   */
  static ProblemPtr CreateFromDialog();

  /**
   * @brief Create a problem from a JSON file
   * @param p Path to JSON file containing problem specification
   * @return ProblemPtr Parsed problem instance
   * @throws std::runtime_error If file I/O error occurs
   * @throws nlohmann::json::exception If JSON parsing fails
   */
  static ProblemPtr CreateFromJson(const Path& p);

  /**
   * @brief Direct problem creation with full parameters
   * @param t Problem type from ProblemType enum
   * @param name Problem name (moved)
   * @param description Problem description (moved)
   * @param spv String collection (moved)
   * @param map Constraint mapping (moved)
   * @return ProblemPtr Created problem instance
   * @throws std::runtime_error If ProblemType::LCS_Base is requested
   *
   * @warning The ProblemType::LCS_Base cannot be instantiated directly
   */
  static ProblemPtr Create(ProblemType t,
                           std::string&& name,
                           std::string&& description,
                           StrPtrVec&& spv,
                           ConstraintMap&& map);

  /// @}

  /// @name Constraint-based Creation
  /// @{

  /**
   * @brief Create a problem with automatic (relaxed) constraint generation
   * @param t Problem type
   * @param s Span of input strings
   * @param at Main algorithm type
   * @param calc_llcs Whether to compute LLCS or the LCS
   * @return ProblemPtr Created problem instance
   *
   * Automatically generates appropriate constraints based on:
   * - Problem type
   * - Input strings
   * - Algorithm capabilities
   */
  static ProblemPtr Create(ProblemType t,
                           std::span<const String> s,
                           AlgoType at,
                           bool calc_llcs = true);

  /**
   * @brief Create a problem with explicit gap constraints
   * @param t Problem type
   * @param s Span of input strings
   * @param gap_tuples Predefined gap constraints
   * @param at Main algorithm type
   * @param calc_llcs Whether to compute LLCS or the LCS
   * @return ProblemPtr Created problem instance
   */
  static ProblemPtr Create(ProblemType t,
                           std::span<const String> s,
                           GapSpan gap_tuples,
                           AlgoType at,
                           bool calc_llcs = true);

  /**
   * @brief Create a problem with explicit gap constraints
   * @param t Problem type
   * @param s Span of input strings
   * @param alph Alphabet used for the strings
   * @param gc Predefined gap constraints with respect to the ordering in alph
   * @param at Main algorithm type
   * @param calc_llcs Whether to compute LLCS or the LCS
   * @return ProblemPtr Created problem instance
   *
   * @note alph and gc must have the same length
   */
  static ProblemPtr Create(ProblemType t,
                           std::span<const String> s,
                           SymbolSpan alph,
                           GapSpan gc,
                           AlgoType at,
                           bool calc_llcs = true);

  /**
   * @brief Create a problem with explicit gap constraints
   * @param t Problem type
   * @param s Span of input strings
   * @param alph Alphabet used for the strings
   * @param l Predefined gap constraints with respect to the ordering in alph
   * @param r Predefined gap constraints with respect to the ordering in alph
   * @param at Main algorithm type
   * @param calc_llcs Whether to compute LLCS or the LCS
   * @return ProblemPtr Created problem instance
   *
   * @note alph, l and r specify together the maps
   * - left[c]  = (l,u) => [l ≤ length (gap before Symbol c) ≤ u]
   * - right[c] = (l,u) => [l ≤ length (gap after Symbol c) ≤ u]
   */
  static ProblemPtr Create(ProblemType t,
                           std::span<const String> s,
                           SymbolSpan alph,
                           GapSpan l,
                           GapSpan r,
                           AlgoType at,
                           bool calc_llcs = true);
  /**
   * @brief Create a problem with explicit gap constraints
   * @details If the constraint data is not provided to allow the creation of a
   *          problem with type t, then a relaxed constraint is created by
   *          calling `Create(t,s,at,calc_llcs)`
   * @param t Problem type
   * @param s Span of input strings
   * @param gc Predefined mc gap constraints
   * @param alph Alphabet used for the strings
   * @param l Predefined gap constraints with respect to the ordering in alph
   * @param r Predefined gap constraints with respect to the ordering in alph
   * @param at Main algorithm type
   * @param calc_llcs Whether to compute LLCS or the LCS
   * @return ProblemPtr Created problem instance
   *
   * @note alph, l and r specify together the maps
   * - left[c]  = (l,u) => [l ≤ length (gap before Symbol c) ≤ u]
   * - right[c] = (l,u) => [l ≤ length (gap after Symbol c) ≤ u]
   */
  static ProblemPtr Create(ProblemType t,
                           std::span<const String> s,
                           GapSpan gc,
                           SymbolSpan alph,
                           GapSpan l,
                           GapSpan r,
                           AlgoType at,
                           bool calc_llcs = true);
  /// @}

  /// @name Metadata and Utilities
  /// @{

  /**
   * @brief Interactive problem type selection
   * @param in Input stream (default: cin)
   * @param os Output stream (default: cout)
   * @return Selected ProblemType
   * @throws std::runtime_error If invalid input provided
   */
  static ProblemType ReadProblemType(std::istream& in = std::cin,
                                     std::ostream& os = std::cout);

  /**
   * Interactively, prompts for the problem's name and description
   * @param t Problem type
   * @param in Input stream (default: cin)
   * @param os Output stream (default: cout)
   * @return Pair of strings containing the name and description
   */
  static std::pair<std::string, std::string> ReadDescription(
      ProblemType t, std::istream& in = std::cin, std::ostream& os = std::cout);

  /**
   * Interactively, prompts for the problem's input sequences
   * @param in Input stream (default: cin)
   * @param os Output stream (default: cout)
   * @return Vector of util::String share pointers
   * @see util::ReadStrPtrVec (in util/StrPtrVector.h)
   */
  static StrPtrVec ReadStrPtrVec(std::istream& in = std::cin,
                                 std::ostream& os = std::cout);

  /**
   * @brief Get a human-readable problem type name
   * @param type ProblemType enumeration value
   * @return String view of a problem type name
   */
  static std::string_view GetName(ProblemType type);

  /**
   * @brief Get problem type description
   * @param type ProblemType enumeration value
   * @return String view of problem description
   */
  static std::string_view GetDescription(ProblemType type);

  /**
   * @brief Get required constraints for a problem type
   * @param type ProblemType enumeration value
   * @return Span of required constraint types
   */
  static ConstraintSpan GetReqConstrains(ProblemType type);

  /**
   * @brief Get optional constraints for a problem type
   * @param type ProblemType enumeration value
   * @return Span of optional constraint types
   */
  static ConstraintSpan GetOptConstrains(ProblemType type);

  /**
   * @brief Get available algorithms for a problem type
   * @param type ProblemType enumeration value
   * @return Span of available algorithm types
   */
  static AlgorithmSpan GetAlgorithms(ProblemType type);

  /**
   * @brief Maps strings back to ProblemTypes
   * @return Constant Reference to Reverse Map
   */
  static const std::map<std::string_view, ProblemType>& GetNameTypeMap();

  /**
   * @brief Save a BaseProblem instance in a json file
   * @param file Filepath to the json file
   * @param ptr Pointer to the BaseProblem instance
   */
  static void SaveToJson(const Path& file, const BaseProblem* ptr);

  /**
   * @brief Old parser for generating Problems from txt files
   * @param file Path to txt file
   * @return ProblemPtr Created problem instance
   */
  static ProblemPtr FromFile(const Path& file);

  /**
   * @brief Old problem save operation
   * @param os An out stream (e.g., std::cout or a ofstream)
   * @param p Pointer to the BaseProblem instance to save
   */
  static void WriteProblem(std::ostream& os, const BaseProblem* p);

  /**
   * Old Wrapper for FromFile
   * @param file Path to txt file
   * @param p Pointer to the BaseProblem instance to save
   * @throws std::runtime_error if the file could not be opened for writing
   */
  static void ToFile(const Path& file, const BaseProblem* p);

  ~ProblemFactory() = default;

  /// @brief Array of all available problem types
  static constexpr std::array<ProblemType, 9> kAvailable = {
      ProblemType::LCS_Base,
      ProblemType::LCS_Classic,
      ProblemType::LCS_MC,
      ProblemType::LCS_MC_INC,
      ProblemType::LCS_MC_1C,
      ProblemType::LCS_MC_O1C_SYNC,
      ProblemType::LCS_Sigma_R,
      ProblemType::LCS_Sigma_L,
      ProblemType::LCS_Sigma,
  };
  /// @}
};

}  // namespace lcs_solver::problems

#endif  // LCS_SOLVER_PROBLEMS_PROBLEMFACTORY_H_