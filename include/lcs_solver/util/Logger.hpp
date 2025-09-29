#ifndef LCS_SOLVER_UTIL_LOGGER_HPP
#define LCS_SOLVER_UTIL_LOGGER_HPP

#include <ostream>
#include <vector>

#include "CommonTypes.h"
namespace lcs_solver::util {
enum class LogLevel {
  ERROR = 0,
  WARNING,
  INFO,
  DETAIL
};

class Logger {
 public:
  Logger(const Logger &) = delete;
  virtual ~Logger();

  static Logger &Get();
  static Logger &SetLogLevel(LogLevel level);
  void setStream(std::ostream *stream, bool owner);
  static std::ostream &Error();
  static std::ostream &Warning();
  static std::ostream &Info();
  static std::ostream &Detail();
  static std::string Matrix(const std::vector<std::vector<int>> &m);
  static std::string Matrix(const std::vector<std::vector<util::uint>> &m);

 private:
  Logger();
  static Logger &InternalSet(LogLevel level);
  std::ostream &InternalError();
  std::ostream &InternalWarning();
  std::ostream &InternalInfo();
  std::ostream &InternalDetail();
  static std::string InternalMatrix(const std::vector<std::vector<int>> &m);
  static std::string InternalMatrix(const std::vector<std::vector<util::uint>> &m);

  std::ostream *m_out;
  std::ostream *m_err;
  bool m_owner;
  LogLevel m_log_level = LogLevel::INFO;

  struct DummyStream : std::ostream {
    template<typename T>
    DummyStream &operator<<(const T &) {
      return *this;
    }
  } m_dummy_stream;
};
}
#endif /* LCS_SOLVER_UTIL_LOGGER_HPP */