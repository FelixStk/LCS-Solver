/******************************************************************************
 * @file Logger.cc
 * @author Steinkopp:Felix
 * @version 1.0
 * @brief Logger is a singleton for debugging with levels and changeable stream
 ******************************************************************************/

#include "util/Logger.hpp"
#include <iostream>
#include <string>
namespace lcs_solver::util {
Logger::Logger()
    : m_out(&std::cout), m_err(&std::cerr), m_owner(false) {
  setStream(&std::cout, false);
  setStream(&std::cerr, false);
}

Logger::~Logger() {
  setStream(nullptr, false);
}

Logger &Logger::Get() {
  static Logger instance;
  return instance;
}

Logger &Logger::SetLogLevel(LogLevel level) {
  return Get().InternalSet(level);
}

Logger &Logger::InternalSet(LogLevel level) {
  Get().m_log_level = level;
  return Get();
}

void Logger::setStream(std::ostream *stream, bool owner) {
  if (m_owner)
    delete m_out;
  m_out = stream;
  m_owner = owner;
}

std::ostream &Logger::Error() {
  return Get().InternalError();
}

std::ostream &Logger::Warning() {
  return Get().InternalWarning();
}

std::ostream &Logger::Info() {
  return Get().InternalInfo();
}

std::ostream &Logger::Detail() {
  return Get().InternalDetail();
}

std::ostream &Logger::InternalError() {
  if (!m_out)
    throw std::runtime_error("No stream set for Logger class");

  if (m_log_level >= LogLevel::ERROR) {
    (*m_err) << "[Error]: ";
    return *m_out;
  }
  return m_dummy_stream;
}

std::ostream &Logger::InternalWarning() {
  if (!m_out)
    throw std::runtime_error("No stream set for Logger class");

  if (m_log_level >= LogLevel::WARNING) {
    (*m_out) << "[Warning]: ";
    return *m_out;
  }
  return m_dummy_stream;
}

std::ostream &Logger::InternalInfo() {
  if (!m_out)
    throw std::runtime_error("No stream set for Logger class");

  if (m_log_level >= LogLevel::INFO) {
    (*m_out) << "[Info]: ";
    return *m_out;
  }
  return m_dummy_stream;
}

std::ostream &Logger::InternalDetail() {
  if (!m_out)
    throw std::runtime_error("No stream set for Logger class");

  if (m_log_level >= LogLevel::DETAIL) {
    (*m_out) << "[Detail]: ";
    return *m_out;
  }
  return m_dummy_stream;
}

std::string Logger::Matrix(const std::vector<std::vector<int>> &matrix) {
  return Get().InternalMatrix(matrix);
}

std::string Logger::Matrix(const std::vector<std::vector<util::uint>> &matrix) {
  return Get().InternalMatrix(matrix);
}

std::string Logger::InternalMatrix(const std::vector<std::vector<int>> &matrix) {
  std::string result;
  if (!matrix.empty()) {
    for (const auto &row : matrix) {
      for (const auto &element : row) {
        result += std::to_string(element) + " ";
      }
      result += "\n";
    }
  } else {
    result = "Matrix is empty\n";
  }
  return result + "\n";
}

std::string Logger::InternalMatrix(const std::vector<std::vector<util::uint>> &matrix) {
  std::string result;
  if (!matrix.empty()) {
    for (const auto &row : matrix) {
      for (const auto &element : row) {
        result += std::to_string(element) + " ";
      }
      result += "\n";
    }
  } else {
    result = "Matrix is empty\n";
  }
  return result + "\n";
}
}