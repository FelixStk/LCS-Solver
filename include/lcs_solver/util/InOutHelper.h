#ifndef LCS_SOLVER_INCLUDE_UTIL_INOUTHELPER_H_
#define LCS_SOLVER_INCLUDE_UTIL_INOUTHELPER_H_

#include "CommonTypes.h"

#include <ostream>
#include <istream>
#include <string>
#include <iostream>

namespace lcs_solver::util{

char ReadChar(std::string_view msg, char default_value, std::istream &in=std::cin, std::ostream &os=std::cout);
std::string ReadStdString(std::string_view msg = "", std::string default_value = "", std::istream &in=std::cin, std::ostream &os=std::cout);
String ReadUtilString(std::string_view msg = "", std::istream &in=std::cin, std::ostream &os=std::cout);
String ReadUtilString(uint n, std::istream &in=std::cin, std::ostream &os=std::cout);
String ReadUtilString(uint n, StringView default_view, std::istream &in=std::cin, std::ostream &os=std::cout);
uint ReadUnsigned(std::string_view msg = "", std::istream &in=std::cin, std::ostream &os=std::cout);
uint ReadUnsigned(std::string_view msg, uint default_value,std::istream &in=std::cin, std::ostream &os=std::cout);
uint ReadUnsigned(std::string_view var_name, uint default_value, const std::pair<uint,uint> &range, std::istream &in=std::cin, std::ostream &os=std::cout);
std::pair<uint,uint> ReadSubinterval(std::string_view var_name, std::pair<uint,uint> default_value, const std::pair<uint, uint> &range, std::istream &in=std::cin, std::ostream &os=std::cout);

bool ReadBool(std::string_view msg, bool default_value, std::istream &in=std::cin, std::ostream &os=std::cout);
bool IsNumber(std::string_view sv);
std::string to_string(const std::pair<uint, uint> &pair);
std::string & Trim(std::string & input);
}
#endif //LCS_SOLVER_INCLUDE_UTIL_INOUTHELPER_H_