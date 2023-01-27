#pragma once

#include <string>
#include <sstream>
#include <fstream>

namespace utils {
	template<typename _Elem> using string_type = std::basic_string<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>>;
	template<typename _Elem> using istream_type = std::basic_istream<_Elem, std::char_traits<_Elem>>;
	template<typename _Elem> using ostream_type = std::basic_ostream<_Elem, std::char_traits<_Elem>>;
	template<typename _Elem> using ifstream_type = std::basic_ifstream<_Elem, std::char_traits<_Elem>>;
	template<typename _Elem> using ofstream_type = std::basic_ofstream<_Elem, std::char_traits<_Elem>>;
	template<typename _Elem> using stringstream_type = std::basic_stringstream<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>>;
	template<typename _Elem> using istringstream_type = std::basic_istringstream<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>>;
	template<typename _Elem> using ostringstream_type = std::basic_ostringstream<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>>;
}