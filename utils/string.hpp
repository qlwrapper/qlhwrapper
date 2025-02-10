/*
 Copyright (C) 2022 QLHWrapper

 This file is part of QLHWrapper, a free-software/open-source library
 for financial quantitative analysts and developers

 QLHWrapper is free software: you can redistribute it and/or modify it
 under the terms of the The 2-Clause BSD License license - https://opensource.org/licenses/BSD-2-Clause.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <utils/string-io-types.hpp>
#include <utils/required-cpp17.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm> 
#include <cctype>
#include <iomanip>
#include <sstream>
#include <functional>

namespace utils {
	// convert between utf8 (char) string and wstring of type char16_t, char32_t, or wchar_t
	// _Elem = wchar_t, char16_t, char32_t
#ifdef REQUIRED_CPP17
	template <typename _Elem> using utf8_wstring_converter = std::for_cpp17::wstring_convert<std::for_cpp17::codecvt_utf8<_Elem>, _Elem>;
#else
	template <typename _Elem> using utf8_wstring_converter = std::wstring_convert<std::codecvt_utf8<_Elem>, _Elem>;
#endif

	template<typename T_Elem>
	inline std::string from_tstring(const string_type<T_Elem>& t_s) {
		utf8_wstring_converter<T_Elem> converter;
		return converter.to_bytes(t_s);
	}
	template<typename T_Elem>
	inline string_type<T_Elem> to_tstring(const std::string& s) {
		utf8_wstring_converter<T_Elem> converter;
		return converter.from_bytes(s);
	}

	// parse string into vector of string by a delimter character
	template <typename _Elem>
	inline std::shared_ptr<std::vector<string_type<_Elem>>> parse_delimited(const string_type<_Elem>& s, unsigned int delimAsci) {
		auto ch = (_Elem)delimAsci;
		auto delim = string_type<_Elem>(1, ch);
		using vector_type = std::vector<string_type<_Elem>>;
		std::shared_ptr<vector_type> ret(new vector_type());
		boost::split(*ret, s, boost::is_any_of(delim));
		return ret;
	}
	// trim from start (in place)
	template <typename _Elem>
	inline void ltrim(string_type<_Elem>& s) {
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
			return !std::isspace(ch);
		}));
	}
	template <typename _Elem>
	inline string_type<_Elem> ltrim_copy(const string_type<_Elem>& s) {
		auto copy = s;
		ltrim<_Elem>(copy);
		return copy;
	}

	// trim from end (in place)
	template <typename _Elem>
	inline void rtrim(string_type<_Elem>& s) {
		s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
			return !std::isspace(ch);
		}).base(), s.end());
	}
	template <typename _Elem>
	inline string_type<_Elem> rtrim_copy(const string_type<_Elem>& s) {
		auto copy = s;
		rtrim<_Elem>(copy);
		return copy;
	}

	// trim from both ends (in place)
	template <typename _Elem>
	inline void trim(string_type<_Elem>& s) {
		ltrim<_Elem>(s);
		rtrim<_Elem>(s);
	}
	template <typename _Elem>
	inline string_type<_Elem> trim_copy(string_type<_Elem>& s) {
		auto copy = s;
		trim<_Elem>(s);
		return copy;
	}

	inline void ucase(std::string& s) {
		std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
	}
	inline std::string ucase_copy(const std::string& s) {
		std::string copy = s;
		ucase(copy);
		return copy;
	}
	inline void lcase(std::string& s) {
		std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
	}
	inline std::string lcase_copy(const std::string& s) {
		std::string copy = s;
		lcase(copy);
		return copy;
	}
	inline std::string JSON_escape_string(
		const std::string& s
	) {
		std::stringstream ss;
		for (size_t i = 0; i < s.length(); ++i) {
			if (unsigned(s[i]) < '\x20') {
				ss << "\\u" << std::setfill('0') << std::setw(4) << std::hex << unsigned(s[i]);
			}
			else if (s[i] == '\\') {
				ss << "\\\\";
			}
			else if (s[i] == '"') {
				ss << "\\\"";
			}
			else {
				ss << s[i];
			}
		}
		return ss.str();
	}
	template<typename VT>
	inline std::string value_vector_join(
		const std::vector<VT>& v,
		char sep = ' ',
		size_t precision = 6
	) {
		std::string s;
		if (!v.empty()) {
			auto n = v.size();
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(precision);
			for (decltype(n) i = 0; i < n; ++i) {
				if (i > 0) {
					oss << " ";
				}
				oss << v[i];
			}
			s = oss.str();
		}
		return s;
	}

	template<typename VT>
	inline std::string value_vector_join_json(
		const std::vector<VT>& v,
		char sep = ' ',
		size_t precision = 6
	) {
		auto s = value_vector_join<VT>(v, sep, precision);
		if (!s.empty()) {
			std::ostringstream oss;
			oss << "\"" << JSON_escape_string(s) << "\"";
			s = oss.str();
		}
		else {
			s = "null";
		}
		return s;
	}
}