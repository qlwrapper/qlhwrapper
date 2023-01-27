#pragma once

#include <vector>
#include <string>
#include <memory>
#include <locale>
#include <codecvt>
#include <utils/types.h>
#include <boost/algorithm/string.hpp>

namespace utils {
	// convert between utf8 (char) string and wstring of type char16_t, char32_t, or wchar_t
	template <typename _Elem> using utf8_wstring_converter = std::wstring_convert<std::codecvt_utf8<_Elem>, _Elem>;

	// parse string into vector of string by a delimter character
	template <typename _Elem>
	std::shared_ptr<std::vector<string_type<_Elem>>> parse_delimited(const string_type<_Elem>& s, unsigned int delimAsci) {
		static string_type<_Elem> delim;
		if (delim.length() == 0) {
			auto ch = (_Elem)delimAsci;
			delim = string_type<_Elem>(1, ch);
		}
		using vector_type = std::vector<string_type<_Elem>>;
		std::shared_ptr<vector_type> ret(new vector_type());
		boost::split(*ret, s, boost::is_any_of(delim));
		return ret;
	}
};