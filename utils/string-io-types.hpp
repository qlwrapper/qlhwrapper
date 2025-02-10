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