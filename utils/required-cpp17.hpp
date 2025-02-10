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

#ifdef _MSC_VER
#if defined(_MSVC_LANG) && _MSVC_LANG >= 201703L
#define REQUIRED_CPP17 
#endif
#else
#if __cplusplus >= 201703L
#define REQUIRED_CPP17 
#endif
#endif

#if defined(REQUIRED_CPP17)
#include <string>
#include <clocale>
#include <vector>
#include <exception>
namespace std {
	namespace for_cpp17 {
		template<
			typename _Elem
		>
		struct codecvt_utf8 {
			typedef basic_string<_Elem, char_traits<_Elem>, allocator<_Elem>> string_type;
			string to_utf8(const string_type&) const {
				throw std::exception("Not implemented");
			}
			string_type from_utf8(const string&) const {
				throw std::exception("Not implemented");
			}
		};
		template <>
		inline string codecvt_utf8<char>::to_utf8(const typename codecvt_utf8<char>::string_type& s) const {
			return s;
		}
		template <>
		inline typename codecvt_utf8<char>::string_type codecvt_utf8<char>::from_utf8(const string& s) const {
			return s;
		}

		template <>
		inline string codecvt_utf8<wchar_t>::to_utf8(const wstring& wstr) const {
			// Set the locale to UTF-8
			setlocale(LC_ALL, "en_US.utf8");
			// Determine the required buffer size
			size_t bufferSize;
			wcstombs_s(&bufferSize, nullptr, 0, wstr.c_str(), 0);
			// Allocate buffer
			vector<char> buffer(bufferSize);
			// Perform the conversion
			wcstombs_s(&bufferSize, buffer.data(), bufferSize, wstr.c_str(), bufferSize);
			return string(buffer.data());
		}
		template <>
		inline wstring codecvt_utf8<wchar_t>::from_utf8(const string& utf8str) const {
			// Set the locale to UTF-8
			setlocale(LC_ALL, "en_US.utf8");
			// Determine the required buffer size
			size_t bufferSize;
			mbstowcs_s(&bufferSize, nullptr, 0, utf8str.c_str(), 0);
			// Allocate buffer
			vector<wchar_t> buffer(bufferSize);
			// Perform the conversion
			mbstowcs_s(&bufferSize, buffer.data(), bufferSize, utf8str.c_str(), bufferSize);
			return wstring(buffer.data());
		}
		template <
			typename CODECVD,
			typename _Elem = wchar_t,
			class Wide_alloc = allocator<_Elem>,
			class Byte_alloc = allocator<char>
		>
		struct wstring_convert {
			typedef basic_string<char, char_traits<char>, Byte_alloc> utf8_string;
			typedef basic_string<_Elem, char_traits<_Elem>, Wide_alloc> string_type;
			utf8_string to_bytes(const string_type& s) const {
				return CODECVD().to_utf8(s);
			}
			string_type from_bytes(const utf8_string& s) const {
				return CODECVD().from_utf8(s);
			}
		};
	}
}
#endif