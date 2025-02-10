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

#include <utils/string-io-types.hpp>
#include <iostream>

namespace utils {
    // counsole io streams
    template <typename _Elem>
    struct console {
        ostream_type<_Elem>& cout() const {}
        ostream_type<_Elem>& cerr() const {}
        istream_type<_Elem>& cin() const {}
        ostream_type<_Elem>& clog() const {}
    };

    template <>
    struct console<char> {
        ostream_type<char>& cout() const { return std::cout; }
        ostream_type<char>& cerr() const { return std::cerr; }
        istream_type<char>& cin() const { return std::cin; }
        ostream_type<char>& clog() const { return std::clog; }
    };

    template <>
    struct console<wchar_t> {
        ostream_type<wchar_t>& cout() const { return std::wcout; }
        ostream_type<wchar_t>& cerr() const { return std::wcerr; }
        istream_type<wchar_t>& cin() const { return std::wcin; }
        ostream_type<wchar_t>& clog() const { return std::wclog; }
    };
}