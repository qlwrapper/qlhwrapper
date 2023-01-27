#pragma once

#include <utils/types.h>
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