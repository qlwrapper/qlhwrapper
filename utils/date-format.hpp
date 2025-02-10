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

#include <ql/quantlib.hpp>
#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>

namespace utils {
    template<typename _Elem>
    struct DateFormat {
        typedef std::basic_string<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>> string;
        typedef std::basic_ostringstream<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>> ostringstream;
        static QuantLib::Date from_yyyymmdd(const string& yyyymmdd, bool hasHyphen = false) {
            auto yyyy = yyyymmdd.substr(0, 4);
            auto mm = yyyymmdd.substr(4 + (std::size_t)(hasHyphen ? 1 : 0), 2);
            auto dd = yyyymmdd.substr(6 + (std::size_t)(hasHyphen ? 2 : 0));
            auto year = boost::lexical_cast<QuantLib::Year>(yyyy);
            auto month = QuantLib::Month(boost::lexical_cast<int>(mm));
            auto day = boost::lexical_cast<QuantLib::Day>(dd);
            return QuantLib::Date(day, month, year);
        }
        static string to_yyyymmdd(const QuantLib::Date& d, bool hyphen = false) {
            auto d_i = d.year() * 10000 + int(d.month()) * 100 + d.dayOfMonth();
            auto s = boost::lexical_cast<string>(d_i);
            if (hyphen) {
                auto yyyy = s.substr(0, 4);
                auto mm = s.substr(4, 2);
                auto dd = s.substr(6);
                ostringstream os;
                os << yyyy << "-" << mm << "-" << dd;
                return os.str();
            }
            else {
                return s;
            }
        }
    };
}
