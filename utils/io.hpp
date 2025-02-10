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
#include <ql_utils/all.hpp>
#include <utils/quant_type.hpp>
#include <string>
#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>

namespace boost {
    template<>
    inline QuantLib::Period lexical_cast<QuantLib::Period, std::string>(
        const std::string& s
    ) {
        QL_REQUIRE(!s.empty(), "not a valid period string: string is empty");
        auto s_n = s.substr(0, s.length() - 1);
        auto s_units = s.substr(s.length() - 1, 1);
        QuantLib::Integer n = (s_n.empty() ? 0 : boost::lexical_cast<QuantLib::Integer>(s_n));
        auto unit = boost::lexical_cast<QuantLib::TimeUnit>(s_units);
        return QuantLib::Period(n, unit);
    }
    template<>
    inline std::string lexical_cast<std::string, QuantLib::Period>(
        const QuantLib::Period& p
    ) {
        std::ostringstream oss;
        oss << p;
        return oss.str();
    }
}

std::ostream& operator << (std::ostream&, const utils::OptionAttribs&);
std::ostream& operator << (std::ostream&, const utils::BlackVolData&);
std::ostream& operator << (std::ostream&, const utils::QuotedVol&);
utils::QuotedVol& operator << (utils::QuotedVol&, const boost::property_tree::ptree&);
QLUtils::TermStructureNode<>& operator << (QLUtils::TermStructureNode<>&, const boost::property_tree::ptree&);
QLUtils::TermStructureNodes<>& operator << (QLUtils::TermStructureNodes<>&, const boost::property_tree::ptree&);