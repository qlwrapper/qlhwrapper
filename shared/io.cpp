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

#include <pch.h>
#include <iostream>
#include <ql/quantlib.hpp>
#include <ql_utils/all.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <utils/io.hpp>

using namespace std;
using namespace utils;
using namespace QuantLib;
using namespace QLUtils;
using namespace boost;

#include <common/json-obj-serialization.hpp>

ostream& operator << (
    ostream& os,
    const OptionAttribs& rhs
) {
    return (os << rhs.expiry << "x" << rhs.tenor);
}

ostream& operator << (
    ostream& os,
    const BlackVolData& rhs
) {
    if (rhs.volType == QuantLib::VolatilityType::ShiftedLognormal) {
        os << QuantLib::io::volatility(rhs.data);
    }
    else if (rhs.volType == QuantLib::VolatilityType::Normal) {
        os << QuantLib::io::normal_volatility(rhs.data);
    }
    else {
        QL_FAIL("unknown/unsupported volality type (" << (int)(rhs.volType) << ")");
    }
    return os;
}

ostream& operator << (
    ostream& os,
    const QuotedVol& rhs
) {
    return (os << ((const OptionAttribs&)rhs) << " = " << ((const BlackVolData&)rhs));
}

QuotedVol& operator << (
    QuotedVol& dest,
    const property_tree::ptree& ptree
) {
    QuantLib::Real volConversionFactor = 1.;
    auto op = ptree.get_child_optional("expiry");
    if (op) {
        auto s = op.value().get_value<std::string>();
        dest.expiry = lexical_cast<QuantLib::Period>(s);
    }
    op = ptree.get_child_optional("tenor");
    if (op) {
        auto s = op.value().get_value<std::string>();
        dest.tenor = lexical_cast<QuantLib::Period>(s);
    }
    op = ptree.get_child_optional("isNormalVol");
    if (op) {
        auto isNormalVol = op.value().get_value<bool>();
        dest.volType = (isNormalVol ? QuantLib::VolatilityType::Normal : QuantLib::VolatilityType::ShiftedLognormal);
        volConversionFactor = (isNormalVol ? 0.0001 : 0.01);   // conversion factor to convert volality to decimal unit
    }
    op = ptree.get_child_optional("value");
    if (op) {
        dest.data = op.value().get_value<QuantLib::Volatility>();
        dest.data *= volConversionFactor;   // convert to decimal
    }
    return dest;
}

TermStructureNode<>& operator << (
    TermStructureNode<>& dest,
    const property_tree::ptree& ptree
) {
    auto op = ptree.get_child_optional("term");
    if (op) {
        dest.term = op.value().get_value<QuantLib::Time>();
    }
    else {
        QL_FAIL("term structure node object is missing the field 'term'");
    }
    op = ptree.get_child_optional("rate");
    if (op) {
        dest.rate = op.value().get_value<QuantLib::Rate>();
        dest.rate *= 0.01;  // convert to decimal
    }
    else {
        QL_FAIL("term structure node object is missing the field 'rate'");
    }
    return dest;
}

TermStructureNodes<>& operator << (
    TermStructureNodes<>& dest,
    const property_tree::ptree& ptree
) {
    dest.pTSNodes() << ptree;
    return dest;
}