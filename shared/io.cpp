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
    const SwaptionAttribs& rhs
) {
    return (os << rhs.toString());
}

ostream& operator << (
    ostream& os,
    const VolatilityData& rhs
) {
    if (rhs.volType == VolatilityType::ShiftedLognormal) {
        os << io::volatility(rhs.vol);
    }
    else if (rhs.volType == VolatilityType::Normal) {
        os << io::normal_volatility(rhs.vol);
    }
    else {
        QL_FAIL("unknown/unsupported volality type (" << (int)(rhs.volType) << ")");
    }
    auto hasLognormalShift = rhs.hasLognormalShift();
    auto isSkewed = rhs.isSkewed();
    auto hasExtraInfo = (hasLognormalShift || isSkewed);
    if (hasExtraInfo) {
        os << "(";
        if (hasLognormalShift) {
            os << "shift=" << io::basis_point(rhs.shift);
        }
        if (isSkewed) {
            if (hasLognormalShift) {
                os << ", ";
            }
            os << "skew=" << io::basis_point(rhs.skew);
        }
        os << ")";
    }
    return os;
}

ostream& operator << (
    ostream& os,
    const QuotedSwaptionVol& rhs
) {
    return (os << ((const SwaptionAttribs&)rhs) << " = " << ((const VolatilityData&)rhs));
}

QuotedSwaptionVol& operator << (
    QuotedSwaptionVol& dest,
    const property_tree::ptree& ptree
) {
    Real volConversionFactor = 1.;
    auto op = ptree.get_child_optional("expiry");
    if (op) {
        auto s = op.value().get_value<string>();
        dest.expiry = lexical_cast<Period>(s);
    }
    op = ptree.get_child_optional("tenor");
    if (op) {
        auto s = op.value().get_value<string>();
        dest.tenor = lexical_cast<Period>(s);
    }
    op = ptree.get_child_optional("isNormalVol");
    if (op) {
        auto isNormalVol = op.value().get_value<bool>();
        dest.volType = (isNormalVol ? VolatilityType::Normal : VolatilityType::ShiftedLognormal);
        volConversionFactor = (isNormalVol ? 0.0001 : 0.01);   // conversion factor to convert volality to decimal unit
    }
    op = ptree.get_child_optional("value");
    if (op) {
        dest.vol = op.value().get_value<Volatility>();
        dest.vol *= volConversionFactor;   // convert to decimal
    }
    if (dest.volType == VolatilityType::ShiftedLognormal) {
        op = ptree.get_child_optional("shift");
        if (op) {
            dest.shift = op.value().get_value<Real>();
            dest.shift /= 10000.;   // convert from basis point to decimal
        }
    }
    op = ptree.get_child_optional("skew");
    if (op) {
        dest.skew = op.value().get_value<Real>();
        dest.skew /= 10000.;   // convert from basis point to decimal
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