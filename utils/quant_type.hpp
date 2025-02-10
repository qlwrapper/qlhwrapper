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
#include <vector>
#include <memory>
#include <exception>

namespace utils {
    struct OptionAttribs {
        QuantLib::Period expiry;
        QuantLib::Period tenor;
        OptionAttribs(
            const QuantLib::Period& expiry = QuantLib::Period(),
            const QuantLib::Period& tenor = QuantLib::Period()
        ): expiry(expiry), tenor(tenor) {}
    };

    template <
        typename DT = QuantLib::Volatility
    >
    struct BlackVolData_T {
        typedef DT DataType;
        DataType data;                      // volatility/sigma
        QuantLib::VolatilityType volType;   // volatility type
        BlackVolData_T(
            DataType data = 0,
            QuantLib::VolatilityType volType = QuantLib::VolatilityType::ShiftedLognormal
        ) : data(data), volType(volType) {}
    };

    typedef BlackVolData_T<QuantLib::Volatility> BlackVolData;

    template <
        typename DT = QuantLib::Volatility
    >
    struct OptionVolData :
        public OptionAttribs,
        public BlackVolData_T<DT> {
        OptionVolData(
            typename BlackVolData_T<DT>::DataType data = 0,
            QuantLib::VolatilityType volType = QuantLib::VolatilityType::ShiftedLognormal
        ) : BlackVolData_T<DT>(data, volType) {}
    };

    typedef OptionVolData<QuantLib::Volatility> QuotedVol;
    typedef std::shared_ptr<QuotedVol> pQuotedVol;
    typedef std::vector<pQuotedVol> QuotedVols;
    typedef std::shared_ptr<QuotedVols> pQuotedVols;
}