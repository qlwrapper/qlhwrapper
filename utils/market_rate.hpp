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
#include <vector>
#include <memory>
#include <ql/quantlib.hpp>
#include <ql_utils/all.hpp>
#include <utils/quant_type.hpp>

namespace utils {
    class MarketRate {
    public:
        typedef std::shared_ptr<QLUtils::YieldTSNodes<QuantLib::Time, QuantLib::Rate>> pTSNodes;
    public:
        static QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> getActual365FixedZeroTermStructure(
            const QuantLib::Date& curveRefDate,
            const std::vector<QuantLib::Time>& act365FixedTimes,
            const std::vector<QuantLib::Rate>& zeroRates,
            bool enableExtrapolation = true
        );

        static QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> getThirty360ZeroTermStructure(
            const QuantLib::Date& curveRefDate,
            const std::vector<QuantLib::Time>& act365FixedTimes,
            const std::vector<QuantLib::Rate>& zeroRates,
            bool enableExtrapolation = true
        );

        static pTSNodes readAct365FixedZeroCurveTS(
            const std::string& inputArg
        );

        static QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> readAsThirty360ZeroCurve(
            const std::string& inputArg,
            const QuantLib::Date& curveRefDate,
            bool enableExtrapolation = true
        );

        static QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> readAsActual365FixedZeroCurve(
            const std::string& inputArg,
            const QuantLib::Date& curveRefDate,
            bool enableExtrapolation = true
        );

        static pQuotedSwaptionVols readQuotedSwaptionVols(
            const std::string& inputArg
        );
    };
}