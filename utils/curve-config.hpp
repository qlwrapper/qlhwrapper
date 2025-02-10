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
#include <functional>
#include <map>

namespace utils {
    struct CurveConfig {
        typedef QuantLib::ext::shared_ptr<QuantLib::IborIndex> pIborIndex;
        typedef std::function<pIborIndex()> IndexFactory;
        QuantLib::Frequency parRateCouponFrequency;
        double parRateDayCountMultiplier;
        IndexFactory cashRateIndexFactory;
        CurveConfig(
            QuantLib::Frequency parRateCouponFrequency = QuantLib::Frequency::Semiannual,
            double parRateDayCountMultiplier = 1.,
            const IndexFactory& cashRateIndexFactory = nullptr
        ) :
            parRateCouponFrequency(parRateCouponFrequency),
            parRateDayCountMultiplier(parRateDayCountMultiplier),
            cashRateIndexFactory(cashRateIndexFactory)
        {}
        pIborIndex createCashRateIndex() const {
            QL_ASSERT(cashRateIndexFactory != nullptr, "index factory is not set");
            return cashRateIndexFactory();
        }
    };
    typedef std::map<bool, CurveConfig> CurveConfigs;
    const CurveConfigs& getCurveConfigs();
}
