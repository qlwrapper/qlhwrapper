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
#include <utils/curve-config.hpp>

using namespace std;
using namespace utils;
using namespace QuantLib;

static map<bool, CurveConfig> curveConfigs_{
    {
        false,  // isSofr = false
        {
            Frequency::Semiannual,
            1.,
            []() {return ::ext::shared_ptr<IborIndex>(new USDLibor3M()); }
        }
    },
    {
        true,  // isSofr = true
        {
            Frequency::Annual,
            360. / 365.,
            []() {return ::ext::shared_ptr<IborIndex>(new UsdOvernightCompoundedAverageIndex<Sofr>()); }
        }
    }
};

namespace utils {
    const CurveConfigs& getCurveConfigs() {
        return curveConfigs_;
    }
}