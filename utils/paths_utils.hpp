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
#include <memory>
#include <cmath>

namespace utils {
    class PathsUtils {
    public:
        static std::shared_ptr<QuantLib::TimeGrid> getEvenlySpacingTimeGrid(
            QuantLib::Time maturity,
            QuantLib::Natural freq
        ) {
            auto nSteps = (QuantLib::Size)std::round(maturity * (QuantLib::Time)freq);
            return std::shared_ptr<QuantLib::TimeGrid>(new QuantLib::TimeGrid(maturity, nSteps));
        }
    };
}