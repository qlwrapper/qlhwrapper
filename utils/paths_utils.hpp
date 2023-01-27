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