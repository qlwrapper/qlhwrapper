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
#include <utils/string-io-types.hpp>
#include <utils/string.hpp>
#include <utils/paths.hpp>
#include <vector>
#include <cmath>

namespace utils {
    struct IdentityStatsConverter {
        QuantLib::Real operator() (QuantLib::Real value, QuantLib::Time, QuantLib::Time) const {
            return value;
        }
    };

    template <QuantLib::Compounding COMPOUNDING = QuantLib::Continuous, QuantLib::Frequency FREQ = QuantLib::NoFrequency>
    struct DiscountFactorToZeroRateConverter {
        QuantLib::Rate operator() (QuantLib::DiscountFactor df, QuantLib::Time t, QuantLib::Time = 0.0) const {
            if (t == 0.0) t = 0.0001;
            return QuantLib::InterestRate::impliedRate(1.0 / df, QuantLib::DayCounter(), COMPOUNDING, FREQ, t).rate();
        }
    };

    template <typename _Elem = char>
    class ZVComparison {
    protected:
        ostream_type<_Elem>& ostream_;
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> zvCurve_;
        std::vector<QuantLib::Natural> importantMaturities_;
    public:
        ZVComparison(
            ostream_type<_Elem>& ostream,
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve,
            const std::vector<QuantLib::Natural>& importantMaturities
        ) : ostream_(ostream), zvCurve_(zvCurve), importantMaturities_(importantMaturities)
        {}
        ostream_type<_Elem>& ostream() {
            return ostream_;
        }
        const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() const {
            return zvCurve_;
        }
        const std::vector<QuantLib::Natural>& importantMaturities() const {
            return importantMaturities_;
        }
        template<
            typename ZV_CALCULATOR,
            typename STATS_CONVERTER = IdentityStatsConverter
        >
        void report(
            const ZV_CALCULATOR& zvCalculator,
            const Paths& paths
        ) {
            auto& os = ostream();
            utf8_wstring_converter<_Elem> utf8_converter;
            auto err = 0.0;
            STATS_CONVERTER statsConverter;
            auto const& timeGrid = *(paths.timeGrid());
            auto statistics = paths.statistics();
            auto grid_delta_t = timeGrid.dt(0);
            auto comma = utf8_converter.from_bytes(",");
            for (auto p = importantMaturities().begin(); p != importantMaturities().end(); ++p) {
                auto const& maturity = *p;
                QuantLib::Time t = maturity;
                auto index = timeGrid.closestIndex(t);
                QL_ASSERT(QuantLib::close_enough(timeGrid.at(index), t), "timeGrid.at(index) (" << timeGrid.at(index) << ") and t (" << t << ") is not close enough");
                auto const& stats = statistics->at(index);
                auto zvValue = zvCalculator(t, grid_delta_t);
                auto mean = statsConverter(stats.first, t, grid_delta_t);
                auto stdev = statsConverter(stats.second, t, grid_delta_t);
                auto diff = mean - zvValue;
                os << t;
                os << comma << utf8_converter.from_bytes("mean=") << QuantLib::io::percent(mean);
                os << comma << utf8_converter.from_bytes("stdev=") << QuantLib::io::percent(stdev);
                os << comma << utf8_converter.from_bytes("zvValue=") << QuantLib::io::percent(zvValue);
                os << comma << utf8_converter.from_bytes("diff=") << QuantLib::io::basis_point(diff);
                os << std::endl;
                err += std::pow(diff, 2.0);
            }
            err = std::sqrt(err);
            os << "err=" << QuantLib::io::basis_point(err) << std::endl;
        }
        template <
            typename ZV_CALCULATOR,
            typename STATS_CONVERTER = IdentityStatsConverter
        >
        void report(
            const Paths& paths
        ) {
            ZV_CALCULATOR zvCalculator(zvCurve());
            report<ZV_CALCULATOR, STATS_CONVERTER>(zvCalculator, paths);
        }
    };
}