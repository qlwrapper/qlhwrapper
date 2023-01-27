#pragma once

#include <ql/quantlib.hpp>
#include <utils/types.h>
#include <utils/string.h>
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
        template <
            typename ZV_CALCULATOR,
            typename STATS_CONVERTER = IdentityStatsConverter
        >
        void report(
            const Paths& paths
        ) {
            auto& os = ostream();
            utf8_wstring_converter<_Elem> utf8_converter;
            auto err = 0.0;
            ZV_CALCULATOR zvCalculator(zvCurve());
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
                os << comma << utf8_converter.from_bytes("mean=") << mean * 100.0;
                os << comma << utf8_converter.from_bytes("stdev=") << stdev * 100.0;
                os << comma << utf8_converter.from_bytes("zvValue=") << zvValue * 100.0;
                os << comma << utf8_converter.from_bytes("diff=") << diff * 10000.0 << " bp";
                os << std::endl;
                err += std::pow(diff, 2.0);
            }
            err = std::sqrt(err);
            os << "err=" << err * 10000.0 << " bp" << std::endl;
        }
    };
}