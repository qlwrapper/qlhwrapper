#pragma once

#include <ql/quantlib.hpp>
#include <utils/interest_rate.hpp>

namespace utils {
    template <typename VT = QuantLib::Rate>
    class ZVCalculator {
    private:
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> zvCurve_;
    public:
        ZVCalculator(
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve
        ) : zvCurve_(zvCurve) {}
        const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() const {
            return zvCurve_;
        }
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() {
            return zvCurve_;
        }
        virtual VT operator() (
            QuantLib::Time t,
            QuantLib::Time = 0.0
            ) const = 0;
        virtual QuantLib::Time tenor() const {
            return 0.0;
        }
    };

    template <
        QuantLib::Compounding COMPOUNDING = QuantLib::Continuous,
        QuantLib::Frequency FREQ = QuantLib::NoFrequency
    >
        class ZeroRateZVCalculator : public ZVCalculator<QuantLib::Rate> {
        public:
            ZeroRateZVCalculator(
                const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve
            ) : ZVCalculator<QuantLib::Rate>(zvCurve) {}
            QuantLib::Rate operator() (
                QuantLib::Time t,
                QuantLib::Time = 0.0
                ) const {
                if (t == 0.0) t = 0.0001;
                return zvCurve()->zeroRate(t, COMPOUNDING, FREQ, true).rate();
            }
    };

    class DiscountFactorZVCalculator : public ZVCalculator<QuantLib::DiscountFactor> {
    public:
        DiscountFactorZVCalculator(
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve
        ) : ZVCalculator<QuantLib::DiscountFactor>(zvCurve) {}
        QuantLib::DiscountFactor operator() (
            QuantLib::Time t,
            QuantLib::Time = 0.0
            ) const {
            return zvCurve()->discount(t, true);
        }
    };

    template <
        QuantLib::Natural TENOR_MONTHS,
        QuantLib::Compounding COMPOUNDING = QuantLib::Continuous,
        QuantLib::Frequency FREQ = QuantLib::NoFrequency
    >
        class ForwardRateZVCalculator : public ZVCalculator<QuantLib::Rate> {
        public:
            ForwardRateZVCalculator(
                const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve
            ) : ZVCalculator<QuantLib::Rate>(zvCurve) {}
            QuantLib::Time tenor() const {
                return (QuantLib::Time)TENOR_MONTHS / 12.0;
            }
            QuantLib::Rate operator() (
                QuantLib::Time t,
                QuantLib::Time = 0.0
                ) const {
                auto T = t + tenor();
                return zvCurve()->forwardRate(t, T, COMPOUNDING, FREQ, true).rate();
            }
    };

    class ShortRateZVCalculator : public ZVCalculator<QuantLib::Rate> {
    public:
        ShortRateZVCalculator(
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve
        ) : ZVCalculator<QuantLib::Rate>(zvCurve) {}
        QuantLib::Rate operator() (
            QuantLib::Time t,
            QuantLib::Time delta_t = 0.0  // delta_t should be very small (instantaneous) 
            ) const {
            auto T = t + delta_t;
            return zvCurve()->forwardRate(t, T, QuantLib::Continuous, QuantLib::NoFrequency, true).rate();
        }
    };

    template <QuantLib::Natural TENOR, QuantLib::Frequency COUPON_FREQ = QuantLib::Semiannual>
    class ParRateZVCalculator : public ZVCalculator<QuantLib::Rate> {
    public:
        ParRateZVCalculator(
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve
        ) : ZVCalculator<QuantLib::Rate>(zvCurve) {}
        QuantLib::Time tenor() const {
            return (QuantLib::Time)TENOR;
        }
        QuantLib::Rate operator() (
            QuantLib::Time t,
            QuantLib::Time = 0.0
            ) const {
            ParRateCalculator rateCalculator(COUPON_FREQ);
            auto numCoupons = rateCalculator.numCoupons(tenor());
            auto dfStart = zvCurve()->discount(t);
            std::vector<QuantLib::DiscountFactor> dfs(numCoupons);
            for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                auto couponTime = rateCalculator.couponTime(t, k);
                auto df = zvCurve()->discount(couponTime) / dfStart;
                dfs[k] = df;
            }
            return rateCalculator(dfs.begin(), dfs.end());
        }
    };
}