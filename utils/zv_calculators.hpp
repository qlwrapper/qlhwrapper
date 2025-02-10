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
        // calculate zero rate @t
        QuantLib::Rate operator() (
            QuantLib::Time t,
            QuantLib::Time = 0.0    // un-used
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
        // calculate discount factor @t
        QuantLib::DiscountFactor operator() (
            QuantLib::Time t,
            QuantLib::Time = 0.0    // un-used
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
        // calculate continuously compounded short rate @t with maturity of very short delta_t
        QuantLib::Rate operator() (
            QuantLib::Time t,
            QuantLib::Time delta_t = 0.0  // delta_t should be very small (instantaneous) 
        ) const {
            auto T = t + delta_t;
            return zvCurve()->forwardRate(t, T, QuantLib::Continuous, QuantLib::NoFrequency, true).rate();
        }
    };

    class ParRateZVCalculator : public ZVCalculator<QuantLib::Rate> {
    protected:
        QuantLib::Time tenor_;
        QuantLib::Frequency couponFreq_;
        double dayCountMultiplier_;
    public:
        ParRateZVCalculator(
            QuantLib::Time tenor,    // tenor in years
            QuantLib::Frequency couponFreq,
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve,
            double dayCountMultiplier = 1.
        ) :
            ZVCalculator<QuantLib::Rate>(zvCurve),
            tenor_(tenor),
            couponFreq_(couponFreq),
            dayCountMultiplier_(dayCountMultiplier)
        {}
        QuantLib::Time tenor() const {
            return tenor_;
        }
        QuantLib::Frequency couponFreq() const {
            return couponFreq_;
        }
        double dayCountMultiplier() const {
            return dayCountMultiplier_;
        }
        // calculate par rate @t
        QuantLib::Rate operator() (
            QuantLib::Time t,
            QuantLib::Time = 0.0    // un-used
        ) const {
            ParRateCalculator rateCalculator(couponFreq());
            auto numCoupons = rateCalculator.numCoupons(tenor());
            auto dfStart = zvCurve()->discount(t);
            std::vector<QuantLib::DiscountFactor> dfs(numCoupons);
            for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                auto couponTime = rateCalculator.couponTime(t, k);
                auto df = zvCurve()->discount(couponTime) / dfStart;
                dfs[k] = df;
            }
            auto parRate = rateCalculator(dfs.begin(), dfs.end());
            parRate *= dayCountMultiplier();
            return parRate;
        }
    };
}