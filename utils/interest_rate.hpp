#pragma once

#include <ql/quantlib.hpp>

namespace utils {

    class InterestRateCalculator {
    private:
        QuantLib::Compounding compounding_;
        QuantLib::Frequency freq_;
    public:
        InterestRateCalculator(
            QuantLib::Compounding compounding = QuantLib::Continuous,
            QuantLib::Frequency freq = QuantLib::NoFrequency
        ) : compounding_(compounding), freq_(freq) {}
        const QuantLib::Compounding& compounding() const {
            return compounding_;
        }
        QuantLib::Compounding compounding() {
            return compounding_;
        }
        const QuantLib::Frequency& freq() const {
            return freq_;
        }
        QuantLib::Frequency& freq() {
            return freq_;
        }
        // calculate rate from discount factor
        QuantLib::Rate rateFromDiscountFactor(
            QuantLib::DiscountFactor df,
            QuantLib::Time t
        ) const {
            return QuantLib::InterestRate::impliedRate(1.0 / df, QuantLib::DayCounter(), compounding(), freq(), t).rate();
        }
         // calculate discount factor from rate
        QuantLib::DiscountFactor discountFactorFromRate(
            QuantLib::Rate rate,
            QuantLib::Time t
        ) const {
            QuantLib::InterestRate ir(rate, QuantLib::DayCounter(), compounding(), freq());
            return ir.discountFactor(t);
        }
        QuantLib::Rate operator() (QuantLib::DiscountFactor df, QuantLib::Time t) const {
            return rateFromDiscountFactor(df, t);
        }
    };

    class ParRateCalculator {
    private:
        QuantLib::Frequency couponFreq_;
        QuantLib::Time dt_; // time between coupon payments
    public:
        ParRateCalculator(QuantLib::Frequency couponFreq = QuantLib::Semiannual) : couponFreq_(couponFreq), dt_(1.0 / (QuantLib::Real)couponFreq) {}
        const QuantLib::Time& dt() const {
            return dt_;
        }
        const QuantLib::Frequency& couponFreq() const {
            return couponFreq_;
        }
        // return the number of coupon payments for the tenor and coupon frequency
        QuantLib::Size numCoupons(QuantLib::Time tenor) const {
            return (QuantLib::Size)std::round(tenor * (QuantLib::Real)couponFreq_);
        }
        QuantLib::Time couponTime(QuantLib::Time tStart, QuantLib::Size couponIndex) const {
            return (tStart + (QuantLib::Time)(couponIndex + 1) * dt());
        }
        template<typename ITER>
        QuantLib::Rate operator() (const ITER& dfBegin, const ITER& dfEnd) const {
            auto dfLast = *(dfEnd - 1);
            auto a = dt() * std::accumulate(dfBegin, dfEnd, 0.0);
            return (1.0 - dfLast) / a;
        }
    };
}