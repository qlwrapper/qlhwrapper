#pragma once

#include <ql/quantlib.hpp>
#include <memory>

namespace utils {
    template <typename DT = QuantLib::Rate>
    class FlatRateShocker {
    private:
        QuantLib::Rate shock_;
    public:
        FlatRateShocker(QuantLib::Rate shock = 0.0) : shock_(shock) {}
        virtual DT operator() (const QuantLib::Time&, const DT&) const {    // return the amount of shock to apply to the value of type DT
            return (DT)shock_;
        }
        const QuantLib::Rate& shock() const {
            return shock_;
        }
        QuantLib::Rate& shock() {
            return shock_;
        }
    };

    template<QuantLib::Compounding COMP = QuantLib::Compounded, QuantLib::Frequency FREQ = QuantLib::Semiannual>
    class DiscountFactorFlatRateShocker : public FlatRateShocker<QuantLib::DiscountFactor> {
    public:
        DiscountFactorFlatRateShocker(QuantLib::Rate shock = 0.0) : FlatRateShocker<QuantLib::Rate>(shock) {}
        QuantLib::DiscountFactor operator() (
            const QuantLib::Time& t,
            const QuantLib::DiscountFactor& df
            ) const {
            if (t != 0.0) {
                auto rate = QuantLib::InterestRate::impliedRate(1.0 / df, QuantLib::DayCounter(), COMP, FREQ, t).rate();
                rate += this->shock();
                QuantLib::InterestRate rateShocked(rate, QuantLib::DayCounter(), COMP, FREQ);
                auto dfShocked = rateShocked.discountFactor(t);
                return dfShocked - df;
            }
            else {
                return 0.0;
            }
        }
    };
}