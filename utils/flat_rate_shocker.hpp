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