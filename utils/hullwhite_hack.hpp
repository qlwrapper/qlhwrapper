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
#include <vector>
#include <functional>

namespace utils {
    // class to "hack" A(t, T), B(t, T) out of the class HullWhite
    class HullWhiteHack : public QuantLib::HullWhite {
    public:
        HullWhiteHack(
            const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
            QuantLib::Real a = 0.1,
            QuantLib::Real sigma = 0.01
        ) : QuantLib::HullWhite(termStructure, a, sigma) {}
        QuantLib::Real A_(QuantLib::Time t, QuantLib::Time T) const {
            return A(t, T);
        }
        QuantLib::Real B_(QuantLib::Time t, QuantLib::Time T) const {
            return B(t, T);
        }
    };

    // class to "hack" A(t, T), B(t, T), V(t, T) out of the class GeneralizedHullWhite
    class GeneralizedHullWhiteHack : public QuantLib::GeneralizedHullWhite {
    public:
        GeneralizedHullWhiteHack(
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure,
            const std::vector<QuantLib::Date>& speedstructure,
            const std::vector<QuantLib::Date>& volstructure,
            const std::vector<QuantLib::Real>& speed,
            const std::vector<QuantLib::Real>& vol,
            const std::function<QuantLib::Real(QuantLib::Real)>& f =
            std::function<QuantLib::Real(QuantLib::Real)>(),
            const std::function<QuantLib::Real(QuantLib::Real)>& fInverse =
            std::function<QuantLib::Real(QuantLib::Real)>())
            : QuantLib::GeneralizedHullWhite(
                yieldtermStructure,
                speedstructure,
                volstructure,
                speed,
                vol,
                f,
                fInverse
            ) {}
        QuantLib::Real A_(QuantLib::Time t, QuantLib::Time T) const {
            return A(t, T);
        }
        QuantLib::Real B_(QuantLib::Time t, QuantLib::Time T) const {
            return B(t, T);
        }
        QuantLib::Real V_(QuantLib::Time t, QuantLib::Time T) const {
            return V(t, T);
        }
    };
}