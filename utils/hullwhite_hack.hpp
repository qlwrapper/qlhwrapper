#pragma once

#include <ql/quantlib.hpp>
#include <vector>

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
            const QuantLib::ext::function<QuantLib::Real(QuantLib::Real)>& f =
            QuantLib::ext::function<QuantLib::Real(QuantLib::Real)>(),
            const QuantLib::ext::function<QuantLib::Real(QuantLib::Real)>& fInverse =
            QuantLib::ext::function<QuantLib::Real(QuantLib::Real)>())
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