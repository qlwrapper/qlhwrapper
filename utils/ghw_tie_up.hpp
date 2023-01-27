#pragma once

#include <ql/quantlib.hpp>
#include <vector>
#include <cmath>
#include <utils/hull_white_params.hpp>

namespace utils {
    template <typename VolTraits = NormalVolTraits<QuantLib::Real>>
    class GeneralizedHullWhiteTieUp {
    public:
        typedef QuantLib::TermStructureFittingParameter::NumericalImpl NumericalImpl;
    private:
        GeneralizedHullWhiteModelParams<VolTraits> modelParams_;
        QuantLib::Handle<QuantLib::YieldTermStructure> termStructure_;
        QuantLib::Interpolation sigmaInterpolation_;
        std::vector<QuantLib::Time> muTimes_;
        std::vector<QuantLib::Real> muValues_;
        QuantLib::Interpolation muInterpolation_;

        class Helper {
        public:
            Helper(const QuantLib::Size i,
                const QuantLib::Real xMin,
                const QuantLib::Real dx,
                const QuantLib::Real discountBondPrice,
                const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
                QuantLib::ext::function<QuantLib::Real(QuantLib::Real)> fInv)
                : size_(tree->size(i)), dt_(tree->timeGrid().dt(i)), xMin_(xMin), dx_(dx),
                statePrices_(tree->statePrices(i)), discountBondPrice_(discountBondPrice),
                fInverse_(std::move(fInv)) {}

            QuantLib::Real operator()(const QuantLib::Real theta) const {
                QuantLib::Real value = discountBondPrice_;
                QuantLib::Real x = xMin_;
                for (QuantLib::Size j = 0; j < size_; j++) {
                    QuantLib::Real discount = std::exp(-fInverse_(theta + x) * dt_);
                    value -= statePrices_[j] * discount;
                    x += dx_;
                }
                return value;
            };

        private:
            QuantLib::Size size_;
            QuantLib::Time dt_;
            QuantLib::Real xMin_, dx_;
            const QuantLib::Array& statePrices_;
            QuantLib::Real discountBondPrice_;
            QuantLib::ext::function<QuantLib::Real(QuantLib::Real)> fInverse_;
        };

    public:
        GeneralizedHullWhiteTieUp(
            const GeneralizedHullWhiteModelParams<VolTraits>& modelParams,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& shortRateTree
        ) : modelParams_(modelParams), termStructure_(termStructure) {
            QuantLib::LinearFlat traits;
            sigmaInterpolation_ = traits.interpolate(modelParams.sigmaTimes().begin(), modelParams.sigmaTimes().end(), modelParams.sigma().begin());
            sigmaInterpolation_.enableExtrapolation();

            auto fInverse = VolTraits::f_inv();
            QuantLib::TermStructureFittingParameter phi(termStructure);
            auto impl = QuantLib::ext::dynamic_pointer_cast<NumericalImpl>(phi.implementation());
            auto const& grid = shortRateTree->timeGrid();
            impl->reset();
            QuantLib::Real value = 1.0;
            QuantLib::Real vMin = -50.0;
            QuantLib::Real vMax = 50.0;

            for (QuantLib::Size i = 0; i < (grid.size() - 1); i++) {
                QuantLib::Real discountBond = termStructure->discount(grid[i + 1]);
                QuantLib::Real xMin = shortRateTree->underlying(i, 0);
                QuantLib::Real dx = (i == 0 ? 0.0 : shortRateTree->underlying(i, 1) - xMin);
                Helper finder(i, xMin, dx, discountBond, shortRateTree, fInverse);
                QuantLib::Brent s1d;
                s1d.setMaxEvaluations(2000);
                value = s1d.solve(finder, 1e-8, value, vMin, vMax);
                // !!! solved value is mu not theta !!!
                // QuantLib's fitting parameters are for mu not theta
                impl->set(grid[i], value);
                muTimes_.push_back(grid[i]);
                muValues_.push_back(value);
            }
            muInterpolation_ = traits.interpolate(muTimes_.begin(), muTimes_.end(), muValues_.begin());
            muInterpolation_.enableExtrapolation();
        }
        const QuantLib::Real& a() const {
            return modelParams_.a();
        }
        const std::vector<QuantLib::Time>& sigmaTimes() const {
            return modelParams_.sigmaTimes();
        }
        const std::vector<QuantLib::Volatility>& sigma() const {
            return modelParams_.sigma();
        }
        const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure() const {
            return termStructure_;
        }
        const QuantLib::Interpolation& sigmaInterpolation() const {
            return sigmaInterpolation_;
        }
        const QuantLib::Interpolation& muInterpolation() const {
            return muInterpolation_;
        }
        QuantLib::Volatility sigmaAt(QuantLib::Time t) const {
            return sigmaInterpolation()(t);
        }
        QuantLib::Real B(QuantLib::Time t, QuantLib::Time T) const {
            return (1.0 - std::exp(-a() * (T - t))) / a();
        }
        QuantLib::Real mu(QuantLib::Time t) const {
            return muInterpolation()(t);
        }
        QuantLib::Real k(QuantLib::Time t, QuantLib::Time T) const {
            auto const& me = *this;
            auto integrationFunctor = [&me, &T](QuantLib::Time s) -> QuantLib::Real {
                auto sigma_s = me.sigmaAt(s);
                auto a = me.a();
                return std::pow(sigma_s * (1.0 - std::exp(-a * (T - s))), 2.0);
            };
            QuantLib::SimpsonIntegral integrator(1e-5, 1000);
            auto ret = integrator(integrationFunctor, 0.0, t);
            return ret;
        }
        QuantLib::Real v(QuantLib::Time t, QuantLib::Time T) const {
            auto u = k(t, t) - k(t, T);
            auto v = 1.0 / (2.0 * std::pow(a(), 2.0)) * u;
            return v;
        }
        // calculate DT / Dt
        QuantLib::DiscountFactor forwardDF(QuantLib::Time t, QuantLib::Time T) const {
            auto Dt = termStructure()->discount(t);
            auto DT = termStructure()->discount(T);
            return DT / Dt;
        }
        QuantLib::Real A(QuantLib::Time t, QuantLib::Time T) const {
            auto A = forwardDF(t, T) * std::exp(B(t, T) * mu(t) + v(t, T));
            return A;
        }
        QuantLib::Real discountBond(QuantLib::Time now, QuantLib::Time maturity, QuantLib::Rate rate) const {
            return A(now, maturity) * std::exp(-B(now, maturity) * rate);
        }
    };
}
