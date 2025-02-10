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
#include <memory>
#include <vector>

namespace utils {
    typedef std::vector<QuantLib::Real> BackwardInductionValues;

    // calculate P(t, T) for each node on the lattice at time slice t
    class BackwardInductionDiscountBondCalculator {
    public:
        // calculation operator
        virtual QuantLib::Array operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Time t,
            QuantLib::Time T
            ) const = 0;
    };

    class BackwardInductionCalculator {
    private:
        std::shared_ptr<BackwardInductionDiscountBondCalculator> discountBondCalculator_;
    protected:
        BackwardInductionCalculator(
            const std::shared_ptr<BackwardInductionDiscountBondCalculator>& discountBondCalculator
        ) : discountBondCalculator_(discountBondCalculator){}
    public:
        const std::shared_ptr<BackwardInductionDiscountBondCalculator>& discountBondCalculator() const {
            return discountBondCalculator_;
        }
        // calculation operator
        virtual std::shared_ptr<BackwardInductionValues> operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Size expectedNumValues,
            QuantLib::Time t,
            QuantLib::Time T
            ) const = 0;
    };

    // calculate backward induction zero coupon bond calculator P(t, T)
    class BackwardInductionZCBCalculator : public BackwardInductionCalculator {
    public:
        BackwardInductionZCBCalculator(
            const std::shared_ptr<BackwardInductionDiscountBondCalculator>& discountBondCalculator
        ) : BackwardInductionCalculator(discountBondCalculator) {}
        // calculation operator
        std::shared_ptr<BackwardInductionValues> operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Size expectedNumValues,
            QuantLib::Time t,
            QuantLib::Time T
            ) const {
            auto const& discountBondCalculator = *(this->discountBondCalculator());
            auto dfs = discountBondCalculator(tree, t, T);
            QL_ASSERT(dfs.size() == expectedNumValues, "dfs.size() != expectedNumValues");
            std::shared_ptr<BackwardInductionValues> values(new BackwardInductionValues(expectedNumValues));
            std::copy(dfs.begin(), dfs.end(), values->begin());
            return values;
        }
    };

    // calculate backward induction forward interest rates
    class BackwardInductionFwdRateCalculator: public BackwardInductionCalculator {
    private:
        InterestRateCalculator rateCalculator;
    public:
        BackwardInductionFwdRateCalculator(
            const std::shared_ptr<BackwardInductionDiscountBondCalculator>& discountBondCalculator,
            QuantLib::Compounding compounding = QuantLib::Continuous,
            QuantLib::Frequency freq = QuantLib::NoFrequency
            
        ) : BackwardInductionCalculator(discountBondCalculator), rateCalculator(compounding, freq) {}
        // calculation operator
        std::shared_ptr<BackwardInductionValues> operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Size expectedNumValues,
            QuantLib::Time t,
            QuantLib::Time T
            ) const {
            auto const& discountBondCalculator = *(this->discountBondCalculator());
            auto dfs = discountBondCalculator(tree, t, T);
            QL_ASSERT(dfs.size() == expectedNumValues, "dfs.size() != expectedNumValues");
            std::shared_ptr<BackwardInductionValues> values(new BackwardInductionValues(expectedNumValues));
            auto tenor = T - t; // tenor of the forward rate
            for (decltype(expectedNumValues) i = 0; i < expectedNumValues; ++i) {
                const auto& df = dfs.at(i);
                auto rate = rateCalculator(df, tenor);
                values->at(i) = rate;
            }
            return values;
        }
    };

    // calculate backward induction forward par rates
    class BackwardInductionParRateCalculator: public BackwardInductionCalculator {
    private:
        ParRateCalculator rateCalculator;
    public:
        BackwardInductionParRateCalculator(
            const std::shared_ptr<BackwardInductionDiscountBondCalculator>& discountBondCalculator,
            QuantLib::Frequency couponFreq = QuantLib::Semiannual
        ) : BackwardInductionCalculator(discountBondCalculator), rateCalculator(couponFreq) {}

        // calculation operator
        std::shared_ptr<BackwardInductionValues> operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Size expectedNumValues,
            QuantLib::Time t,
            QuantLib::Time T
            ) const {
            auto const& discountBondCalculator = *(this->discountBondCalculator());
            auto tenor = T - t; // tenor of the swap
            auto numCoupons = rateCalculator.numCoupons(tenor);
            QuantLib::Matrix dfMatrix(expectedNumValues, numCoupons);
            for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                auto couponTime = rateCalculator.couponTime(t, k);
                auto dfs = discountBondCalculator(tree, t, couponTime);
                QL_ASSERT(dfs.size() == expectedNumValues, "dfs.size() != expectedNumValues");
                std::copy(dfs.begin(), dfs.end(), dfMatrix.column_begin(k));
            }
            std::shared_ptr<BackwardInductionValues> values(new BackwardInductionValues(expectedNumValues));
            for (decltype(expectedNumValues) i = 0; i < expectedNumValues; ++i) {   // for each row/node
                auto parRate = rateCalculator(dfMatrix.row_begin(i), dfMatrix.row_end(i));
                values->at(i) = parRate;
            }
            return values;
        }
    };

    // calculate backward induction for CMS break even rate
    class BackwardInductionCMSBreakevenCalculator : public BackwardInductionCalculator {
    private:
        ParRateCalculator rateCalculator;
        QuantLib::Time dt_CMS_;
    public:
        BackwardInductionCMSBreakevenCalculator(
            const std::shared_ptr<BackwardInductionDiscountBondCalculator>& discountBondCalculator,
            QuantLib::Frequency couponFreq = QuantLib::Semiannual,
            QuantLib::Frequency cmsFreq = QuantLib::Semiannual
        ) : BackwardInductionCalculator(discountBondCalculator), rateCalculator(couponFreq), dt_CMS_(1.0/ (QuantLib::Real)cmsFreq){}

        const QuantLib::Time& dt_CMS() const {
            return dt_CMS_;
        }

        // calculate par rate for each node @ t
        std::shared_ptr<std::vector<QuantLib::Rate>> calculateParRates(
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Size expectedNumValues,
            QuantLib::Time t,
            QuantLib::Time T
        ) const {
            auto const& discountBondCalculator = *(this->discountBondCalculator());
            auto tenor = T - t; // tenor of the swap
            auto numCoupons = rateCalculator.numCoupons(tenor);
            QuantLib::Matrix dfMatrix(expectedNumValues, numCoupons);
            for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                auto couponTime = rateCalculator.couponTime(t, k);
                auto dfs = discountBondCalculator(tree, t, couponTime); // calculate each node's P(t, couponTime) @ t
                QL_ASSERT(dfs.size() == expectedNumValues, "dfs.size() != expectedNumValues");
                std::copy(dfs.begin(), dfs.end(), dfMatrix.column_begin(k));
            }
            std::shared_ptr<std::vector<QuantLib::Rate>> parRates(new std::vector<QuantLib::Rate>(expectedNumValues));
            for (decltype(expectedNumValues) index = 0; index < expectedNumValues; ++index) {   // for each row/node
                auto parRate = rateCalculator(dfMatrix.row_begin(index), dfMatrix.row_end(index));
                parRates->at(index) = parRate;
            }
            return parRates;
        }

        QuantLib::Rate calculateCMSBreakevenRate(
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Time t,
            const std::vector<QuantLib::Rate>& parRates
        ) const {
            auto numNodes = parRates.size();
            auto i = tree->timeGrid().closestIndex(t);
            auto statePrices = tree->statePrices(i);
            QL_ASSERT(statePrices.size() == numNodes, "statePrices.size() != numNodes");
            auto const& discountBondCalculator = *(this->discountBondCalculator());
            auto dt = dt_CMS();
            auto T = t + dt;    // CMS coupon time
            auto dfs = discountBondCalculator(tree, t, T);  // calculate each node's P(t, T) @ t
            QL_ASSERT(dfs.size() == numNodes, "dfs.size() != numNodes");
            auto singleCouponPVCalculator = [&dt, &dfs, &statePrices](const std::vector<QuantLib::Rate>& coupons) -> QuantLib::Real {
                auto numNodes = dfs.size();
                QuantLib::Array cfPVs(numNodes);
                for (decltype(numNodes) index = 0; index < numNodes; ++index) { // for each node at t
                    auto const& df = dfs.at(index); // P(t, T) = discount factor
                    auto const& coupon = coupons.at(index);
                    auto cashflowPV = coupon * dt * df;
                    cfPVs.at(index) = cashflowPV;
                }
                auto PV = QuantLib::DotProduct(cfPVs, statePrices); // PV @ t=0
                return PV;
            };
            auto targetPV = singleCouponPVCalculator(parRates);
            auto solverFunctor = [&numNodes, &singleCouponPVCalculator, &targetPV](QuantLib::Rate breakevenCoupon) -> QuantLib::Real {
                std::vector<QuantLib::Rate> coupons(numNodes, breakevenCoupon);
                auto calcPV = singleCouponPVCalculator(coupons);
                return calcPV - targetPV;
            };
            QuantLib::Brent s1d;
            s1d.setMaxEvaluations(2000);
            auto guess = std::accumulate(parRates.begin(), parRates.end(), 0.0) / numNodes; // guess = average of par rates
            auto breakevenCoupon = s1d.solve(solverFunctor, 1e-12, guess, 0.01);
            return breakevenCoupon;
        }
        // calculation operator
        std::shared_ptr<BackwardInductionValues> operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Size expectedNumValues,
            QuantLib::Time t,
            QuantLib::Time T
            ) const {
            auto parRates = calculateParRates(tree, expectedNumValues, t, T);
            auto breakEvenRate = calculateCMSBreakevenRate(tree, t, *parRates);
            std::shared_ptr<BackwardInductionValues> values(new BackwardInductionValues(expectedNumValues, breakEvenRate));
            return values;
        }
    };

    class DiscretizedDiscountBondCalculator : public BackwardInductionDiscountBondCalculator {
    public:
        QuantLib::Array operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Time t,
            QuantLib::Time T
            ) const {
            QuantLib::DiscretizedDiscountBond zeroCouponBond;
            zeroCouponBond.initialize(tree, T);
            zeroCouponBond.rollback(t);
            return zeroCouponBond.values();
        }
    };
}