#pragma once

#include <ql/quantlib.hpp>
#include <utils/backward-induction.hpp>
#include <utils/short_rate_tree_adaptor.hpp>
#include <cmath>

namespace utils {
    // OneFactorAffineLike is an abtraction for an entity that can calculate A(t, T) and B(t, T) for a one factor affine like model
    struct OneFactorAffineLike {
        // calculate both A(t, T) and B(t, T) simultaneously
        virtual std::pair<QuantLib::Real, QuantLib::Real> calculateAB(
            QuantLib::Time t,
            QuantLib::Time T
        ) const = 0;
    };
    // calculate P(t, T) given a OneFactorAffineLike entity
    class OneFactorAffineDiscountBondCalculator : public BackwardInductionDiscountBondCalculator {
    private:
        std::shared_ptr<OneFactorAffineLike> oneFactorAffineLike_;
    public:
        OneFactorAffineDiscountBondCalculator(
            const std::shared_ptr<OneFactorAffineLike>& oneFactorAffineLike
        ) : oneFactorAffineLike_(oneFactorAffineLike) {}
        const std::shared_ptr<OneFactorAffineLike>& oneFactorAffineLike() const {
            return oneFactorAffineLike_;
        }
        static QuantLib::Real discountBond(
            QuantLib::Real A,
            QuantLib::Real B,
            QuantLib::Rate shortRate
        ) {
            return A * std::exp(-B * shortRate);
        }
        QuantLib::Array operator() (
            const QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree>& tree,
            QuantLib::Time t,
            QuantLib::Time T
            ) const {
            auto AB = oneFactorAffineLike()->calculateAB(t, T);
            auto const& A = AB.first;
            auto const& B = AB.second;
            ShorRateTreeAdaptor treeAdaptor(tree);
            auto i = tree->timeGrid().closestIndex(t);
            auto numNodes = tree->size(i);
            QuantLib::Array values(numNodes);
            for (decltype(numNodes) index = 0; index < numNodes; ++index) { // for each node at the t time slice
                auto shortRate = treeAdaptor.shortRateAtNode(i, index); // get the short rate at this node
                values.at(index) = discountBond(A, B, shortRate);   // calculate P(t, T) for this node
            }
            return values;
        }
    };
}