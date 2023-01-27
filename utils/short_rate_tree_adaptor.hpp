#pragma once

#include <ql/quantlib.hpp>
#include <memory>

namespace utils {
    template <typename RNG>
    class RandomTreeWalker {
    public:
        typedef QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree> pTree;
    private:
        pTree shortRateTree_;
        const RNG& randomNumberGenerator_;
    public:
        RandomTreeWalker(
            const pTree& shortRateTree,
            const RNG& randomNumberGenerator
        ) : shortRateTree_(shortRateTree), randomNumberGenerator_(randomNumberGenerator) {}
        const pTree& tree() const {
            return shortRateTree_;
        }
        const RNG& randomNumberGenerator() const {
            return randomNumberGenerator_;
        }
        // from the specified node on the tree, randomly walk to the next node
        QuantLib::Size walkToNextNode(QuantLib::Size i, QuantLib::Size index) {
            auto const& shortRateTree = tree();
            auto const& rngGenerator = randomNumberGenerator();
            auto p1 = shortRateTree->probability(i, index, 0);  // probability of walk up
            auto p2 = shortRateTree->probability(i, index, 1);  // probability of walk flat
            auto p3 = shortRateTree->probability(i, index, 2);  // probability of walk down
            // p1 + p2 + p3 = 1
            auto randomNumber = rngGenerator.nextReal();    // generate a random number
            QuantLib::Size nextNodeIndex = 0;
            if (randomNumber <= p1) {
                nextNodeIndex = shortRateTree->descendant(i, index, 0);    // walk up
            }
            else if (randomNumber > 1 - p3) {
                nextNodeIndex = shortRateTree->descendant(i, index, 2);    // walk flat
            }
            else {
                nextNodeIndex = shortRateTree->descendant(i, index, 1);    // walk down
            }
            return nextNodeIndex;
        }
    };
    class ShorRateTreeAdaptor {
    public:
        typedef QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree> pTree;
    private:
        pTree shortRateTree_;
    public:
        ShorRateTreeAdaptor(
            const pTree& shortRateTree
        ) : shortRateTree_(shortRateTree) {}
        const pTree& shortRateTree() const {
            return shortRateTree_;
        }
        // returns the short rate at a node
        QuantLib::Rate shortRateAtNode(
            QuantLib::Size i,
            QuantLib::Size index
        ) const {
            auto df = shortRateTree()->discount(i, index);
            auto dt = shortRateTree()->timeGrid().dt(i);
            QuantLib::Rate shortRate = -log(df) / dt;
            return shortRate;
        }
        // returns a random walker for the tree given a random number generator
        template <typename RNG>
        std::shared_ptr<RandomTreeWalker<RNG>> getRandomTreeWalker(
            const RNG& randomNumberGenerator
        ) const {
            return std::shared_ptr<RandomTreeWalker<RNG>>(new RandomTreeWalker<RNG>(shortRateTree(), randomNumberGenerator));
        }
    };
}