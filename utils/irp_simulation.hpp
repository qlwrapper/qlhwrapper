#pragma once

#include <ql/quantlib.hpp>
#include <utils/paths.hpp>
#include <utils/paths_utils.hpp>
#include <utils/short_rate_tree_adaptor.hpp>
#include <utils/backward-induction.hpp>
#include <memory>
#include <vector>
#include <algorithm>

namespace utils {
    class OneFactorShortRateTreeSampling {
    private:
        const QuantLib::OneFactorModel& model_;
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> zvCurve_;
        bool applyDiscountFactorDriftAdj_;
    public:
        struct RandomWalksResult {
            std::shared_ptr<QuantLib::TimeGrid> timeGrid;   // length=nTimes
            QuantLib::ext::shared_ptr<QuantLib::Lattice> lattice;
            QuantLib::ext::shared_ptr<QuantLib::OneFactorModel::ShortRateTree> shortRateTree;
            std::shared_ptr<std::vector<QuantLib::Size>> treeSizes; // tree sizes by time (length=nTimes)
            std::shared_ptr<std::vector<std::vector<QuantLib::Size>>> walkPaths;   // nPaths x nTimes
            std::shared_ptr<Paths> shortRatePaths;   // row x cols = nTimes x nPaths
            std::shared_ptr<Paths> discountFactorPaths;   // row x cols = nTimes x nPaths
            std::shared_ptr<std::vector<std::vector<QuantLib::Size>>> nodeVisitStats;   // length=nTimes
            std::shared_ptr<std::vector<std::vector<QuantLib::Real>>> nodeVisitProbabilities;   // length=nTimes
            std::shared_ptr<std::vector<std::vector<QuantLib::Real>>> nodeVisitCumProbabilities;   // length=nTimes
            std::shared_ptr<std::vector<std::vector<QuantLib::Real>>> statePrices;   // length=nTimes
            std::shared_ptr<std::vector<std::vector<QuantLib::Rate>>> nodeShortRates;   // length=nTimes

            QuantLib::Size nPaths() const {
                return this->walkPaths->size();
            }
            QuantLib::Size nTimes() const {
                return this->timeGrid->size();
            }
            QuantLib::Time maturity() const {
                return this->timeGrid->back();
            }
            QuantLib::Size nTreeNodes() const {
                auto const& treeSizes = *(this->treeSizes);
                return std::accumulate(treeSizes.begin(), treeSizes.end(), 0);
            }
        };

        struct NullProgressObserver {
            void operator() (
                const QuantLib::Size&,
                const QuantLib::Size&
            ) const {}
        };

        struct MonthlyResult {
            std::shared_ptr<QuantLib::TimeGrid> timeGrid;   // length=nTimes
            std::shared_ptr<Paths> shortRatePaths;   // row x cols = nTimes x nPaths
            std::shared_ptr<Paths> discountFactorPaths;   // row x cols = nTimes x nPaths
         };

        OneFactorShortRateTreeSampling(
            const QuantLib::OneFactorModel& model,
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve,
            bool applyDiscountFactorDriftAdj = true
        ) : model_(model), zvCurve_(zvCurve), applyDiscountFactorDriftAdj_(applyDiscountFactorDriftAdj) {}
        const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() const {
            return zvCurve_;
        }
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() {
            return zvCurve_;
        }
        const bool& applyDiscountFactorDriftAdj() const {
            return applyDiscountFactorDriftAdj_;
        }
        bool& applyDiscountFactorDriftAdj() {
            return applyDiscountFactorDriftAdj_;
        }
        static std::shared_ptr<Paths> calculateDiscountFactorPaths(
            const Paths& shortRatePaths
        ) {
            auto nTimes = shortRatePaths.nTimes();
            auto nPaths = shortRatePaths.nPaths();
            std::shared_ptr<Paths> ret(new Paths(nTimes, nPaths, shortRatePaths.timeGrid()));
            auto const& timeGrid = *(shortRatePaths.timeGrid());
            const auto& shortRateMatrix = *(shortRatePaths.matrix());
            auto& discountFactorMatrix = *(ret->matrix());
            for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                auto areaUnderShortRateCurve = 0.0;
                for (decltype(nTimes) i = 0; i < nTimes; ++i) {  // for each row/ gird time
                    auto discountFactor = std::exp(-areaUnderShortRateCurve);
                    discountFactorMatrix[i][path] = discountFactor;
                    if (i < nTimes - 1) {
                        areaUnderShortRateCurve += shortRateMatrix[i][path] * timeGrid.dt(i);
                    }
                }
            }
            return ret;
        }
        // perform random walks to sample short rates and discount factors
        template <
            typename RandomNumberGenerator,
            typename ProgressObserver = NullProgressObserver
        >
        RandomWalksResult randomWalks(
            QuantLib::Time maturity,
            QuantLib::Natural freq,
            QuantLib::Size nPaths,
            const RandomNumberGenerator& rngGenerator,
            const ProgressObserver& progressObserver = ProgressObserver()
        ) const {
            RandomWalksResult ret;
            ret.timeGrid = PathsUtils::getEvenlySpacingTimeGrid(maturity, freq);    // create the time grid
            auto nTimes = ret.timeGrid->size();
            ret.treeSizes.reset(new std::vector<QuantLib::Size>(nTimes));
            auto& treeSizes = *(ret.treeSizes);
            // based on the current short rate sampling function (shortRateAtNode(), specifically the line: auto dt = shortRateTree->timeGrid().dt(i)),
            // the sampling time grid need to have one extra time interval more than the true time grid
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            auto dt = (QuantLib::Time)(1.0 / (QuantLib::Real)freq);
            auto maturitySampling = maturity + dt;  // one extra time interval
            auto nStepsSampling = (QuantLib::Size)std::round(maturitySampling * (QuantLib::Time)freq);
            QuantLib::TimeGrid gridSampling(maturitySampling, nStepsSampling);
            QL_ASSERT(gridSampling.size() == nTimes + 1, "gridSampling.size() != nTimes + 1");
            ////////////////////////////////////////////////////////////////////////////////////////////////////
            auto lattice = model_.tree(gridSampling);   // create the tree lattice
            ret.lattice = lattice;
            auto shortRateTree = QuantLib::ext::dynamic_pointer_cast<QuantLib::OneFactorModel::ShortRateTree>(lattice);
            ret.shortRateTree = shortRateTree;
            ShorRateTreeAdaptor treeAdaptor(shortRateTree);
            auto randomTreeWalker = treeAdaptor.getRandomTreeWalker(rngGenerator);
            for (decltype(nTimes) i = 0; i < nTimes; i++) {    // nTimes iterations
                treeSizes[i] = shortRateTree->size(i);
            }
            ret.walkPaths.reset(new std::vector<std::vector<QuantLib::Size>>(nPaths));
            ret.shortRatePaths.reset(new Paths(nTimes, nPaths, ret.timeGrid));
            // perform random walks to sample short rates
            /////////////////////////////////////////////////////////////////////////////////
            auto& shortRateMatrix = *(ret.shortRatePaths->matrix());
            for (decltype(nPaths) path = 0; path < nPaths; path++) { // for each path
                auto& visitedNodesForPath = ret.walkPaths->at(path);
                visitedNodesForPath.resize(nTimes);
                // there should be nTimes # of short rate sampling calls and nTimes-1 (nSteps) # of randomly picking of the next branch
                QuantLib::Size nodeIndexVisited = 0;
                for (decltype(nTimes) i = 0; i < nTimes; ++i) {  // for each time step (nTimes iterations)
                    visitedNodesForPath[i] = nodeIndexVisited;
                    auto shortRateSampled = treeAdaptor.shortRateAtNode(i, nodeIndexVisited);  // sample the short rate
                    shortRateMatrix[i][path] = shortRateSampled; // put the sampled short rate on the paths matrix
                    if (i < nTimes - 1) {
                        nodeIndexVisited = randomTreeWalker->walkToNextNode(i, nodeIndexVisited);   // nodeIndexVisited is now the visited node index for i+1 time index
                    }
                }
                progressObserver(path + 1, nPaths);
            }
            /////////////////////////////////////////////////////////////////////////////////

            ret.discountFactorPaths = calculateDiscountFactorPaths(*(ret.shortRatePaths));
            if (applyDiscountFactorDriftAdj()) {
                auto& discountFactorMatrix = *(ret.discountFactorPaths->matrix());
                auto stats = ret.discountFactorPaths->statistics();
                auto const& timeGrid = *(ret.discountFactorPaths->timeGrid());
                for (decltype(nTimes) i = 0; i < nTimes; i++) {    // nTimes iterations
                    auto const& t = timeGrid.at(i);
                    auto const& meanDF = stats->at(i).first;
                    auto zvDF = zvCurve()->discount(t);
                    auto adjMultiplier = zvDF/meanDF;
                    for (decltype(nPaths) path = 0; path < nPaths; ++path) { // for each path
                        discountFactorMatrix[i][path] *= adjMultiplier;
                    }
                }
            }

            // calculate random walk node visit statistics
            /////////////////////////////////////////////////////////////////////////////////
            ret.nodeVisitStats.reset(new std::vector<std::vector<QuantLib::Size>>(nTimes));
            auto& nodeVisitStats = *(ret.nodeVisitStats);
            for (decltype(nTimes) i = 0; i < nTimes; i++) {    // nTimes iterations
                nodeVisitStats.at(i).resize(treeSizes.at(i), 0);
            }
            for (decltype(nPaths) path = 0; path < nPaths; path++) { // for each path
                const auto& visitedNodesForPath = ret.walkPaths->at(path);
                for (decltype(nTimes) i = 0; i < nTimes; ++i) {  // for each time (nTimes iterations)
                    auto const& visitedNodeIndex = visitedNodesForPath.at(i);
                    auto& stats = nodeVisitStats[i];
                    QL_ASSERT(visitedNodeIndex >= 0 && visitedNodeIndex < treeSizes.at(i), "bad node index");
                    auto& visitedCount = stats[visitedNodeIndex];
                    visitedCount++;
                }
            }
            for (decltype(nTimes) i = 0; i < nTimes; ++i) {
                auto const& stats = nodeVisitStats.at(i);
                auto totalVisitCount = std::accumulate(stats.begin(), stats.end(), 0);
                QL_ASSERT(totalVisitCount == nPaths, "total visit count (" << totalVisitCount << ") at time step " << i << " is not equal to the number of paths (" << nPaths  << ") for time index " << i);
            }
            /////////////////////////////////////////////////////////////////////////////////

            // get state prices from the tree lattice
            /////////////////////////////////////////////////////////////////////////////////
            ret.statePrices.reset(new std::vector<std::vector<QuantLib::Real>>(nTimes));
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time step
                auto const& src = shortRateTree->statePrices(i);
                auto& dest = ret.statePrices->at(i);
                dest.resize(src.size());
                std::copy(src.begin(), src.end(), dest.begin());
            }
            /////////////////////////////////////////////////////////////////////////////////

            // calculate the probabilities of nodes being visited by random walk
            /////////////////////////////////////////////////////////////////////////////////
            ret.nodeVisitProbabilities.reset(new std::vector<std::vector<QuantLib::Real>>(nTimes));
            auto& nodeVisitProbabilities = *(ret.nodeVisitProbabilities);
            nodeVisitProbabilities[0] = std::vector<QuantLib::Real>(1, 1.0);
            for (decltype(nTimes) i = 1; i < nTimes; ++i) {
                auto prevTimeStepIndex = i - 1;
                auto const& prevTimeStepProbs = nodeVisitProbabilities.at(prevTimeStepIndex);
                auto& currTimeStepProbs = nodeVisitProbabilities[i];
                currTimeStepProbs.resize(treeSizes.at(i), 0.0);
                for (QuantLib::Size n = 0; n < prevTimeStepProbs.size(); n++) {   // for each node of prev time step
                    auto const& probVisitPrevTimeStepNode = prevTimeStepProbs.at(n);
                    for (QuantLib::Size d = 0; d < 3; ++d) {  // for each walk direction (up, flat, down)
                        auto currTimeStepNodeIndex = shortRateTree->descendant(prevTimeStepIndex, n, d);    // node index @ time i
                        auto& probVisitCurrTimeStepNode = currTimeStepProbs[currTimeStepNodeIndex];
                        //cout << "prob (" << prevTimeStepIndex << "," << n << "," << d << ")=" << shortRateTree->probability(prevTimeStepIndex, n, d) << " added to node @(" << i << "," << currNodeIndex << ")" << endl;
                        probVisitCurrTimeStepNode += probVisitPrevTimeStepNode * shortRateTree->probability(prevTimeStepIndex, n, d);
                    }
                }
            }
            for (decltype(nTimes) i = 0; i < nTimes; ++i) {
                auto const& probs = nodeVisitProbabilities.at(i);
                auto cumProb = std::accumulate(probs.begin(), probs.end(), 0.0);
                //QL_ASSERT(QuantLib::close_enough(cumProb, 1.0), "cum. of probability visitng each node (" << cumProb << ") at time step " << i << " is not close to 1.0");
            }
            /////////////////////////////////////////////////////////////////////////////////

            // calculate the cum. probabilities of nodes being visited by random walk
            /////////////////////////////////////////////////////////////////////////////////
            ret.nodeVisitCumProbabilities.reset(new std::vector<std::vector<QuantLib::Real>>(nTimes));
            auto& nodeVisitCumProbabilities = *(ret.nodeVisitCumProbabilities);
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time step
                auto const& probs = nodeVisitProbabilities.at(i);
                auto n = probs.size();
                auto& cumProbs = nodeVisitCumProbabilities[i];
                cumProbs.resize(n);
                decltype(n) j = 0;
                QuantLib::Real cumProb = 0.0;
                for (auto p = probs.begin(); p != probs.end(); ++p, ++j) {  // for each node
                    cumProb += *p;
                    cumProbs[j] = cumProb;
                }
                //QL_ASSERT(QuantLib::close_enough(cumProb, 1.0), "cum. of probability visitng each node (" << cumProb << ") at time step " << i << " is not close to 1.0");
            }
            /////////////////////////////////////////////////////////////////////////////////

            // get all the short rates for the nodes
            /////////////////////////////////////////////////////////////////////////////////
            ret.nodeShortRates.reset(new std::vector<std::vector<QuantLib::Rate>>(nTimes));
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time step
                auto numNodes = treeSizes.at(i);
                auto& rates = ret.nodeShortRates->at(i);
                rates.resize(numNodes);
                for (decltype(numNodes) j = 0; j < numNodes; ++j) {
                    rates.at(j) = treeAdaptor.shortRateAtNode(i, j);
                }
            }
            /////////////////////////////////////////////////////////////////////////////////
            return ret;
        }

        static MonthlyResult toMonthlyResult(
            const RandomWalksResult& randomWalksResult
        ) {
            auto const& srcShortRateMatrix = *(randomWalksResult.shortRatePaths->matrix());
            auto const& srcDiscountFactorMatrix = *(randomWalksResult.discountFactorPaths->matrix());
            auto const& shortRateGrid = *(randomWalksResult.timeGrid);
            MonthlyResult ret;
            ret.timeGrid = PathsUtils::getEvenlySpacingTimeGrid(randomWalksResult.maturity(), (QuantLib::Natural)QuantLib::Monthly);    // create the time grid
            auto const& monthlyTimeGrid = *(ret.timeGrid);
            auto nTimes = monthlyTimeGrid.size();
            auto nPaths = randomWalksResult.nPaths();
            ret.shortRatePaths.reset(new Paths(nTimes, nPaths, ret.timeGrid));
            ret.discountFactorPaths.reset(new Paths(nTimes, nPaths, ret.timeGrid));
            auto& shortRatePaths = *(ret.shortRatePaths->matrix());
            auto& discountFactorPaths = *(ret.discountFactorPaths->matrix());
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time slice on the monthly time grid
                auto t = monthlyTimeGrid.at(i);
                auto srTimeIndex = shortRateGrid.closestIndex(t); // the corresponding time index on the short rate tree
                auto srTime = shortRateGrid.at(srTimeIndex);
                QL_ASSERT(QuantLib::close_enough(t, srTime), "short rate grid and monthly grid are mis-aligned. t=" << t << ", srTime=" << srTime);
                std::copy(srcShortRateMatrix.row_begin(srTimeIndex), srcShortRateMatrix.row_end(srTimeIndex), shortRatePaths.row_begin(i));
                std::copy(srcDiscountFactorMatrix.row_begin(srTimeIndex), srcDiscountFactorMatrix.row_end(srTimeIndex), discountFactorPaths.row_begin(i));
            }
            return ret;
        }
    };

    class SamplingBackwardInductionValues {
    private:
        const OneFactorShortRateTreeSampling::RandomWalksResult& randomWalksResult_;
    public:
        struct Result {
            std::shared_ptr<std::vector<std::shared_ptr<BackwardInductionValues>>> valuesTree; // length=nTimes
            std::shared_ptr<Paths> valuePaths;   // row x cols = nTimes x nPaths
            QuantLib::Size nTimes() const {
                return this->valuePaths->nTimes();
            }
            QuantLib::Size nPaths() const {
                return this->valuePaths->nPaths();
            }
            std::shared_ptr<QuantLib::TimeGrid> timeGrid() const {
                return this->valuePaths->timeGrid();
            }
        };
        SamplingBackwardInductionValues(
            const OneFactorShortRateTreeSampling::RandomWalksResult& randomWalksResult
        ) : randomWalksResult_(randomWalksResult) {}
        const OneFactorShortRateTreeSampling::RandomWalksResult& randomWalksResult() const {
            return randomWalksResult_;
        }

        struct NullProgressObserver {
            void operator() (
                const QuantLib::Size&,
                const QuantLib::Size&,
                const std::shared_ptr<BackwardInductionValues>&
            ) const {}
        };

        template <
            typename BackwardInductionCalculator,
            typename ProgressObserver = NullProgressObserver
        >
        Result sample(
            QuantLib::Size tenorMonths,
            const BackwardInductionCalculator& backwardInductionCalculator,
            const ProgressObserver& progressObserver = ProgressObserver()
        ) const {
            auto freq = QuantLib::Monthly;
            Result ret;
            auto tenor = (QuantLib::Time)((QuantLib::Real)tenorMonths / (QuantLib::Real)freq);
            auto const& randomWalksGrid = *(randomWalksResult().timeGrid);
            auto const& shortRateTreeSizes = *(randomWalksResult().treeSizes);
            auto const& shortRateTree = randomWalksResult().shortRateTree;
            auto randomWalkMaturity = randomWalksResult().maturity();
            QL_REQUIRE(randomWalkMaturity >= tenor, "tenor has to be less or equal to the random walks grid maturity");
            auto maturity = randomWalkMaturity - tenor;  // maturity need to be shorten by the tenor so the discretized values calculation won't go out of bound on the random walk grid
            auto pTimeGrid = PathsUtils::getEvenlySpacingTimeGrid(maturity, freq);
            auto const& monthlyGrid = *(pTimeGrid);
            auto nTimes = monthlyGrid.size();
            // compute the backward induction values for each time slice
            ///////////////////////////////////////////////////////////////////////////////////////
            ret.valuesTree.reset(new std::vector<std::shared_ptr<BackwardInductionValues>>(nTimes));
            auto& valuesTree = *(ret.valuesTree);
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time in the monthly grid
                auto const& t_start = monthlyGrid.at(i);
                auto t_end = t_start + tenor;
                auto indexStart = randomWalksGrid.closestIndex(t_start);  // start time index on the random walks grid
                auto indexEnd = randomWalksGrid.closestIndex(t_end);  // end time index on the random walks grid
                QL_ASSERT(QuantLib::close_enough(randomWalksGrid[indexStart], t_start) && QuantLib::close_enough(randomWalksGrid[indexEnd], t_end), "random walks grid and monthly grid are mis-aligned: randomWalksGrid[indexStart]=" << randomWalksGrid[indexStart] << ",t_start=" << t_start << ",randomWalksGrid[indexEnd]=" << randomWalksGrid[indexEnd] << ",t_end=" << t_end);
                auto t = randomWalksGrid.at(indexStart); // start time on the random walks grid
                auto T = randomWalksGrid.at(indexEnd);   // end time on the random walks grid
                auto const& expectedNumValues = shortRateTreeSizes.at(indexStart); // the values calculator must return this number of backward induction values
                auto values = backwardInductionCalculator(shortRateTree, expectedNumValues, t, T);   // call the backward induction calculator to get the backward induction values
                QL_ASSERT(values->size() == expectedNumValues, "incorrect number of backward induction values calculated");
                valuesTree[i] = values;   // store the calculated backward induction values on the values tree
                progressObserver(i + 1, nTimes, values);    // notify the progress observer
            }
            ///////////////////////////////////////////////////////////////////////////////////////
            auto const& randomWalkPaths = *(randomWalksResult().walkPaths);
            auto nPaths = randomWalksResult().nPaths(); // number of simulated short rate paths
            ret.valuePaths.reset(new Paths(nTimes, nPaths, pTimeGrid));
            auto& valuesMatrix = *(ret.valuePaths->matrix());
            for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                auto const& visitedNodesForPath = randomWalkPaths.at(path);
                for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time in the monthly grid
                    auto const& t = monthlyGrid.at(i);
                    auto rwTimeIndex = randomWalksGrid.closestIndex(t); // get the time index on the random walks grid
                    auto const& visitedNodeIndex = visitedNodesForPath.at(rwTimeIndex);   // node visited for this time index on the random walks grid
                    auto const& values = *(valuesTree.at(i)); // get the backward induction values from the tree storage at this time step
                    auto const& visitedValue = values.at(visitedNodeIndex);   // backward induction value for the visited node
                    valuesMatrix[i][path] = visitedValue;   // put the visit node value on the paths matrix
                }
            }
            return ret;
        }
    };

    class SamplingDriftAdjustZCBPrices {
    private:
        const SamplingBackwardInductionValues& biValuesSampler_;
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> zvCurve_;
        std::shared_ptr<utils::Paths> discountFactorPaths_;
    public:
        enum AdjustmentType {
            FrontEndAdj = 0,    // adjustment at t
            BackEndAdj = 1      // adjustment at T
        };
        struct Result {
            std::shared_ptr<Paths> B;   // un-adjusted P(t, T), tenor = T - t, size = nTimes x n Paths
            std::shared_ptr<Paths> adj_B;   // drift-adjusted P(t, T), tenor = T - t, size = nTimes x n Paths
            std::shared_ptr<std::vector<QuantLib::Real>> adj;   // adjustment factor, size = nTimes
            QuantLib::Size nTimes() const {
                return this->B->nTimes();
            }
            QuantLib::Size nPaths() const {
                return this->B->nPaths();
            }
            std::shared_ptr<QuantLib::TimeGrid> timeGrid() const {
                return this->B->timeGrid();
            }
        };
        SamplingDriftAdjustZCBPrices(
            const SamplingBackwardInductionValues& biValuesSampler,
            const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve,
            const std::shared_ptr<utils::Paths>& discountFactorPaths
        ) :
            biValuesSampler_(biValuesSampler),
            zvCurve_(zvCurve),
            discountFactorPaths_(discountFactorPaths)
        {}
        const SamplingBackwardInductionValues& biValuesSampler() const {
            return biValuesSampler_;
        }
        const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() const {
            return zvCurve_;
        }
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve() {
            return zvCurve_;
        }
        const std::shared_ptr<utils::Paths>& discountFactorPaths() const {
            return discountFactorPaths_;
        }
        std::shared_ptr<utils::Paths>& discountFactorPaths() {
            return discountFactorPaths_;
        }
        template <
            typename ProgressObserver = SamplingBackwardInductionValues::NullProgressObserver
        >
        Result calculate(
            QuantLib::Size tenorMonths,
            const BackwardInductionZCBCalculator& biZCBCalculator,
            AdjustmentType adjustmentType,
            const ProgressObserver& progressObserver = ProgressObserver()
        ) const {
            auto result = biValuesSampler().sample(
                tenorMonths,
                biZCBCalculator,
                progressObserver
            );
            auto tenorTime = (QuantLib::Time)tenorMonths / 12.0;
            Result ret;
            ret.B = result.valuePaths;
            auto const& timeGrid = *(result.timeGrid());
            auto nTimes = result.nTimes();
            auto nPaths = result.nPaths();
            auto const& B = *(result.valuePaths->matrix());
            ret.adj_B.reset(new Paths(nTimes, nPaths, result.timeGrid()));
            auto& adj_B = *(ret.adj_B->matrix());
            ret.adj.reset(new std::vector<QuantLib::Real>(nTimes));
            auto& adj = *(ret.adj);
            auto const& dfMatrix = *(discountFactorPaths()->matrix());
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time slice
                auto t = timeGrid.at(i);
                auto T = t + tenorTime;
                auto adj_t_T = 1.0;
                if (adjustmentType == AdjustmentType::FrontEndAdj) {
                    auto j = discountFactorPaths()->timeGrid()->closestIndex(T);
                    auto x_t = 0.0; // = Avg[P(0, T)/B(t, T)] = mean of P(0, t)
                    for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                        x_t += (dfMatrix[j][path] / B[i][path]) / nPaths;
                    }
                    auto adj_t_T = zvCurve()->discount(t) / x_t; // ratio of ZV P(0, t) to mean of P(0, t)
                    for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                        adj_B[i][path] = B[i][path] / adj_t_T;
                    }
                }
                else if (adjustmentType == AdjustmentType::BackEndAdj) {
                    auto j = discountFactorPaths()->timeGrid()->closestIndex(t);
                    auto x_T = 0.0; // = Avg[P(0, t)*B(t, T)] = mean of P(0, T)
                    for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                        x_T += (dfMatrix[j][path] * B[i][path]) / nPaths;
                    }
                    auto adj_t_T = zvCurve()->discount(T) / x_T; // ratio of ZV P(0, T) to mean of P(0, T)
                    for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                        adj_B[i][path] = B[i][path] * adj_t_T;
                    }
                }
                else {
                    QL_FAIL("unknown/unsupported adjustment type");
                }
                adj.at(i) = adj_t_T;
            }
            return ret;
        }
    };
}