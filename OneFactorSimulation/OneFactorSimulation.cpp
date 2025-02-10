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

#define BOOST_LIB_DIAGNOSTIC
#include <iostream>
#include <boost/config.hpp>
#include <ql/quantlib.hpp>
#include <exception>
#include <utils/utilities.h>
#include <map>
#include <chrono>
#include <output.hpp>

#include <boost/program_options.hpp>
// auto link with the correct boost libraries
#ifdef BOOST_MSVC
#define BOOST_LIB_NAME boost_program_options
#include <boost/config/auto_link.hpp>
#undef BOOST_LIB_NAME
#endif

using namespace std;
using namespace QuantLib;
using namespace utils;

static std::shared_ptr<ShortRateModelParams> instantiateModelParams(
    ModelType modelType,
    const std::string& paramsJSON
) {
    std::shared_ptr<ShortRateModelParams> ret;
    switch (modelType) {
    case ModelType::HullWhite1F:
        ret.reset(new HullWhiteModelParams());
        break;
    case ModelType::GeneralizedHullWhite1F:
        ret.reset(new GeneralizedHullWhiteModelParams<NormalVolTraits<QuantLib::Real>>());
        break;
    case ModelType::GeneralizedHullWhiteLogNormal1F:
        ret.reset(new GeneralizedHullWhiteModelParams<LogNormalVolTraits<QuantLib::Real>>());
        break;
    default:
        QL_FAIL("unknown model type");
    }
    ret->JSON_parse(paramsJSON);
    return ret;
}

static std::pair<std::shared_ptr<utils::Paths>, std::shared_ptr<utils::Paths>> getSimpleRatePaths(
    Natural tenorMonths,
    const DayCounter& cashRateDayCounter,
    const utils::Paths& adj_B,
    const Date& evalDate = Date()
) {
    auto referenceDate = (evalDate == Date() ? Settings::instance().evaluationDate() : evalDate);
    auto tenorTime = (Time)tenorMonths / 12.0;
    auto nTimes = adj_B.nTimes();
    auto nPaths = adj_B.nPaths();
    auto const& timeGrid = *(adj_B.timeGrid());
    std::shared_ptr<utils::Paths> pSimpleFwdPaths(new utils::Paths(nTimes, nPaths, adj_B.timeGrid()));
    std::shared_ptr<utils::Paths> pCashRateFwdPaths(new utils::Paths(nTimes, nPaths, adj_B.timeGrid()));
    auto const& B_adj = *(adj_B.matrix());
    auto& simpleFwdMatrix = *(pSimpleFwdPaths->matrix());
    auto& cashRateFwdMatrix = *(pCashRateFwdPaths->matrix());
    for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each month
        auto t = timeGrid.at(i);
        auto monthStart = (Natural)std::round(t * 12.0);
        auto monthEnd = monthStart + tenorMonths;
        auto dtStart = referenceDate + monthStart * Months;
        auto dtEnd = referenceDate + monthEnd * Months;
        auto tenorCashRate = cashRateDayCounter.yearFraction(dtStart, dtEnd);
        for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
            auto compounding = 1.0 / B_adj[i][path];
            simpleFwdMatrix[i][path] = (compounding - 1.0) / tenorTime;
            cashRateFwdMatrix[i][path] = (compounding - 1.0) / tenorCashRate;
        }
    }
    std::pair<std::shared_ptr<utils::Paths>, std::shared_ptr<utils::Paths>> ret(pSimpleFwdPaths, pCashRateFwdPaths);
    return ret;
}

std::shared_ptr<Paths> calculateParRatePaths(
    std::ostream& os,
    Natural tenor,
    Frequency couponFrequency,
    double dayCountMultiplier,
    const ext::shared_ptr<YieldTermStructure>& zvCurve,
    Size nPaths,
    const SamplingDriftAdjustZCBPrices& driftAdjSampling,
    const BackwardInductionZCBCalculator& biZCBCalculator,
    ZVComparison<char>& zvComparison
) {
    os << endl;
    auto freq = (Natural)couponFrequency;
    Natural couponTenorMonths = 12 / freq;  // number of months between coupon payments
    os << tenor << " yr Par Rate Calculation" << endl;
    os << "calculating drift Adj B matrices..." << endl;
    Natural numCoupons = tenor * freq;
    vector<std::shared_ptr<Matrix>> driftAdjBs(numCoupons);  // numCoupons drift-adj B matrices, one for each coupon => P(t, t+0.5), P(t, t+1.0),... P(t, t+10.0)
    std::shared_ptr<TimeGrid> parTimeGrid;
    for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
        auto driftAdjResult = driftAdjSampling.calculate(
            (k + 1) * couponTenorMonths,
            biZCBCalculator,
            SamplingDriftAdjustZCBPrices::BackEndAdj
        );
        driftAdjBs[k] = driftAdjResult.adj_B->matrix();
        parTimeGrid = driftAdjResult.adj_B->timeGrid();
        os << "coupon " << k + 1 << "/" << numCoupons;
        auto now = std::chrono::system_clock::now();
        os << ", secSinceEpoch=" << std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
        os << std::endl;
    }
    os << "done calculating drift Adj B matrices. calculating par rate paths..." << endl;
    auto nTimes = parTimeGrid->size();
    ParRateCalculator parRateClaculator(couponFrequency);
    std::shared_ptr<Paths> parRatePaths(new Paths(nTimes, nPaths, parTimeGrid));
    auto& parRateMatrix = *(parRatePaths->matrix());
    for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time slice
        for (decltype(nPaths) path = 0; path < nPaths; ++path) { // for each path
            vector<DiscountFactor> dfs(numCoupons);
            for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                auto const& adj_B = *(driftAdjBs.at(k));
                dfs[k] = adj_B[i][path];
            }
            auto parRate = parRateClaculator(dfs.begin(), dfs.end());
            parRate *= dayCountMultiplier;
            parRateMatrix[i][path] = parRate;
        }
    }
    os << tenor << " yr par rate paths calculation completed" << endl;
    os << endl;
    os << tenor << " yr Par Rate ZV Comparison Report" << endl;
    ParRateZVCalculator zvCalculator((Time)tenor, couponFrequency, zvCurve, dayCountMultiplier);
    zvComparison.report<>(zvCalculator, *parRatePaths);
    return parRatePaths;
}

#define DEFAULT_RANDOM_NUMBER_GENERATOR_SEED (28749)
#define DEFAULT_N_PATHS (500)
#define DEFAULT_NO_BACKWARD_INDUCTION_VALUES (false)

static string SHORT_RATE_FILENAME_SUFFIX = "short_rate";
static string DISCOUNT_FACTOR_FILENAME_SUFFIX = "discount_factor";
static string FWD_CASH_RATE_FILENAME_SUFFIX = "fwd_cash_rate";

static string CMS_BREAKEVEN_RATE_TREE_FILENAME_SUFFIX = "cms_breakeven_rate_tree";
static string CMS_BREAKEVEN_RATE_SIMULATION_FILENAME_SUFFIX = "cms_breakeven_rate_simulation";

namespace po = boost::program_options;

int main(int argc, char* argv[], char* envp[]) {
    try {
        string asofdate;
        string s_zeroCurve;
        string flatForwardRate;
        string model = DEFAULT_MODEL;
        string params;
        unsigned long rngSeed = DEFAULT_RANDOM_NUMBER_GENERATOR_SEED;
        Real maturity = DEFAULT_MATURITY;
        Natural gridFreq = DEFAULT_GRID_FREQ;
        Size nPaths = DEFAULT_N_PATHS;
        bool isSofr = DEFAULT_IS_SOFR;
        bool doNotApplyDiscountFactorDriftAdj = false;
        bool noBackwardInductionValues = DEFAULT_NO_BACKWARD_INDUCTION_VALUES;
        string outputFolder;
        bool silent = DEFAULT_SILENT;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("asofdate,d", po::value<string>(&asofdate)->default_value(""), "(optional) evaluation date in yyyymmdd format. default to today's date")
            ("zero-curve,z", po::value<string>(&s_zeroCurve)->default_value(""), "(required) zero/spot curve input json file. \"@-\" to piped in from stdin")
            ("flat-forward", po::value<string>(&flatForwardRate)->default_value(""), "(optional) flat forward rate in % unit")
            ("model,m", po::value<string>(&model)->default_value(DEFAULT_MODEL), "(optional) model type (hw|ghw|ghwln). default value is hw")
            ("params", po::value<string>(&params)->default_value(""), "(optional) model parameters json. formats: JSON string literal, \"@{filePath}\", or \"@-\" to piped in from stdin")
            ("rng-seed", po::value<unsigned long>(&rngSeed)->default_value(DEFAULT_RANDOM_NUMBER_GENERATOR_SEED), "(optional) seed for the random number generator. the default value is 28749")
            ("maturity", po::value<Real>(&maturity)->default_value(DEFAULT_MATURITY), "(optional) simulation duration in years. The default value is 50")
            ("grid-freq,f", po::value<Natural>(&gridFreq)->default_value(DEFAULT_GRID_FREQ), "grid frequency per year: 12=monthly|24=bi-monthly|36=10d|60=6d|72=5d|120=3d|180=2d|360=1d. The default value is 12 (monthly)")
            ("paths,p", po::value<Size>(&nPaths)->default_value(DEFAULT_N_PATHS), "(optional) number of interest rate paths to be generated. The default value is 500")
            ("is-sofr", po::value<bool>(&isSofr)->default_value(DEFAULT_IS_SOFR), "(optional) the input zero curve is a SOFR zero curve. true or false. The default value is true")
            ("irp-output-dir,o", po::value<string>(&outputFolder)->default_value(""), "(optional) interest rate paths output folder path. The default is the current working directory")
            ("no-adj-df-drift", po::value<bool>(&doNotApplyDiscountFactorDriftAdj)->default_value(false), "(optional) do not apply adjustment to discount factor paths to match zv discount factor. true or false. The default value is false")
            ("no-bi-vals", po::value<bool>(&noBackwardInductionValues)->default_value(DEFAULT_NO_BACKWARD_INDUCTION_VALUES), "(optional) do not generate and output backward induction values. true or false. The default value is false")
            ("silent,s", po::value<bool>(&silent)->default_value(DEFAULT_SILENT), "(optional) silent running mode. true or false. The default is false")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        QL_REQUIRE(!s_zeroCurve.empty(), "zero/spot curve input json file is not optional");
        if (model.empty()) {
            model = DEFAULT_MODEL;
        }
        QL_REQUIRE(!params.empty(), "model params cannot be empty");
        QL_REQUIRE(maturity > 0., "maturity (" << maturity << ") must be a positive");
        QL_REQUIRE(gridFreq > 0, "grid frequency (" << gridFreq << ") must be a positive integer");
        QL_REQUIRE(nPaths > 0, "number of interest rate paths (" << nPaths << ") must be a positive integer");
        auto evalDate = (asofdate.empty() ? Date::todaysDate() : DateFormat<char>::from_yyyymmdd(asofdate, false));
        auto modelType = boost::lexical_cast<ModelType>(model);
        auto paramsJSON = get_input_content(params);

        auto& os = get_nullable_cout<char>(silent);

        os << endl;
        os << "One-Factor Model Simulation" << endl;
        os << "asofdate=" << DateFormat<char>::to_yyyymmdd(evalDate, true) << endl;
        os << "zero-curve=" << s_zeroCurve << endl;
        os << "flat-forward=" << (flatForwardRate.empty() ? "(not set - zv curve will be used)" : flatForwardRate + "%") << endl;
        os << "model=" << boost::lexical_cast<string>(modelType) << endl;
        os << "params=" << paramsJSON << endl;
        os << "rng-seed=" << rngSeed << endl;
        os << "maturity=" << maturity << endl;
        os << "grid-freq=" << gridFreq << endl;
        os << "paths=" << nPaths << endl;
        os << "is-sofr=" << (isSofr ? "true" : "false") << endl;
        os << "irp-output-dir=" << (!outputFolder.empty() ? outputFolder : "(current working directory)") << endl;
        os << "no-adj-df-drift=" << (doNotApplyDiscountFactorDriftAdj ? "true" : "false") << endl;
        os << "no-bi-vals=" << (noBackwardInductionValues ? "true" : "false") << endl;
        os << "silent=" << (silent ? "true" : "false") << endl;
        os << endl;

        Settings::instance().evaluationDate() = evalDate;

        auto pParams = instantiateModelParams(modelType, paramsJSON);
        os << endl;
        os << boost::lexical_cast<string>(modelType) << " model params" << endl;
        os << *pParams;
        const auto& curveConfig = getCurveConfigs().at(isSofr);
        const auto& parRateCouponFrequency = curveConfig.parRateCouponFrequency;
        const auto& parRateDayCountMultiplier = curveConfig.parRateDayCountMultiplier;
        auto pCashRateIndex = curveConfig.createCashRateIndex();
        auto cashRateIndexTenor = pCashRateIndex->tenor();
        auto cashRateIndexTenorUnits = cashRateIndexTenor.units();
        QL_REQUIRE(cashRateIndexTenorUnits == Years || cashRateIndexTenorUnits == Months, "cash rate index's time unit must be either Months or Years");
        auto cashRateTenorMonths = (cashRateIndexTenorUnits == Years ? cashRateIndexTenor.length() * 12 : cashRateIndexTenor.length());
        auto cashRateDayCounter = pCashRateIndex->dayCounter();
        auto yieldCurve = (flatForwardRate.empty()
            ? MarketRate::readAsThirty360ZeroCurve(s_zeroCurve, evalDate, true)
            : ext::shared_ptr<YieldTermStructure>(new FlatForward(evalDate, boost::lexical_cast<Rate>(flatForwardRate)/100.0, Thirty360(Thirty360::BondBasis))));
        Handle<YieldTermStructure> zvCurve(yieldCurve);
        auto shortRateModel = pParams->createModel(zvCurve);
        auto oneFactorModel = ext::dynamic_pointer_cast<OneFactorModel>(shortRateModel); // the short rate model must be a one factor model
        QL_ASSERT(oneFactorModel != nullptr, "the short rate model is not a one factor model");
        auto oneFactorAffineModel = ext::dynamic_pointer_cast<OneFactorAffineModel>(shortRateModel); // the short rate model must also be a one factor affine model
        QL_ASSERT(oneFactorModel != nullptr, "the short rate model is not a one factor affine model");
        std::shared_ptr<OneFactorAffineLike> oneFactorAffineLike(new OneFactorAffineLikeHack(modelType, oneFactorAffineModel));
        OneFactorShortRateTreeSampling shortRateSampling(*oneFactorModel, yieldCurve, !doNotApplyDiscountFactorDriftAdj);
        // create mersenne twister uniform random generator
        MersenneTwisterUniformRng rngGenerator(rngSeed);
        // start the random walks
        os << endl;
        os << "Random Walking..." << endl;
        auto randomWalksResult = shortRateSampling.randomWalks(
            maturity,
            gridFreq,
            nPaths,
            rngGenerator,
            [&os](Size completed, Size total) -> void {
                os << "path: " << completed << "/" << total << endl;
            }
        );
        os << "Done with random walks" << endl;

        os << endl;
        os << "nTimes=" << randomWalksResult.nTimes() << endl;
        os << "nTreeNodes=" << randomWalksResult.nTreeNodes() << endl;

        std::vector<Natural> importantMaturities{ 0, 1, 2, 3, 5, 7, 10, 12, 15, 20, 25, 30, 40 };
        ZVComparison<char> zvComparison(os, yieldCurve, importantMaturities);

        os << endl;
        os << "Short Rate ZV Comparison Report" << endl;
        zvComparison.report<ShortRateZVCalculator>(*(randomWalksResult.shortRatePaths));

        os << endl;
        os << "Discount Factor ZV Comparison Report" << endl;
        zvComparison.report<ZeroRateZVCalculator<Continuous>, DiscountFactorToZeroRateConverter<Continuous>>(*(randomWalksResult.discountFactorPaths));

        if (modelType == ModelType::HullWhite1F) {
            auto pHWModelParams = std::dynamic_pointer_cast<HullWhiteModelParams>(pParams);
            auto const& a = pHWModelParams->a();
            auto const& sigma = pHWModelParams->sigma();
            auto hwVariance = [&a, &sigma](Time t) -> Real {
                return sigma * sigma / (2.0 * a) * (1.0 - std::exp(-2.0 * a * t));
            };
            auto hwAlpha = [&a, &sigma](Rate forward, Time t) -> Real {
                return forward + 0.5 * std::pow(sigma / a * (1.0 - exp(-a * t)), 2.0);
            };
            auto hwStdev = [&hwVariance](Time t) -> Real {
                return std::sqrt(hwVariance(t));
            };
            auto hwMean = hwAlpha;
            os << std::endl;
            os << modelType << " Short Rate Random Walk Actual vs Theoretical Comparison Report" << std::endl;
            auto timeGrid = randomWalksResult.shortRatePaths->timeGrid();
            auto statistics = randomWalksResult.shortRatePaths->statistics();
            for (const auto& year : importantMaturities) {
                auto t = (Time)year;
                auto i = timeGrid->closestIndex(t);
                auto const& stats = statistics->at(i);
                auto forward = zvCurve->forwardRate(t, t, Continuous, NoFrequency);
                os << t;
                os << "," << "forward=" << forward * 100.0;
                os << "," << "mean=" << stats.first * 100.0;
                os << "," << "theoMean=" << hwMean(forward, t) * 100.0;
                os << "," << "stdev=" << stats.second * 100.0;
                os << "," << "theoStdev=" << hwStdev(t) * 100.0;
                os << std::endl;
            }
        }

        auto monthlyResult = OneFactorShortRateTreeSampling::toMonthlyResult(randomWalksResult);
        auto const& monthlyTimeGrid = *(monthlyResult.timeGrid);

        auto get_Model_OutputFilePath = [&outputFolder, &model](const string& fnSuffix, const string& extension) -> string {
            std::ostringstream os;
            os << model << "_" << fnSuffix << extension;
            auto filename = os.str();
            return (outputFolder.length() > 0 ? outputFolder + "\\" + filename : filename);
        };
        auto get_OutputFilePath = [&outputFolder](const std::string& fnSuffix, const string& extension) -> std::string {
            std::ostringstream oss;
            oss << fnSuffix << extension;
            auto filename = oss.str();
            return (outputFolder.length() > 0 ? outputFolder + "\\" + filename : filename);
        };
        OneFactorSimulationOutput output(
            os,
            get_Model_OutputFilePath,
            get_OutputFilePath
        );

        std::shared_ptr<Paths> pDFPaths = monthlyResult.discountFactorPaths;
        std::shared_ptr<Paths> pShortRatePaths = monthlyResult.shortRatePaths;
        
        output.dumpPaths(
            *(monthlyResult.shortRatePaths),
            OutputTraits("short rate", SHORT_RATE_FILENAME_SUFFIX),
            true
        );
        output.dumpPaths(
            *(monthlyResult.discountFactorPaths),
            OutputTraits("discount factor", DISCOUNT_FACTOR_FILENAME_SUFFIX),
            true
        );

        bool generateBackwardInductionValues = !noBackwardInductionValues;

        if (generateBackwardInductionValues) {
            bool doDiscretizedDiscountBondCalc = true;
            std::shared_ptr<BackwardInductionDiscountBondCalculator> biDiscountBondCalcator;    // calculator that calculat P(t, T)
            if (doDiscretizedDiscountBondCalc) {
                biDiscountBondCalcator.reset(new DiscretizedDiscountBondCalculator());
            }
            else {
                biDiscountBondCalcator.reset(new OneFactorAffineDiscountBondCalculator(oneFactorAffineLike));
            }

            BackwardInductionZCBCalculator biZCBCalculator(biDiscountBondCalcator); // backward induction calculator that calculat P(t, T)
            //BackwardInductionFwdRateCalculator biFwdRateCalculator(biDiscountBondCalcator, QuantLib::Continuous);
            BackwardInductionParRateCalculator biParRateCalculator(biDiscountBondCalcator, QuantLib::Semiannual); // backward induction calculator that calculate par rate

            auto backwardInductionProgressMonitor = [&os](
                const Size& completed,
                const Size& total,
                const std::shared_ptr<BackwardInductionValues>& values
                ) -> void {
                auto now = std::chrono::system_clock::now();
                os << "month " << completed;
                os << "/" << total;
                os << ", numNodes=" << values->size();
                os << ", secSinceEpoch=" << std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
                os << endl;
            };

            SamplingBackwardInductionValues backwardInductionValuesSampling(randomWalksResult);
            /*
            {
                QuantLib::Time timeSlice = 5.;
                auto timeIndex = randomWalksResult.timeGrid->closestIndex(timeSlice);
                auto result = backwardInductionValuesSampling.sample(120, biParRateCalculator);
                auto pCalculateValues = result.valuesTree->at(timeIndex);
                auto probabilityVector = randomWalksResult.nodeVisitProbabilities->at(timeIndex);
                QL_ASSERT(pCalculateValues != nullptr, "");
                const auto& valueVector = *pCalculateValues;
                QL_ASSERT(probabilityVector.size() == valueVector.size(), "");
                auto n_nodes = probabilityVector.size();
                auto weightedAvg = 0.;
                for (decltype(n_nodes) i = 0; i < n_nodes; ++i) {
                    const auto& valueAtTreeNode = valueVector[i];
                    const auto& probabilityAtTreeNode = probabilityVector[i];
                    weightedAvg += valueAtTreeNode * probabilityAtTreeNode;
                }
                cout << fixed << "weighted-average backward-induction 10yr par rate @" << timeSlice << "yr = " << weightedAvg * 100.0 << "%" << endl;
            }
            {
                QuantLib::Time timeSlice = 10.;
                auto timeIndex = randomWalksResult.timeGrid->closestIndex(timeSlice);
                auto result = backwardInductionValuesSampling.sample(120, biParRateCalculator);
                auto pCalculateValues = result.valuesTree->at(timeIndex);
                auto probabilityVector = randomWalksResult.nodeVisitProbabilities->at(timeIndex);
                QL_ASSERT(pCalculateValues != nullptr, "");
                const auto& valueVector = *pCalculateValues;
                QL_ASSERT(probabilityVector.size() == valueVector.size(), "");
                auto n_nodes = probabilityVector.size();
                auto weightedAvg = 0.;
                for (decltype(n_nodes) i = 0; i < n_nodes; ++i) {
                    const auto& valueAtTreeNode = valueVector[i];
                    const auto& probabilityAtTreeNode = probabilityVector[i];
                    weightedAvg += valueAtTreeNode * probabilityAtTreeNode;
                }
                cout << fixed << "weighted-average backward-induction 10yr par rate @" << timeSlice << "yr = " << weightedAvg * 100.0 << "%" << endl;
            }
            */

            SamplingDriftAdjustZCBPrices driftAdjSampling(backwardInductionValuesSampling, yieldCurve, monthlyResult.discountFactorPaths);
            // get simple-copounded forward cash rates
            std::shared_ptr<Paths> pCashRatePaths;
            {
                os << endl;
                auto tenorMonths = cashRateTenorMonths;
                auto driftAdjResult = driftAdjSampling.calculate(
                    tenorMonths,
                    biZCBCalculator,
                    SamplingDriftAdjustZCBPrices::FrontEndAdj,
                    backwardInductionProgressMonitor
                );
                auto pr = getSimpleRatePaths(tenorMonths, cashRateDayCounter, *(driftAdjResult.adj_B), evalDate);
                pCashRatePaths = pr.second;
                output.dumpPaths(
                    *(pr.second),
                    OutputTraits("fwd cash rate", FWD_CASH_RATE_FILENAME_SUFFIX),
                    false
                );
            }

            vector<Natural> parTenors{ 5, 10 };
            map<Natural, std::shared_ptr<utils::Paths>> parRatePathsByTenor;
            for (const auto& parTenor : parTenors) {    // for each par tenor
                // get the backend adjusted par rates
                auto pParRatePaths = calculateParRatePaths(
                    os,
                    parTenor,
                    parRateCouponFrequency,
                    parRateDayCountMultiplier,
                    yieldCurve,
                    nPaths,
                    driftAdjSampling,
                    biZCBCalculator,
                    zvComparison
                );
                ostringstream oss;
                oss << parTenor << " yr par rate";
                auto type = oss.str();
                oss.str("");
                oss.clear();
                oss << "swap_rate_" << parTenor << "y";
                auto fnSuffix = oss.str();
                output.dumpPaths(
                    *pParRatePaths,
                    OutputTraits(type, fnSuffix),
                    false
                );
                parRatePathsByTenor[parTenor] = pParRatePaths;
            }

            {
                os << endl;
                IRPSimulationOutput simulationOutput;
                simulationOutput.pDFs = pDFPaths;
                simulationOutput.pShortRates = pShortRatePaths;
                simulationOutput.pBenchmarkCashRates = pCashRatePaths;
                simulationOutput.pBenchmark5yrParRates = parRatePathsByTenor.at(5);
                simulationOutput.pBenchmark10yrParRates = parRatePathsByTenor.at(10);
                output.dumpIRPSimulationOutput(
                    simulationOutput,
                    OutputTraits("interest rate path simulation json output", "irp_simulation_output", ".json")
                );
            }

            /*
            // CMS breakeven rate calculations
            /////////////////////////////////////////////////////////////////////////////////////////////
            QuantLib::Frequency cmsFloatingLegFreq = QuantLib::Quarterly;   // one CMS par coupon payment @ 3mo time
            // backward induction calculator that calculate the tree CMS breakeven rates
            BackwardInductionCMSBreakevenCalculator biCMSBreakevenCalculator(
                biDiscountBondCalcator,
                PAR_COUPON_FREQ,
                cmsFloatingLegFreq
            );
            os << endl;
            os << "Tree CMS Breakeven Rate Calculation" << endl;
            auto cmsBreakevenTreeResult = backwardInductionValuesSampling.sample(
                PAR_TENOR * 12,
                biCMSBreakevenCalculator,
                backwardInductionProgressMonitor
            );  // calculate the tree breakevent CMS rates. there should be one tree CMS breakevent rate per time slice
            auto const& cmdBreakevenTreeValues = *(cmsBreakevenTreeResult.valuesTree);

            os << endl;
            os << "Tree vs Simulation CMS Breakeven Rate Comparison" << endl;
            os << "t";
            os << "," << "Tree Rate";
            os << "," << "Simulation Rate";
            os << "," << "Diff (bp)";
            os << endl;
            os << std::fixed << std::setprecision(16);
            std::vector<QuantLib::Real> cmsBreakevenRatesTree(nTimes);
            std::vector<QuantLib::Real> cmsBreakevenRatesSimulation(nTimes);
            auto discountFactorPaths = randomWalksResult.discountFactorPaths;
            auto const& discountFactorTimeGrid = *(discountFactorPaths->timeGrid());
            auto const& discountFactorMatrix = *(discountFactorPaths->matrix());
            QuantLib::Time dt = 1.0 / (QuantLib::Real)cmsFloatingLegFreq;   // first coupon time
            for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time slice
                auto t = parTimeGrid->at(i);
                auto T = t + dt;
                auto j = discountFactorTimeGrid.closestIndex(T);
                auto avgPV = 0.0;   // average PV @ t=0 of the single par rate coupon paid @ T
                for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
                    auto const& parRate = parRateMatrix[i][path];   // par rate is the coupon rate
                    auto const& df = discountFactorMatrix[j][path]; // discount factor at T (P(0, T))
                    auto pathPV = parRate * dt * df;    // PV of the coupon casflow for this path
                    avgPV += pathPV / nPaths;
                }
                auto zvDF_T = zvCurve->discount(T);
                auto simulationBreakEventCMSRate = avgPV / dt / zvDF_T;   // avgPV = simulationBreakEventCMSRate * dt * zvDF_T
                auto const& treeBreakEventCMSRate = cmdBreakevenTreeValues.at(i)->at(0);
                cmsBreakevenRatesTree[i] = treeBreakEventCMSRate;
                cmsBreakevenRatesSimulation[i] = simulationBreakEventCMSRate;
                os << t;
                os << "," << treeBreakEventCMSRate * 100.0;
                os << "," << simulationBreakEventCMSRate * 100.0;
                os << "," << (simulationBreakEventCMSRate - treeBreakEventCMSRate) * 10000.0;
                os << endl;
            }
            output.dumpTreeCMSBreakevenRates(
                cmsBreakevenRatesTree,
                *parTimeGrid,
                OutputTraits("tree CMS breakeven rate", CMS_BREAKEVEN_RATE_TREE_FILENAME_SUFFIX),
                true,
                100.0
            );
            output.dumpSimulationCMSBreakevenRates(
                cmsBreakevenRatesSimulation,
                *parTimeGrid,
                OutputTraits("simulation CMS breakeven rate", CMS_BREAKEVEN_RATE_SIMULATION_FILENAME_SUFFIX),
                true,
                100.0
            );
            /////////////////////////////////////////////////////////////////////////////////////////////
            */
        }

        os << endl;
        os << "Done" << endl;

        return 0;
    }
    catch (const std::exception& e) {
        cerr << "!!! Error: " << e.what() << endl;
        return 1;
    }
}

/*
{
    auto pGHWParams = std::dynamic_pointer_cast<GeneralizedHullWhiteModelParams<NormalVolTraits<QuantLib::Real>>>(pParams);
    QuantLib::ext::shared_ptr<GeneralizedHullWhiteTieUp<>> pNGHWWCMR(new GeneralizedHullWhiteTieUp<>(*pGHWParams, zvCurve, randomWalksResult.shortRateTree));
    typedef OneFactorAffineDiscountBondCalculator<GeneralizedHullWhiteTieUp<>> OneFactorAffineCalculator;
    std::shared_ptr<OneFactorAffineCalculator> oneFactorAffineDiscountBondCalc(new OneFactorAffineCalculator(pNGHWWCMR));
            
    BackwardInductionZCBCalculator oneFactorAffineZCBCalc(oneFactorAffineDiscountBondCalc);

    {
        Natural tenorMonths = 120;
        Time timeSlice = 10.0;

        os << std::fixed << std::setprecision(12);
        auto tenor = (QuantLib::Time)((QuantLib::Real)tenorMonths / 12.0);
        auto const& randomWalksGrid = *(randomWalksResult.timeGrid);
        auto const& shortRateTreeSizes = *(randomWalksResult.treeSizes);
        auto const& shortRateTree = randomWalksResult.shortRateTree;
        auto randomWalkMaturity = randomWalksResult.maturity();
        QL_REQUIRE(randomWalkMaturity >= tenor, "tenor has to be less or equal to the random walks grid maturity");
        auto indexStart = randomWalksGrid.closestIndex(timeSlice);
        auto indexEnd = randomWalksGrid.closestIndex(timeSlice + tenor);
        auto t = randomWalksGrid.at(indexStart); // start time on the random walks grid
        auto T = randomWalksGrid.at(indexEnd);   // end time on the random walks grid
        auto numNodes = shortRateTreeSizes.at(indexStart);
        ShorRateTreeAdaptor treeAdaptor(shortRateTree);
        GeneralizedHullWhiteTieUp<> nghwcmr(*pGHWParams, zvCurve, randomWalksResult.shortRateTree);
        auto pGHWModel = QuantLib::ext::dynamic_pointer_cast<QuantLib::GeneralizedHullWhite>(oneFactorModel);
        auto const& ghwModel = *(pGHWModel);
        auto const& ghwModelHack = reinterpret_cast<const GeneralizedHullWhiteHack&>(ghwModel);
        os << endl;
        for (decltype(numNodes) index = 0; index < numNodes; ++index) { // for each node
            auto shortRate = treeAdaptor.shortRateAtNode(indexStart, index);
            auto PtT = nghwcmr.discountBond(t, T, shortRate);
            auto AtT = nghwcmr.A(t, T);
            auto BtT = nghwcmr.B(t, T);
            auto mu_t = nghwcmr.mu(t);
            auto vtT = nghwcmr.v(t, T);
            auto DFtT = nghwcmr.forwardDF(t, T);
            auto ql_PtT = ghwModel.discountBond(t, T, shortRate);
            auto ql_AtT = ghwModelHack.A_(t, T);
            auto ql_BtT = ghwModelHack.B_(t, T);
            os << "nodeIndex=" << index;
            os << "," << "shortRate=" << shortRate*100.0;
            os << "," << "PtT=" << PtT;
            os << "," << "AtT=" << AtT;
            os << "," << "BtT=" << BtT;
            os << "," << "mu_t=" << mu_t;
            os << "," << "vtT=" << vtT;
            os << "," << "DFtT=" << DFtT;
            os << "," << "ql_PtT=" << ql_PtT;
            os << "," << "ql_AtT=" << ql_AtT;
            os << "," << "ql_BtT=" << ql_BtT;
            os << endl;
        }
    }
}
*/