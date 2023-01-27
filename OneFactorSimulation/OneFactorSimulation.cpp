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

#include <iostream>
#include <boost/config.hpp>
#include <ql/quantlib.hpp>
#include <exception>
#include <utils/utilities.h>
#include <map>
#include <chrono>
#include <output.hpp>

// auto link with the correct QuantLib libraries
#ifdef BOOST_MSVC
#include <ql/auto_link.hpp>
#endif

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
    QuantLib::Natural tenorMonths,
    const utils::Paths& adj_B,
    const QuantLib::Date& evalDate = QuantLib::Date()
) {
    auto referenceDate = (evalDate == QuantLib::Date() ? QuantLib::Settings::instance().evaluationDate() : evalDate);
    auto tenorTime = (QuantLib::Time)tenorMonths / 12.0;
    auto nTimes = adj_B.nTimes();
    auto nPaths = adj_B.nPaths();
    auto const& timeGrid = *(adj_B.timeGrid());
    QuantLib::DayCounter act360 = QuantLib::Actual360();
    std::shared_ptr<utils::Paths> pSimpleFwdPaths(new utils::Paths(nTimes, nPaths, adj_B.timeGrid()));
    std::shared_ptr<utils::Paths> pAct360FwdPaths(new utils::Paths(nTimes, nPaths, adj_B.timeGrid()));
    auto const& B_adj = *(adj_B.matrix());
    auto& simpleFwdMatrix = *(pSimpleFwdPaths->matrix());
    auto& act360FwdMatrix = *(pAct360FwdPaths->matrix());
    for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each month
        auto t = timeGrid.at(i);
        auto monthStart = (QuantLib::Natural)std::round(t * 12.0);
        auto monthEnd = monthStart + tenorMonths;
        auto dtStart = referenceDate + monthStart * Months;
        auto dtEnd = referenceDate + monthEnd * Months;
        auto tenorAct360 = act360.yearFraction(dtStart, dtEnd);
        for (decltype(nPaths) path = 0; path < nPaths; ++path) {    // for each path
            auto compounding = 1.0 / B_adj[i][path];
            simpleFwdMatrix[i][path] = (compounding - 1.0) / tenorTime;
            act360FwdMatrix[i][path] = (compounding - 1.0) / tenorAct360;
        }
    }
    std::pair<std::shared_ptr<utils::Paths>, std::shared_ptr<utils::Paths>> ret(pSimpleFwdPaths, pAct360FwdPaths);
    return ret;
}

using DH = utils::DateHelper<char>;

#define DEFAULT_ASOFDATE "T-1"
#define DEFAULT_MR_XML ""
#define DEFAULT_FLAT_FORWARD_RATE ""
#define DEFAULT_MODEL "hw"
#define DEFAULT_MODEL_PARAMS "{\"a\":\"0.03\",\"sigma\":\"0.0001\"}"
#define DEFAULT_MATURITY (50.0) // 50 years
#define DEFAULT_GRID_FREQ (12)  // monthly
#define DEFAULT_N_PATHS (1000)

static string SHORT_RATE_FILENAME_SUFFIX = "short_rate";
static string DISCOUNT_FACTOR_FILENAME_SUFFIX = "discount_factor";
static string PAR_RATE_10_YR_FILENAME_SUFFIX = "swap_rate_10y";
static string FWD_1_MO_FILENAME_SUFFIX = "fwd_1mo";
static string FWD_3_MO_FILENAME_SUFFIX = "fwd_3mo";
static string FWD_3_MO_ACT_360_FILENAME_SUFFIX = "fwd_3mo_act_360";
static string FWD_1_MO_ACT_360_FILENAME_SUFFIX = "fwd_1mo_act_360";
static string FWD_3_MO_AFFINE_FILENAME_SUFFIX = "fwd_3mo_affine";

static string ZV_DISCOUNT_FACTOR_FILENAME_SUFFIX = "zv_discount_factor";
static string ZV_PAR_RATE_10_YR_FILENAME_SUFFIX = "zv_par_rate_10y";

static string CMS_BREAKEVEN_RATE_TREE_FILENAME_SUFFIX = "cms_breakeven_rate_tree";
static string CMS_BREAKEVEN_RATE_SIMULATION_FILENAME_SUFFIX = "cms_breakeven_rate_simulation";

namespace po = boost::program_options;

int main(int argc, char* argv[], char* envp[]) {
    try {
        string asofdate = DEFAULT_ASOFDATE;
        string mrXML = DEFAULT_MR_XML;
        string flatForwardRate = DEFAULT_FLAT_FORWARD_RATE;
        string model = DEFAULT_MODEL;
        string params = DEFAULT_MODEL_PARAMS;
        Real maturity = DEFAULT_MATURITY;
        Natural gridFreq = DEFAULT_GRID_FREQ;
        Size nPaths = DEFAULT_N_PATHS;
        string outputFolder;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("asofdate,d", po::value<string>(&asofdate)->default_value(DEFAULT_ASOFDATE), "evaluation date (T-0|T-1|yyyymmdd)")
            ("mr-xml", po::value<string>(&mrXML)->default_value(DEFAULT_MR_XML), "market rate xml file")
            ("flat-forward", po::value<string>(&flatForwardRate)->default_value(DEFAULT_FLAT_FORWARD_RATE), "flat forward rate in % unit")
            ("model,m", po::value<string>(&model)->default_value(DEFAULT_MODEL), "model type (hw|ghw|ghwln)")
            ("params", po::value<string>(&params)->default_value(DEFAULT_MODEL_PARAMS), "model parameters")
            ("maturity", po::value<Real>(&maturity)->default_value(DEFAULT_MATURITY), "simulation duration in years")
            ("grid-freq,f", po::value<Natural>(&gridFreq)->default_value(DEFAULT_GRID_FREQ), "grid frequency per year: 12=monthly|24=bi-monthly|36=10d|60=6d|72=5d|120=3d|180=2d|360=1d")
            ("paths,p", po::value<Size>(&nPaths)->default_value(DEFAULT_N_PATHS), "number of interest rate paths to be generated")
            ("irp-output-dir,o", po::value<string>(&outputFolder)->default_value(""), "interest rate paths output folder path")
            ("no-adj-df-drift", "do not apply adjustment to discount factor paths to match zv discount factor")
            ("no-bi-vals", "do not generate and output backward induction values")
            ("silent,s", "silent running mode")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        bool silent = (vm.count("silent") > 0);
        auto doNotApplyDiscountFactorDriftAdj = (vm.count("no-adj-df-drift") > 0);
        auto noBackwardInductionValues = (vm.count("no-bi-vals") > 0);
        auto evalDate = DH::getEvaluationDate(asofdate);
        auto modelType = ModelTypeLookup::from_string(model);
        auto paramsJSON = get_input_content(params);

        auto& os = get_nullable_cout<char>(silent);

        os << "asofdate=" << (asofdate.length() > 0 ? asofdate : "(empty)") << endl;
        os << "evalDate=" << DH::to_yyyymmdd(evalDate, true) << endl;
        os << "mrXML=" << mrXML << endl;
        os << "flatForwardRate=" << (flatForwardRate == "" ? "(empty - zv curve will be used)" : flatForwardRate + "%") << endl;
        os << "modelType=" << modelType << endl;
        os << "params=" << paramsJSON << endl;
        os << "maturity=" << maturity << endl;
        os << "gridFreq=" << gridFreq << endl;
        os << "paths=" << nPaths << endl;
        os << "irp-output-dir=" << (outputFolder.length() > 0 ? outputFolder : "(current working directory)") << endl;
        os << "doNotApplyDiscountFactorDriftAdj=" << (doNotApplyDiscountFactorDriftAdj ? "true" : "false") << endl;
        os << "noBackwardInductionValues=" << (noBackwardInductionValues ? "true" : "false") << endl;
        os << "silent=" << (silent ? "true" : "false") << endl;

        if (nPaths == 0) {
            throw std::exception("paths has to be positive (> 0)");
        }

        Settings::instance().evaluationDate() = evalDate;
        DayCounter dayCounter = Thirty360(Thirty360::BondBasis);
        Natural settlementDays = 2;
        Calendar calendar = TARGET();
        auto d = calendar.adjust(evalDate);
        auto settleDate = calendar.advance(d, settlementDays * Days);
        auto const& referenceDate = settleDate;
        auto pParams = instantiateModelParams(modelType, paramsJSON);
        os << endl;
        os << modelType << " model params" << endl;
        os << *pParams;
        auto yieldCurve = (flatForwardRate == "" ? MarketRate(evalDate).load(mrXML).getThirty360ZVDiscountCurve(referenceDate) : ext::shared_ptr<YieldTermStructure>(new FlatForward(referenceDate, boost::lexical_cast<Rate>(flatForwardRate)/100.0, dayCounter)));
        Handle<YieldTermStructure> zvCurve(yieldCurve);
        auto shortRateModel = pParams->createModel(zvCurve);
        auto oneFactorModel = ext::dynamic_pointer_cast<OneFactorModel>(shortRateModel); // the short rate model must be a one factor model
        QL_ASSERT(oneFactorModel != nullptr, "the short rate model is not a one factor model");
        auto oneFactorAffineModel = ext::dynamic_pointer_cast<OneFactorAffineModel>(shortRateModel); // the short rate model must also be a one factor affine model
        QL_ASSERT(oneFactorModel != nullptr, "the short rate model is not a one factor affine model");
        std::shared_ptr<OneFactorAffineLike> oneFactorAffineLike(new OneFactorAffineLikeHack(modelType, oneFactorAffineModel));
        OneFactorShortRateTreeSampling shortRateSampling(*oneFactorModel, yieldCurve, !doNotApplyDiscountFactorDriftAdj);
        // create mersenne twister uniform random generator
        unsigned long seed = 28749;
        MersenneTwisterUniformRng rngGenerator(seed);
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
        using Scenario = std::string;
        auto get_Paths_Scenario_Model_OutputFilePath = [&outputFolder, &evalDate, &model, &nPaths](const string& fnSuffix, const Scenario& scenario) -> string {
            std::ostringstream os;
            os << DH::to_yyyymmdd(evalDate, false);
            os << "_" << nPaths;
            os << "_" << scenario;
            os << "_" << model;
            os << "_" << fnSuffix;
            os << ".txt";
            auto filename = os.str();
            return (outputFolder.length() > 0 ? outputFolder + "\\" + filename : filename);
        };
        auto get_Paths_Model_OutputFilePath = [&outputFolder, &evalDate, &model, &nPaths](const string& fnSuffix) -> string {
            std::ostringstream os;
            os << DH::to_yyyymmdd(evalDate, false);
            os << "_" << nPaths;
            os << "_" << model;
            os << "_" << fnSuffix;
            os << ".txt";
            auto filename = os.str();
            return (outputFolder.length() > 0 ? outputFolder + "\\" + filename : filename);
        };
        auto get_Model_OutputFilePath = [&outputFolder, &evalDate, &model](const string& fnSuffix) -> string {
            std::ostringstream os;
            os << DH::to_yyyymmdd(evalDate, false);
            os << "_" << model;
            os << "_" << fnSuffix;
            os << ".txt";
            auto filename = os.str();
            return (outputFolder.length() > 0 ? outputFolder + "\\" + filename : filename);
        };
        auto get_OutputFilePath = [&outputFolder, &evalDate](const std::string& fnSuffix) -> std::string {
            std::ostringstream oss;
            oss << DH::to_yyyymmdd(evalDate, false);
            oss << "_" << fnSuffix << ".txt";
            auto filename = oss.str();
            return (outputFolder.length() > 0 ? outputFolder + "\\" + filename : filename);
        };
        OneFactorSimulationOutput output(
            os,
            get_Paths_Scenario_Model_OutputFilePath,
            get_Paths_Model_OutputFilePath,
            get_Model_OutputFilePath,
            get_OutputFilePath
        );

        const Rate UP_25 = 25.0 / 10000.0;
        const Rate DN_25 = -UP_25;
        ShockTraits<FlatRateShocker<Rate>> up(UP_25, "up");
        ShockTraits<FlatRateShocker<Rate>> dn(DN_25, "dn");
        ShockTraits<DiscountFactorFlatRateShocker<QuantLib::Compounded, QuantLib::Semiannual>> up_df(UP_25, "up");
        ShockTraits<DiscountFactorFlatRateShocker<QuantLib::Compounded, QuantLib::Semiannual>> dn_df(DN_25, "dn");
        
        output.dumpPathsWithShocks(
            *(monthlyResult.shortRatePaths),
            up,
            dn,
            OutputTraits("short rate", SHORT_RATE_FILENAME_SUFFIX),
            true
        );
        output.dumpPathsWithShocks(
            *(monthlyResult.discountFactorPaths),
            up_df,
            dn_df,
            OutputTraits("discount factor", DISCOUNT_FACTOR_FILENAME_SUFFIX),
            true
        );

        output.dumpZVValues<DiscountFactorZVCalculator>(
            yieldCurve,
            monthlyTimeGrid,
            OutputTraits("zv discount factor", ZV_DISCOUNT_FACTOR_FILENAME_SUFFIX),
            true
        );
        output.dumpZVValues<ParRateZVCalculator<10, Semiannual>>(
            yieldCurve,
            monthlyTimeGrid,
            OutputTraits("zv 10 yr par rate", ZV_PAR_RATE_10_YR_FILENAME_SUFFIX),
            false,
            100.0
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

            SamplingDriftAdjustZCBPrices driftAdjSampling(backwardInductionValuesSampling, yieldCurve, monthlyResult.discountFactorPaths);
            // get 3 month simple forward rates
            {
                os << endl;
                Natural tenorMonths = 3;
                auto driftAdjResult = driftAdjSampling.calculate(
                    tenorMonths,
                    biZCBCalculator,
                    SamplingDriftAdjustZCBPrices::FrontEndAdj,
                    backwardInductionProgressMonitor
                );
                auto pr = getSimpleRatePaths(tenorMonths, *(driftAdjResult.adj_B), evalDate);
                output.dumpPathsWithShocks(
                    *(pr.first),
                    up,
                    dn,
                    OutputTraits("fwd 3mo rate", FWD_3_MO_FILENAME_SUFFIX),
                    false
                );
                output.dumpPathsWithShocks(
                    *(pr.second),
                    up,
                    dn,
                    OutputTraits("fwd 3mo rate (act/360)", FWD_3_MO_ACT_360_FILENAME_SUFFIX),
                    false
                );
                /*
                auto tenorTime = (QuantLib::Time)tenorMonths / 12.0;
                auto const& discountFactorMatrix = *(monthlyResult.discountFactorPaths->matrix());
                auto const& simpleFwdMatrix = *(pr.first->matrix());
                {
                    {
                        os << endl;
                        auto avgPV = 0.0;
                        for (decltype(nPaths) path = 0; path < nPaths; ++path) {
                            auto pv = simpleFwdMatrix[0][path] * tenorTime * discountFactorMatrix[3][path];
                            pv += discountFactorMatrix[3][path];
                            avgPV += pv / nPaths;
                        }
                        os << "avgPV=" << avgPV << endl;
                    }
                    {
                        os << endl;
                        auto avgPV = 0.0;
                        for (decltype(nPaths) path = 0; path < nPaths; ++path) {
                            auto pv = simpleFwdMatrix[0][path] * tenorTime * discountFactorMatrix[3][path];
                            pv += simpleFwdMatrix[3][path] * tenorTime * discountFactorMatrix[6][path];
                            pv += discountFactorMatrix[6][path];
                            avgPV += pv / nPaths;
                        }
                        os << "avgPV=" << avgPV << endl;
                    }
                    {
                        os << endl;
                        auto avgPV = 0.0;
                        for (decltype(nPaths) path = 0; path < nPaths; ++path) {
                            auto pv = simpleFwdMatrix[0][path] * tenorTime * discountFactorMatrix[3][path];
                            pv += simpleFwdMatrix[3][path] * tenorTime * discountFactorMatrix[6][path];
                            pv += simpleFwdMatrix[6][path] * tenorTime * discountFactorMatrix[9][path];
                            pv += discountFactorMatrix[9][path];
                            avgPV += pv / nPaths;
                        }
                        os << "avgPV=" << avgPV << endl;
                    }
                }
                */
            }

#define PAR_TENOR (10)
#define PAR_COUPON_FREQ (QuantLib::Semiannual)
            // get the backend adjusted 10 yr par rates
            {
                os << endl;
                QuantLib::Natural couponTenorMonths = 12 / (QuantLib::Natural)PAR_COUPON_FREQ;
                os << PAR_TENOR << " yr Par Rate Calculation" << endl;
                QuantLib::Natural numCoupons = PAR_TENOR * (QuantLib::Natural)PAR_COUPON_FREQ;
                std::vector<std::shared_ptr<QuantLib::Matrix>> driftAdjBs(numCoupons);  // numCoupons drift-adj B matrices, one for each coupon => P(t, t+0.5), P(t, t+1.0),... P(t, t+10.0)
                std::shared_ptr<QuantLib::TimeGrid> parTimeGrid;
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
                auto nTimes = parTimeGrid->size();
                ParRateCalculator parRateClaculator(PAR_COUPON_FREQ);
                std::shared_ptr<Paths> parRatePaths(new Paths(nTimes, nPaths, parTimeGrid));
                auto& parRateMatrix = *(parRatePaths->matrix());
                for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each time slice
                    for (decltype(nPaths) path = 0; path < nPaths; ++path) { // for each path
                        std::vector<QuantLib::DiscountFactor> dfs(numCoupons);
                        for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                            auto const& adj_B = *(driftAdjBs.at(k));
                            dfs[k] = adj_B[i][path];
                        }
                        auto parRate = parRateClaculator(dfs.begin(), dfs.end());
                        parRateMatrix[i][path] = parRate;
                    }
                }

                os << endl;
                os << PAR_TENOR << " yr Par Rate ZV Comparison Report" << endl;
                zvComparison.report<ParRateZVCalculator<PAR_TENOR, PAR_COUPON_FREQ>>(*parRatePaths);
                output.dumpPathsWithShocks(
                    *parRatePaths,
                    up,
                    dn,
                    OutputTraits("10 yr par rate", PAR_RATE_10_YR_FILENAME_SUFFIX),
                    false
                );

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
            }

            /*
            // get unadjusted 10 yr par rates
            {
                os << endl;
                auto result = backwardInductionValuesSampling.sample(PAR_TENOR*12, biParRateCalculator, backwardInductionProgressMonitor);
                os << endl;
                os << PAR_TENOR << " yr Par Rate ZV Comparison Report" << endl;
                zvComparison.report<ParRateZVCalculator<PAR_TENOR, PAR_COUPON_FREQ>>(*(result.valuePaths));
                output.dumpPathsWithShocks(
                    *(result.valuePaths),
                    up,
                    dn,
                    OutputTraits("10 yr par rate", PAR_RATE_10_YR_FILENAME_SUFFIX),
                    false
                );
            }
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