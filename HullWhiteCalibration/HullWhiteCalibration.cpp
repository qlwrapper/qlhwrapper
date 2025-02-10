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
#include <sstream>
#include <vector>
#include <utils/utilities.h>
#include <boost/lexical_cast.hpp>
#include <common/json-obj-serialization.hpp>

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

class ShortRateCalibrationProcessBase : public ShortRateCalibrationProcess {
public:
    void showCalibrationInput(
        std::ostream& os,
        const ShortRateModelParams& modelParams,
        const std::vector<bool>& fixedParameters = std::vector<bool>()
    ) const {
        os << std::endl;
        os << "< Calibration Input >" << endl;
        os << modelParams;
        if (!fixedParameters.empty()) {
            os << "fixedParameters=" << fixedParameters << endl;
        }
    };
    void showCalibrationResult(
        std::ostream& os,
        BlackModelCalibrator& blackCalibrator,
        const QuantLib::Array& calibratedParams
    ) {
        os << std::endl;
        os << "< Calibration Result >" << endl;
        os << "functionEvaluation=" << blackCalibrator.functionEvaluation << std::endl;
        os << "endCriteriaType=" << blackCalibrator.endCriteriaType << std::endl;
        os << "durationSeconds=" << blackCalibrator.durationSeconds << std::endl;
        os << "calibrated params=" << calibratedParams << endl;
    };
};

// base class for all HullWhite based calibration process
class HullWhiteCalibrationProcessBase : public ShortRateCalibrationProcessBase {
protected:
    QuantLib::Real meanReversionFixed_;
public:
    HullWhiteCalibrationProcessBase(
        QuantLib::Real meanReversionFixed = QuantLib::Null<QuantLib::Real>()
    ) : meanReversionFixed_(meanReversionFixed) {}
    const QuantLib::Real& meanReversionFixed() const {
        return meanReversionFixed_;
    }
    QuantLib::Real& meanReversionFixed() {
        return meanReversionFixed_;
    }
    bool meanReversionToBeFixed() const {
        return (meanReversionFixed() != QuantLib::Null<QuantLib::Real>());
    }
};

class HullWhiteCalibrationProcess : public HullWhiteCalibrationProcessBase {
public:
    HullWhiteCalibrationProcess(QuantLib::Real meanReversionFixed = QuantLib::Null<QuantLib::Real>()) : HullWhiteCalibrationProcessBase(meanReversionFixed) {}
    std::shared_ptr<ShortRateModelParams> operator() (ShortRateBlackCalibration& calibration, BlackModelCalibrator& blackCalibrator) {
        auto& os = calibration.ostream();
        auto const& yieldtermStructure = calibration.yieldtermStructure();
        auto const& grid = calibration.grid();
        HullWhiteModelParams modelParams(meanReversionToBeFixed() ? meanReversionFixed() : 0.1, 0.01);
        auto model = modelParams.createModel(yieldtermStructure);
        ext::shared_ptr<QuantLib::PricingEngine> engine(new QuantLib::TreeSwaptionEngine(model, grid));
        std::vector<bool> fixedParameters{ meanReversionToBeFixed(), false};
        calibration.ensureValidFixedParams(fixedParameters);
        showCalibrationInput(os, modelParams, fixedParameters);
        blackCalibrator.calibrate(model, engine, fixedParameters);
        auto params = model->params();
        showCalibrationResult(os, blackCalibrator, params);
        calibration.verifyCalibration(engine);
        auto a = params[0];
        auto sigma = params[1];
        std::shared_ptr<ShortRateModelParams> ret(new HullWhiteModelParams(a, sigma));
        return ret;
    }
};

template <typename VolTraits = utils::NormalVolTraits<QuantLib::Real>>
class GeneralizedHullWhiteCalibrationProcess : public HullWhiteCalibrationProcessBase {
public:
    GeneralizedHullWhiteCalibrationProcess(QuantLib::Real meanReversionFixed = QuantLib::Null<QuantLib::Real>()) : HullWhiteCalibrationProcessBase(meanReversionFixed) {}
    std::shared_ptr<ShortRateModelParams> operator() (ShortRateBlackCalibration& calibration, BlackModelCalibrator& blackCalibrator) {
        auto& os = calibration.ostream();
        auto const& yieldtermStructure = calibration.yieldtermStructure();
        auto const& grid = calibration.grid();
        // calibrate one node of mean reversion and one node of sigma simultaneously
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        Real a = (meanReversionToBeFixed() ? meanReversionFixed() : 0.1);
        Real vol = 0.01;
        {
            os << endl << "Generalized Hull-White calibration - Step 1" << endl;
            vector<Time> sigmaTimes{ 0.0 };
            vector<Real> sigma_(sigmaTimes.size(), vol);
            GeneralizedHullWhiteModelParams<VolTraits> modelParams(a, sigmaTimes, sigma_);
            auto model = modelParams.createModel(yieldtermStructure);
            ext::shared_ptr<PricingEngine> engine(new TreeSwaptionEngine(model, grid));
            std::vector<bool> fixedParameters{ meanReversionToBeFixed(), false };
            calibration.ensureValidFixedParams(fixedParameters);
            showCalibrationInput(os, modelParams, fixedParameters);
            blackCalibrator.calibrate(model, engine, fixedParameters);
            auto params = model->params();
            showCalibrationResult(os, blackCalibrator, params);
            calibration.verifyCalibration(engine);

            a = params[0];
            vol = params[1];
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        vector<Real> sigma;
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        {
            os << endl << "Generalized Hull-White calibration - Step 2" << endl;
            vector<Time> sigmaTimes{ 0.0, 1.0, 5.0, 10.0, 15.0 };
            vector<Real> sigma_(sigmaTimes.size(), vol);
            GeneralizedHullWhiteModelParams<VolTraits> modelParams(a, sigmaTimes, sigma_);
            auto model = modelParams.createModel(yieldtermStructure);
            ext::shared_ptr<PricingEngine> engine(new TreeSwaptionEngine(model, grid));
            std::vector<bool> fixedParameters{ true, false, false, false, false, false };   // fixed the mean reversion
            calibration.ensureValidFixedParams(fixedParameters);
            showCalibrationInput(os, modelParams, fixedParameters);
            blackCalibrator.calibrate(model, engine, fixedParameters);
            auto params = model->params();
            showCalibrationResult(os, blackCalibrator, params);
            calibration.verifyCalibration(engine);

            sigma.insert(sigma.end(), std::next(params.begin(), 1), params.end());
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        vector<Time> sigmaTimes{ 0.0, 1.0, 5.0, 10.0, 15.0 };
        std::shared_ptr<ShortRateModelParams> ret(new GeneralizedHullWhiteModelParams<VolTraits>(a, sigmaTimes, sigma));
        return ret;
    }
};

class G2CalibrationProcess : public ShortRateCalibrationProcessBase {
public:
    std::shared_ptr<ShortRateModelParams> operator() (
        ShortRateBlackCalibration& calibration,
        BlackModelCalibrator& blackCalibrator
    ) {
        auto& os = calibration.ostream();
        auto const& yieldtermStructure = calibration.yieldtermStructure();
        auto const& grid = calibration.grid();
        QuantLib::Real a_guess = 0.04;
        QuantLib::Real sigma_guess = 0.004;
        QuantLib::Real b_guess = a_guess;
        QuantLib::Real eta_guess = sigma_guess;
        QuantLib::Real rho_guess = 0.05;
        G2ModelParams modelParams(a_guess, sigma_guess, b_guess, eta_guess, rho_guess);
        auto model = modelParams.createModel(yieldtermStructure);
        auto g2Model = ::ext::dynamic_pointer_cast<QuantLib::G2>(model);
        ext::shared_ptr<QuantLib::PricingEngine> engine(new QuantLib::G2SwaptionEngine(g2Model, 10, 1000));
        //ext::shared_ptr<QuantLib::PricingEngine> engine(new QuantLib::TreeSwaptionEngine(model, grid));
        showCalibrationInput(os, modelParams);
        blackCalibrator.calibrate(model, engine);
        auto params = model->params();
        showCalibrationResult(os, blackCalibrator, params);
        calibration.verifyCalibration(engine);
        auto a = params[0];
        auto sigma = params[1];
        auto b = params[2];
        auto eta = params[3];
        auto rho = params[4];
        std::shared_ptr<ShortRateModelParams> ret(new G2ModelParams(a, sigma, b, eta, rho));
        return ret;
    }
};

static std::shared_ptr<ShortRateCalibrationProcess> createShortRateCalibrationProcess(
    ModelType modelType,
    const string& a_fixed
) {
    auto meanReversionFixed = (!a_fixed.empty() ? boost::lexical_cast<QuantLib::Real>(a_fixed) : QuantLib::Null<QuantLib::Real>());
    std::shared_ptr<ShortRateCalibrationProcess> process;
    switch (modelType) {
    case ModelType::HullWhite1F:
        process.reset(new HullWhiteCalibrationProcess(meanReversionFixed));
        break;
    case ModelType::GeneralizedHullWhite1F:
        process.reset(new GeneralizedHullWhiteCalibrationProcess<NormalVolTraits<QuantLib::Real>>(meanReversionFixed));
        break;
    case ModelType::GeneralizedHullWhiteLogNormal1F:
        process.reset(new GeneralizedHullWhiteCalibrationProcess<LogNormalVolTraits<QuantLib::Real>>(meanReversionFixed));
        break;
    case ModelType::G2:
        process.reset(new G2CalibrationProcess());
        break;
    default:
        QL_FAIL("unknown/unsupported short rate model type (" << modelType << ")");
    }
    return process;
}

static void verifyATMStrikeForQuotedVols(
    ostream& os,
    const Handle<YieldTermStructure>& hTS,
    const Handle<YieldTermStructure>& hTSAct365F,
    bool isSofr,
    const QuotedVols& quotedVols
) {
    QL_ASSERT(!quotedVols.empty(), "quoted vols object cannot be empty");
    vector<OptionAttribs> v_optAttibutes;
    v_optAttibutes.push_back(OptionAttribs(0 * Months, quotedVols[0]->tenor));    // added spot a swap
    auto n = quotedVols.size();
    for (decltype(n) i = 0; i < n; ++i) {
        const OptionAttribs& optionAttribs = *(quotedVols[i]);
        v_optAttibutes.push_back(optionAttribs);
    }
    n = v_optAttibutes.size();
    os << "ATM strikes (30/360 vs Act/365F):" << endl;
    os << "=============================================================================================" << endl;
    for (decltype(n) i = 0; i < n; ++i) {
        const auto& optionAttribs = v_optAttibutes[i];
        const auto& forward = optionAttribs.expiry;
        const auto& tenor = optionAttribs.tenor;
        auto swap = ShortRateBlackCalibration::makeFwdSwap(isSofr, forward, tenor, hTS);
        auto strike = swap->fairRate();
        auto swapAct365F = ShortRateBlackCalibration::makeFwdSwap(isSofr, forward, tenor, hTSAct365F);
        auto strikeAct365F = swapAct365F->fairRate();
        os << optionAttribs;
        os << " [" << DateFormat<char>::to_yyyymmdd(swap->startDate()) << "," << DateFormat<char>::to_yyyymmdd(swap->maturityDate()) << "]";
        os << " ==> " << "30/360=" << io::percent(strike) << " vs Act/365F=" << io::percent(strikeAct365F);
        os << ", diff=" << io::basis_point(strike - strikeAct365F);
        os << endl;
    }
    os << "=============================================================================================" << endl;
}

namespace po = boost::program_options;

int main(int argc, char* argv[], char* envp[]) {
    try {
        string asofdate;
        string s_vols;
        string s_zeroCurve;
        string model = DEFAULT_MODEL;
        string meanReversionFixed;
        Natural gridFreq = DEFAULT_GRID_FREQ;
        bool isSofr = DEFAULT_IS_SOFR;
        bool silent = DEFAULT_SILENT;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("asofdate,d", po::value<string>(&asofdate)->default_value(""), "(optional) evaluation date in yyyymmdd format. default to today's date")
            ("vols,v", po::value<string>(&s_vols)->default_value(""), "(required) ATM volatility calibration input json file. \"@-\" to piped in from stdin")
            ("zero-curve,z", po::value<string>(&s_zeroCurve)->default_value(""), "(required) zero/spot curve input json file. \"@-\" to piped in from stdin")
            ("model,m", po::value<string>(&model)->default_value(DEFAULT_MODEL), "(optional) model type (hw|ghw|ghwln|g2). The default hw")
            ("mean-reversion,a", po::value<string>(&meanReversionFixed)->default_value(""), "(optional) locked mean reversion speed (a). if not specified, it will be calibrated")
            ("grid-freq,f", po::value<Natural>(&gridFreq)->default_value(DEFAULT_GRID_FREQ), "(optional) grid frequency per year: 12=monthly|24=bi-monthly|36=10d|60=6d|72=5d|120=3d|180=2d|360=1d. The default value is 12 (monthly)")
            ("is-sofr", po::value<bool>(&isSofr)->default_value(DEFAULT_IS_SOFR), "(optional) calibrating the model using SOFR swaptions. true or false. The default value is true")
            ("silent,s", po::value<bool>(&silent)->default_value(DEFAULT_SILENT), "(optional) silent running mode. true or false. The default value is false")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }
        QL_REQUIRE(!s_vols.empty(), "ATM volatility calibration input json file is not optional");
        QL_REQUIRE(!s_zeroCurve.empty(), "zero/spot curve input json file is not optional");
        if (model.empty()) {
            model = DEFAULT_MODEL;
        }
        QL_REQUIRE(gridFreq > 0, "grid frequency (" << gridFreq << ") must be a positive integer");
        auto evalDate = (asofdate.empty() ? Date::todaysDate() : DateFormat<char>::from_yyyymmdd(asofdate, false));
        auto modelType = boost::lexical_cast<ModelType>(model);

        auto& os = get_nullable_cout<char>(silent);

        os << endl;
        os << "Hull-White Model Calibrator" << endl;
        os << "asofdate=" << DateFormat<char>::to_yyyymmdd(evalDate, true) << endl;
        os << "vols=" << s_vols << endl;
        os << "zero-curve=" << s_zeroCurve << endl;
        os << "model=" << boost::lexical_cast<string>(modelType) << endl;
        os << "mean-reversion=" << (meanReversionFixed.empty() ? "(to be calibrated)" : meanReversionFixed) << endl;
        os << "grid-freq=" << gridFreq << endl;
        os << "is-sofr=" << (isSofr ? "true" : "false") << endl;
        os << "silent=" << (silent ? "true" : "false") << endl;
        
        Settings::instance().evaluationDate() = evalDate;

        auto pQuotedVols = MarketRate::readQuotedVols(s_vols);
        const auto& quotedVols = *pQuotedVols;
        os << endl;
        os << "quoted vols input:" << endl;
        os << "==============================" << endl;
        os << quotedVols << endl;
        os << "==============================" << endl;

        auto pYieldCurve = MarketRate::readAsThirty360ZeroCurve(s_zeroCurve, evalDate, true);
        auto pYieldCurveAct365F = MarketRate::readAsActual365FixedZeroCurve(s_zeroCurve, evalDate, true);
        Handle<YieldTermStructure> hTS(pYieldCurve);
        Handle<YieldTermStructure> hTSAct365F(pYieldCurveAct365F);

        os << endl;
        verifyATMStrikeForQuotedVols(os, hTS, hTSAct365F, isSofr, quotedVols);

        ShortRateBlackCalibration shortRateCalibration(os, pQuotedVols, hTS);
        auto calibrationProcess = createShortRateCalibrationProcess(modelType, meanReversionFixed);
        auto calibratedParams = shortRateCalibration.calibrate(isSofr, *calibrationProcess, FixedFrequencyCalibrationGridFactory(gridFreq));

        os << endl;
        os << (*calibratedParams);
        os << endl;
        os << "calibration grid min=" << shortRateCalibration.grid().front() << endl;
        os << "calibration grid max=" << shortRateCalibration.grid().back() << endl;
        os << "calibration grid size=" << shortRateCalibration.grid().size() << endl;
        os << endl;
        cout << calibratedParams->JSON_stringify<char>(true);
        os << endl;
        os << "Done" << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "!!! Error: " << e.what() << endl;
        return 1;
    }
}