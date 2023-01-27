#include <iostream>
#include <boost/config.hpp>
#include <ql/quantlib.hpp>
#include <exception>
#include <utils/utilities.h>
#include <boost/lexical_cast.hpp>

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

static vector<OptionAttribs> toBeCalibratedSwaptions{
    {3 * Months, 10 * Years}
    ,{1 * Years, 10 * Years}
    ,{5 * Years, 10 * Years}
    ,{10 * Years, 10 * Years}
    ,{15 * Years, 10 * Years}
};

class HullWhiteCalibrationProcessBase : public ShortRateCalibrationProcess {
protected:
    QuantLib::Real meanReversionFixed_;
public:
    HullWhiteCalibrationProcessBase(QuantLib::Real meanReversionFixed = QuantLib::Null<QuantLib::Real>()) : meanReversionFixed_(meanReversionFixed) {}
    const QuantLib::Real& meanReversionFixed() const {
        return meanReversionFixed_;
    }
    QuantLib::Real& meanReversionFixed() {
        return meanReversionFixed_;
    }
    bool meanReversionToBeFixed() const {
        return (meanReversionFixed() != QuantLib::Null<QuantLib::Real>());
    }
    void showCalibrationInput(std::ostream& os, const ShortRateModelParams& modelParams, const std::vector<bool>& fixedParameters) const {
        os << std::endl;
        os << "< Calibration Input >" << endl;
        os << modelParams;
        os << "fixedParameters=" << fixedParameters << endl;
    };
    void showCalibrationResult(std::ostream& os, BlackModelCalibrator& blackCalibrator, const QuantLib::Array& calibratedParams) {
        os << std::endl;
        os << "< Calibration Result >" << endl;
        os << "functionEvaluation=" << blackCalibrator.functionEvaluation << std::endl;
        os << "endCriteriaType=" << blackCalibrator.endCriteriaType << std::endl;
        os << "durationSeconds=" << blackCalibrator.durationSeconds << std::endl;
        os << "calibrated params=" << calibratedParams << endl;
    };
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

/*
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

        //////////////////////////////////////////////////////////////////////////////////////////////////////
        vector<Real> sigmaPrev;
        {
            os << endl << "Generalized Hull-White calibration - Step 2" << endl;
            vector<Time> sigmaTimes{ 0.0, 1.0, 5.0};
            vector<Real> sigma(sigmaTimes.size(), vol);
            GeneralizedHullWhiteModelParams<VolTraits> modelParams(a, sigmaTimes, sigma);
            auto model = modelParams.createModel(yieldtermStructure);
            ext::shared_ptr<PricingEngine> engine(new TreeSwaptionEngine(model, grid));
            std::vector<bool> fixedParameters{ true, false, false, false};
            calibration.ensureValidFixedParams(fixedParameters);
            showCalibrationInput(os, modelParams, fixedParameters);
            blackCalibrator.calibrate(model, engine, fixedParameters);
            auto params = model->params();
            showCalibrationResult(os, blackCalibrator, params);
            calibration.verifyCalibration(engine);

            sigmaPrev.insert(sigmaPrev.end(), std::next(params.begin(), 1), params.end());
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////

        vector<Real> sigma;
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        {
            os << endl << "Generalized Hull-White calibration - Step 3" << endl;
            vector<Time> sigmaTimes{ 0.0, 1.0, 5.0, 10.0, 15.0 };
            vector<Real> sigma_(sigmaPrev.begin(), sigmaPrev.end());
            sigma_.resize(sigmaTimes.size(), *std::prev(sigmaPrev.end()));
            GeneralizedHullWhiteModelParams<VolTraits> modelParams(a, sigmaTimes, sigma_);
            auto model = modelParams.createModel(yieldtermStructure);
            ext::shared_ptr<PricingEngine> engine(new TreeSwaptionEngine(model, grid));
            std::vector<bool> fixedParameters{ true, true, true, false, false, false };
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
*/

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

static std::shared_ptr<HullWhiteCalibrationProcessBase> createHWCalibrationProcess(ModelType modelType, const string& a_fixed) {
    auto meanReversionFixed = (a_fixed != "" ? boost::lexical_cast<QuantLib::Real>(a_fixed) : QuantLib::Null<QuantLib::Real>());
    std::shared_ptr<HullWhiteCalibrationProcessBase> process;
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
    default:
        QL_FAIL("unknown model type");
    }
    return process;
}

using DH = utils::DateHelper<char>;

#define DEFAULT_ASOFDATE "T-1"
#define DEFAULT_MR_XML ""
#define DEFAULT_MODEL "hw"
#define DEFAULT_MEAN_VERSION_FIXED ""
#define DEFAULT_GRID_FREQ (12)  // monthly

namespace po = boost::program_options;

int main(int argc, char* argv[], char* envp[]) {
    try {
        string asofdate = DEFAULT_ASOFDATE;
        string mrXML = DEFAULT_MR_XML;
        string model = DEFAULT_MODEL;
        string meanReversionFixed = DEFAULT_MEAN_VERSION_FIXED;
        Natural gridFreq = DEFAULT_GRID_FREQ;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("asofdate,d", po::value<string>(&asofdate)->default_value(DEFAULT_ASOFDATE), "evaluation date (T-0|T-1|yyyymmdd)")
            ("mr-xml", po::value<string>(&mrXML)->default_value(DEFAULT_MR_XML), "market rate xml file")
            ("model,m", po::value<string>(&model)->default_value(DEFAULT_MODEL), "model type (hw|ghw|ghwln)")
            ("mean-reversion,a", po::value<string>(&meanReversionFixed)->default_value(DEFAULT_MEAN_VERSION_FIXED), "mean reversion speed. if not specified, it will be calibrated")
            ("grid-freq,f", po::value<Natural>(&gridFreq)->default_value(DEFAULT_GRID_FREQ), "grid frequency per year: 12=monthly|24=bi-monthly|36=10d|60=6d|72=5d|120=3d|180=2d|360=1d")
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
        auto evalDate = DH::getEvaluationDate(asofdate);
        auto modelType = ModelTypeLookup::from_string(model);

        auto& os = get_nullable_cout<char>(silent);

        os << "asofdate=" << (asofdate.length() > 0 ? asofdate : "(empty)") << endl;
        os << "mrXML=" << mrXML << endl;
        os << "evalDate=" << DH::to_yyyymmdd(evalDate, true) << endl;
        os << "modelType=" << modelType << endl;
        os << "mean-reversion=" << (meanReversionFixed == "" ? "(to be calibrated)" : meanReversionFixed) << endl;
        os << "gridFreq=" << gridFreq << endl;
        os << "silent=" << (silent ? "true" : "false") << endl;

        // general parameters
        Settings::instance().evaluationDate() = evalDate;
        DayCounter dayCounter = Thirty360(Thirty360::BondBasis);
        Natural settlementDays = 2;
        Calendar calendar = TARGET();
        auto d = calendar.adjust(evalDate);
        auto settleDate = calendar.advance(d, settlementDays * Days);
        auto const& referenceDate = settleDate;

        MarketRate mr(evalDate);
        mr.load(mrXML);

        auto yieldCurve = mr.getThirty360ZVDiscountCurve(referenceDate);
        Handle<YieldTermStructure> yieldtermStructure(yieldCurve);
        ext::shared_ptr<IborIndex> iborIndex(new USDLibor(Period(3, Months), yieldtermStructure));
        auto logNormalVols = mr.getATMSwaptionVolData<LogNormalVol>();
        auto volData = mr.filterOptionVolData(*logNormalVols, toBeCalibratedSwaptions);
        auto const& vols = *volData;
        ShortRateBlackCalibration shortRateCalibration(os, vols, iborIndex, yieldtermStructure);
        auto calibrationProcess = createHWCalibrationProcess(modelType, meanReversionFixed);
        auto calibratedParams = shortRateCalibration.calibrate<USDSwaptionHelperFactory>(*calibrationProcess, FixedFrequencyCalibrationGridFactory(gridFreq));

        os << endl;
        os << (*calibratedParams);
        os << endl;
        os << "calibration grid min=" << shortRateCalibration.grid().front() << endl;
        os << "calibration grid max=" << shortRateCalibration.grid().back() << endl;
        os << "calibration grid size=" << shortRateCalibration.grid().size() << endl;
        os << endl;
        cout << calibratedParams->JSON_stringify<char>(true);

        {
            auto pVols = mr.filterOptionVolData(*logNormalVols, std::vector<OptionAttribs>{ {10 * Years, 10 * Years}});
            auto log_normal_vol_10y_10y = pVols->at(0).data;
            // calculate forward par rate from t to T given yield term structure and coupon frequency
            auto zvParRateCalc = [&yieldCurve](Time t, Time T, QuantLib::Frequency couponFreq = QuantLib::Semiannual) -> Rate {
                auto tenor = T - t;
                ParRateCalculator rateCalculator(couponFreq);
                auto numCoupons = rateCalculator.numCoupons(tenor);
                auto dfStart = yieldCurve->discount(t);
                std::vector<QuantLib::DiscountFactor> dfs(numCoupons);
                for (decltype(numCoupons) k = 0; k < numCoupons; ++k) { // for each coupon
                    auto couponTime = rateCalculator.couponTime(t, k);
                    auto df = yieldCurve->discount(couponTime) / dfStart;
                    dfs[k] = df;
                }
                return rateCalculator(dfs.begin(), dfs.end());
            };
            auto parRate_10y_10y = zvParRateCalc(10, 20);
            auto normal_vol_10y_10y = log_normal_vol_10y_10y * parRate_10y_10y;

            auto normalVols = mr.getATMSwaptionVolData<NormalVol>();
            auto ret = mr.filterOptionVolData(*normalVols, std::vector<OptionAttribs>{ {10 * Years, 10 * Years}});
            auto mr_normal_vol_10y_10y = ret->at(0).data;

            os << endl;
            os << "parRate_10y_10y=" << parRate_10y_10y << endl;
            os << "normal_vol_10y_10y=" << normal_vol_10y_10y << endl;
            os << "mr_normal_vol_10y_10y=" << mr_normal_vol_10y_10y << endl;
        }
               
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "!!! Error: " << e.what() << endl;
        return 1;
    }
}