#pragma once

#include <ql/quantlib.hpp>
#include <iostream> 
#include <memory>
#include <algorithm>
#include <vector>
#include <list>
#include <utils/json_io.h>
#include <utils/quant_type.hpp>
#include <utils/black_model_calibrator.hpp>

namespace utils {
    class ShortRateModelParams : public json_io::json_serializable {
    public:
        virtual QuantLib::ext::shared_ptr<QuantLib::ShortRateModel> createModel(const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure) const = 0;
        virtual void dump(std::ostream& os) const = 0;
    };

    class ShortRateBlackCalibration;

    struct ShortRateCalibrationProcess {
        virtual std::shared_ptr<ShortRateModelParams> operator() (ShortRateBlackCalibration& calibration, BlackModelCalibrator& blackCalibrator) = 0;
    };

    class CalibrationGridFactory {
    public:
        virtual std::shared_ptr<QuantLib::TimeGrid> operator() (const std::list<QuantLib::Time>& mandatoryTimes) const = 0;
    };

    class MandatoryTimesWithStepsCalibrationGridFactory : public CalibrationGridFactory {
    protected:
        QuantLib::Size steps_;
    public:
        MandatoryTimesWithStepsCalibrationGridFactory(QuantLib::Size steps = 0): steps_(steps) {}
        const QuantLib::Size& step() const { return steps_; }
        QuantLib::Size& step() { return steps_; }
        std::shared_ptr<QuantLib::TimeGrid> operator() (const std::list<QuantLib::Time>& mandatoryTimes) const {
            std::shared_ptr<QuantLib::TimeGrid> ret;
            ret.reset(step() != 0 ? new QuantLib::TimeGrid(mandatoryTimes.begin(), mandatoryTimes.end(), step()) : new QuantLib::TimeGrid(mandatoryTimes.begin(), mandatoryTimes.end()));
            return ret;
        }
    };

    class FixedFrequencyCalibrationGridFactory : public CalibrationGridFactory {
    protected:
        QuantLib::Natural frequency_;
    public:
        FixedFrequencyCalibrationGridFactory(QuantLib::Natural frequency = 12) : frequency_(frequency) {
            QL_REQUIRE(frequency_ > 0, "frequency (" << frequency_ << ") must be greater than 0");
        }
        const QuantLib::Natural& frequency() const { return frequency_; }
        QuantLib::Natural& frequency() { return frequency_; }
        std::shared_ptr<QuantLib::TimeGrid> operator() (const std::list<QuantLib::Time>& mandatoryTimes) const {
            std::vector<QuantLib::Time> times(mandatoryTimes.begin(), mandatoryTimes.end());
            std::sort(times.begin(), times.end(), std::less<QuantLib::Time>());
            auto maxTime = times.back();
            decltype(maxTime) t = 0.0;
            while (t < maxTime) {
                times.push_back(t);
                t += 1.0 / (QuantLib::Time)frequency();
            }
            std::sort(times.begin(), times.end(), std::less<QuantLib::Time>());
            return std::shared_ptr<QuantLib::TimeGrid>(new QuantLib::TimeGrid(times.begin(), times.end()));
        }
    };

    // short rate calibration class with Black option vols
    // aslo support calibration verification
    class ShortRateBlackCalibration {
    public:
        using OptionsData = std::vector<OptionVolData<QuantLib::Volatility>>;
        using CalibrationHelpers = std::vector<QuantLib::ext::shared_ptr<QuantLib::BlackCalibrationHelper>>;

    protected:
        std::ostream& ostream_;
        // input of the calibration
        const OptionsData& vols_;   // Black vol data
        QuantLib::ext::shared_ptr<QuantLib::IborIndex> iborIndex_;
        QuantLib::Handle<QuantLib::YieldTermStructure> yieldtermStructure_;
        // output of the calibration
        std::shared_ptr<CalibrationHelpers> helpers_;   // Black calibration helpers
        std::shared_ptr<QuantLib::TimeGrid> grid_;      // calibration time grid
        std::shared_ptr<ShortRateModelParams> calibratedParams_; // calibrated params;
    public:
        ShortRateBlackCalibration(
            std::ostream& ostream,
            const OptionsData& vols,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& iborIndex,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure
        ) : ostream_(ostream), vols_(vols), iborIndex_(iborIndex), yieldtermStructure_(yieldtermStructure) {}
        std::ostream& ostream() {
            return ostream_;
        }
        const OptionsData& vols() const {
            return vols_;
        }
        const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& iborIndex() const {
            return iborIndex_;
        }
        const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure() const {
            return yieldtermStructure_;
        }
        const CalibrationHelpers& helpers() const {
            return *helpers_;
        }
        const QuantLib::TimeGrid& grid() const {
            return *grid_;
        }
        const std::shared_ptr<ShortRateModelParams>& calibratedParams() const {
            return calibratedParams_;
        }
    protected:
        template <typename HELPER_FACTORY>
        static std::shared_ptr<CalibrationHelpers> createHelpers(
            const OptionsData& vols,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& iborIndex,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure
        ) {
            std::shared_ptr<CalibrationHelpers> helpers(new CalibrationHelpers());
            HELPER_FACTORY helperFactory;
            for (auto const& i : vols) {    // for each option vol data 
                QuantLib::Handle<QuantLib::Quote> quoteHandle(QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(i.data)));
                auto helper = helperFactory(quoteHandle, i.expiry, i.tenor, iborIndex, yieldtermStructure);
                helpers->push_back(helper);
            }
            return helpers;
        }
        static std::shared_ptr<QuantLib::TimeGrid> createCalibrationGrid(const CalibrationHelpers& helpers, const CalibrationGridFactory& calibrationGridFactory) {
            std::list<QuantLib::Time> times;
            for (auto const& helper : helpers) {    // for each helper
                helper->addTimesTo(times);
            }
            // create the calibration time grid by calling the factory object
            return calibrationGridFactory(times);
        }
    public:
        // calibration verification report
        static void verifyCalibration(
            std::ostream& os,
            const QuantLib::ext::shared_ptr<QuantLib::PricingEngine>& engine,
            const CalibrationHelpers& calibrationHelpers,
            const OptionsData& vols
        ) {
            QL_ASSERT(calibrationHelpers.size() == vols.size(), "calibration verification data size mis-matched");
            auto n = calibrationHelpers.size();
            os << std::endl << "* Model Calibration Verification *" << std::endl;
            os << "******************************************************************************************" << std::endl;
            for (decltype(n) i = 0; i < n; i++) {
                auto const& volData = vols[i];
                auto const& helper = calibrationHelpers[i];
                auto const& marketVol = volData.data;
                //auto marketBlackNVP = helper->blackPrice(marketVol);
                //os << marketBlackNVP << endl;
                helper->setPricingEngine(engine);
                auto modelNVP = helper->modelValue();
                //os << modelNVP << endl;
                auto modelImpliedVol = helper->impliedVolatility(modelNVP, 1e-16, 10000, 0.00001, 1.0);
                auto diffVol = modelImpliedVol - marketVol;
                os << volData.expiry << "x" << volData.tenor;
                os << std::setprecision(5) << std::noshowpos;
                os << ", " << "model implied vol=" << std::setw(7) << QuantLib::io::volatility(modelImpliedVol);
                os << ", " << "market vol=" << std::setw(7) << QuantLib::io::volatility(marketVol);
                os << ", " << "diff=" << std::setw(7) << std::showpos << QuantLib::io::volatility(diffVol) << std::noshowpos << ")";
                os << std::endl;
            }
            os << "******************************************************************************************" << std::endl;
        }
        // calibration verification report
        void verifyCalibration(const QuantLib::ext::shared_ptr<QuantLib::PricingEngine>& engine) {
            verifyCalibration(ostream(), engine, helpers(), vols());
        }
        void ensureValidFixedParams(const std::vector<bool>& fixedParameters) const {
            auto nHelpers = helpers().size();
            decltype(nHelpers) nFreeParams = std::count_if(fixedParameters.begin(), fixedParameters.end(), [](const bool& fixed) -> bool {
                return !fixed;
            });
            QL_ENSURE(nFreeParams <= nHelpers, "the number of free calibration parameters (" << nFreeParams << ") can not be more than the number of calibration helpers (" << nHelpers << ")");
        };
        // run the calibration
        template <typename HELPER_FACTORY>
        std::shared_ptr<ShortRateModelParams> calibrate(ShortRateCalibrationProcess& calibrationProcess, const CalibrationGridFactory& calibrationGridFactory) {
            helpers_ = createHelpers<HELPER_FACTORY>(vols(), iborIndex(), yieldtermStructure());
            grid_ = createCalibrationGrid(helpers(), calibrationGridFactory);
            QuantLib::EndCriteria endCriteria(10000, 100, 0.000001, 0.00000001, 0.00000001);
            BlackModelCalibrator blackCalibrator(endCriteria, helpers_);
            calibratedParams_ = calibrationProcess(*this, blackCalibrator);
            return calibratedParams_;
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const utils::ShortRateModelParams& params) {
    params.dump(os);
    return os;
}