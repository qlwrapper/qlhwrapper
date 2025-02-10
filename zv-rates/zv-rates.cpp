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
#include <string>
#include <boost/config.hpp>
#include <ql/quantlib.hpp>
#include <exception>
#include <utils/utilities.h>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <sstream>
#include <iomanip>
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

static void verifyZVZeroCurve(
    ostream& os,
    bool isSofr,
    ::ext::shared_ptr<YieldTermStructure>& pYieldCurve,
    ::ext::shared_ptr<YieldTermStructure>& pYieldCurveAct365F
) {
    Handle<YieldTermStructure> hTS(pYieldCurve);
    Handle<YieldTermStructure> hTSAct365F(pYieldCurveAct365F);
    vector<Integer> swapTenorYears{ 2,3,4,5,7,10,12,15,20,25,30,40,50 };
    os << endl;
    os << "swap rates (30/360 vs Act/365F):" << endl;
    os << "========================================================================================" << endl;
    for (const auto& tenorYear : swapTenorYears) {
        auto tenor = tenorYear * Years;
        auto swap = ShortRateBlackCalibration::makeFwdSwap(isSofr, 0 * Months, tenor, hTS);
        auto swapRate = swap->fairRate();
        auto swapAct365F = ShortRateBlackCalibration::makeFwdSwap(isSofr, 0 * Months, tenor, hTSAct365F);
        auto swapRateAct365F = swapAct365F->fairRate();
        os << tenor;
        os << " [" << DateFormat<char>::to_yyyymmdd(swap->startDate()) << "," << DateFormat<char>::to_yyyymmdd(swap->maturityDate()) << "]";
        os << " ==> " << "30/360=" << io::percent(swapRate) << " vs Act/365F=" << io::percent(swapRateAct365F);
        os << ", diff=" << io::basis_point(swapRate - swapRateAct365F);
        os << endl;
    }
    os << "========================================================================================" << endl;
}

namespace po = boost::program_options;

int main(int argc, char* argv[], char* envp[]) {
    try {
        string asofdate;
        string s_zeroCurve;
        Real maturity = DEFAULT_MATURITY;
        Natural gridFreq = DEFAULT_GRID_FREQ;
        bool isSofr = DEFAULT_IS_SOFR;
        bool silent = DEFAULT_SILENT;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("asofdate,d", po::value<string>(&asofdate)->default_value(""), "(optional) evaluation date in yyyymmdd format. default to today's date")
            ("zero-curve,z", po::value<string>(&s_zeroCurve)->default_value(""), "(required) zero/spot curve input json file. \"@-\" to piped in from stdin")
            ("maturity", po::value<Real>(&maturity)->default_value(DEFAULT_MATURITY), "(optional) simulation duration in years. default to 50")
            ("grid-freq,f", po::value<Natural>(&gridFreq)->default_value(DEFAULT_GRID_FREQ), "(optional) grid frequency per year: 12=monthly|24=bi-monthly|36=10d|60=6d|72=5d|120=3d|180=2d|360=1d. The default value is 12 (monthly)")
            ("is-sofr", po::value<bool>(&isSofr)->default_value(DEFAULT_IS_SOFR), "(optional) the input zero curve is a SOFR zero curve. true or false. The default value is true")
            ("silent,s", po::value<bool>(&silent)->default_value(DEFAULT_SILENT), "(optional) silent running mode. true or false. the default is false")
            ;
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        QL_REQUIRE(!s_zeroCurve.empty(), "zero/spot curve input json file is not optional");
        QL_REQUIRE(maturity > 0., "maturity (" << maturity << ") must be a positive");
        QL_REQUIRE(gridFreq > 0, "grid frequency (" << gridFreq << ") must be a positive integer");
        auto evalDate = (asofdate.empty() ? Date::todaysDate() : DateFormat<char>::from_yyyymmdd(asofdate, false));

        auto& os = get_nullable_cout<char>(silent);

        os << endl;
        os << "ZV Rates Generator" << endl;
        os << "asofdate=" << DateFormat<char>::to_yyyymmdd(evalDate, true) << endl;
        os << "zero-curve=" << s_zeroCurve << endl;
        os << "maturity=" << maturity << endl;
        os << "grid-freq=" << gridFreq << endl;
        os << "is-sofr=" << (isSofr ? "true" : "false") << endl;
        os << "silent=" << (silent ? "true" : "false") << endl;

        Settings::instance().evaluationDate() = evalDate;
        auto pYieldCurve = MarketRate::readAsThirty360ZeroCurve(s_zeroCurve, evalDate, true);
        //auto pYieldCurveAct365F = MarketRate::readAsActual365FixedZeroCurve(s_zeroCurve, evalDate, true);
        //verifyZVZeroCurve(os, isSofr, pYieldCurve, pYieldCurveAct365F);

        const auto& curveConfig = getCurveConfigs().at(isSofr);
        const auto& parRateCouponFrequency = curveConfig.parRateCouponFrequency;
        const auto& parRateDayCountMultiplier = curveConfig.parRateDayCountMultiplier;
        auto pCashRateIndex = curveConfig.createCashRateIndex();
        auto cashRateIndexTenor = pCashRateIndex->tenor();
        auto cashRateIndexTenorUnits = cashRateIndexTenor.units();
        QL_REQUIRE(cashRateIndexTenorUnits == Years || cashRateIndexTenorUnits == Months, "cash rate index's time unit must be either Months or Years");
        auto cashRateTenorMonths = (cashRateIndexTenorUnits == Years ? cashRateIndexTenor.length() * 12 : cashRateIndexTenor.length());
        auto cashRateDayCounter = pCashRateIndex->dayCounter();
        auto n = boost::lexical_cast<size_t>(std::round(maturity * 12.));
        Time grid_delta_t = 1. / ((Real)gridFreq); // this is required for calculating short rate that is comparable to the irp process
        ShortRateZVCalculator shortRateCalculator(pYieldCurve);
        DiscountFactorZVCalculator dfCalculator(pYieldCurve);
        ParRateZVCalculator parRate5yrCalculator(5, parRateCouponFrequency, pYieldCurve, parRateDayCountMultiplier);
        ParRateZVCalculator parRate10yrCalculator(10, parRateCouponFrequency, pYieldCurve, parRateDayCountMultiplier);
        PathOutput output;
        auto& dfs = output.dfs;
        auto& shortRates = output.shortRates;
        auto& benchmarkCashRates = output.benchmarkCashRates;
        auto& benchmark5yrParRates = output.benchmark5yrParRates;
        auto& benchmark10yrParRates = output.benchmark10yrParRates;
        for (decltype(n) month = 0; month <= n; month++) {  // [0, n] => n+1 iterations/months projection
            auto t = (Time)month / 12.0;
            auto df = dfCalculator(t);
            dfs.push_back(df);
            auto shortRate = shortRateCalculator(t, grid_delta_t);
            shortRates.push_back(shortRate * 100.);
            if (month <= n - cashRateTenorMonths) {
                auto startDt = evalDate + Period(month, Months);
                auto endDt = evalDate + Period(month + cashRateTenorMonths, Months);
                auto dfStart = pYieldCurve->discount(startDt);
                auto dfEnd = pYieldCurve->discount(endDt);
                auto dt = cashRateDayCounter.yearFraction(startDt, endDt);
                auto df = dfEnd / dfStart;
                auto cashRate = ((1. / df) - 1.) / dt;
                benchmarkCashRates.push_back(cashRate * 100.);
            }
            if (month <= n - 5 * 12) {
                auto parRate = parRate5yrCalculator(t);
                benchmark5yrParRates.push_back(parRate * 100.);
            }
            if (month <= n - 10 * 12) {
                auto parRate = parRate10yrCalculator(t);
                benchmark10yrParRates.push_back(parRate * 100.);
            }
        }
        os << endl;
        cout << output << endl;
        os << endl;
        os << "Done" << endl;
        return 0;
    }
    catch (const std::exception& e) {
        cerr << "!!! Error: " << e.what() << endl;
        return 1;
    }
}
