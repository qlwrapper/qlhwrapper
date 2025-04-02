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

#include <pch.h>
#include <utils/market_rate.hpp>
#include <utils/misc.h>
#include <utils/json_io.h>
#include <utils/io.hpp>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <exception>
#include <common/json-obj-serialization.hpp>

using namespace std;
using namespace utils;
using namespace QuantLib;
using namespace QLUtils;

namespace {
    boost::property_tree::ptree loadTreeFromJSONInput(
        const string& inputArg
    ) {
        auto json = get_input_text<char>(inputArg);
        auto ptree = json_io::parse_json<char>(json);
        return ptree;
    }
    pair<vector<Date>, vector<Rate>> convertToThirty360TermVectors(
        const Date& baseDate,
        const vector<Time>& act365FixedTimes,
        const vector<Rate>& zeroRates
    ) {
        auto n = act365FixedTimes.size();
        QL_REQUIRE(zeroRates.size() == n, "number of zero rates (" << zeroRates.size() << ") is not what's expected (" << n << ")");
        QL_REQUIRE(n > 1, "too few term structure nodes (" << n << "). minimum is 2");
        QL_REQUIRE(act365FixedTimes[0] == 0., "first term structure node must have a term (" << act365FixedTimes[0] << ") of 0.");
        DayCounter dcAct365F = Actual365Fixed();
        DayCounter dc30360 = Thirty360(Thirty360::BondBasis);
        using Thirty360Time = Time;
        using Thirty360Rate = Rate;
        map<Thirty360Time, pair<Date, Thirty360Rate>> m;    // map from 30/360 term to date/zero rate pair. this to ensure the result 30/360 curve nodes have unique term
        LinearInterpolation interp(act365FixedTimes.begin(), act365FixedTimes.end(), zeroRates.begin());
        auto dtStart = baseDate + 1;
        auto dtLast = baseDate + (Integer)(std::round(act365FixedTimes[n-1] * 365.));
        for (decltype(dtStart) date = dtStart; date <= dtLast; ++date) {
            auto t_Act365F = dcAct365F.yearFraction(baseDate, date);
            auto zeroRate_Act365F = interp(t_Act365F, true);
            auto df = std::exp(-zeroRate_Act365F * t_Act365F);
            Thirty360Time t_30360 = dc30360.yearFraction(baseDate, date);
            if (t_30360 != 0.) {
                Thirty360Rate zeroRate_30360 = -std::log(df) / t_30360;
                m[t_30360] = pair<Date, Thirty360Rate>(date, zeroRate_30360);
            }
        }
        n = m.size();
        vector<pair<Date, Rate>> v(n);
        decltype(n) i = 0;
        for (const auto& pr : m) {
            v[i++] = pr.second;
        }
        std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {return a.first < b.first; });
        auto firstRate = v[0].second;
        v.insert(v.begin(), pair<Date, Rate> {baseDate, firstRate});
        n = v.size();
        pair<vector<Date>, vector<Rate>> ret;
        auto& dates = ret.first;
        auto& rates = ret.second;
        dates.resize(n);
        rates.resize(n);
        for (decltype(n) i = 0; i < n; ++i) {
            dates[i] = v[i].first;
            rates[i] = v[i].second;
        }
        return ret;
    }

    ::ext::shared_ptr<YieldTermStructure> constuctZeroCurve(
        const Date& curveRefDate,
        const vector<Date>& dates,
        const vector<Rate>& rates,
        const DayCounter& dayCounter,
        bool enableExtrapolation
    ) {
        ext::shared_ptr<YieldTermStructure> pTS(new InterpolatedZeroCurve<Linear>(
            dates,
            rates,
            dayCounter
        ));
        QL_ASSERT(pTS->referenceDate() == curveRefDate, "curve ref date (" << pTS->referenceDate() << ") is not what's expected (" << curveRefDate << ")");
        if (enableExtrapolation) {
            pTS->enableExtrapolation(enableExtrapolation);
        }
        return pTS;
    }
    void assertValidQuotedSwaptionVols(
        const pQuotedSwaptionVols& pVols
    ) {
        QL_REQUIRE(pVols != nullptr && !pVols->empty(), "quoted vols object is null or empty");
        const auto& vols = *pVols;
        auto n = vols.size();
        for (decltype(n) i = 0; i < n; ++i) {
            try {
                const auto& pVol = vols[i];
                QL_REQUIRE(pVol != nullptr, "quoted vol object is null");
            }
            catch (const std::exception& e) {
                QL_FAIL("vols[" << i << "]: " << e.what());
            }
        }
    }
}

namespace utils {
    ::ext::shared_ptr<YieldTermStructure> MarketRate::getActual365FixedZeroTermStructure(
        const Date& baseDate,
        const vector<Time>& act365FixedTimes,
        const vector<Rate>& zeroRates,
        bool enableExtrapolation
    ) {
        auto n = act365FixedTimes.size();
        QL_REQUIRE(zeroRates.size() == n, "number of zero rates (" << zeroRates.size() << ") is not what's expected (" << n << ")");
        vector<Date> dates(n);
        for (decltype(n) i = 0; i < n; ++i) {
            const auto& t_Act365Fixed = act365FixedTimes[i];
            auto days = (Integer)(std::round(t_Act365Fixed * 365.));
            auto date = baseDate + days * QuantLib::Days;
            dates[i] = date;
        }
        return constuctZeroCurve(baseDate, dates, zeroRates, Actual365Fixed(), enableExtrapolation);
    }

    ext::shared_ptr<YieldTermStructure> MarketRate::getThirty360ZeroTermStructure(
        const Date& baseDate,
        const vector<Time>& act365FixedTimes,
        const vector<Rate>& zeroRates,
        bool enableExtrapolation
    ) {
        auto pr = convertToThirty360TermVectors(baseDate, act365FixedTimes, zeroRates);
        const auto& dates = pr.first;
        const auto& rates = pr.second;
        return constuctZeroCurve(baseDate, dates, rates, Thirty360(Thirty360::BondBasis), enableExtrapolation);
    }

    MarketRate::pTSNodes MarketRate::readAct365FixedZeroCurveTS(
        const string& inputArg
    ) {
        auto ptree = loadTreeFromJSONInput(inputArg);
        pTermStructureNodes<QuantLib::Time, QuantLib::Rate> pNodes;
        pNodes << ptree;
        MarketRate::pTSNodes pTSNodes = *(pNodes);
        return pTSNodes;
    }

    ext::shared_ptr<YieldTermStructure> MarketRate::readAsThirty360ZeroCurve(
        const string& inputArg,
        const Date& curveRefDate,
        bool enableExtrapolation
    ) {
        auto pTSNodes = readAct365FixedZeroCurveTS(inputArg);
        const auto& act365FixedTimes = pTSNodes->maturities;
        const auto& zeroRates = pTSNodes->rates;
        auto pYieldCurve = getThirty360ZeroTermStructure(curveRefDate, act365FixedTimes, zeroRates, enableExtrapolation);
        return pYieldCurve;
    }

    ext::shared_ptr<YieldTermStructure> MarketRate::readAsActual365FixedZeroCurve(
        const string& inputArg,
        const Date& curveRefDate,
        bool enableExtrapolation
    ) {
        auto pTSNodes = readAct365FixedZeroCurveTS(inputArg);
        const auto& act365FixedTimes = pTSNodes->maturities;
        const auto& zeroRates = pTSNodes->rates;
        auto pYieldCurve = getActual365FixedZeroTermStructure(curveRefDate, act365FixedTimes, zeroRates, enableExtrapolation);
        return pYieldCurve;
    }

    pQuotedSwaptionVols MarketRate::readQuotedSwaptionVols(
        const string& inputArg
    ) {
        auto ptree = loadTreeFromJSONInput(inputArg);
        pQuotedSwaptionVols pVols;
        pVols << ptree;
        assertValidQuotedSwaptionVols(pVols);
        return pVols;
    }
}