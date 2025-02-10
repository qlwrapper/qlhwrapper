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

#pragma once

#include <ql/quantlib.hpp>
#include <ql_utils/all.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <cmath>
#include <exception>
#include <utils/paths.hpp>

namespace utils {
    // per-path/zv output data
    class PathOutput {
    public:
        typedef QuantLib::Real Percent;
        typedef std::vector<Percent> RateVector;
        typedef std::size_t TenorMonth;
        typedef QLUtils::YieldTSNodes<TenorMonth, Percent> YieldTermStructureNodes;
        typedef QLUtils::SimpleParRateCalculator<QLUtils::RateUnit::Percent, QuantLib::Frequency::Semiannual> ParRateCalculator;
    public:
        std::vector<QuantLib::DiscountFactor> dfs;      // monthly discount factors
        RateVector shortRates;                          // monthly short rates (% unit)
        RateVector benchmarkCashRates;                  // monthly benchmark cash rates (% unit)
        RateVector benchmark5yrParRates;                // monthly benchmark 5yr par rates (% unit)
        RateVector benchmark10yrParRates;               // monthly benchmark 10yr par rates (% unit)
    private:
        mutable RateVector zeroRates_;                  // lazily calculated zero rates (% unit)
        mutable YieldTermStructureNodes parYieldCurve_; // lazily calculated par yield nodes (% unit)
    public:
        // returns monthly semi-annual zero rates (% unit)
        const RateVector& zeroRates() const {
            if (zeroRates_.empty()) {
                auto n = dfs.size();
                zeroRates_.resize(n);
                for (decltype(n) month = 1; month < n; ++month) {
                    auto const& df = dfs[month];
                    auto& zeroRate = zeroRates_[month];
                    auto t = ((QuantLib::Time)month) / 12.;
                    zeroRate = (std::pow(1. / df, 1. / (t * 2.)) - 1.) * 2. * 100.; // 1/df = (1 + zeroRate/2)^(t * 2)
                }
                zeroRates_[0] = zeroRates_[1];
            }
            return zeroRates_;
        }
        // returns (par) yield curve nodes (% unit)
        const YieldTermStructureNodes& parYieldCurve() const {
            if (parYieldCurve_.empty()) {
                const auto& monthlyZeroRates = zeroRates();
                auto n = monthlyZeroRates.size();
                parYieldCurve_.resize(n - 1);
                ParRateCalculator parRateCalculator(monthlyZeroRates);
                for (decltype(n) month = 1; month < n; ++month) {
                    auto& tenorMonth = parYieldCurve_.maturities[month-1];
                    auto& parYield = parYieldCurve_.rates[month-1];
                    tenorMonth = month;
                    parYield = parRateCalculator(tenorMonth, 0);
                }
            }
            return parYieldCurve_;
        }
        template <
            typename MT
        >
        static std::string maturitiesToString(
            const std::vector<MT>& v
        ) {
            std::ostringstream oss;
            auto n = v.size();
            for (decltype(n) i = 0; i < n; ++i) {
                if (i > 0) {
                    oss << " ";
                }
                oss << v[i];
            }
            return oss.str();
        }
        template<
            typename VT
        >
        static std::string vectorToString(
            const std::vector<VT>& v,
            std::streamsize precision = 6
        ) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(precision);
            auto n = v.size();
            for (decltype(n) i = 0; i < n; ++i) {
                if (i > 0) {
                    oss << " ";
                }
                oss << v[i];
            }
            return oss.str();
        }

        friend std::ostream& operator << (
            std::ostream& os,
            const PathOutput& rhs
        ) {
            auto dfsRamp = rhs.vectorToString(rhs.dfs, 10);
            auto shortRatesRamp = rhs.vectorToString(rhs.shortRates, 6);
            auto zeroRatesRamp = rhs.vectorToString(rhs.zeroRates(), 6);
            const auto& parYieldCurve = rhs.parYieldCurve();
            auto maturitiesRamp = rhs.maturitiesToString(parYieldCurve.maturities);
            auto parYieldsRamp = rhs.vectorToString(parYieldCurve.rates, 6);
            auto benchmarkCashRatesRamp = rhs.vectorToString(rhs.benchmarkCashRates, 6);
            auto benchmark5yrParRatesRamp = rhs.vectorToString(rhs.benchmark5yrParRates, 6);
            auto benchmark10yrParRatesRamp = rhs.vectorToString(rhs.benchmark10yrParRates, 6);
            std::ostringstream oss;
            oss << "{" << std::endl;
            oss << "\"dfs\":" << "\"" << dfsRamp << "\"" << std::endl;
            oss << "," << "\"shortRates\":" << "\"" << shortRatesRamp << "\"" << std::endl;
            oss << "," << "\"zeroRates\":" << "\"" << zeroRatesRamp << "\"" << std::endl;
            oss << "," << "\"yieldCurve\":" << "{" << std::endl;
            oss << "\"maturities\":" << "\"" << maturitiesRamp << "\"" << std::endl;
            oss << "," << "\"vals\":" << "\"" << parYieldsRamp << "\"" << std::endl;
            oss << "}" << std::endl;
            oss << "," << "\"benchmarkCashRates\":" << "\"" << benchmarkCashRatesRamp << "\"" << std::endl;
            oss << "," << "\"benchmark5yrParRates\":" << "\"" << benchmark5yrParRatesRamp << "\"" << std::endl;
            oss << "," << "\"benchmark10yrParRates\":" << "\"" << benchmark10yrParRatesRamp << "\"" << std::endl;
            oss << "}";
            return os << oss.str();
        }
    };

    struct IRPSimulationOutput {
        std::shared_ptr<Paths> pDFs;
        std::shared_ptr<Paths> pShortRates;
        std::shared_ptr<Paths> pBenchmarkCashRates;
        std::shared_ptr<Paths> pBenchmark5yrParRates;
        std::shared_ptr<Paths> pBenchmark10yrParRates;
        typedef size_t PathIndex;
        void verify() const {
            QL_REQUIRE(pDFs != nullptr, "discount factor paths is null");
            QL_REQUIRE(pShortRates != nullptr, "short rate paths is null");
            QL_REQUIRE(pBenchmarkCashRates != nullptr, "benchmark cash rate paths is null");
            QL_REQUIRE(pBenchmark5yrParRates != nullptr, "benchmark 5yr par rate paths is null");
            QL_REQUIRE(pBenchmark10yrParRates != nullptr, "benchmark 10yr par rate paths is null");
            auto nPaths = this->nPaths();
            QL_REQUIRE(nPaths > 0, "no irp path");
            QL_REQUIRE(pShortRates->nPaths() == nPaths, "incorrect number of short rate paths (" << pShortRates->nPaths() << "). expecting " << nPaths << ".");
            QL_REQUIRE(pBenchmarkCashRates->nPaths() == nPaths, "incorrect number of benchmark cash rate paths (" << pBenchmarkCashRates->nPaths() << "). expecting " << nPaths << ".");
            QL_REQUIRE(pBenchmark5yrParRates->nPaths() == nPaths, "incorrect number of benchmark 5yr par rate paths (" << pBenchmark5yrParRates->nPaths() << "). expecting " << nPaths << ".");
            QL_REQUIRE(pBenchmark10yrParRates->nPaths() == nPaths, "incorrect number of benchmark 10yr par rate paths (" << pBenchmark10yrParRates->nPaths() << "). expecting " << nPaths << ".");
        }
        size_t nPaths() const {
            return pDFs->nPaths();
        }
        static void copyVectorForPath(
            std::vector<QuantLib::Real>& dest,
            const Paths& scrPaths,
            const PathIndex& pathIndex,
            QuantLib::Real scalingMultiplier = 1.   // 100. = convert from decimal to % unit
        ) {
            const auto& srcMatrix = *(scrPaths.matrix());
            dest.resize(scrPaths.nTimes());
            std::copy(srcMatrix.column_begin(pathIndex), srcMatrix.column_end(pathIndex), dest.begin());
            if (scalingMultiplier != 1.) {
                auto scaler = [&scalingMultiplier](QuantLib::Real& value) {value *= scalingMultiplier; };
                std::for_each(dest.begin(), dest.end(), scaler);
            }
        }
        // path operator
        std::shared_ptr<PathOutput> operator [] (
            const PathIndex& pathIndex
        ) const {
            QL_ASSERT(pathIndex >= 0 && pathIndex < nPaths(), "path index (" << pathIndex << ") is out of range. range=[0, " << nPaths() << ")");
            std::shared_ptr<PathOutput> ret(new PathOutput());
            auto& pathOutput = *ret;
            // dfs
            copyVectorForPath(pathOutput.dfs, *pDFs, pathIndex, 1.);
            // short rates
            copyVectorForPath(pathOutput.shortRates, *pShortRates, pathIndex, 100.);
            // benchmark cash rates
            copyVectorForPath(pathOutput.benchmarkCashRates, *pBenchmarkCashRates, pathIndex, 100.);
            // benchmark 5yr par rates
            copyVectorForPath(pathOutput.benchmark5yrParRates, *pBenchmark5yrParRates, pathIndex, 100.);
            // benchmark 10yr par rates
            copyVectorForPath(pathOutput.benchmark10yrParRates, *pBenchmark10yrParRates, pathIndex, 100.);
            return ret;
        }
        friend std::ostream& operator<< (
            std::ostream& os,
            const IRPSimulationOutput& rhs
        ) {
            rhs.verify();
            auto nPaths = rhs.nPaths();
            typedef size_t PathIndex;
            typedef std::pair<PathIndex, const PathOutput&> PathItem;
            auto getPathItemJSON = [](const PathItem& item) {
                std::ostringstream oss;
                oss << "{" << std::endl;
                oss << "\"pathIndex\":" << std::fixed << item.first << std::endl;
                oss << "," << "\"pathOutput\":" << item.second << std::endl;
                oss << "}";
                return oss.str();
            };
            std::ostringstream oss;
            oss << "[" << std::endl;
            for (decltype(nPaths) path = 0; path < nPaths; ++path) {  // for each path
                try {
                    auto pPathOutput = rhs[path];
                    QL_ASSERT(pPathOutput != nullptr, "path output is null");
                    PathItem pathItem(path, *pPathOutput);
                    if (path > 0) {
                        oss << ",";
                    }
                    auto pathItemJSON = getPathItemJSON(pathItem);
                    oss << pathItemJSON << std::endl;
                }
                catch (const std::exception& e) {
                    QL_FAIL("irp path " << path << " error: " << e.what());
                }
            }
            oss << "]";
            return os << oss.str();
        }
    };
}