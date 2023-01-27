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

#define ATTR_SET ".<xmlattr>"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <cstdlib>
#include <sstream>
#include <memory>
#include <exception>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <ql/quantlib.hpp>
#include <utils/quant_type.hpp>
#include <utils/misc.h>

namespace utils {

    enum SwaptionVolDataType {
        LogNormalVol = 0,
        NormalVol = 1,
        ForwardSwapRate = 2
    };

    template <SwaptionVolDataType SVDT>
    struct ATMSwaptionVolDataStorageSpec {
        typedef QuantLib::Real data_type;
        std::string volNodePath() const {
            throw std::exception("volNodePath() must be template specialized");
        }
        QuantLib::Real valueScale() const {
            throw std::exception("valueScale() must be template specialized");
        }
    };

    template<>
    struct ATMSwaptionVolDataStorageSpec<LogNormalVol> {
        typedef QuantLib::Volatility data_type;
        std::string volNodePath() const {
            return "MarketRates.SwaptionVolCube.Slice";
        }
        QuantLib::Real valueScale() const {
            return 100.0;
        }
    };
   
    template<>
    struct ATMSwaptionVolDataStorageSpec<NormalVol> {
        typedef QuantLib::Volatility data_type;
        std::string volNodePath() const {
            return "MarketRates.SwaptionVolCubeNormal.Slice";
        }
        QuantLib::Real valueScale() const {
            return 10000.0;
        }
    };

    template<>
    struct ATMSwaptionVolDataStorageSpec<ForwardSwapRate> {
        typedef QuantLib::Rate data_type;
        std::string volNodePath() const {
            return "MarketRates.SwaptionATMSwap";
        }
        QuantLib::Real valueScale() const {
            return 100.0;
        }
    };

    class MarketRate {
    private:
        QuantLib::Date date_;
        boost::property_tree::ptree pt_root_;

        static const boost::property_tree::ptree& empty_ptree() {
            static boost::property_tree::ptree t;
            return t;
        }
        static std::string defaultMRXMLFolder() {
            return "\\\\cerbdata09\\polypath_data$\\Data\\daily\\alib";
        }
        static std::string defaultMRXMLFilename(const QuantLib::Date& date) {
            std::ostringstream os;
            os << utils::DateHelper<char>::to_yyyymmdd(date) << "_unch_1.xml";
            return os.str();
        }
    public:
        static std::string defaultMRXMLFilePath(const QuantLib::Date& date) {
            std::ostringstream os;
            os << defaultMRXMLFolder() << "\\" << defaultMRXMLFilename(date);
            return os.str();
        }
        MarketRate(const QuantLib::Date& date) : date_(date) {}
        const QuantLib::Date& date() const {
            return date_;
        }
        const boost::property_tree::ptree& root() const {
            return pt_root_;
        }
        boost::property_tree::ptree& root() {
            return pt_root_;
        }
        MarketRate& load(const std::string& mr_xml = "") {
            auto filename = (mr_xml == "" ? defaultMRXMLFilePath(date()) : mr_xml);
            boost::property_tree::read_xml(filename, pt_root_);
            return *this;
        }
    private:
        void getZVDiscountFactors(std::vector<QuantLib::Date>& dates, std::vector<QuantLib::DiscountFactor>& dfs) const {
            BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, root().get_child("MarketRates")) {
                boost::property_tree::ptree subtree = node.second;
                if (node.first == "TermStructure" && node.second.get_child("<xmlattr>.Name").data() == "SWAP_CURVE") {
                    BOOST_FOREACH(boost::property_tree::ptree::value_type const& v, subtree.get_child("")) {
                        if (v.first == "TS") {
                            const boost::property_tree::ptree& attributes = v.second.get_child("<xmlattr>", empty_ptree());
                            QuantLib::Time term_act_365f = 0.0;
                            QuantLib::Rate zeroRateContinuous = 0.0;
                            BOOST_FOREACH(boost::property_tree::ptree::value_type const& elm, attributes) {
                                if (elm.first == "Term") {
                                    term_act_365f = stod(elm.second.data());
                                }
                                else if (elm.first == "Rate") {
                                    zeroRateContinuous = stod(elm.second.data()) / 100.0;
                                }
                            }
                            auto days = (QuantLib::Natural)std::round(term_act_365f * 365.0);
                            auto dt = date() + days * QuantLib::Days;
                            dates.push_back(dt);
                            QuantLib::DiscountFactor df = std::exp(-zeroRateContinuous * term_act_365f);
                            dfs.push_back(df);
                        }
                    }
                }
            }
        }
        static QuantLib::Period term_30_360_to_period(QuantLib::Time term_30_360) {
            auto days = (QuantLib::Natural)std::round(term_30_360 * 360.0);
            if (days == 0) {
                return QuantLib::Period(0, QuantLib::Days);
            }
            else if (days % 360 == 0) {
                return QuantLib::Period(days / 360, QuantLib::Years);
            }
            else if (days % 30 == 0) {
                return QuantLib::Period(days / 30, QuantLib::Months);
            }
            else if (days % 7 == 0) {
                return QuantLib::Period(days / 7, QuantLib::Weeks);
            }
            else {
                return QuantLib::Period(days, QuantLib::Days);
            }
        }
        static std::tuple<std::vector<QuantLib::Time>, std::vector<QuantLib::Date>, std::vector<QuantLib::DiscountFactor>> ensureUnique30360Times(const QuantLib::Date& baseDate, const std::vector<QuantLib::Date>& dates, std::vector<QuantLib::DiscountFactor> dfs) {
            QuantLib::DayCounter dayCounter30360 = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis);
            using TupleItem = std::tuple<QuantLib::Time, QuantLib::Date, QuantLib::DiscountFactor>;
            std::map<QuantLib::Time, TupleItem> map_unique_30_360_time;
            auto itDF = dfs.begin();
            for (auto p = dates.begin(); p != dates.end(); ++p, ++itDF) {
                auto const& date = *p;
                auto time_30_360 = dayCounter30360.yearFraction(baseDate, date);
                if (map_unique_30_360_time.find(time_30_360) == map_unique_30_360_time.end()) {
                    map_unique_30_360_time[time_30_360] = TupleItem(time_30_360, date, *itDF);
                }
            }
            std::vector<TupleItem> v;
            for (auto p = map_unique_30_360_time.begin(); p != map_unique_30_360_time.end(); ++p) {
                auto tp = p->second;
                v.push_back(tp);
            }
            std::sort(v.begin(), v.end(), [](const TupleItem& a, const TupleItem& b) -> bool {
                auto const& time_a = std::get<0>(a);
                auto const& time_b = std::get<0>(b);
                return (time_a < time_b);
            });
            if (std::get<0>(v.front()) != 0.0) {
                v.insert(v.begin(), TupleItem(0.0, baseDate, 1.0));
            }
            auto n = v.size();
            std::tuple<std::vector<QuantLib::Time>, std::vector<QuantLib::Date>, std::vector<QuantLib::DiscountFactor>> ret;
            auto& times = std::get<0>(ret);
            auto& dts = std::get<1>(ret);
            auto& discountFactors = std::get<2>(ret);
            times.resize(n);
            dts.resize(n);
            discountFactors.resize(n);
            for (decltype(n) i = 0; i < n; i++) {
                auto const& tp = v.at(i);
                times[i] = std::get<0>(tp);
                dts[i] = std::get<1>(tp);
                discountFactors[i] = std::get<2>(tp);
            }
            return ret;
        }
    public:
        template<typename DT>
        static std::shared_ptr<std::vector<OptionVolData<DT>>> filterOptionVolData(const std::vector<OptionVolData<DT>>& data, const std::vector<OptionAttribs>& filter) {
            std::shared_ptr<std::vector<OptionVolData<DT>>> ret(new std::vector<OptionVolData<DT>>());
            for (auto p = filter.begin(); p != filter.end(); ++p) {
                auto found = false;
                for (auto it = data.begin(); it != data.end(); ++it) {
                    if (p->expiry == it->expiry && p->tenor == it->tenor) {
                        ret->push_back(*it);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    QL_FAIL("cannot find data for option " << p->expiry << "x" << p->tenor);
                }
            }
            return ret;
        }
        template<SwaptionVolDataType SVDT> using SWAPTION_DT = typename ATMSwaptionVolDataStorageSpec<SVDT>::data_type;
        template<SwaptionVolDataType SVDT>
        std::shared_ptr<std::vector<OptionVolData<SWAPTION_DT<SVDT>>>> getATMSwaptionVolData() const {
            using DT = SWAPTION_DT<SVDT>;
            ATMSwaptionVolDataStorageSpec<SVDT> storageSpec;
            std::shared_ptr<std::vector<OptionVolData<SWAPTION_DT<SVDT>>>> ret(new std::vector<OptionVolData<SWAPTION_DT<SVDT>>>());
            BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, root().get_child(storageSpec.volNodePath())) {
                if (node.first == "Row") {
                    boost::property_tree::ptree subtree = node.second;
                    BOOST_FOREACH(boost::property_tree::ptree::value_type const& v, subtree.get_child("")) {
                        bool loaded = true;
                        if (v.first == "Col" && loaded) {
                            OptionVolData<DT> data;
                            data.expiry = term_30_360_to_period((QuantLib::Time)stod(node.second.get_child("<xmlattr>.Term").data()));
                            const boost::property_tree::ptree& attributes = v.second.get_child("<xmlattr>", empty_ptree());
                            BOOST_FOREACH(boost::property_tree::ptree::value_type const& elm, attributes) {
                                if (elm.first == "Term") {
                                    data.tenor = term_30_360_to_period((QuantLib::Time)stod(elm.second.data()));
                                }
                                else if (elm.first == "Rate") {
                                    data.data = (DT)(stod(elm.second.data()) / storageSpec.valueScale());
                                    loaded = false;
                                    ret->push_back(data);
                                }
                            }
                        }
                    }
                }
            }
            return ret;
        }
        QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> getThirty360ZVDiscountCurve(const QuantLib::Date& curveRefDate) const {
            auto baseDate = date();
            QL_ASSERT(curveRefDate >= baseDate, "curveRefDate must be greater of equalt to the base date");
            std::vector<QuantLib::Date> dates_;
            std::vector<QuantLib::DiscountFactor> dfs_;
            getZVDiscountFactors(dates_, dfs_);
            QL_ASSERT(dates_.size() == dfs_.size(), "bad curve data");
            auto res = ensureUnique30360Times(baseDate, dates_, dfs_);
            // re-base the dates and discount factors to base on curveRefDate
            /////////////////////////////////////////////////////////////////////////////////////////////////
            auto const& times = std::get<0>(res);
            auto const& dts = std::get<1>(res);
            auto const& discountFactors = std::get<2>(res);
            QuantLib::DayCounter dayCounter30360 = QuantLib::Thirty360(QuantLib::Thirty360::BondBasis);
            auto t_curveRefDate = dayCounter30360.yearFraction(baseDate, curveRefDate);
            QuantLib::Linear linear;
            auto interp = linear.interpolate(times.begin(), times.end(), discountFactors.begin());
            auto df_curveRefDate = interp(t_curveRefDate, true);
            std::vector<QuantLib::Date> datesRebased;
            std::vector<QuantLib::DiscountFactor> discountFactorsRebased;
            auto it = discountFactors.begin();
            for (auto p = dts.begin(); p != dts.end(); ++p, ++it) {
                auto const& date = *p;
                auto const& df = *it;
                if (date >= curveRefDate) {
                    datesRebased.push_back(date);
                    discountFactorsRebased.push_back(df / df_curveRefDate);
                }
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////
            res = ensureUnique30360Times(curveRefDate, datesRebased, discountFactorsRebased);
            auto const& dates = std::get<1>(res);
            auto const& dfs = std::get<2>(res);
            QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure> yc(new QuantLib::InterpolatedDiscountCurve<QuantLib::Linear>(dates, dfs, dayCounter30360));
            QL_ASSERT(yc->referenceDate() == curveRefDate, "curve's reference is not what's expected");
            yc->enableExtrapolation();
            return yc;
        }

        static QuantLib::Rate usdSwapRate(const QuantLib::Period& swapTenor, const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& iborIndex, const QuantLib::Handle<QuantLib::YieldTermStructure>& discountCurve) {
            QuantLib::ext::shared_ptr<QuantLib::PricingEngine> swapPricingEngine(new QuantLib::DiscountingSwapEngine(discountCurve));
            QuantLib::ext::shared_ptr<QuantLib::VanillaSwap> swap = QuantLib::MakeVanillaSwap(swapTenor, iborIndex)
                .withSettlementDays(2)
                .withFixedLegCalendar(QuantLib::UnitedStates(QuantLib::UnitedStates::GovernmentBond))
                .withFixedLegTenor(QuantLib::Period(6, QuantLib::Months))
                .withFixedLegConvention(QuantLib::ModifiedFollowing)
                .withFixedLegDayCount(QuantLib::Thirty360(QuantLib::Thirty360::BondBasis))
                .withFixedLegEndOfMonth(false)
                .withFloatingLegCalendar(QuantLib::UnitedStates(QuantLib::UnitedStates::GovernmentBond))
                .withFloatingLegEndOfMonth(false)
                .withDiscountingTermStructure(discountCurve);
            swap->setPricingEngine(swapPricingEngine);
            auto swapRate = swap->fairRate();
            return swapRate;
        }
    };
}