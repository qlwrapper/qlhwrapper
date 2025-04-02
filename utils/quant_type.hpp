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
#include <vector>
#include <memory>
#include <exception>
#include <sstream>
#include <string>

namespace utils {
    struct OptionAttribs {
        QuantLib::Period expiry;
        OptionAttribs(
            const QuantLib::Period& expiry = QuantLib::Period()
        ) : expiry(expiry)
        {}
        virtual std::string toString() const {
            return (std::ostringstream() << expiry).str();
        }
    };

    struct SwaptionAttribs : public OptionAttribs {
        QuantLib::Period tenor; // swap tenor
        SwaptionAttribs(
            const QuantLib::Period& expiry = QuantLib::Period(),
            const QuantLib::Period& tenor = QuantLib::Period()
        ): OptionAttribs(expiry), tenor(tenor)
        {}
        std::string toString() const {
            return (std::ostringstream() << expiry << "x" << tenor).str();
        }
    };

    typedef QuantLib::Real Skew;    // strike - forward (K - F)

    struct VolatilityData {
        QuantLib::Volatility vol;          // volatility/sigma
        QuantLib::VolatilityType volType;   // volatility type
        QuantLib::Real shift;
        Skew skew;                          // strike - forward (K - F)
        VolatilityData(
            QuantLib::Volatility vol = 0.,
            QuantLib::VolatilityType volType = QuantLib::VolatilityType::ShiftedLognormal,
            QuantLib::Real shift = 0.,
            Skew skew = 0.
        ) : vol(vol), volType(volType), shift(shift), skew(skew)
        {}
        bool hasLognormalShift() const {
            return (volType == QuantLib::VolatilityType::ShiftedLognormal && shift != 0.);
        }
        bool isATM() const {
            return (skew == 0.);
        }
        bool isSkewed() const {
            return (skew != 0.);
        }
    };

    struct QuotedSwaptionVol :
        public SwaptionAttribs,
        public VolatilityData {
        QuotedSwaptionVol(
            QuantLib::Volatility vol = 0.,
            QuantLib::VolatilityType volType = QuantLib::VolatilityType::ShiftedLognormal,
            QuantLib::Real shift = 0.,
            Skew skew = 0.
        ) : VolatilityData(vol, volType, shift, skew)
        {}
    };

    typedef std::shared_ptr<QuotedSwaptionVol> pQuotedSwaptionVol;
    typedef std::vector<pQuotedSwaptionVol> QuotedSwaptionVols;
    typedef std::shared_ptr<QuotedSwaptionVols> pQuotedSwaptionVols;
}