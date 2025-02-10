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

#include <utils/enum-lookup.hpp>
#include <utils/model_type.hpp>
#include <ql/quantlib.hpp>

namespace utils {
    // ModelType_Lookup
    BEGIN_DECLARE_ENUM_LOOKUP(ModelType)
        ENUM_LOOKUP_ITEM(ModelType::HullWhite1F, "Hull-White 1F", true);
        ENUM_LOOKUP_ITEM(ModelType::HullWhite1F, "hw", false);
        ENUM_LOOKUP_ITEM(ModelType::GeneralizedHullWhite1F, "Generalized Hull-White 1F", true);
        ENUM_LOOKUP_ITEM(ModelType::GeneralizedHullWhite1F, "ghw", false);
        ENUM_LOOKUP_ITEM(ModelType::GeneralizedHullWhiteLogNormal1F, "Generalized Hull-White Log Normal 1F", true);
        ENUM_LOOKUP_ITEM(ModelType::GeneralizedHullWhiteLogNormal1F, "ghwln", false);
        ENUM_LOOKUP_ITEM(ModelType::G2, "Two-factor Additive Gaussian", true);
        ENUM_LOOKUP_ITEM(ModelType::G2, "g2", false);
    END_DECLARE_ENUM_LOOKUP()

    // TimeUnit_Lookup
    BEGIN_DECLARE_ENUM_LOOKUP_NS(QuantLib, TimeUnit)
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Days, "D", true);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Days, "Days", false);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Weeks, "W", true);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Weeks, "Weeks", false);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Months, "M", true);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Months, "Months", false);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Years, "Y", true);
        ENUM_LOOKUP_ITEM(QuantLib::TimeUnit::Years, "Years", false);
    END_DECLARE_ENUM_LOOKUP()

    // VolatilityType_Lookup
    BEGIN_DECLARE_ENUM_LOOKUP_NS(QuantLib, VolatilityType)
        ENUM_LOOKUP_ITEM(QuantLib::VolatilityType::Normal, "Normal", true);
        ENUM_LOOKUP_ITEM(QuantLib::VolatilityType::ShiftedLognormal, "ShiftedLognormal", true);
    END_DECLARE_ENUM_LOOKUP()
}

namespace boost {
    ENUM_BOOST_LEXICAL_CAST_SUPPORT(utils, ModelType, utils)
    ENUM_BOOST_LEXICAL_CAST_SUPPORT(QuantLib, TimeUnit, utils)
    ENUM_BOOST_LEXICAL_CAST_SUPPORT(QuantLib, VolatilityType, utils)
}