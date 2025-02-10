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

namespace utils {
    // swap traits for ATMSwaptionHelperFactory
    template <
        typename SWAP_INDEX,
        typename IBOR_INDEX
    >
    struct SwaptionHelperSwapTraits {
        typedef typename IBOR_INDEX IborIndexType;
        QuantLib::Period fixedLegTenor(
            const QuantLib::Period& tenor
        ) const {
            SWAP_INDEX swapIndex(tenor);
            return swapIndex.fixedLegTenor();
        }
        QuantLib::DayCounter fixedLegDayCounter(
            const QuantLib::Period& tenor
        ) const {
            SWAP_INDEX swapIndex(tenor);
            return swapIndex.dayCounter();
        }
        QuantLib::Natural settlementDays(
            const QuantLib::Period& tenor
        ) const {
            SWAP_INDEX swapIndex(tenor);
            return swapIndex.fixingDays();
        }
    };

    // ATM SwaptionHelper factory object for a given helper swap traits object 
    template <
        typename HelperSwapTraits
    >
    struct ATMSwaptionHelperFactory {
        typedef typename HelperSwapTraits::IborIndexType IborIndexType;
        // swaption helper factoring operator
        QuantLib::ext::shared_ptr<QuantLib::SwaptionHelper> operator ()(
            QuantLib::Handle<QuantLib::Quote>& quoteVol,    // quoted ATM volality
            const QuantLib::Period& expiry, // expiry of the quoted ATM volality
            const QuantLib::Period& tenor,  // tenor of the quoted ATM volality
            const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,    // term structure used to estimate ibor index and discounting cashflows
            QuantLib::VolatilityType quoteVolType = QuantLib::VolatilityType::ShiftedLognormal  // volality type of the quoted ATM volality (ShiftedLognormal or Normal)
        ) const {
            HelperSwapTraits swapTraits;
            QuantLib::ext::shared_ptr<QuantLib::IborIndex> index(new IborIndexType(termStructure)); // use the same term structure to estimate the ibor index
            QuantLib::ext::shared_ptr<QuantLib::SwaptionHelper> helper(new QuantLib::SwaptionHelper(
                expiry, // maturity
                tenor,  // length
                quoteVol,   // volatility
                index,   // index of type IborIndex
                swapTraits.fixedLegTenor(tenor),    // fixedLegTenor
                swapTraits.fixedLegDayCounter(tenor),   // fixedLegDayCounter
                index->dayCounter(),  // floatingLegDayCounter
                termStructure,  // termStructure
                QuantLib::BlackCalibrationHelper::RelativePriceError,  // errorType
                QuantLib::Null<QuantLib::Real>(),   // strike, Null<Real>() mean ATM, and quoteVol is ATM vol
                1.,  // nominal
                quoteVolType,    //  type
                0.,  // shift
                swapTraits.settlementDays(tenor)    // settlementDays
            ));
            return helper;
        }
        // fwd swap factoring operator, can be used to calculate ATM strike
        QuantLib::ext::shared_ptr<QuantLib::FixedVsFloatingSwap> operator () (
            const QuantLib::Period& forward, // forward period of the swap
            const QuantLib::Period& tenor,  // tenor of the swap
            const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure    // term structure used to estimate ibor index and discounting cashflows
        ) const {
            QuantLib::Handle<QuantLib::Quote> quoteHandle(QuantLib::ext::shared_ptr<QuantLib::Quote>(new QuantLib::SimpleQuote(0.)));
            auto pSwaptionHelper = this->operator()(quoteHandle, forward, tenor, termStructure, QuantLib::VolatilityType::ShiftedLognormal);
            auto pSwap = pSwaptionHelper->underlying();
            return pSwap;
        }
    };

    typedef SwaptionHelperSwapTraits<QuantLib::UsdLiborSwapIsdaFixAm, QuantLib::USDLibor3M> USDLibor3MSwapTraits;
    typedef ATMSwaptionHelperFactory<USDLibor3MSwapTraits> USDLibor3MATMSwaptionHelperFactory;

    typedef SwaptionHelperSwapTraits<QuantLib::UsdFwdOISVanillaSwapIndex<QuantLib::Sofr>, QuantLib::UsdOvernightCompoundedAverageIndex<QuantLib::Sofr>> USDSofrSwapTraits;
    //typedef SwaptionHelperSwapTraits<QuantLib::UsdOvernightIndexedSwapIsdaFix<QuantLib::Sofr>, QuantLib::Sofr> USDSofrSwapTraits;
    typedef ATMSwaptionHelperFactory<USDSofrSwapTraits> USDSofrATMSwaptionHelperFactory;
}