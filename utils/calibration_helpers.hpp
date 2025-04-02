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
        typename SWAP_INDEX,    // swap index
        typename IBOR_INDEX     // floating leg ibor index
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
    struct SwaptionHelperFactory {
        typedef typename HelperSwapTraits::IborIndexType IborIndexType;
        // swaption helper factoring operator
        QuantLib::ext::shared_ptr<QuantLib::SwaptionHelper> operator ()(
            QuantLib::Handle<QuantLib::Quote>& quoteVol,    // quoted volality
            const QuantLib::Period& expiry, // expiry of the quoted volality
            const QuantLib::Period& tenor,  // tenor of the quoted volality
            const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure,    // term structure used to estimate ibor index and discounting cashflows
            QuantLib::VolatilityType quoteVolType = QuantLib::VolatilityType::ShiftedLognormal,  // volality type of the quoted volality (ShiftedLognormal or Normal)
            QuantLib::Real shift = 0.,
            QuantLib::Real strike = QuantLib::Null<QuantLib::Real>()
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
                strike,   // strike, Null<Real>() mean ATM, and quoteVol is ATM vol
                1.,  // nominal
                quoteVolType,    //  type
                shift,  // shift
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
            auto pSwaptionHelper = this->operator()(
                quoteHandle,
                forward,
                tenor,
                termStructure,
                QuantLib::VolatilityType::ShiftedLognormal, // volType
                0., // shift
                QuantLib::Null<QuantLib::Real>()    // strike @ATM
            );
            auto pSwap = pSwaptionHelper->underlying();
            return pSwap;
        }
        QuantLib::Rate atmRate(
            const QuantLib::Period& forward, // forward period of the swap
            const QuantLib::Period& tenor,  // tenor of the swap
            const QuantLib::Handle<QuantLib::YieldTermStructure>& termStructure    // term structure used to estimate ibor index and discounting cashflows
        ) const {
            auto swap = this->operator()(
                forward,
                tenor,
                termStructure
            );
            auto fairRate = swap->fairRate();
            return fairRate;
        }
    };

    // legacy usd libor 3m swaption
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    typedef SwaptionHelperSwapTraits<
        QuantLib::UsdLiborSwapIsdaFixAm,    // SWAP_INDEX
        QuantLib::USDLibor3M                // IBOR_INDEX
    > USDLibor3MSwapTraits;
    typedef SwaptionHelperFactory<USDLibor3MSwapTraits> USDLibor3MSwaptionHelperFactory;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // usd sofr swaption with annual cashflow exchange
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    typedef SwaptionHelperSwapTraits<
        QuantLib::UsdFwdOISVanillaSwapIndex<QuantLib::Sofr, QuantLib::Frequency::Annual>,           // SWAP_INDEX
        QuantLib::UsdOvernightCompoundedAverageIndex<QuantLib::Sofr, QuantLib::Frequency::Annual>   // IBOR_INDEX
    > USDSofrSwapTraits;
    typedef SwaptionHelperFactory<USDSofrSwapTraits> USDSofrSwaptionHelperFactory;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}