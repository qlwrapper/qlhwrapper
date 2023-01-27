#pragma once

#include <ql/quantlib.hpp>

namespace utils {
    struct USDSwaptionHelperFactory {
        QuantLib::ext::shared_ptr<QuantLib::SwaptionHelper> operator() (
            QuantLib::Handle<QuantLib::Quote>& quote,
            const QuantLib::Period& expiry,
            const QuantLib::Period& tenor,
            const QuantLib::ext::shared_ptr<QuantLib::IborIndex>& iborIndex,
            const QuantLib::Handle<QuantLib::YieldTermStructure>& curveHandle
            ) const {
            QuantLib::ext::shared_ptr<QuantLib::SwaptionHelper> helper(new QuantLib::SwaptionHelper(
                expiry,
                tenor,
                quote,
                iborIndex,
                QuantLib::Period(6, QuantLib::Months),
                QuantLib::Thirty360(QuantLib::Thirty360::BondBasis),
                QuantLib::Actual360(),
                curveHandle
            ));
            return helper;
        }
    };
}