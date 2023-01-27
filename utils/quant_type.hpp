#pragma once

#include <ql/quantlib.hpp>

namespace utils {
    struct OptionAttribs {
        QuantLib::Period expiry;
        QuantLib::Period tenor;
    };

    template <typename DT>
    struct OptionVolData : public OptionAttribs {
        DT data;    // "data" (vol or rate)
    };
}