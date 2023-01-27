#pragma once

#include <map>
#include <string>
#include <utils/types.h>
#include <utils/string.h>
#include <utils/enum_helper.h>

namespace utils {
    enum ModelType {
        HullWhite1F = 0,
        GeneralizedHullWhite1F = 1,
        GeneralizedHullWhiteLogNormal1F = 2,
    };

    struct ModelTypeMapLoader {
        void operator() (std::map<ModelType, std::string>& map) const {
#define INSERT_ELEMENT(p, s) map[p] = s
            INSERT_ELEMENT(ModelType::HullWhite1F, "hw");
            INSERT_ELEMENT(ModelType::GeneralizedHullWhite1F, "ghw");
            INSERT_ELEMENT(ModelType::GeneralizedHullWhiteLogNormal1F, "ghwln");
#undef INSERT_ELEMENT
        }
    };
    typedef EnumHelper<ModelType, ModelTypeMapLoader> ModelTypeLookup;
}

template<typename _Elem>
inline utils::ostream_type<_Elem>& operator << (utils::ostream_type<_Elem>& os, const utils::ModelType& modelType) {
    utils::utf8_wstring_converter<_Elem> utf8_converter;
    switch (modelType) {
    case utils::ModelType::HullWhite1F:
        os << utf8_converter.from_bytes("Hull-White 1F");
        break;
    case utils::ModelType::GeneralizedHullWhite1F:
        os << utf8_converter.from_bytes("Generalized Hull-White 1F");
        break;
    case utils::ModelType::GeneralizedHullWhiteLogNormal1F:
        os << utf8_converter.from_bytes("Generalized Hull-White Log Normal 1F");
        break;
    }
    return os;
}