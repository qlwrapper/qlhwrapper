#pragma once

#include <map>
#include <string>
#include <string>
#include <exception>
#include <sstream>

namespace utils {
    // helper class to translate enum type to string and vice versa
    template <typename EnumType, typename MapLoader>
    struct EnumHelper {
    protected:
        static const std::map<EnumType, std::string>& get_lookup_map() {
            static std::map<EnumType, std::string> map;
            if (map.size() == 0) {
                MapLoader loader;
                loader(map);
            }
            return map;
        }
    public:
        static std::string to_string(const EnumType& value) {
            auto const& map_ = get_lookup_map();
            auto iter = map_.find(value);
            if (iter != map_.end()) {
                return iter->second;
            }
            throw std::exception("unknown/unsupported enum");
        }
        static EnumType from_string(const std::string& value) {
            auto const& map_ = get_lookup_map();
            for (auto const& i : map_) {
                if (i.second == value) {
                    return i.first;
                }
            }
            std::ostringstream os;
            os << "unknown/unsupported enum string '" << value << "'";
            throw std::exception(os.str().c_str());
        }
    };
};