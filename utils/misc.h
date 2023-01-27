#pragma once

#include <ql/quantlib.hpp>
#include <utils/types.h>
#include <boost/lexical_cast.hpp>
#include <utils/console.hpp>
#include <utils/string.h>
#include <fstream>
#include <sstream>
#include <ostream>
#include <vector>
#include <exception>

namespace utils {

    template<typename _Elem>
    struct DateHelper {
        typedef string_type<_Elem> string;
        static QuantLib::Date from_yyyymmdd(const string& yyyymmdd, bool hasHyphen = false) {
            auto yyyy = yyyymmdd.substr(0, 4);
            auto mm = yyyymmdd.substr(4 + (std::size_t)(hasHyphen ? 1 : 0), 2);
            auto dd = yyyymmdd.substr(6 + (std::size_t)(hasHyphen ? 2 : 0));
            auto year = boost::lexical_cast<QuantLib::Year>(yyyy);
            auto month = QuantLib::Month(boost::lexical_cast<int>(mm));
            auto day = boost::lexical_cast<QuantLib::Day>(dd);
            return QuantLib::Date(day, month, year);
        }
        static string to_yyyymmdd(const QuantLib::Date& d, bool hyphen = false) {
            auto d_i = d.year() * 10000 + int(d.month()) * 100 + d.dayOfMonth();
            auto s = boost::lexical_cast<string>(d_i);
            if (hyphen) {
                auto yyyy = s.substr(0, 4);
                auto mm = s.substr(4, 2);
                auto dd = s.substr(6);
                ostringstream_type<_Elem> os;
                os << yyyy << "-" << mm << "-" << dd;
                return os.str();
            }
            else {
                return s;
            }
        }

        static QuantLib::Date getEvaluationDate(const string& dateString) {
            auto T_MINUS_ONE = boost::lexical_cast<string>("T-1");
            auto T_MINUS_ZERO = boost::lexical_cast<string>("T-0");
            auto ds = (dateString.length() == 0 ? T_MINUS_ONE : dateString); // T-1 as default
            if (ds == T_MINUS_ZERO || ds == T_MINUS_ONE) {
                QuantLib::Date dt = QuantLib::Settings::instance().evaluationDate();
                if (ds == T_MINUS_ONE) {
                    dt -= 1 * QuantLib::TimeUnit::Days;
                }
                // avoid weekend by move to Friday prior if necessary
                if (dt.weekday() == QuantLib::Weekday::Saturday)
                    return dt - 1 * QuantLib::TimeUnit::Days;
                else if (dt.weekday() == QuantLib::Weekday::Sunday)
                    return dt - 2 * QuantLib::TimeUnit::Days;
                else
                    return dt;
            }
            else {  // dateString = yyyymmdd
                const char* message = "expecting date format to be yyyymmdd";
                QL_REQUIRE(ds.length() == 8, message);
                auto yyyy = ds.substr(0, 4);
                auto mm = ds.substr(4, 2);
                auto dd = ds.substr(6, 2);
                try {
                    auto year = boost::lexical_cast<QuantLib::Year>(yyyy);
                    auto month = (QuantLib::Month)(boost::lexical_cast<QuantLib::Integer>(mm));
                    auto day = boost::lexical_cast<QuantLib::Day>(dd);
                    return QuantLib::Date(day, month, year);
                }
                catch (const std::exception&) {
                    QL_FAIL(message);
                }
            }
        }
    };

    namespace impl_ {
        template <typename _Elem>
        class null_buffer : public std::basic_streambuf<_Elem, std::char_traits<_Elem>>
        {
        public:
            int overflow(int c) { return c; }
        };

        template <typename _Elem>
        struct get_nullable_cout {
            ostream_type<_Elem>& get_null_ostream() const {
                static null_buffer<_Elem> null_buffer;
                static ostream_type<_Elem> null_stream(&null_buffer);
                return null_stream;
            }
            ostream_type<_Elem>& operator() (bool silent) const {
                return (silent ? get_null_ostream() : console<_Elem>().cout());
            }
        };
    }

    // get cout with silent option
    template <typename _Elem>
    ostream_type<_Elem>& get_nullable_cout(bool silent = false) {
        impl_::get_nullable_cout<_Elem> go;
        return go(silent);
    }

    template<typename _Elem>
    string_type<_Elem> get_input_content(const string_type<_Elem>& input_arg) {
        utf8_wstring_converter<_Elem> utf8_converter;
        auto ret = input_arg;
        auto readFromTextStream = [&utf8_converter](istream_type<_Elem>& is) -> string_type<_Elem> {
            auto newline_s = utf8_converter.from_bytes("\n");
            string_type<_Elem> str;
            ostringstream_type<_Elem> os;
            while (std::getline(is, str)) {
                os << str;
                os << newline_s;
            }
            return os.str();
        };
        auto at_s = utf8_converter.from_bytes("@");
        auto hyphen_s = utf8_converter.from_bytes("-");
        if (input_arg.length() > 0) {
            if (input_arg.substr(0, 1) == at_s) {
                if (input_arg.length() >= 2) {
                    if (input_arg.substr(1, 1) == hyphen_s) {
                        ret = readFromTextStream(console<_Elem>().cin());
                    }
                    else {
                        auto filePath = input_arg.substr(1);
                        ifstream_type<_Elem> file(filePath);
                        if (file.is_open()) {
                            ret = readFromTextStream(file);
                            file.close();
                        }
                        else {
                            auto s = utf8_converter.to_bytes(filePath);
                            std::ostringstream os;
                            os << "unable to open file " << s;
                            throw std::exception(os.str().c_str());
                        }
                    }
                }
            }
        }
        return ret;
    }
}

// dump vector to ostream
template<typename _Elem, typename T>
inline std::basic_ostream<_Elem, std::char_traits<_Elem>>& operator<<(std::basic_ostream<_Elem, std::char_traits<_Elem>>& out, const std::vector<T>& v) {
    utils::utf8_wstring_converter<_Elem> utf8_converter;
    out << utf8_converter.from_bytes("[ ");
    if (!v.empty()) {
        for (size_t n = 0; n < v.size() - 1; ++n)
            out << v.at(n) << utf8_converter.from_bytes("; ");
        out << v.back();
    }
    out << utf8_converter.from_bytes(" ]");
    return out;
}