#pragma once

#include <ql/quantlib.hpp>
#include <utils/json_io.h>
#include <utils/string.h>
#include <utils/short_rate_calibration.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <vector>

namespace utils {
    class HullWhiteModelParams : public ShortRateModelParams {
    protected:
        QuantLib::Real a_;
        QuantLib::Volatility sigma_;
    public:
        HullWhiteModelParams(QuantLib::Real a = 0.0, QuantLib::Volatility sigma = 0.0) : a_(a), sigma_(sigma) {}

        const QuantLib::Real& a() const { return a_; }
        QuantLib::Real& a() { return a_; }
        const QuantLib::Volatility& sigma() const { return sigma_; }
        QuantLib::Volatility& sigma() { return sigma_; }

        QuantLib::ext::shared_ptr<QuantLib::ShortRateModel> createModel(const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure) const {
            return QuantLib::ext::shared_ptr<QuantLib::HullWhite>(new QuantLib::HullWhite(yieldtermStructure, a(), sigma()));
        }
        void dump(std::ostream& os) const {
            os << std::fixed;
            os << "a=" << a() << std::endl;
            os << "sigma=" << sigma() << std::endl;
        }
    public:
        prop_tree_type serialize() const {
            prop_tree_type pt;
            pt.put("a", a());
            pt.put("sigma", sigma());
            return pt;
        }
        void deserialize(const prop_tree_type& pt) {
            a() = pt.get<QuantLib::Real>("a");
            sigma() = pt.get<QuantLib::Volatility>("sigma");
        }
    };

    template <typename T = QuantLib::Real>
    struct LogNormalVolTraits {
        static T log_r(T x) {
            return std::log(x);
        }
        static T exp_r(T x) {
            return std::exp(x);
        }
        static QuantLib::ext::function<T(T)> f() { return &log_r; }
        static QuantLib::ext::function<T(T)> f_inv() { return &exp_r; }
    };

    template <typename T = QuantLib::Real>
    struct NormalVolTraits {
        static T identity_r(T x) {
            return x;
        }
        static QuantLib::ext::function<T(T)> f() { return &identity_r; }
        static QuantLib::ext::function<T(T)> f_inv() { return &identity_r; }
    };

    template <typename VolTraits = NormalVolTraits<QuantLib::Real>>
    class GeneralizedHullWhiteModelParams : public ShortRateModelParams {
    protected:
        QuantLib::Real a_;
        std::vector<QuantLib::Time> sigmaTimes_;
        std::vector<QuantLib::Volatility> sigma_;
    private:
        void verifySigma() const {
            QL_ASSERT(sigmaTimes_.size() == sigma_.size(), "sigma and sigma times must be the same size");
        }
        void loadSigmaFromString(std::string s) {
            auto ret = parse_delimited<char>(s, 44);
            auto n = ret->size();
            std::vector<std::pair<QuantLib::Time, QuantLib::Volatility>> tmp(n);
            for (decltype(n) i = 0; i < n; ++i) {
                auto const& s = ret->at(i);
                auto x = s.find('@');
                if (x != s.npos) {
                    auto s_vol = s.substr(0, x);
                    auto s_time = s.substr(x + 1);
                    tmp[i].first = boost::lexical_cast<QuantLib::Time>(s_time);
                    tmp[i].second = boost::lexical_cast<QuantLib::Volatility>(s_vol);
                }
                else {
                    QL_FAIL("bad sigma format: " << s);
                }
            }
            std::sort(tmp.begin(), tmp.end(), [](const std::pair<QuantLib::Time, QuantLib::Volatility>& a, const std::pair<QuantLib::Time, QuantLib::Volatility>& b) -> bool {
                return (a.first < b.first);
            });
            sigma().resize(n, 0.0);
            sigmaTimes().resize(n, 0.0);
            for (decltype(n) i = 0; i < n; ++i) {
                auto const& t = tmp.at(i).first;
                auto const& vol = tmp.at(i).second;
                sigma()[i] = vol;
                sigmaTimes()[i] = t;
            }
        }
    public:
        const QuantLib::Real& a() const { return a_; }
        QuantLib::Real& a() { return a_; }
        const std::vector<QuantLib::Time>& sigmaTimes() const { return sigmaTimes_; }
        std::vector<QuantLib::Time>& sigmaTimes() { return sigmaTimes_; }
        const std::vector<QuantLib::Volatility>& sigma() const { return sigma_; }
        std::vector<QuantLib::Volatility>& sigma() { return sigma_; }

        GeneralizedHullWhiteModelParams() : a_(0.0) {}
        GeneralizedHullWhiteModelParams(QuantLib::Real a, const std::vector<QuantLib::Time>& sigmaTimes, const std::vector<QuantLib::Volatility>& sigma)
            :a_(a)
            , sigmaTimes_(sigmaTimes)
            , sigma_(sigma)
        {}
        GeneralizedHullWhiteModelParams(QuantLib::Real a, const std::string& sigma) :a_(a) {
            loadSigmaFromString(sigma);
        }

        // from the calibrated model parameters, create the GHW model
        QuantLib::ext::shared_ptr<QuantLib::ShortRateModel> createModel(const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure) const {
            verifySigma();
            auto f = VolTraits::f();
            auto f_inv = VolTraits::f_inv();
            auto baseDate = yieldtermStructure->referenceDate();
            std::vector<QuantLib::Date> speedstructure, volstructure;
            std::vector<QuantLib::Real> a_;

            speedstructure.push_back(baseDate);
            a_.push_back(a());

            auto timeToPeriod = [](const QuantLib::Time& t) -> QuantLib::Period {
                auto days = (QuantLib::Natural)std::round(t * 360.0);
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
            };

            for (auto p = sigmaTimes().begin(); p < sigmaTimes().end(); ++p) {
                auto const& t = *p;
                auto period = timeToPeriod(t);
                auto dt = baseDate + period;
                volstructure.push_back(dt);
            }
            QuantLib::ext::shared_ptr<QuantLib::ShortRateModel> model(new QuantLib::GeneralizedHullWhite(yieldtermStructure, speedstructure, volstructure, a_, sigma(), f, f_inv));
            return model;
        }
        std::string sigmaString(bool fullPrecision = true) const {
            verifySigma();
            auto n = sigma().size();
            std::ostringstream os;
            os << std::fixed;
            if (fullPrecision) {
                os << std::setprecision(16);
            }
            for (decltype(n) i = 0; i < n; ++i) {
                if (i > 0) {
                    os << ",";
                }
                os << sigma()[i] << "@" << sigmaTimes()[i];
            }
            return os.str();
        }
        void dump(std::ostream& os) const {
            os << std::fixed;
            os << "a=" << a() << std::endl;
            os << "sigma=" << sigmaString(false) << std::endl;
        }
        prop_tree_type serialize() const {
            prop_tree_type pt;
            pt.put("a", a());
            pt.put("sigma", sigmaString(true));
            return pt;
        }
        void deserialize(const prop_tree_type& pt) {
            a() = pt.get<QuantLib::Real>("a");
            loadSigmaFromString(pt.get<std::string>("sigma"));
        }
    };
}
