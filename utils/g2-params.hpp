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
#include <utils/json_io.h>
#include <utils/short_rate_calibration.hpp>
#include <iostream>

namespace utils {
    class G2ModelParams : public ShortRateModelParams {
    protected:
        QuantLib::Real a_;
        QuantLib::Real sigma_;
        QuantLib::Real b_;
        QuantLib::Real eta_;
        QuantLib::Real rho_;
    public:
        G2ModelParams(
            QuantLib::Real a = 0.1,
            QuantLib::Real sigma = 0.01,
            QuantLib::Real b = 0.1,
            QuantLib::Real eta = 0.01,
            QuantLib::Real rho = -0.75
        ) :
            a_(a),
            sigma_(sigma),
            b_(b),
            eta_(eta),
            rho_(rho)
        {}
        const QuantLib::Real& a() const { return a_; }
        QuantLib::Real& a() { return a_; }
        const QuantLib::Real& sigma() const { return sigma_; }
        QuantLib::Real& sigma() { return sigma_; }
        const QuantLib::Real& b() const { return b_; }
        QuantLib::Real& b() { return b_; }
        const QuantLib::Real& eta() const { return eta_; }
        QuantLib::Real& eta() { return eta_; }
        const QuantLib::Real& rho() const { return rho_; }
        QuantLib::Real& rho() { return rho_; }

        QuantLib::ext::shared_ptr<QuantLib::ShortRateModel> createModel(
            const QuantLib::Handle<QuantLib::YieldTermStructure>& yieldtermStructure
        ) const {
            return QuantLib::ext::shared_ptr<QuantLib::G2>(
                new QuantLib::G2(
                    yieldtermStructure,
                    a(),
                    sigma(),
                    b(),
                    eta(),
                    rho()
                )
            );
        }
        void dump(
            std::ostream& os
        ) const {
            os << std::fixed;
            os << "a=" << a() << std::endl;
            os << "sigma=" << sigma() << std::endl;
            os << "b=" << b() << std::endl;
            os << "eta=" << eta() << std::endl;
            os << "rho=" << rho() << std::endl;
        }
    public:
        prop_tree_type serialize() const {
            prop_tree_type pt;
            pt.put("a", a());
            pt.put("sigma", sigma());
            pt.put("b", b());
            pt.put("eta", eta());
            pt.put("rho", rho());
            return pt;
        }
        void deserialize(const prop_tree_type& pt) {
            a() = pt.get<QuantLib::Real>("a");
            sigma() = pt.get<QuantLib::Real>("sigma");
            b() = pt.get<QuantLib::Real>("b");
            eta() = pt.get<QuantLib::Real>("eta");
            rho() = pt.get<QuantLib::Real>("rho");
        }
    };
}