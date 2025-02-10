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
#include <memory>
#include <chrono>
#include <map>
#include <string>

namespace utils {
    class BlackModelCalibrator {
    public:
        // calibration input
        const QuantLib::EndCriteria& endCriteria;
        std::shared_ptr<std::vector<QuantLib::ext::shared_ptr<QuantLib::BlackCalibrationHelper>>> pHelpers;
        // calibration output:
        QuantLib::Integer functionEvaluation;
        QuantLib::EndCriteria::Type endCriteriaType;
        long long durationSeconds;
    private:
        void resetOutput() {
            functionEvaluation = 0;
            endCriteriaType = QuantLib::EndCriteria::None;
            durationSeconds = 0;
        }
    public:
        BlackModelCalibrator(
            const QuantLib::EndCriteria& ec,
            const std::shared_ptr<std::vector<QuantLib::ext::shared_ptr<QuantLib::BlackCalibrationHelper>>>& helpers
        ) : endCriteria(ec), pHelpers(helpers) {}
        // calibrate the model
        void calibrate(
            const QuantLib::ext::shared_ptr<QuantLib::CalibratedModel>& model,
            const QuantLib::ext::shared_ptr<QuantLib::PricingEngine>& engine,
            const std::vector<bool>& fixParameters = std::vector<bool>()
        ) {
            QL_ASSERT(pHelpers != nullptr && pHelpers->size() > 0, "helpers cannot be null or empty");
            resetOutput();
            auto const& helpers = *pHelpers;
            // assign pricing engine to all calibration helpers
            for (auto it = helpers.begin(); it != helpers.end(); it++) {
                auto const& helper = *it;
                helper->setPricingEngine(engine);
            }
            auto n = helpers.size();
            // convert helpers (BlackCalibrationHelper) to vector of CalibrationHelper toa void deprecation warning
            std::vector<QuantLib::ext::shared_ptr<QuantLib::CalibrationHelper>> tmp(helpers.size());
            for (decltype(n) i = 0; i <n; ++i)
                tmp[i] = QuantLib::ext::static_pointer_cast<QuantLib::CalibrationHelper>(helpers[i]);

            QuantLib::LevenbergMarquardt method;

            auto start = std::chrono::high_resolution_clock::now();
            const std::vector<QuantLib::Time> weights;
            model->calibrate(tmp, method, endCriteria, QuantLib::NoConstraint(), weights, fixParameters);
            auto stop = std::chrono::high_resolution_clock::now();
            functionEvaluation = model->functionEvaluation();
            endCriteriaType = model->endCriteria();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            durationSeconds = duration.count();
        }
    };
}