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
#include <vector>
#include <algorithm>
#include <iterator>

namespace utils {
	class DiscretizedFixedRateBond : public QuantLib::DiscretizedAsset {
	public:
		DiscretizedFixedRateBond(
			const std::vector<QuantLib::Time>& paymentTimes,
			const std::vector<QuantLib::Rate>& coupons,
			QuantLib::Real redemption = 1.0
		) : paymentTimes_(paymentTimes), coupons_(coupons), redemption_(redemption) {}
		std::vector<QuantLib::Time> madatoryTimes() const {
			return paymentTimes_;
		}
		void reset(QuantLib::Size size) {
			values_ = QuantLib::Array(size, redemption_);
			adjustValues();
		}
	private:
		void postAdjustValuesImpl() {
			auto nPayments = paymentTimes_.size();
			for (decltype(nPayments) i = 0; i < nPayments; ++i) {	// for each payment
				auto t = paymentTimes_[i];
				if (t >= 0.0 && isOnTime(t)) {
					addCoupon(i);
				}
			}
		}
		void addCoupon(QuantLib::Size i) {
			values_ += coupons_[i];
		}

		std::vector<QuantLib::Time> paymentTimes_;
		std::vector<QuantLib::Rate> coupons_;
		QuantLib::Real redemption_;
	};

	class DiscretizedFloatingRateBond : public QuantLib::DiscretizedAsset {
	public:
		DiscretizedFloatingRateBond(
			const std::vector<QuantLib::Time>& paymentTimes,
			const std::vector<QuantLib::Time>& fixingTimes,
			const std::vector<QuantLib::Rate>& coupons = std::vector<QuantLib::Rate>(),
			QuantLib::Real notional = 1.0
		) : paymentTimes_(paymentTimes), fixingTimes_(fixingTimes), coupons_(coupons), notional_(notional) {}
		std::vector<QuantLib::Time> madatoryTimes() const {
			std::vector<QuantLib::Time> times = paymentTimes_;
			std::copy(fixingTimes_.begin(), fixingTimes_.end(), std::back_inserter(times));
			return times;
		}
		void reset(QuantLib::Size size) {
			values_ = QuantLib::Array(size, notional_);
			adjustValues();
		}
	private:
		void preAdjustValuesImpl() {
			auto nFixings = fixingTimes_.size();
			for (decltype(nFixings) i = 0; i < nFixings; ++i) {	// for each fixing time
				auto t = fixingTimes_[i];
				if (t >= 0.0 && isOnTime(t)) {
					addCoupon(i);
				}
			}
		}
		void addCoupon(QuantLib::Size i) {
			QuantLib::DiscretizedDiscountBond bond;
			auto T = paymentTimes_[i];
			auto t = time_;	// == fixingTimes_[i];
			auto dt = T - t;
			bond.initialize(method(), T);
			bond.rollback(t);
			auto dfs = bond.values();
			auto nNodes = values_.size();
			QL_ASSERT(dfs.size() == nNodes, "dfs.size() != values_.size()");
			for (decltype(nNodes) index = 0; index < nNodes; ++index) {	// for each node @ time_/fixingTimes_[i]
				auto const& df = dfs.at(index);	// discount factor at this node
				auto coupon = (coupons_.size() == 0 ? (1.0/df - 1.0) / dt: coupons_[i]);
				auto cashflow = notional_ * coupon * dt;
				values_[index] += cashflow * df;
			}
		}

		std::vector<QuantLib::Time> paymentTimes_;
		std::vector<QuantLib::Time> fixingTimes_;
		std::vector<QuantLib::Rate> coupons_;
		QuantLib::Real notional_;
	};
}