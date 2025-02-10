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
#include <utils/model_type.hpp>
#include <utils/hullwhite_hack.hpp>
#include <utils/one_factor_affine_calc.hpp>

namespace utils {
	class OneFactorAffineLikeHack : public OneFactorAffineLike {
	private:
		ModelType modelType_;
		QuantLib::ext::shared_ptr<QuantLib::OneFactorAffineModel> oneFactorAffineModel_;
	public:
		OneFactorAffineLikeHack(
			ModelType modelType,
			const QuantLib::ext::shared_ptr<QuantLib::OneFactorAffineModel>& oneFactorAffineModel
		) : modelType_(modelType), oneFactorAffineModel_(oneFactorAffineModel){}
		const ModelType& modelType() const {
			return modelType_;
		}
		const QuantLib::ext::shared_ptr<QuantLib::OneFactorAffineModel>& oneFactorAffineModel() const {
			return oneFactorAffineModel_;
		}
	private:
		template<typename Model, typename ModelHack>
		std::pair<QuantLib::Real, QuantLib::Real> calculateABImpl(
			QuantLib::Time t,
			QuantLib::Time T
		) const {
			auto pModel = QuantLib::ext::dynamic_pointer_cast<Model>(oneFactorAffineModel());
			auto const& model = *pModel;
			auto const& modelHack = reinterpret_cast<const ModelHack&>(model);
			auto A = modelHack.A_(t, T);
			auto B = modelHack.B_(t, T);
			return std::pair<QuantLib::Real, QuantLib::Real>(A, B);
		}
	public:
		std::pair<QuantLib::Real, QuantLib::Real> calculateAB(
			QuantLib::Time t,
			QuantLib::Time T
		) const {
			switch (modelType()) {
			case ModelType::HullWhite1F:
				return calculateABImpl<QuantLib::HullWhite, HullWhiteHack>(t, T);
			case ModelType::GeneralizedHullWhite1F:
			case ModelType::GeneralizedHullWhiteLogNormal1F:
				return calculateABImpl<QuantLib::GeneralizedHullWhite, GeneralizedHullWhiteHack>(t, T);
			default:
				QL_FAIL("unknown model type");
			}
		}
	};
}