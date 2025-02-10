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

#ifdef _MSC_VER
	#include <utils/auto-link.hpp>
#endif
#include <utils/string-io-types.hpp>
#include <utils/string.hpp>
#include <utils/json_io.h>
#include <utils/text-file-reader.hpp>
#include <utils/date-format.hpp>
#include <utils/enum-lookup.hpp>
#include <utils/misc.h>
#include <utils/enum-lookups.hpp>
#include <utils/io.hpp>
#include <utils/quant_type.hpp>
#include <utils/interest_rate.hpp>
#include <utils/black_model_calibrator.hpp>
#include <utils/calibration_helpers.hpp>
#include <utils/short_rate_calibration.hpp>
#include <utils/hull_white_params.hpp>
#include <utils/g2-params.hpp>
#include <utils/model_type.hpp>
#include <utils/short_rate_tree_adaptor.hpp>
#include <utils/paths.hpp>
#include <utils/irp_simulation.hpp>
#include <utils/backward-induction.hpp>
#include <utils/discretizedbonds.hpp>
#include <utils/market_rate.hpp>
#include <utils/flat_rate_shocker.hpp>
#include <utils/one_factor_affine_calc.hpp>
#include <utils/paths_utils.hpp>
#include <utils/zv_calculators.hpp>
#include <utils/simulation_verification.hpp>
#include <utils/path-output.hpp>
#include <utils/curve-config.hpp>
#include <utils/default-defs.hpp>
// hack and debug helpers
///////////////////////////////////////////////////////////////
#include <utils/hullwhite_hack.hpp>
#include <utils/one_factor_affine_like_hack.hpp>
#include <utils/ghw_tie_up.hpp>
///////////////////////////////////////////////////////////////