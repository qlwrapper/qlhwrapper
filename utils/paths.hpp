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
#include <cmath>

namespace utils {
	class Paths {
	private:
		std::shared_ptr<QuantLib::Matrix> matrix_;
		std::shared_ptr<QuantLib::TimeGrid> timeGrid_;
	public:
		Paths(
			const std::shared_ptr<QuantLib::Matrix>& matrix = std::shared_ptr<QuantLib::Matrix>(),
			const std::shared_ptr<QuantLib::TimeGrid>& timeGrid = std::shared_ptr<QuantLib::TimeGrid>()
		) : matrix_(matrix), timeGrid_(timeGrid) {}
		Paths(
			QuantLib::Size rows,
			QuantLib::Size columns,
			const std::shared_ptr<QuantLib::TimeGrid>& timeGrid = std::shared_ptr<QuantLib::TimeGrid>()
		) : matrix_(new QuantLib::Matrix(rows, columns)), timeGrid_(timeGrid) {}
		const std::shared_ptr<QuantLib::Matrix>& matrix() const {
			return matrix_;
		}
		std::shared_ptr<QuantLib::Matrix>& matrix() {
			return matrix_;
		}
		const std::shared_ptr<QuantLib::TimeGrid>& timeGrid() const {
			return timeGrid_;
		}
		std::shared_ptr<QuantLib::TimeGrid>& timeGrid() {
			return timeGrid_;
		}
		QuantLib::Size nTimes() const {
			return matrix()->rows();
		}
		QuantLib::Size nPaths() const {
			return matrix()->columns();
		}
		template <typename ITER>
		static std::pair<double, double> calculateMeanSD(
			ITER iter_begin,
			ITER iter_end
		) {
			double sum = 0.0, mean, standardDeviation = 0.0;
			size_t i = 0;
			for (decltype(iter_begin) p = iter_begin; p != iter_end; ++p, ++i) {
				sum += *p;
			}
			mean = sum / (double)i;
			for (decltype(iter_begin) p = iter_begin; p != iter_end; ++p) {
				standardDeviation += std::pow(*p - mean, 2);
			}
			standardDeviation = std::sqrt(standardDeviation / (double)i);
			return std::pair<double, double>{mean, standardDeviation};
		}
		static std::shared_ptr<std::vector<std::pair<double, double>>> calculatePathsStatistics(
			const QuantLib::Matrix& paths
		) {
			auto rows = paths.rows();
			std::shared_ptr<std::vector<std::pair<double, double>>> ret(new std::vector<std::pair<double, double>>(rows));
			auto& statisitcs = *ret;
			for (decltype(rows) i = 0; i < rows; ++i) { // for each grid time/matrix row
				auto pr = calculateMeanSD(paths.row_begin(i), paths.row_end(i));
				statisitcs[i] = pr;
			}
			return ret;
		}
		std::shared_ptr<std::vector<std::pair<double, double>>> statistics() const {
			return calculatePathsStatistics(*matrix());
		}
		// merge the paths with the time grid at the first column
		static std::shared_ptr<QuantLib::Matrix> mergeWithTimeGrid(
			const QuantLib::Matrix& paths,
			const QuantLib::TimeGrid& timeGrid
		) {
			auto nTimes = timeGrid.size();
			QL_REQUIRE(paths.rows() == timeGrid.size(), "paths.rows() != timeGrid.size()");
			auto nPaths = paths.columns();
			std::shared_ptr<QuantLib::Matrix> ret(new QuantLib::Matrix(nTimes, nPaths + 1));
			auto& m = *ret;
			for (decltype(nTimes) i = 0; i < nTimes; ++i) { // for each row
				m[i][0] = timeGrid.at(i);
				std::copy(paths.row_begin(i), paths.row_end(i), m.row_begin(i) + 1);
			}
			return ret;
		}
		std::shared_ptr<QuantLib::Matrix> mergeWithTime() const {
			return mergeWithTimeGrid(*matrix(), *timeGrid());
		}
		// shock the paths with a Shocker
		// the Shocker must implement a function call operator
		// Quant::Real operator() (const QuantLib::Time& t, const Quant::Real& value) const;
		template <
			typename Shocker
		>
		static std::shared_ptr<QuantLib::Matrix> shockPaths(
			const QuantLib::Matrix& paths,
			const QuantLib::TimeGrid& timeGrid,
			const Shocker& shocker
		) {
			auto rows = paths.rows();
			auto columns = paths.columns();
			QL_ASSERT(paths.rows() == timeGrid.size(), "paths.rows() != timeGrid.size()");
			std::shared_ptr<QuantLib::Matrix> ret(new QuantLib::Matrix(paths));
			auto& m = *ret;
			for (decltype(rows) row = 0; row < rows; ++row) {   // for each row
				auto const& t = timeGrid.at(row);
				for (decltype(columns) col = 0; col < columns; ++col) { // for ech column
					m[row][col] += shocker(t, paths[row][col]);
				}
			}
			return ret;
		}
		template <
			typename Shocker
		>
		std::shared_ptr<Paths> shock(
			const Shocker& shocker
		) const {
			auto shockedMatrix = shockPaths(*matrix(), *timeGrid(), shocker);
			return std::shared_ptr<Paths>(new Paths(shockedMatrix, timeGrid()));
		}
	};
}