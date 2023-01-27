#pragma once

#include <ql/quantlib.hpp>
#include <ostream>
#include <memory>
#include <functional>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <utils/utilities.h>

class OutputTraits {
private:
    std::string type_;
    std::string fnSuffix_;
public:
    OutputTraits(
        const std::string& type,
        const std::string& fnSuffix
    ) : type_(type), fnSuffix_(fnSuffix) {}
    const std::string& type() const {
        return type_;
    }
    std::string& type() {
        return type_;
    }
    const std::string& fnSuffix() const {
        return fnSuffix_;
    }
    std::string& fnSuffix() {
        return fnSuffix_;
    }
};

template <typename RateShocker>
class ShockTraits {
private:
    RateShocker rateShocker_;
    std::string scenarioName_;
public:
    ShockTraits(
        QuantLib::Rate rateShock,
        const std::string& scenarioName
    ) : rateShocker_(rateShock), scenarioName_(scenarioName) {}
    const RateShocker& rateShocker() const {
        return rateShocker_;
    }
    const std::string& scenarioName() const {
        return scenarioName_;
    }
    std::string& scenarioName() {
        return scenarioName_;
    }
};

class OneFactorSimulationOutput {
private:
	std::ostream& ostream_;
	std::function<std::string(const std::string&, const std::string&)> get_Paths_Scenario_Model_OutputFilePath_;
    std::function<std::string(const std::string&)> get_Paths_Model_OutputFilePath_;
    std::function<std::string(const std::string&)> get_Model_OutputFilePath_;
	std::function<std::string(const std::string&)> get_OutputFilePath_;
public:
	OneFactorSimulationOutput(
		std::ostream& ostream,
		const std::function<std::string(const std::string&, const std::string&)>& get_Paths_Scenario_Model_OutputFilePath,
        const std::function<std::string(const std::string&)>& get_Paths_Model_OutputFilePath,
        const std::function<std::string(const std::string&)>& get_Model_OutputFilePath,
		const std::function<std::string(const std::string&)>& get_OutputFilePath
	) :
        ostream_(ostream),
        get_Paths_Scenario_Model_OutputFilePath_(get_Paths_Scenario_Model_OutputFilePath),
        get_Paths_Model_OutputFilePath_(get_Paths_Model_OutputFilePath),
        get_Model_OutputFilePath_(get_Model_OutputFilePath),
        get_OutputFilePath_(get_OutputFilePath)
    {}
	std::ostream& os() {
		return ostream_;
	}
	const std::function<std::string(const std::string&, const std::string&)>& get_Paths_Scenario_Model_OutputFilePath() const {
		return get_Paths_Scenario_Model_OutputFilePath_;
	}
    const std::function<std::string(const std::string&)>& get_Paths_Model_OutputFilePath() const {
        return get_Paths_Model_OutputFilePath_;
    }
    const std::function<std::string(const std::string&)>& get_Model_OutputFilePath() const {
        return get_Model_OutputFilePath_;
    }
	const std::function<std::string(const std::string&)>& get_OutputFilePath() const {
		return get_OutputFilePath_;
	}

    template<typename ITER_TERM, typename ITER_VALUE>
    static void dumpTermStructure(
        const ITER_TERM& itTermStart,
        const ITER_TERM& itTermEnd,
        const ITER_VALUE& itValueStart,
        const std::string& filename,
        bool includeTerm = false,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        std::ofstream file(filename);
        file << std::fixed << std::setprecision(precision);
        auto itValue = itValueStart;
        for (auto itTerm = itTermStart; itTerm != itTermEnd; ++itTerm, ++itValue) {
            auto const& t = *itTerm;
            auto const& value = *itValue;
            if (includeTerm) {
                file << t << "\t";
            }
            file << (value * scaling) << std::endl;
        }
        file.close();
    }
    template <
        typename ZVCalculator
    >
    void dumpZVValues(
        const QuantLib::ext::shared_ptr<QuantLib::YieldTermStructure>& zvCurve,
        const QuantLib::TimeGrid& timeGrid,
        const OutputTraits& outputTraits,
        bool mergeWithGrid = false,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        auto const& getOutputFilePath = get_OutputFilePath();
        auto filename = getOutputFilePath(outputTraits.fnSuffix());
        ZVCalculator zvCalculator(zvCurve);
        auto tMax = timeGrid.back() - zvCalculator.tenor();
        QuantLib::TimeGrid::const_iterator itEnd = timeGrid.end();
        std::vector<QuantLib::Real> values;
        for (auto p = timeGrid.begin(); p != timeGrid.end(); ++p) {
            auto const& t = *p;
            if (t > tMax) {
                itEnd = p;
                break;
            }
            else {
                auto value = zvCalculator(t);
                values.push_back(value);
            }
        }
        dumpTermStructure(timeGrid.begin(), itEnd, values.begin(), filename, mergeWithGrid, scaling, precision);
        os() << std::endl;
        os() << outputTraits.type() << " file -> " << filename << std::endl;
    }

    void dumpTreeCMSBreakevenRates(
        const std::vector<QuantLib::Rate>& rates,
        const QuantLib::TimeGrid& timeGrid,
        const OutputTraits& outputTraits,
        bool mergeWithGrid = false,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        auto const& getOutputFilePath = get_Model_OutputFilePath();
        auto filename = getOutputFilePath(outputTraits.fnSuffix());
        dumpTermStructure(timeGrid.begin(), timeGrid.end(), rates.begin(), filename, mergeWithGrid, scaling, precision);
        os() << std::endl;
        os() << outputTraits.type() << " file -> " << filename << std::endl;
    }

    void dumpSimulationCMSBreakevenRates(
        const std::vector<QuantLib::Rate>& rates,
        const QuantLib::TimeGrid& timeGrid,
        const OutputTraits& outputTraits,
        bool mergeWithGrid = false,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        auto const& getOutputFilePath = get_Paths_Model_OutputFilePath();
        auto filename = getOutputFilePath(outputTraits.fnSuffix());
        dumpTermStructure(timeGrid.begin(), timeGrid.end(), rates.begin(), filename, mergeWithGrid, scaling, precision);
        os() << std::endl;
        os() << outputTraits.type() << " file -> " << filename << std::endl;
    }

    static void printMatrix(
        const QuantLib::Matrix& matrix,
        std::string filePathName,
        std::string sep,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        auto rows = matrix.rows();
        auto columns = matrix.columns();
        // open text file for input, loop through matrix rows
        std::ofstream file(filePathName);
        for (decltype(rows) i = 0; i < rows; ++i) { // for each row
            std::ostringstream oss;
            oss << std::fixed;
            oss << std::setprecision(precision);
            for (decltype(columns) j = 0; j < columns; j++) {   // for each column
                oss << (matrix[i][j] * scaling);
                if (j != columns - 1) {
                    oss << sep;
                }
                else {
                    oss << std::endl;
                }
            }
            file << oss.str();
        }
        // close text file
        file.close();
    }

    template <
        typename UpShockTraits,
        typename DownShockTraits
    >
    void dumpPathsWithShocks(
        const utils::Paths& basePaths,
        const UpShockTraits& upShockTraits,
        const DownShockTraits& downShockTraits,
        const OutputTraits& outputTraits,
        bool mergeWithGrid = false,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        auto p_up = basePaths.shock(upShockTraits.rateShocker());
        auto p_dn = basePaths.shock(downShockTraits.rateShocker());
        auto m_base = basePaths.matrix();
        auto m_up = p_up->matrix();
        auto m_dn = p_dn->matrix();
        if (mergeWithGrid) {
            m_base = basePaths.mergeWithTime();
            m_up = p_up->mergeWithTime();
            m_dn = p_dn->mergeWithTime();
        }
        std::vector<std::pair<std::string, std::shared_ptr<QuantLib::Matrix>>> scenarios{
            {"base", m_base},
            {upShockTraits.scenarioName(), m_up},
            {downShockTraits.scenarioName(), m_dn}
        };
        os() << std::endl;
        auto const& getOutputFilePath = get_Paths_Scenario_Model_OutputFilePath();
        for (auto const& scenario : scenarios) {   // for each scenario
            auto const& scenarioName = scenario.first;
            auto const& matrix = *(scenario.second);
            auto filename = getOutputFilePath(outputTraits.fnSuffix(), scenarioName);
            printMatrix(matrix, filename, "\t", scaling, precision);
            os() << scenarioName << " " << outputTraits.type() << " file -> " << filename << std::endl;
        }
    }
};
