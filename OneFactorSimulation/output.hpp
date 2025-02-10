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
    std::string extension_;
public:
    OutputTraits(
        const std::string& type,
        const std::string& fnSuffix,
        const std::string& extension = ".txt"
    ) : type_(type), fnSuffix_(fnSuffix), extension_(extension){}
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
    const std::string& extension() const {
        return extension_;
    }
    std::string& extension() {
        return extension_;
    }
};

class OneFactorSimulationOutput {
public:
    typedef std::function<std::string(const std::string&, const std::string&)> FileNameFunction;
private:
	std::ostream& ostream_;
    FileNameFunction get_Model_OutputFilePath_;
    FileNameFunction get_OutputFilePath_;
public:
	OneFactorSimulationOutput(
		std::ostream& ostream,
        const FileNameFunction& get_Model_OutputFilePath,
		const FileNameFunction& get_OutputFilePath
	) :
        ostream_(ostream),
        get_Model_OutputFilePath_(get_Model_OutputFilePath),
        get_OutputFilePath_(get_OutputFilePath)
    {}
	std::ostream& os() {
		return ostream_;
	}
    const FileNameFunction& get_Model_OutputFilePath() const {
        return get_Model_OutputFilePath_;
    }
	const FileNameFunction& get_OutputFilePath() const {
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
        auto filename = getOutputFilePath(outputTraits.fnSuffix(), outputTraits.extension());
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
        auto filename = getOutputFilePath(outputTraits.fnSuffix(), outputTraits.extension());
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
        auto const& getOutputFilePath = get_Model_OutputFilePath();
        auto filename = getOutputFilePath(outputTraits.fnSuffix(), outputTraits.extension());
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
        std::ostringstream oss;
        oss << std::fixed;
        oss << std::setprecision(precision);
        for (decltype(rows) i = 0; i < rows; ++i) { // for each row
            for (decltype(columns) j = 0; j < columns; j++) {   // for each column
                oss << (matrix[i][j] * scaling);
                if (j != columns - 1) {
                    oss << sep;
                }
                else {
                    oss << std::endl;   // last column of the row => newline to a new row
                }
            }
        }
        // open text file for input, loop through matrix rows
        std::ofstream file(filePathName);
        file << oss.str();
        file.close();   // close text file
    }

    void dumpPaths(
        const utils::Paths& basePaths,
        const OutputTraits& outputTraits,
        bool mergeWithGrid = false,
        QuantLib::Real scaling = 1.0,
        std::streamsize precision = 6
    ) {
        auto m_base = basePaths.matrix();
        if (mergeWithGrid) {
            m_base = basePaths.mergeWithTime();
        }
        os() << std::endl;
        const auto& getOutputFilePath = get_Model_OutputFilePath();
        const auto& matrix = *m_base;
        auto filename = getOutputFilePath(outputTraits.fnSuffix(), outputTraits.extension());
        printMatrix(matrix, filename, "\t", scaling, precision);
        os() << outputTraits.type() << " file -> " << filename << std::endl;
    }

    void dumpIRPSimulationOutput(
        const utils::IRPSimulationOutput& simulationOutput,
        const OutputTraits& outputTraits
    ) {
        const auto& filePathFunc = get_Model_OutputFilePath();
        auto filename = filePathFunc(outputTraits.fnSuffix(), outputTraits.extension());
        std::ofstream file(filename);
        file << simulationOutput << std::endl;
        file.close();
        os() << outputTraits.type() << " file -> " << filename << std::endl;
    }
};
