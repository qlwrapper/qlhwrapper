#pragma once

#include <exception>
#include <memory>
#include <vector>
#include <functional>
#include <utils/string-io-types.hpp>
#include <utils/string.hpp>

namespace utils {

	template <typename _Elem>
	struct PassthroughLineProcessor {
		typedef string_type<_Elem> ret_type;
		ret_type operator()(const string_type<_Elem>& line) const {
			return line;
		}
	};

	template <typename _Elem, unsigned int delimAsci>
	struct ParseDelimitedLineProcessor {
		typedef std::vector<string_type<_Elem>> vector_type;
		typedef std::shared_ptr<vector_type> ret_type;
		ret_type operator()(const string_type<_Elem>& line) const {
			return parse_delimited<_Elem>(line, delimAsci);
		}
	};

	template <
		typename _Elem
	>
	class TextFileReader {
	protected:
		using istream_type = istream_type<_Elem>;
		template<typename RT>
		static RT readFile(
			const string_type<_Elem>& filePath,
			const std::function<RT(istream_type&)>& readProcessor
		) {
			ifstream_type<_Elem> f(filePath.c_str());
			if (!f.good()) {
				std::stringstream ss;
				ss << "error reading file " << filePath;
				throw std::exception(ss.str().c_str());
			}
			RT ret = readProcessor(f);
			f.close();
			return ret;
		}
	public:
		static auto getAllAsStringProcessor() {
			return (
				[](istream_type& is) -> string_type<_Elem> {
					stringstream_type<_Elem> buffer;
					buffer << is.rdbuf();
					return buffer.str();
				}
			);
		}
		static string_type<_Elem> readAsString(
			const string_type<_Elem>& filePath
		) {
			return readFile<string_type<_Elem>>(filePath, getAllAsStringProcessor());
		}
		static string_type<_Elem> readAsString(
			istream_type& is
		) {
			return getAllAsStringProcessor()(is);
		}

		template <
			typename LineProcessorTraits = PassthroughLineProcessor<_Elem>
		>
		static auto getLineByLineProcessor(
			bool ignoreEmptyLine = true
		) {
			using vector_type = std::vector<typename LineProcessorTraits::ret_type>;
			using ret_type = std::shared_ptr<vector_type>;
			return [ignoreEmptyLine](istream_type& is) -> ret_type {
				string_type<_Elem> line;
				LineProcessorTraits lineProcessor;
				ret_type ret(new vector_type());
				while (std::getline(is, line)) {
					if (!ignoreEmptyLine || line.length() != 0) {
						auto retLine = lineProcessor(line);
						ret->push_back(retLine);
					}
				}
				return ret;
			};
		}
		template <
			typename LineProcessorTraits = PassthroughLineProcessor<_Elem>
		>
			static auto readLineByLine(
				const string_type<_Elem>& filePath,
				bool ignoreEmptyLine = true
			) {
			using vector_type = std::vector<typename LineProcessorTraits::ret_type>;
			using ret_type = std::shared_ptr<vector_type>;
			return readFile<ret_type>(filePath, getLineByLineProcessor<LineProcessorTraits>(ignoreEmptyLine));
		}
		template <
			typename LineProcessorTraits = PassthroughLineProcessor<_Elem>
		>
			static auto readLineByLine(
				istream_type& is,
				bool ignoreEmptyLine = true
			) {
			auto proc = getLineByLineProcessor<LineProcessorTraits>(ignoreEmptyLine);
			return proc(is);
		}
	};
};