#pragma once

#include <exception>
#include <memory>
#include <vector>
#include <functional>
#include <utils/types.h>
#include <utils/string.h>

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

	template <typename _Elem>
	class TextFileReader {
	protected:
		template<typename RT>
		static RT readFile(const string_type<_Elem>& filePath, const std::function<RT(ifstream_type<_Elem>&)>& readProcessor) {
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
		static string_type<_Elem> readAsString(const string_type<_Elem>& filePath) {
			return readFile<string_type<_Elem>>(filePath, [](ifstream_type<_Elem>& f) -> string_type<_Elem> {
				stringstream_type<_Elem> buffer;
				buffer << f.rdbuf();
				return buffer.str();
			});
		}

		template <typename LineProcessorTraits = PassthroughLineProcessor<_Elem>>
		static std::shared_ptr<std::vector<typename LineProcessorTraits::ret_type>> readLineByLine(const string_type<_Elem>& filePath, bool ignoreEmptyLine = true) {
			using vector_type = std::vector<typename LineProcessorTraits::ret_type>;
			using ret_type = std::shared_ptr<vector_type>;
			return readFile<ret_type>(filePath, [ignoreEmptyLine](ifstream_type<_Elem>& f) -> ret_type {
				string_type<_Elem> line;
				LineProcessorTraits lineProcessor;
				ret_type ret(new vector_type());
				while (std::getline(f, line)) {
					if (!ignoreEmptyLine || line.length() != 0) {
						auto retLine = lineProcessor(line);
						ret->push_back(retLine);
					}
				}
				return ret;
			});
		}
	};
};