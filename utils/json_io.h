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

#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>

namespace utils {
	namespace json_io {
		template <typename T>
		struct ValueTypeSerializerTraits {
			static boost::property_tree::ptree serialize(const T& value) {
				boost::property_tree::ptree pt;
				pt.put("", value);
				return pt;
			}
			static T deserialize(const boost::property_tree::ptree& pt) {
				return pt.get_value<T>();
			}
		};

		template <typename T, typename TSerializerTraits = ValueTypeSerializerTraits<T>>
		struct vector_serializer {
			static boost::property_tree::ptree serialize(const std::vector<T>& v) {
				boost::property_tree::ptree pt;
				for (auto const& item : v) {
					auto child = TSerializerTraits::serialize(item);
					pt.push_back(std::make_pair("", child));
				}
				return pt;
			}
			static std::vector<T> deserialize(const boost::property_tree::ptree& pt) {
				std::vector<T> v;
				for (auto const& pair : pt) {
					auto const& child = pair.second;
					auto item = TSerializerTraits::deserialize(child);
					v.push_back(item);
				}
				return v;
			}
		};

		template <class _Elem>
		std::basic_string<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>> stringify_ptree(const boost::property_tree::ptree& pt, bool pretty = false) {
			std::basic_ostringstream<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>> buf;
			boost::property_tree::write_json(buf, pt, pretty);
			return buf.str();
		}

		template <class _Elem>
		boost::property_tree::ptree parse_json(const std::basic_string<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>>& json) {
			std::basic_istringstream<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>> is(json);
			boost::property_tree::ptree pt;
			boost::property_tree::read_json(is, pt);
			return pt;
		}

		struct json_serializable {
			typedef boost::property_tree::ptree prop_tree_type;
			template<typename _Elem> using string_type = std::basic_string<_Elem, std::char_traits<_Elem>, std::allocator<_Elem>>;
			virtual prop_tree_type serialize() const = 0;
			virtual void deserialize(const prop_tree_type& pt) = 0;

			template <class _Elem>
			inline string_type<_Elem> JSON_stringify(bool pretty = false) const {
				auto pt = serialize();
				return stringify_ptree<_Elem>(pt, pretty);
			}
			template <class _Elem>
			inline void JSON_parse(const string_type<_Elem>& json) {
				auto pt = parse_json<_Elem>(json);
				deserialize(pt);
			}
		};
	};
};