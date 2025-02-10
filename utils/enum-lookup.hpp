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

#include <utils/string.hpp>
#include <string>
#include <map>
#include <sstream>
#include <iostream>
#include <ql/quantlib.hpp>
#include <boost/lexical_cast.hpp>

namespace utils {
	// helper class to lookup from string to an enum type
	// lookup is case-insensitive
	template<
		typename ENUM_TYPE
	>
	class EnumLookupBase {
	public:
		typedef ENUM_TYPE EnumType;
	private:
		std::map<std::string, EnumType> lookup_;
		std::map<EnumType, std::string> lookupReverse_;
	protected:
		EnumLookupBase() {}
		// add enum alias
		EnumLookupBase& add(
			const EnumType& e,
			const std::string& alias,
			bool isMainAlias = false
		) {
			QL_REQUIRE(!alias.empty(), "enum alias cannot be empty");
			if (isMainAlias) {
				lookupReverse_[e] = alias;
			}
			auto s = lcase_copy(alias);
			lookup_[s] = e;
			// lookup via numeric value string
			std::ostringstream oss;
			oss << (int)e;
			s = oss.str();
			lookup_[s] = e;
			return *this;
		}
	public:
		virtual std::string getEnumName() const = 0;
		// string to enum lookup operator
		EnumType operator () (
			const std::string& str
		) const {
			auto s = lcase_copy(str);
			try {
				return lookup_.at(s);
			}
			catch (...) {
				QL_FAIL("bad/unsupported enum string '" << str << "' for enum type '" << getEnumName() << "'");
			}
		}
		std::string operator () (
			const EnumType& en
		) const {
			try {
				return lookupReverse_.at(en);
			}
			catch (...) {
				QL_FAIL("bad/unsupported enum (" << (int)en << ") for enum type '" << getEnumName() << "'");
			}
		}
	};
}

#define BEGIN_DECLARE_ENUM_LOOKUP_NS(ENUM_NAMESPACE, ENUM) \
class ENUM##_Lookup: public utils::EnumLookupBase<ENUM_NAMESPACE##::ENUM> {	\
public: std::string getEnumName() const {return #ENUM;}	\
private: static const ENUM##_Lookup& getSingleton() {static ENUM##_Lookup lookup; return lookup;}	\
public: static ENUM_NAMESPACE##::ENUM fromString(const std::string& s) {return getSingleton()(s);} \
public: static std::string toString(const ENUM_NAMESPACE##::ENUM& en) {return getSingleton()(en);}	\
private: ENUM##_Lookup() {

#define BEGIN_DECLARE_ENUM_LOOKUP(ENUM) \
class ENUM##_Lookup: public utils::EnumLookupBase<ENUM> {	\
public: std::string getEnumName() const {return #ENUM;}	\
private: static const ENUM##_Lookup& getSingleton() {static ENUM##_Lookup lookup; return lookup;}	\
public: static ENUM fromString(const std::string& s) {return getSingleton()(s);} \
public: static std::string toString(const ENUM& en) {return getSingleton()(en);}	\
private: ENUM##_Lookup() {

#define ENUM_LOOKUP_ITEM(item, alias, isMainAlias)  add(item, alias, isMainAlias)
#define END_DECLARE_ENUM_LOOKUP() }};

#define ENUM_BOOST_LEXICAL_CAST_SUPPORT(ENUM_NAMESPACE, ENUM, LOOKUP_NAMESPACE) \
template<>	\
inline ENUM_NAMESPACE##::ENUM lexical_cast<ENUM_NAMESPACE##::ENUM, std::string>(const std::string& s) {	\
	using lookup_type = LOOKUP_NAMESPACE##::ENUM##_Lookup;	\
	return lookup_type::fromString(s);	\
}	\
template<>	\
inline std::string lexical_cast<std::string, ENUM_NAMESPACE##::ENUM>(const ENUM_NAMESPACE##::ENUM& en) {	\
	using lookup_type = LOOKUP_NAMESPACE##::ENUM##_Lookup;	\
	return lookup_type::toString(en);	\
}	\
inline std::ostream& operator<<(std::ostream& os, const ENUM_NAMESPACE##::ENUM& rhs) {	\
	using lookup_type = LOOKUP_NAMESPACE##::ENUM##_Lookup;	\
	return os << lookup_type::toString(rhs);	\
}