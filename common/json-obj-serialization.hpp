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

#include <iostream>
#include <vector>
#include <memory>
#include <boost/property_tree/ptree.hpp>

// user must define std::ostream& operator << (std::ostream&, const T&)
template <typename T>
inline std::ostream& operator << (
    std::ostream& os,
    const std::shared_ptr<T>& rhs
) {
    if (rhs == nullptr) {
        return (os << "null");
    }
    else {
        return (os << *rhs);
    }
}

template <typename T>
inline std::ostream& operator << (
    std::ostream& os,
    const std::vector<T>& rhs
) {
    os << "[" << std::endl;
    auto n = rhs.size();
    for (decltype(n) i = 0; i < n; ++i) {
        if (i > 0) {
            os << ",";
        }
        os << rhs[i] << std::endl;
    }
    os << "]";
    return os;
}

template <typename T>
inline std::ostream& operator << (
    std::ostream& os,
    const std::shared_ptr<std::vector<T>>& rhs
) {
    if (rhs == nullptr) {
        return (os << "null");
    }
    else {
        return (os << *rhs);
    }
}

// user must define std::ostream& operator << (std::ostream&, std::pair<const T&, const CONTEXT&>&)
template <typename T, typename CONTEXT>
inline std::ostream& operator << (
    std::ostream& os,
    const std::pair<const std::shared_ptr<T>&, const CONTEXT&>& rhs
) {
    const auto& p = rhs.first;
    const auto& context = rhs.second;
    if (p == nullptr) {
        return (os << "null");
    }
    else {
        std::pair<const T&, const CONTEXT&> pr(*p, context);
        return (os << pr);
    }
}

template <typename T, typename CONTEXT>
inline std::ostream& operator << (
    std::ostream& os,
    const std::pair<const std::vector<T>&, const CONTEXT&>& rhs
) {
    os << "[" << std::endl;
    const auto& v = rhs.first;
    const auto& context = rhs.second;
    auto n = v.size();
    for (decltype(n) i = 0; i < n; ++i) {
        if (i > 0) {
            os << ",";
        }
        os << std::pair<const T&, const CONTEXT&>(v[i], context) << std::endl;
    }
    os << "]";
    return os;
}

// user must define T& operator << (T&, const boost::property_tree::ptree& ptree)
template <typename T>
inline std::shared_ptr<T>& operator << (
    std::shared_ptr<T>& dest,
    const boost::property_tree::ptree& ptree
) {
    if (!ptree.empty() && ptree.data().empty()) {
        dest.reset(new T());
        (*dest) << ptree;
    }
    else {
        dest = nullptr;
    }
    return dest;
}

template <typename T>
inline std::shared_ptr<std::vector<std::shared_ptr<T>>>& operator << (
    std::shared_ptr<std::vector<std::shared_ptr<T>>>& dest,
    const boost::property_tree::ptree& ptree
) {
    for (const auto& node : ptree) {
        const auto& o = node.second;
        if (dest == nullptr) {
            dest.reset(new std::vector<std::shared_ptr<T>>());
        }
        std::shared_ptr<T> pT;
        pT << o;
        dest->push_back(pT);
    }
    return dest;
}

// user must define T& operator << (T&, const std::pair<const boost::property_tree::ptree&, const CONTEXT&>&)
template <typename T, typename CONTEXT>
inline std::shared_ptr<T>& operator << (
    std::shared_ptr<T>& dest,
    const std::pair<const boost::property_tree::ptree&, const CONTEXT&>& rhs
) {
    const auto& ptree = rhs.first;
    if (!ptree.empty() && ptree.data().empty()) {
        dest.reset(new T());
        (*dest) << rhs;
    }
    else {
        dest = nullptr;
    }
    return dest;
}

template <typename T, typename CONTEXT>
inline std::shared_ptr<std::vector<std::shared_ptr<T>>>& operator << (
    std::shared_ptr<std::vector<std::shared_ptr<T>>>& dest,
    const std::pair<const boost::property_tree::ptree&, const CONTEXT&>& rhs
) {
    const auto& ptree = rhs.first;
    const auto& context = rhs.second;
    for (const auto& node : ptree) {
        const auto& o = node.second;
        if (dest == nullptr) {
            dest.reset(new std::vector<std::shared_ptr<T>>());
        }
        std::shared_ptr<T> pT;
        pT << std::pair<const boost::property_tree::ptree&, const CONTEXT&>(o, context);
        dest->push_back(pT);
    }
    return dest;
}
