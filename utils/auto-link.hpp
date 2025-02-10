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

#ifdef _M_X64
#  define SHARED_LIB_PLATFORM "-x64"
#else
#  define SHARED_LIB_PLATFORM "-x86"
#endif

#ifdef _DEBUG
#  define SHARED_LIB_CONFIGURATION "-Debug"
#else
#  define SHARED_LIB_CONFIGURATION "-Release"
#endif

#define SHARED_LIB_NAME "shared" SHARED_LIB_PLATFORM SHARED_LIB_CONFIGURATION ".lib"

#pragma comment(lib, SHARED_LIB_NAME)
#ifdef BOOST_LIB_DIAGNOSTIC
#  pragma message("Will (need to) link to lib file: " SHARED_LIB_NAME)
#endif

#undef SHARED_LIB_PLATFORM
#undef SHARED_LIB_CONFIGURATION
#undef SHARED_LIB_NAME