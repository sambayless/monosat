/*****************************************************************************************[Main.cc]
 The MIT License (MIT)

 Copyright (c) 2017, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/
#include "monosat/Version.h"

#define QUOTE(name) #name
#define STR_VALUE(macro) QUOTE(macro)

//define MONOSAT_VERSION and MONOSAT_BUILD when building to include more precise version information
#define VERSION_NUM "1.6.0"

#ifdef MONOSAT_BUILD
#define BUILD_TYPE ", " STR_VALUE(MONOSAT_BUILD)
#else
#define BUILD_TYPE ""
#endif


#ifdef MONOSAT_VERSION
#define VERSION_STR ", (" STR_VALUE(MONOSAT_VERSION) ")"
#else
#define VERSION_STR ""
#endif

const char* Monosat::MONOSAT_VERSION_STR = VERSION_NUM BUILD_TYPE VERSION_STR;