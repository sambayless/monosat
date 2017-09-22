/*****************************************************************************************[Main.cc]
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless
 Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 Copyright (c) 2007-2010, Niklas Sorensson

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
#include <cstddef>
#include <gmpxx.h>
#include <fstream>
#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <signal.h>
#include <zlib.h>
#include <sstream>
#include <string>

#include <unistd.h>
#include <sys/time.h>
#include <algorithm>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "monosat/api/Monosat.h"

using namespace std;
//=================================================================================================



int main(int argc, char** argv) {
	try {
		setUsageHelp(
				"USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");

		SolverPtr solver =  newSolver_args(argc,argv);
        //BVTheoryPtr bvtheory = initBVTheory(solver);


        lbool ret =  toLbool(readGNF(solver,argv[1]));
		fflush(stdout);

		return (ret == l_True ? 10 : ret == l_False ? 20 : 0);

	} catch (OutOfMemoryException&) {
		printf("===============================================================================\n");
		printf("INDETERMINATE\n");
		exit(0);
	}
}
