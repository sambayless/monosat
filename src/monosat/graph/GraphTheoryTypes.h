/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

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

#ifndef GRAPHTHEORYTYPES_H_
#define GRAPHTHEORYTYPES_H_

#include <cstddef>
#include "monosat/core/SolverTypes.h"
#include "monosat/mtl/Rnd.h"

namespace Monosat {


struct Edge {
    Var v;
    Var outerVar;
    int from;
    int to;
    int edgeID;
    int bvID;

    Edge(Var v, Var outerVar, int from, int to, int edgeID, int bvID = -1) :
            v(v), outerVar(outerVar), from(from), to(to), edgeID(edgeID), bvID(bvID){ //,weight(weight){

    }

    Edge() :
            v(var_Undef), outerVar(var_Undef), from(-1), to(-1), edgeID(-1), bvID(-1){ //,weight(0){

    }
};
};

#endif /* GRAPHTHEORYTYPES_H_ */
