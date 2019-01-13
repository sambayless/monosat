/*******************************************************************************************[Alg.h]
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

#ifndef Minisat_Alg_h
#define Minisat_Alg_h

#include "monosat/mtl/Vec.h"
#include <cmath>

namespace Monosat {

/*
 Finite subsequences of the Luby-sequence:

 0: 1
 1: 1 1 2
 2: 1 1 2 1 1 2 4
 3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
 ...


 */

static double luby(double y, int x){
    assert(x >= 0);
    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for(size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1);

    assert(size > x);
    while(size - 1 != x){
        size = (size - 1) >> 1;
        seq--;
        //According to Coverity: size can be zero at this line, leading to a mod by zero...
        //However, since size must be >= x+1 above, and always > x in this loop, and x is positive, this is safe.
        x = x % size;
    }

    return pow(y, seq);
}
//=================================================================================================
// Useful functions on vector-like types:

//=================================================================================================
// Removing and searching for elements:
//

template<class V, class T>
static inline void remove(V& ts, const T& t){
    int j = 0;
    for(; j < ts.size() && ts[j] != t; j++);
    assert(j < ts.size());
    for(; j < ts.size() - 1; j++)
        ts[j] = ts[j + 1];
    ts.pop();
}

template<class V, class T>
static inline bool find(V& ts, const T& t){
    int j = 0;
    for(; j < ts.size() && ts[j] != t; j++);
    return j < ts.size();
}

//=================================================================================================
// Copying vectors with support for nested vector types:
//

// Base case:
template<class T>
static inline void copy(const T& from, T& to){
    to = from;
}

// Recursive case:
template<class T>
static inline void copy(const vec<T>& from, vec<T>& to, bool append = false){
    if(!append)
        to.clear();
    for(int i = 0; i < from.size(); i++){
        to.push();
        copy(from[i], to.last());
    }
}

template<class T>
static inline void append(const vec<T>& from, vec<T>& to){
    copy(from, to, true);
}

template<class T>
static inline vec<T>& reverse(vec<T>& from){
    int i = 0;
    int j = from.size();
    while((i != j) && (i != --j)){
        std::swap(from[i], from[j]);
        ++i;
    }
    return from;
}

//=================================================================================================
}

#endif
