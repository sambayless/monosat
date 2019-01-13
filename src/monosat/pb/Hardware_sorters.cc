/*****************************************************************************[Hardware_sorters.cc]
Copyright (c) 2005-2010, Niklas Een, Niklas Sorensson

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

#include "Hardware.h"

namespace Monosat {
namespace PB {
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

macro Formula operator&&(Formula f, Formula g){
    if(f == _0_ || g == _0_) return _0_;
    else if(f == _1_) return g;
    else if(g == _1_) return f;
    else if(f == g) return f;
    else if(f == ~g) return _0_;

    if(g < f) swp(f, g);
    return Bin_new(op_And, f, g);
}

macro Formula operator||(Formula f, Formula g){
    return ~(~f && ~g);
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


static inline void cmp2(vec<Formula>& fs, int begin){
    Formula a = fs[begin];
    Formula b = fs[begin + 1];
#if 1
    fs[begin] = a | b;
    fs[begin + 1] = a & b;
#else
    fs[begin]     = a || b;
    fs[begin + 1] = a && b;
#endif
}

static void riffle(vec<Formula>& fs){
    vec<Formula> tmp;
    fs.copyTo(tmp);
    for(int i = 0; i < fs.size() / 2; i++){
        fs[i * 2] = tmp[i];
        fs[i * 2 + 1] = tmp[i + fs.size() / 2];
    }
}

static void unriffle(vec<Formula>& fs){
    vec<Formula> tmp;
    fs.copyTo(tmp);
    for(int i = 0; i < fs.size() / 2; i++){
        fs[i] = tmp[i * 2];
        fs[i + fs.size() / 2] = tmp[i * 2 + 1];
    }
}

static void oddEvenMerge(vec<Formula>& fs, int begin, int end){
    assert(end - begin > 1);
    if(end - begin == 2)
        cmp2(fs, begin);
    else{
        int mid = (end - begin) / 2;
        vec<Formula> tmp;
        for(int i = 0; i < end - begin; i++)
            tmp.push(fs[begin + i]);
        unriffle(tmp);
        oddEvenMerge(tmp, 0, mid);
        oddEvenMerge(tmp, mid, tmp.size());
        riffle(tmp);
        for(int i = 1; i < tmp.size() - 1; i += 2)
            cmp2(tmp, i);
        for(int i = 0; i < tmp.size(); i++)
            fs[i + begin] = tmp[i];
    }
}

// Inputs to the circuit is the formulas in fs, which is overwritten
// by the resulting outputs of the circuit.
// NOTE: The number of comparisons is bounded by: n * log n * (log n + 1)
void oddEvenSort(vec<Formula>& fs){
    int orig_sz = fs.size();
    int sz;
    for(sz = 1; sz < fs.size(); sz *= 2);
    fs.growTo(sz, _0_);

    for(int i = 1; i < fs.size(); i *= 2)
        for(int j = 0; j + 2 * i <= fs.size(); j += 2 * i)
            oddEvenMerge(fs, j, j + 2 * i);
    fs.shrink(sz - orig_sz);
}
}
}