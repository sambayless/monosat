/*****************************************************************************************[FEnv.cc]
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

#include "FEnv.h"

namespace Monosat {
namespace PB {
namespace FEnv {
vec<NodeData> nodes;
Map<NodeData, int> uniqueness_table;

vec<int> stack;
}


//=================================================================================================


static
bool eval(Formula f, AMap<char>& values, CMap<char>& memo){
    if(Const_p(f))
        return !sign(f);
    else if(Atom_p(f))
        return sign(f) ? !values.at(f) : values.at(f);
    else{
        int ret = memo.at(f);
        if(ret == -1){
            if(Bin_p(f)){
                bool l = eval(left(f), values, memo);
                bool r = eval(left(f), values, memo);
                if(op(f) == op_And)
                    ret = l & r;
                else{
                    assert(op(f) == op_Equiv);
                    ret = (l ^ r) ^ 1;
                }
            }else if(ITE_p(f)){
                bool sel = eval(cond(f), values, memo);
                bool tru = eval(tt(f), values, memo);
                bool fal = eval(ff(f), values, memo);
                ret = sel ? tru : fal;
            }else{
                assert(FA_p(f));
                bool x = eval(FA_x(f), values, memo);
                bool y = eval(FA_y(f), values, memo);
                bool c = eval(FA_c(f), values, memo);
                if(isCarry(f))
                    ret = ((int) x + (int) y + (int) c) >= 2;
                else
                    ret = x ^ y ^ c;
            }
            memo.set(f, ret);
        }
        return sign(f) ? !ret : ret;
    }
}


bool eval(Formula f, AMap<char>& values){
    CMap<char> memo(-1);;
    return eval(f, values, memo);
}
}
}