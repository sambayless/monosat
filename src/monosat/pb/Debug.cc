/****************************************************************************************[Debug.cc]
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

#include "Debug.h"

namespace Monosat {
namespace PB {

vec<cchar*>* debug_names = nullptr;

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void dump(Int num){
    if(num == Int_MIN) reportf("-oo");
    else if(num == Int_MAX) reportf("+oo");
    else{
        char* tmp = toString(num);
        reportf("%s", tmp);
        xfree(tmp);
    }
}


void dump(Lit p){
    if(debug_names == nullptr)
        reportf("%sx%d", sign(p) ? "~" : "", var(p));
    else{
        if(var(p) < debug_names->size())
            reportf("%s%s", sign(p) ? "~" : "", (*debug_names)[var(p)]);
        else
            reportf("%s@%d", sign(p) ? "~" : "", var(p));
    }
}


void dump(Formula f){
    if(Atom_p(f))
        reportf("%sv%d", sign(f) ? "~" : "", index(f));
    else
        reportf("%s@%d", sign(f) ? "~" : "", index(f));
}


void dump(const vec<Lit>& ps, const vec<Int>& Cs){
    assert(ps.size() == Cs.size());
    for(int i = 0; i < ps.size(); i++){
        dump(Cs[i]);
        reportf("*");
        dump(ps[i]);
        if(i + 1 < ps.size()) reportf("  ");
    }
}

void dump(const vec<Formula>& ps, const vec<Int>& Cs){
    assert(ps.size() == Cs.size());
    for(int i = 0; i < ps.size(); i++){
        dump(Cs[i]);
        reportf("*");
        dump(ps[i]);
        if(i + 1 < ps.size()) reportf("  ");
    }
}


void dump(const vec<Lit>& ps, const vec<Int>& Cs, const vec<int>& assigns){
    assert(ps.size() == Cs.size());
    for(int i = 0; i < ps.size(); i++){
        dump(Cs[i]);
        reportf("*");
        dump(ps[i]);
        if(assigns[var(ps[i])] == toInt(l_Undef))
            reportf(":?");
        else if((assigns[var(ps[i])] == toInt(l_True) && !sign(ps[i])) ||
                (assigns[var(ps[i])] == toInt(l_False) && sign(ps[i])))
            reportf(":1");
        else
            reportf(":0");
        if(i + 1 < ps.size()) reportf("  ");
    }
}


void dump(const Linear& pb){
    for(int i = 0; i < pb.size; i++){
        dump(pb(i));
        reportf("*");
        dump(pb[i]);
        reportf("  ");
    }
    reportf("in [");
    dump(pb.lo);
    reportf(",");
    dump(pb.hi);
    reportf("]");
}


void dump(const Linear& pb, const vec<int>& assigns){
    for(int i = 0; i < pb.size; i++){
        dump(pb(i));
        reportf("*");
        dump(pb[i]);
        if(assigns[var(pb[i])] == toInt(l_Undef))
            reportf(":?");
        else if((assigns[var(pb[i])] == toInt(l_True) && !sign(pb[i])) ||
                (assigns[var(pb[i])] == toInt(l_False) && sign(pb[i])))
            reportf(":1");
        else
            reportf(":0");
        if(i + 1 < pb.size) reportf("  ");
    }
    reportf("in [");
    dump(pb.lo);
    reportf(",");
    dump(pb.hi);
    reportf("]");
}
}
}