/*************************************************************************[PbSolver_convertSort.cc]
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

#include "PbSolver.h"
#include "Hardware.h"
#include "Debug.h"

namespace Monosat {
namespace PB {
//#define pf(format, args...) (reportf(format, ## args), fflush(stdout))
#define pf(format, args...) nothing()

void nothing(void){}


//=================================================================================================


//#define PickSmallest
#define ExpensiveBigConstants
#define AllDigitsImportant

int primes[] = {2, 3, 5, 7, 11, 13, 17};
//int primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
//int primes[] = { 2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };


static
void optimizeBase(vec<Int>& seq, int carry_ins, vec<Int>& rhs, int cost, vec<int>& base, int& cost_bestfound,
                  vec<int>& base_bestfound){
    if(cost >= cost_bestfound)
        return;

    // "Base case" -- don't split further, build sorting network for current sequence:
    int final_cost = 0;
    for(int i = 0; i < seq.size(); i++){
        if(seq[i] > INT_MAX)
            goto TooBig;
#ifdef ExpensiveBigConstants
        final_cost += toint(seq[i]);
#else
        int c; for (c = 1; c*c < seq[i]; c++);
        final_cost += c;
#endif
        if(final_cost < 0)
            goto TooBig;
    }
    if(cost + final_cost < cost_bestfound){
        base.copyTo(base_bestfound);
        cost_bestfound = cost + final_cost;
    }
    TooBig:;

    /**/static int depth = 0;

    // <<== could count 1:s here for efficiency

    vec<Int> new_seq;
    vec<Int> new_rhs;
#ifdef PickSmallest
    int p = -1;
    for (int i = 0; i < seq.size(); i++)
        if (seq[i] > 1){ p = seq[i]; break; }
    if (p != -1){
#else
    //int upper_lim = (seq.size() == 0) ? 1 : seq.last(); // <<== Check that sqRoot is an 'int' (no truncation of 'Int')
    //for (int i = 0; i < (int)elemsof(primes) && primes[i] <= upper_lim; i++){
    for(int i = 0; i < (int) elemsof(primes); i++){
        int p = primes[i];
#endif
        int rest = carry_ins;   // Sum of all the remainders.
        Int div, rem;

        /**/for(int n = depth; n != 0; n--) pf("  ");
        pf("prime=%d   carry_ins=%d\n", p, carry_ins);
        /**/for(int n = depth; n != 0; n--) pf("  ");
        pf("New seq:");
        for(int j = 0; j < seq.size(); j++){
            rest += toint(seq[j] % Int(p));
            div = seq[j] / Int(p);
            if(div > 0)
                //**/pf(" %d", div),
                new_seq.push(div);
        }
        /**/pf("\n");
        /**/for(int n = depth; n != 0; n--) pf("  ");
        pf("rest=%d\n", rest);

        /**/for(int n = depth; n != 0; n--) pf("  ");
        pf("New rhs:");
#ifdef AllDigitsImportant
        bool digit_important = true;
#else
        bool    digit_important = false;
#endif
        for(int j = 0; j < rhs.size(); j++){
            div = rhs[j] / p;
            if(new_rhs.size() == 0 || div > new_rhs.last()){
                rem = rhs[j] % p;
                /**/pf(" %d:%d", div, rem),
                        new_rhs.push(div);
                if(!(rem == 0 && rest < p) && !(rem > rest))
                    digit_important = true;
            }
        }
        pf("\n\n");

        base.push(p);
        /**/depth++;
        optimizeBase(new_seq, rest / p, new_rhs, cost + (digit_important ? rest : 0), base, cost_bestfound,
                     base_bestfound);
        /**/depth--;
        base.pop();

        new_seq.clear();
        new_rhs.clear();
    }
}


static
void optimizeBase(vec<Int>& seq, vec<Int>& rhs, int& cost_bestfound, vec<int>& base_bestfound){
    vec<int> base;
    cost_bestfound = INT_MAX;
    base_bestfound.clear();
    optimizeBase(seq, 0, rhs, 0, base, cost_bestfound, base_bestfound);
}


//=================================================================================================

#define lit2fml(p) id(var(var(p)),sign(p))


static
void buildSorter(vec<Formula>& ps, vec<int>& Cs, vec<Formula>& out_sorter){
    out_sorter.clear();
    for(int i = 0; i < ps.size(); i++)
        for(int j = 0; j < Cs[i]; j++)
            out_sorter.push(ps[i]);
    oddEvenSort(out_sorter); // (overwrites inputs)
}

static
void buildSorter(vec<Formula>& ps, vec<Int>& Cs, vec<Formula>& out_sorter){
    vec<int> Cs_copy;
    for(int i = 0; i < Cs.size(); i++)
        Cs_copy.push(toint(Cs[i]));
    buildSorter(ps, Cs_copy, out_sorter);
}


class Exception_TooBig : public std::exception {
};

static
void buildConstraint(vec<Formula>& ps, vec<Int>& Cs, vec<Formula>& carry, vec<int>& base, int digit_no,
                     vec<vec<Formula>>& out_digits, int max_cost){
    assert(ps.size() == Cs.size());

    if(FEnv::topSize() > max_cost) throw Exception_TooBig();
    /**
    pf("buildConstraint(");
    for (int i = 0; i < ps.size(); i++)
        pf("%d*%s ", Cs[i], (*debug_names)[index(ps[i])]);
    pf("+ %d carry)\n", carry.size());
    **/

    if(digit_no == base.size()){
        // Final digit, build sorter for rest:
        // -- add carry bits:
        for(int i = 0; i < carry.size(); i++)
            ps.push(carry[i]),
                    Cs.push(1);
        out_digits.push();
        buildSorter(ps, Cs, out_digits.last());

    }else{
        vec<Formula> ps_rem;
        vec<int> Cs_rem;
        vec<Formula> ps_div;
        vec<Int> Cs_div;

        // Split sum according to base:
        int B = base[digit_no];
        for(int i = 0; i < Cs.size(); i++){
            Int div = Cs[i] / Int(B);
            int rem = toint(Cs[i] % Int(B));
            if(div > 0){
                ps_div.push(ps[i]);
                Cs_div.push(div);
            }
            if(rem > 0){
                ps_rem.push(ps[i]);
                Cs_rem.push(rem);
            }
        }

        // Add carry bits:
        for(int i = 0; i < carry.size(); i++)
            ps_rem.push(carry[i]),
                    Cs_rem.push(1);

        // Build sorting network:
        vec<Formula> result;
        buildSorter(ps_rem, Cs_rem, result);

        // Get carry bits:
        carry.clear();
        for(int i = B - 1; i < result.size(); i += B)
            carry.push(result[i]);

        out_digits.push();
        for(int i = 0; i < B - 1; i++){
            Formula out = _0_;
            for(int j = 0; j < result.size(); j += B){
                int n = j + B - 1;
                if(j + i < result.size())
                    out |= result[j + i] & ((n >= result.size()) ? _1_ : ~result[n]);
            }
            out_digits.last().push(out);
        }

        buildConstraint(ps_div, Cs_div, carry, base, digit_no + 1, out_digits, max_cost); // <<== change to normal loop
    }
}

/*
Naming:
  - a 'base' is a vector of integers, stating how far you count at that position before you wrap to the next digit (generalize base).
  - A 'dig' is an integer representing a digit in a number of some base.
  - A 'digit' is a vector of formulas, where the number of 1:s represents a digit in a number of some base.
*/


static
void convert(Int num, vec<int>& base, vec<int>& out_digs){
    for(int i = 0; i < base.size(); i++){
        out_digs.push(toint(num % Int(base[i])));
        num /= Int(base[i]);
    }
    out_digs.push(toint(num));
}


// Compare number lexicographically to output digits from sorter networks.
// Formula is TRUE when 'sorter-digits >= num'.
//
static
Formula lexComp(int sz, vec<int>& num, vec<vec<Formula>>& digits){
    if(sz == 0)
        return _1_;
    else{
/**/
        pf("num    :");
        for(int i = 0; i < sz; i++) pf(" %d", num[i]);
        pf("\n");
        pf("#digits:");
        for(int i = 0; i < sz; i++) pf(" %d", digits[i].size());
        pf("\n");
/**/
        sz--;
        vec<Formula>& digit = digits[sz];
        int dig = num[sz];

        Formula gt = (digit.size() > dig) ? digit[dig] : _0_;       // This digit is greater than the "dig" of 'num'.
        Formula ge = (dig == 0) ? _1_ :
                     (digit.size() > dig - 1) ? digit[dig - 1]
                                              : _0_;   // This digit is greater than or equal to the "dig" of 'num'.

        /**/if(sz == 0) return ge;
        return gt | (ge & lexComp(sz, num, digits));
    }
}

static
Formula lexComp(vec<int>& num, vec<vec<Formula>>& digits){
    assert(num.size() == digits.size());
    return lexComp(num.size(), num, digits);
}


static
Formula buildConstraint(vec<Formula>& ps, vec<Int>& Cs, vec<int>& base, Int lo, Int hi, int max_cost){
    vec<Formula> carry;
    vec<vec<Formula>> digits;
    buildConstraint(ps, Cs, carry, base, 0, digits, max_cost);
    if(FEnv::topSize() > max_cost) throw Exception_TooBig();

    vec<int> lo_digs;
    vec<int> hi_digs;
    if(lo != Int_MIN)
        convert(lo, base, lo_digs);
    if(hi != Int_MAX)
        convert(hi + 1, base, hi_digs);   // (+1 because we will change '<= x' to '!(... >= x+1)'


    /*DEBUG
    pf("Networks:");
    for (int i = 0; i < digits.size(); i++)
        pf(" %d", digits[i].size());
    pf("\n");

    if (lo != Int_MIN){
        pf("lo=%d :", lo); for (int i = 0; i < lo_digs.size(); i++) pf(" %d", lo_digs[i]); pf("\n"); }
    if (hi != Int_MAX){
        pf("hi+1=%d :", hi+1); for (int i = 0; i < hi_digs.size(); i++) pf(" %d", hi_digs[i]); pf("\n"); }
    END*/

/*
Base:  (1)    8    24   480
       aaa bbbbbb ccc ddddddd
Num:    2    0     5     6
*/

    Formula ret = ((lo == Int_MIN) ? _1_ : lexComp(lo_digs, digits))
                  & ((hi == Int_MAX) ? _1_ : ~lexComp(hi_digs, digits));
    if(FEnv::topSize() > max_cost) throw Exception_TooBig();
    return ret;
}


/*
a7...a1   b
0001111   001111111  00111
  ^^         ^        ^

a5 | (a4 & (b7 | b6 & (c3)))

a4
~a5 -> b6
~a6 & ~b7 -> c3
...

>= 404
*/



// Will return '_undef_' if 'cost_limit' is exceeded.
//
Formula buildConstraint(const Linear& c, int max_cost){
    vec<Formula> ps;
    vec<Int> Cs;

    for(int j = 0; j < c.size; j++)
        ps.push(lit2fml(c[j])),
                Cs.push(c(j));

    vec<Int> dummy;
    int cost;
    vec<int> base;
    optimizeBase(Cs, dummy, cost, base);
    FEnv::push();

    Formula ret;
    try{
        ret = buildConstraint(ps, Cs, base, c.lo, c.hi, max_cost);
    }catch(Exception_TooBig){
        FEnv::pop();
        return _undef_;
    }

    if(opt_verbosity >= 1){
        reportf("Sorter-cost:%5d     ", FEnv::topSize());
        reportf("Base:");
        for(int i = 0; i < base.size(); i++) reportf(" %d", base[i]);
        reportf("\n");
    }
    FEnv::keep();
    return ret;
}
}
}