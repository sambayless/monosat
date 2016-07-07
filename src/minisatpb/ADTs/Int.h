/*******************************************************************************************[Int.h]
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

#ifndef Int_h
#define Int_h
#include "Global.h"
//=================================================================================================
namespace Monosat {
namespace PB {

struct Exception_IntOverflow {
    char *where;

    Exception_IntOverflow(char *w) : where(w) { }    // (takes ownership of string)
};
}
}

#ifdef NO_GMP
//=================================================================================================
// Fake bignums using 'int64':
//=================================================================================================

namespace Monosat {
namespace PB {
#define Int_Max__   9223372036854775807LL
#define Int_Min__   (-Int_Max__)
#define Int_Undef__ (-Int_Max__ - 1LL)

//-------------------------------------------------------------------------------------------------

#define A1 assert(data != Int_Undef__);
#define A2 assert(data != Int_Undef__), assert(other.data != Int_Undef__);

// NOTE! This is not a proper abstraction of big numbers. It just includes the operators
// used in 'PbSolver'. It should be easy enough to add the missing operators on demand.
//
class Int {
    int64   data;
public:
    Int() : data(Int_Undef__) {}
    Int(int   x) : data(x) {}
    Int(int64 x) : data(x) {}

    // "operator =" and copy-constructor "Int(const Int& src)" are default defined to the right thing.

    uint hash() const {A1 return (uint)data ^ (uint)(data >> 32); }

    bool operator == (Int other) const {A2 return data == other.data; }
    bool operator != (Int other) const {A2 return data != other.data; }
    bool operator <  (Int other) const {A2 return data <  other.data; }
    bool operator >  (Int other) const {A2 return data >  other.data; }
    bool operator <= (Int other) const {A2 return data <= other.data; }
    bool operator >= (Int other) const {A2 return data >= other.data; }

    Int  operator &  (Int other) const {A2 return Int(data & other.data); }
    Int& operator >>= (int n) {A1 data >>= n; return *this; }

    Int  operator -  ()          const {A1 return Int(-data); }
    Int& operator ++ ()                {A1 ++data; return *this; }
    Int& operator -= (Int other)       {A2 data -= other.data; return *this; }
    Int& operator += (Int other)       {A2 data += other.data; return *this; }
    Int& operator *= (Int other)       {A2 data *= other.data; return *this; }
    Int& operator /= (Int other)       {A2 data /= other.data; return *this; }
    Int  operator +  (Int other) const {A2 return Int(data + other.data); }
    Int  operator -  (Int other) const {A2 return Int(data - other.data); }
    Int  operator *  (Int other) const {A2 return Int(data * other.data); }
    Int  operator /  (Int other) const {A2 return Int(data / other.data); }
    Int  operator %  (Int other) const {A2 return Int(data % other.data); }

    friend char* toString(Int num) { char buf[32]; sprintf(buf, "%lld", num.data); return xstrdup(buf); }   // Caller must free string.
    friend int   toint   (Int num) { if (num > INT_MAX || num < INT_MIN) throw Exception_IntOverflow(xstrdup("toint")); return (int)num.data; }
};


#define Int_MAX Int(Int_Max__)
#define Int_MIN Int(Int_Min__)

#undef A1
#undef A2
}
}

#else
//=================================================================================================
// Real bignums using "GNU Multiple Precision Arithmetic Library"
//=================================================================================================
#include <cstddef>
#include "gmp.h"

//=================================================================================================
namespace Monosat {
namespace PB {

#define A1 assert(!small());
#define A2 assert(!small()); assert(!other.small());

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define Int_MIN Int((mpz_t*)-1)
#define Int_MAX Int((mpz_t*)1)


class Int {
    mpz_t *data;       // This pointer is meant to contain small integers when bit 0 is set (for efficiency).
    // Currently the only small integers used are the special values 'Int_MIN' and 'Int_MAX'.
    bool small() const { return ((intp) data & 1) != 0; }

public:
    // Constructors/Destructor (+assignment operator)
    //
    Int(mpz_t *d) : data(d) { }      // Low-level constructor -- don't use!

    Int() {
        data = xmalloc<mpz_t>(1);
        assert(((intp) data & 1) == 0);
        mpz_init(*data);
    }

    Int(int x) {
        data = xmalloc<mpz_t>(1);
        assert(((intp) data & 1) == 0);
        mpz_init_set_si(*data, x);
    }

    Int(const Int &src) {
        if (src.small())
            data = src.data;
        else {
            data = xmalloc<mpz_t>(1);
            assert(((intp) data & 1) == 0);
            mpz_init_set(*data, *src.data);
        }
    }

    ~Int() {
        if (!small()) {
            mpz_clear(*data);
            xfree(data);
        }
        data = 0;
    }

    Int &operator=(const Int &other) {
        if (&other != this) {
            if (other.small()) {
                this->~Int();
                data = other.data;
            } else {
                if (small()) {
                    data = xmalloc<mpz_t>(1);
                    assert(((intp) data & 1) == 0);
                    mpz_init_set(*data, *other.data);
                } else
                    mpz_set(*data, *other.data);
            }
        }
        return *this;
    }

    // Operators:
    //

    // -- Comparison (supports infinity)
    //    '+oo' and '-oo' are treated as two unique points beyond the integers. For instanse '+oo' is not < than itself, but <= than itself.
    bool operator==(const Int &other) const {
        if (small())
            return other.small() ? (data == other.data) : false;
        else
            return other.small() ? false : mpz_cmp(*data, *other.data) == 0;
    }

    bool operator<(const Int &other) const {
        if (small()) {
            if (data == Int_MIN.data)
                return (!other.small() || other.data != Int_MIN.data);
            else {
                assert(data == Int_MAX.data);
                return false;
            }
        } else {
            if (other.small()) {
                if (other.data == Int_MIN.data)
                    return false;
                else {
                    assert(other.data == Int_MAX.data);
                    return true;
                }
            } else
                return mpz_cmp(*data, *other.data) < 0;
        }
    }

    bool operator!=(const Int &other) const { return !(*this == other); }

    bool operator>=(const Int &other) const { return !(*this < other); }

    bool operator>(const Int &other) const { return other < *this; }

    bool operator<=(const Int &other) const { return !(*this > other); }


    // -- Arithmetic (not allowed on infinity except for unary '-')
    Int  operator+(const Int &other) const {
        A2
        Int ret;
        mpz_add(*ret.data, *data, *other.data);
        return ret;
    }

    Int  operator-(const Int &other) const {
        A2
        Int ret;
        mpz_sub(*ret.data, *data, *other.data);
        return ret;
    }

    Int  operator*(const Int &other) const {
        A2
        Int ret;
        mpz_mul(*ret.data, *data, *other.data);
        return ret;
    }

    Int  operator/(const Int &other) const {
        A2
        Int ret;
        mpz_tdiv_q(*ret.data, *data, *other.data);
        return ret;
    }

    Int  operator%(const Int &other) const {
        A2
        Int ret;
        mpz_tdiv_r(*ret.data, *data, *other.data);
        return ret;
    }

    Int &operator+=(const Int &other) {
        A2
        mpz_add(*data, *data, *other.data);
        return *this;
    }

    Int &operator-=(const Int &other) {
        A2
        mpz_sub(*data, *data, *other.data);
        return *this;
    }

    Int &operator*=(const Int &other) {
        A2
        mpz_mul(*data, *data, *other.data);
        return *this;
    }

    Int &operator/=(const Int &other) {
        A2
        mpz_tdiv_q(*data, *data, *other.data);
        return *this;
    }

    Int &operator%=(const Int &other) {
        A2
        mpz_tdiv_r(*data, *data, *other.data);
        return *this;
    }

    Int &operator++() { return *this += Int(1); }

    Int &operator--() { return *this -= Int(1); }

    Int operator-() const {
        if (small())
            return Int((mpz_t *) (-(intp) data));
        else {
            Int ret;
            mpz_neg(*ret.data, *data);
            return ret;
        }
    }

    // -- Bit operators (incomplete; we don't need more at the moment)
    Int  operator&(const Int &other) const {
        A2
        Int ret;
        mpz_and(*ret.data, *data, *other.data);
        return ret;
    }

    Int &operator>>=(int n) {
        A1
        mpz_fdiv_q_2exp(*data, *data, n);
        return *this;
    }

    // Methods:
    //
    friend char *toString(Int num) {
        if (num == Int_MIN) return xstrdup("-oo");
        else if (num == Int_MAX) return xstrdup("+oo");
        assert(!num.small());
        char *tmp = xmalloc<char>(mpz_sizeinbase(*num.data, 10) + 2);
        mpz_get_str(tmp, 10, *num.data);
        return tmp;
    }

    friend int toint(Int num) {
        if (num.small() || !mpz_fits_sint_p(*num.data))
            throw Exception_IntOverflow(xstrdup("toint"));
        return (int) mpz_get_si(*num.data);
    }

    uint hash() const {   // primitive hash function -- not good with bit-shifts
        mp_size_t size = mpz_size(*data);
        mp_limb_t val = 0;
        for (mp_size_t i = 0; i < size; i++) {
            mp_limb_t limb = mpz_getlimbn(*data, i);
            val ^= limb;
        }
#ifdef LP64
        return (uint) val | (uint) (val >> 32);
#else
        return (uint)val;
#endif
    }
};

}
}
//=================================================================================================
#endif

#endif
