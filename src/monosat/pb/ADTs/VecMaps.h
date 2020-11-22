/***************************************************************************************[VecMaps.h]
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

#ifndef VecMaps_h
#define VecMaps_h

//=================================================================================================

// TODO: Adapt constructors of 'BitMap' and 'DeckMap' to those of VecMap.
namespace Monosat {
namespace PB {

class BitMap {
    vec<unsigned> v;
    bool bool_null;

    unsigned word(int index) const{return (unsigned) index / (sizeof(int) * 8);}

    unsigned mask(int index) const{return 1 << ((unsigned) index % (sizeof(int) * 8));}

public:
    typedef int Key;
    typedef bool Datum;

    BitMap(void) : bool_null(){}

    BitMap(bool null) : bool_null(null){}

    BitMap(bool null, int capacity) : v((capacity + (sizeof(int) * 8) - 1) / (sizeof(int) * 8), -(int) null),
                                      bool_null(null){}

    bool at(int index) const{
        if(word(index) >= (unsigned) v.size()) return bool_null;
        else return v[word(index)] & mask(index);
    }

    void set(int index, bool value){
        if(word(index) >= (unsigned) v.size()) assert(index >= 0), v.growTo(word(index) + 1, -(int) bool_null);
        if(value == false)
            v[word(index)] &= ~mask(index);
        else
            v[word(index)] |= mask(index);
    }

    void clear(void){v.clear();}
};


template<class T>
class VecMap {
    vec<T> v;
    T T_null;

public:
    typedef int Key;
    typedef T Datum;

    VecMap(void) : T_null(){}

    VecMap(T null) : T_null(null){}

    VecMap(T null, int capacity) : v(capacity, null), T_null(null){}

    T at(int index) const{
        if((unsigned) index >= (unsigned) v.size()) return T_null;
        else return v[index];
    }

    void set(int index, T value){
        if((unsigned) index >= (unsigned) v.size()) assert(index >= 0), v.growTo(index + 1, T_null);
        v[index] = value;
    }

    void clear(void){v.clear();}
};

template<>
struct VecMap<bool> : BitMap {
    VecMap(void){}

    VecMap(bool null) : BitMap(null){}

    VecMap(bool null, int capacity) : BitMap(null, capacity){}
};


template<class T>
class DeckMap {
    VecMap<T> pos, neg;

public:
    typedef int Key;
    typedef T Datum;

    DeckMap(void){}

    DeckMap(T null) : pos(null), neg(null){}

    DeckMap(T null, int capacity) : pos(null, capacity), neg(null, capacity){}

    T at(int index) const{
        if(index >= 0) return pos.at(index);else return neg.at(~index);
    }

    void set(int index, Datum value){
        if(index >= 0) pos.set(index, value);else neg.set(~index, value);
    }

    void clear(void){
        pos.clear();
        neg.clear();
    }
};
}
}
//=================================================================================================

#endif
