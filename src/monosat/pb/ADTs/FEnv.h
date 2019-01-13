/******************************************************************************************[FEnv.h]
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

#ifndef FEnv_h
#define FEnv_h

//=================================================================================================


#define ENV FEnv
#define FML Formula
#define UNDEF_INDEX -2

#include "Map.h"
#include "VecMaps.h"

namespace Monosat {
namespace PB {
enum Op {
    op_And, op_Equiv
};

//-------------------------------------------------------------------------------------------------

// Layout: 31..4=index, 3=sign, 2=compo, 1..0=unused
class FML {
    unsigned data;
public:
    FML(unsigned d) : data(d){}

    FML(bool compo = false, int index = UNDEF_INDEX, bool sign = false) :
            data(((unsigned) index << 4) | ((unsigned) sign << 3) | ((unsigned) compo << 2)){}

    operator unsigned(void) const{return data;}
};

//-------------------------------------------------------------------------------------------------

namespace ENV {
struct NodeData {
    unsigned data0;
    unsigned data1;
    unsigned data2;

    NodeData(unsigned d0, unsigned d1, unsigned d2) : data0(d0), data1(d1), data2(d2){}

//      unsigned hash(void) const { return data0 ^ data1 ^ data2; }
    unsigned hash(void) const{return data0 ^ ((data1 << 16) | (data1 >> 16)) ^ data2;}

    bool operator==(const NodeData& other) const{
        return data0 == other.data0 && data1 == other.data1 && data2 == other.data2;
    }
};

extern vec<NodeData> nodes;
extern Map<NodeData, int> uniqueness_table;
}

//-------------------------------------------------------------------------------------------------

macro bool operator==(FML f, FML g){return (unsigned) f == (unsigned) g;}

macro bool operator!=(FML f, FML g){return (unsigned) f != (unsigned) g;}

macro bool operator<(FML f, FML g){return (unsigned) f < (unsigned) g;}

template<>
struct Hash<FML> {
    unsigned operator()(FML f) const{return (unsigned) f;}
};

macro FML neg(FML f){return FML((unsigned) f ^ 8);}

macro FML unsign(FML f){return FML((unsigned) f & ~8);}

macro FML id(FML f, bool sign){return FML((unsigned) f ^ ((unsigned) sign << 3));}

macro bool sign(FML f){return ((unsigned) f & 8) != 0;}

macro bool compo(FML f){return ((unsigned) f & 4) != 0;}

macro int index(FML f){return (int) (unsigned) f >> 4;}

macro int sindex(FML f){return (int) (unsigned) f >> 3;}

macro FML FormulaC(int index, bool sign = false){return FML(true, index, sign);}

macro FML FormulaA(int index, bool sign = false){return FML(false, index, sign);}

namespace ENV {
macro FML Comp_new(int index, bool sign = false){return FML(true, index, sign);}

macro FML Atom_new(int index, bool sign = false){return FML(false, index, sign);}
}

//-------------------------------------------------------------------------------------------------

#ifndef tag_Atom
#define tag_Atom (-1)
#endif
#define tag_Bin  0
#define tag_ITE  1
#define tag_FA   2

macro int ctag(FML f){return (int) (ENV::nodes[index(f)].data0 & 3);}

macro int ctag(int index){return (int) (ENV::nodes[index].data0 & 3);}

macro int tag(FML f){return compo(f) ? ctag(f) : tag_Atom;}

macro bool Atom_p(FML f){return !compo(f);}

macro bool Bin_p(FML f){return tag(f) == tag_Bin;}

macro bool ITE_p(FML f){return tag(f) == tag_ITE;}

macro bool FA_p(FML f){return tag(f) == tag_FA;}

//-------------------------------------------------------------------------------------------------

macro Op op(FML f){return (Op) (ENV::nodes[index(f)].data1 & 0x3);}

macro FML left(FML f){return (FML) (ENV::nodes[index(f)].data1 & 0xFFFFFFFC);}

macro FML right(FML f){return (FML) (ENV::nodes[index(f)].data2 & 0xFFFFFFFC);}

macro FML cond(FML f){return (FML) (ENV::nodes[index(f)].data0 & 0xFFFFFFFC);}

macro FML tt(FML f){return (FML) (ENV::nodes[index(f)].data1 & 0xFFFFFFFC);}

macro FML ff(FML f){return (FML) (ENV::nodes[index(f)].data2 & 0xFFFFFFFC);}

macro bool isCarry(FML f){return (bool) (ENV::nodes[index(f)].data1 & 0x1);}

macro FML FA_x(FML f){return (FML) (ENV::nodes[index(f)].data0 & 0xFFFFFFFC);}

macro FML FA_y(FML f){return (FML) (ENV::nodes[index(f)].data1 & 0xFFFFFFFC);}

macro FML FA_c(FML f){return (FML) (ENV::nodes[index(f)].data2 & 0xFFFFFFFC);}

//-------------------------------------------------------------------------------------------------

macro ENV::NodeData Bin_newData(Op op, FML left, FML right){
#ifdef PARANOID
    assert((unsigned)op < 0x4U);
    assert(((unsigned)left & 0x3) == 0);
    assert(((unsigned)right & 0x3) == 0);
#endif
    return ENV::NodeData(tag_Bin, (unsigned) op | (unsigned) left, (unsigned) right);
}

macro ENV::NodeData ITE_newData(FML cond, FML tt, FML ff){
#ifdef PARANOID
    assert(((unsigned)cond & 0x3) == 0);
    assert(((unsigned)tt & 0x3) == 0);
    assert(((unsigned)ff & 0x3) == 0);
#endif
    return ENV::NodeData(tag_ITE | (unsigned) cond, (unsigned) tt, (unsigned) ff);
}

macro ENV::NodeData FA_newData(bool isCarry, FML FA_x, FML FA_y, FML FA_c){
#ifdef PARANOID
    assert((unsigned)isCarry < 0x2U);
    assert(((unsigned)FA_x & 0x3) == 0);
    assert(((unsigned)FA_y & 0x3) == 0);
    assert(((unsigned)FA_c & 0x3) == 0);
#endif
    return ENV::NodeData(tag_FA | (unsigned) FA_x, (unsigned) isCarry | (unsigned) FA_y, (unsigned) FA_c);
}

namespace ENV {
macro FML new_helper(ENV::NodeData node, bool sign){
    int index = ENV::nodes.size();
    ENV::nodes.push(node);
    return ENV::Comp_new(index, sign);
}

macro FML newS_helper(ENV::NodeData node, bool sign){
    int index;
    if(!uniqueness_table.peek(node, index)){
        index = ENV::nodes.size();
        ENV::nodes.push(node);
        uniqueness_table.set(node, index);
    }
    return ENV::Comp_new(index, sign);
}
}

macro FML Bin_new(Op op, FML left, FML right, bool sign = false, bool sharing = false){
    return sharing ? ENV::newS_helper(Bin_newData(op, left, right), sign)
                   : ENV::new_helper(Bin_newData(op, left, right), sign);
}

macro FML Bin_newS(Op op, FML left, FML right, bool sign = false){return Bin_new(op, left, right, sign, true);}

macro FML ITE_new(FML cond, FML tt, FML ff, bool sign = false, bool sharing = false){
    return sharing ? ENV::newS_helper(ITE_newData(cond, tt, ff), sign)
                   : ENV::new_helper(ITE_newData(cond, tt, ff), sign);
}

macro FML ITE_newS(FML cond, FML tt, FML ff, bool sign = false){return ITE_new(cond, tt, ff, sign, true);}

macro FML FA_new(bool isCarry, FML FA_x, FML FA_y, FML FA_c, bool sign = false, bool sharing = false){
    return sharing ? ENV::newS_helper(FA_newData(isCarry, FA_x, FA_y, FA_c), sign)
                   : ENV::new_helper(FA_newData(isCarry, FA_x, FA_y, FA_c), sign);
}

macro FML FA_newS(bool isCarry, FML FA_x, FML FA_y, FML FA_c, bool sign = false){
    return FA_new(isCarry, FA_x, FA_y, FA_c, sign, true);
}

//-------------------------------------------------------------------------------------------------

#define atom_Undef (-2)
#define atom_True  (-1)

const FML _undef_ = ENV::Atom_new(atom_Undef, false);
const FML _error_ = ENV::Atom_new(atom_Undef, true);
const FML _1_ = ENV::Atom_new(atom_True, false);
const FML _0_ = ENV::Atom_new(atom_True, true);

// Must only be used for atoms:
macro bool Exeception_p(FML f){return index(f) == atom_Undef;}

macro bool Exeception_p(int index){return index == atom_Undef;}

macro bool Const_p(FML f){return index(f) == atom_True;}

macro bool Const_p(int index){return index == atom_True;}

macro bool Var_p(FML f){return index(f) >= 0;}

macro bool Var_p(int index){return index >= 0;}

// Use indices >= 0 for variables!
//
macro Formula var(int index){return FormulaA(index);}


//-------------------------------------------------------------------------------------------------

namespace ENV {
template<class T, bool sgn = false>
struct AtomMap : public VecMap<T> {
    typedef FML Key;
    typedef T Datum;

    T at(FML f){return VecMap<T>::at(sgn ? sindex(f) + 4 : index(f) + 2);}

    void set(FML f, const T& value){VecMap<T>::set((sgn ? sindex(f) + 4 : index(f) + 2), value);}
};
}

namespace ENV {
template<class T, bool sgn = false>
class CompMap : DeckMap<T> {
    int offset;
public:
    typedef FML Key;
    typedef T Datum;

    CompMap(void) : DeckMap<T>(), offset(ENV::nodes.size()){}

    CompMap(T null) : DeckMap<T>(null), offset(ENV::nodes.size()){}

    T at(FML f){return DeckMap<T>::at((sgn ? sindex(f) : Monosat::PB::index(f)) - offset);}

    void set(FML f, T value){DeckMap<T>::set((sgn ? sindex(f) : Monosat::PB::index(f)) - offset, value);}
};

template<class T, bool sgn = false>
class FullMap {
    AtomMap<T, sgn> amap;
    CompMap<T, sgn> cmap;
public:
    typedef FML Key;
    typedef T Datum;

    T at(FML f){return compo(f) ? cmap.at(f) : amap.at(f);}

    void set(FML f, T value){compo(f) ? cmap.set(f, value) : amap.set(f, value);}
};
}

#define AMap FEnv::AtomMap
#define CMap FEnv::CompMap
#define FMap FEnv::FullMap

//-------------------------------------------------------------------------------------------------

macro Formula operator~(Formula f){
    return neg(f);
}

macro Formula operator&(Formula f, Formula g){
    if(f == _0_ || g == _0_) return _0_;
    else if(f == _1_) return g;
    else if(g == _1_) return f;
    else if(f == g) return f;
    else if(f == ~g) return _0_;

    if(g < f) swp(f, g);
    return Bin_newS(op_And, f, g);
}

macro Formula operator|(Formula f, Formula g){
    return ~(~f & ~g);
}

macro Formula operator^(Formula f, Formula g){
    if(f == _0_) return g;
    else if(f == _1_) return ~g;
    else if(g == _0_) return f;
    else if(g == _1_) return ~f;
    else if(f == g) return _0_;
    else if(f == ~g) return _1_;

    if(g < f) swp(f, g);
    return Bin_newS(op_Equiv, unsign(f), unsign(g), sign(f) == sign(g));
}

macro Formula operator&=(Formula& f, Formula g){return f = f & g;}

macro Formula operator|=(Formula& f, Formula g){return f = f | g;}

macro Formula operator^=(Formula& f, Formula g){return f = f ^ g;}


macro Formula FAs(Formula x, Formula y, Formula c)      // XOR of 3 arguments: x # y # c
{
    bool sgn = sign(x) ^sign(y) ^sign(c);
    x = unsign(x);
    y = unsign(y);
    c = unsign(c);

    if(x == _1_) return id(y ^ c, !sgn);
    if(y == _1_) return id(x ^ c, !sgn);
    if(c == _1_) return id(x ^ y, !sgn);
    if(x == y) return id(c, sgn);
    if(x == c) return id(y, sgn);
    if(y == c) return id(x, sgn);

    if(c < y) swp(c, y);
    if(y < x) swp(y, x);
    if(c < y) swp(c, y);

    return FA_newS(false, x, y, c, sgn);
}


macro Formula FAc(Formula x, Formula y, Formula c)      // x + y + c >= 2
{
    if(x == _0_) return y & c;
    if(x == _1_) return y | c;
    if(y == _0_) return x & c;
    if(y == _1_) return x | c;
    if(c == _0_) return x & y;
    if(c == _1_) return x | y;

    if(x == y) return x;
    if(x == c) return c;
    if(y == c) return y;
    if(x == ~y) return c;
    if(x == ~c) return y;
    if(y == ~c) return x;

    if(c < y) swp(c, y);
    if(y < x) swp(y, x);
    if(c < y) swp(c, y);

    bool s = sign(c);
    return FA_newS(true, id(x, s), id(y, s), unsign(c), s);
}


macro Formula ITE(Formula c, Formula t, Formula f){
    if(c == _0_) return f;
    if(c == _1_) return t;
    if(t == f) return t;
    if(t == ~f) return c ^ f;
    if(t == _0_ || t == ~c) return ~c & f;
    if(t == _1_ || t == c) return c | f;
    if(f == _0_ || f == c) return c & t;
    if(f == _1_ || f == ~c) return ~c | t;

    if(t < f)
        swp(t, f),
                c = ~c;

    return ITE_newS(c, id(t, sign(f)), unsign(f), sign(f));
}

//-------------------------------------------------------------------------------------------------

bool eval(Formula f, AMap<char>& values);

namespace FEnv {
extern vec<int> stack;

macro void clear(){
    nodes.clear();
    uniqueness_table.clear();
}

macro void push(){stack.push(nodes.size());}

macro void pop(){
    while(nodes.size() > stack.last())
        uniqueness_table.remove(nodes.last()),
                nodes.pop();
}

macro void keep(){stack.pop();}

macro int topSize(){
    return (stack.size() == 0) ? nodes.size() : nodes.size() - stack.last();
}
}

}
}
//-------------------------------------------------------------------------------------------------

#undef ENV
#undef FML
#undef UNDEF_INDEX


//=================================================================================================

#endif
