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
//Barebones bit-vector class.

#ifndef BITSET_H_
#define BITSET_H_

#include "monosat/mtl/Vec.h"
#include <limits>

namespace Monosat {

//=================================================================================================
#define BITSET_ELEMENT_SIZE (sizeof(uint64_t))

class Bitset {
    vec<uint64_t> buf;
    int sz;

public:
    Bitset() : sz(0){}

    Bitset(int size) : sz(size){buf.growTo(size / BITSET_ELEMENT_SIZE + 1);}

    Bitset(int size, bool default_value) : sz(size){
        buf.growTo(size / BITSET_ELEMENT_SIZE + 1);
        if(default_value){
            for(int i = 0; i < buf.size(); i++){
                buf[i] = std::numeric_limits<uint64_t>::max();
            }
        }
    }

    void memset(bool to){
        if(to){
            for(int i = 0; i < buf.size(); i++){
                buf[i] = std::numeric_limits<uint64_t>::max();
            }
        }else{
            for(int i = 0; i < buf.size(); i++){
                buf[i] = 0;
            }
        }
    }

    void growTo(int size){
        buf.growTo(size / BITSET_ELEMENT_SIZE + 1);
        sz = size;
    }

    void clear(bool dealloc = false){
        buf.clear(dealloc);
        sz = 0;
    }

    void zero(){
        for(int i = 0; i < buf.size(); i++)
            buf[i] = 0;
    }

    void ones(){
        for(int i = 0; i < buf.size(); i++)
            buf[i] = std::numeric_limits<uint64_t>::max();
    }

    void invert(){
        for(int i = 0; i < buf.size(); i++){
            buf[i] = ~buf[i];
        }
    }

    int size() const{return sz;}

private:
    static int popcount(uint64_t i){
        //http://stackoverflow.com/a/2709523
        //should replace this with a platform specific popcnt instruction
        i = i - ((i >> 1) & 0x5555555555555555UL);
        i = (i & 0x3333333333333333UL) + ((i >> 2) & 0x3333333333333333UL);
        return (int) ((((i + (i >> 4)) & 0xF0F0F0F0F0F0F0FUL) * 0x101010101010101UL) >> 56);
    }

    // Return index of lowest non-zero bit (or -1 if no bit is set)
    static int bit_pos64(uint64_t value) {
        // this can be done more quickly using intrinsics, etc
        int p = -1;
        while(value){
            value>>=1;
            p++;
        }
        return p;
    }
public:

    void swap(Bitset& b){
        buf.swap(b.buf);
        std::swap(sz, b.sz);
    }

    //population count (number of 1s in the bitset)
    int count() const{
        int n = 0;
        for(int i = 0; i < buf.size(); i++){
            uint64_t b = buf[i];
            n += popcount(b);
        }
        return n;
    }

    void copyFrom(const Bitset& from){
        from.buf.copyTo(buf);
        sz = from.size();
    }

    void copyTo(Bitset& to) const{
        buf.copyTo(to.buf);
        to.sz = sz;

    }

    inline const bool operator[](int index) const{
        int i = index / BITSET_ELEMENT_SIZE;
        int rem = index % BITSET_ELEMENT_SIZE;
        assert(i < buf.size());
        return buf[i] & (1 << rem);
    }

    inline void set(int index){
        assert(index < size());
        int i = index / BITSET_ELEMENT_SIZE;
        int r = index % BITSET_ELEMENT_SIZE;

        buf[i] |= (1 << r);

    }

    inline void clear(int index){
        assert(index < size());
        int i = index / BITSET_ELEMENT_SIZE;
        int r = index % BITSET_ELEMENT_SIZE;

        buf[i] &= ~(1 << r);

    }

    inline void toggle(int index){
        assert(index < size());
        int i = index / BITSET_ELEMENT_SIZE;
        int r = index % BITSET_ELEMENT_SIZE;

        buf[i] ^= (1 << r);

    }

    bool equals(Bitset& c){
        if(c.size() != size())
            return false;
        int max_i = size() / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){
            uint64_t b = c.buf[i];
            if(b != buf[i]){
                return false;
            }
        }
        return true;
    }

    //true if all 1's of c are 1's of this bitset
    bool contains(Bitset& c){
        int max = size();
        if(max > c.size()){
            max = c.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){

            uint64_t b = c.buf[i];
            if(b & ~buf[i]){
                return false;
            }
        }
        return true;
    }

    void Not(Bitset& out){
        out.clear();
        out.growTo(size());

        for(int i = 0; i < buf.size(); i++){
            uint64_t a = buf[i];
            out.buf[i] = ~a;
        }
        out.sz = sz;
    }

    void And(const Bitset& with){

        int max = size();
        if(max > with.size()){
            max = with.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){

            uint64_t b = with.buf[i];
            buf[i] &= b;
        }

    }

    void Or(const Bitset& with){

        int max = size();
        if(max > with.size()){
            max = with.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){

            uint64_t b = with.buf[i];
            buf[i] |= b;
        }

    }

    void And(const Bitset& with, Bitset& out){
        out.clear();
        out.growTo(size());
        int max = size();
        if(max > with.size()){
            max = with.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){
            uint64_t a = buf[i];
            uint64_t b = with.buf[i];
            out.buf[i] = a & b;
        }
        for(int i = max_i; i < buf.size(); i++){
            uint64_t a = buf[i];
            out.buf[i] = a;
        }
        out.sz = sz;
    }

    void Or(const Bitset& with, Bitset& out){
        out.clear();
        out.growTo(size());
        int max = size();
        if(max > with.size()){
            max = with.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){
            uint64_t a = buf[i];
            uint64_t b = with.buf[i];
            out.buf[i] = a | b;
        }
        for(int i = max_i; i < buf.size(); i++){
            uint64_t a = buf[i];
            out.buf[i] = a;
        }
        out.sz = sz;
    }

    void Xor(const Bitset& with, Bitset& out){
        out.clear();
        out.growTo(size());
        int max = size();
        if(max > with.size()){
            max = with.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){
            uint64_t a = buf[i];
            uint64_t b = with.buf[i];
            out.buf[i] = a ^ b;
        }
        for(int i = max_i; i < buf.size(); i++){
            uint64_t a = buf[i];
            out.buf[i] = a;
        }
        out.sz = sz;
    }

    bool GreaterThan(const Bitset& with){
        //Thanks to Tobias for this implementation
        int max = size();
        if(max != with.size()){
            return false;
        }
        if(size() == 0){
            return false;
        }
        int max_i = max / BITSET_ELEMENT_SIZE;
        for(int i = max_i; i >= 0; i--){
            if(buf[i] < with.buf[i]){
                return false;
            }
            if(buf[i] > with.buf[i]){
                return true;
            }
        }
        return false;
    }

    // Return index of an arbitrary bit that is set in this bitset
    int getIndexOfSetBit() const{
        int max = size();
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){
            uint64_t a = buf[i];
            if (a != 0){
                return i*BITSET_ELEMENT_SIZE + bit_pos64(a);
            }
        }
        return -1;
    }
    // Return index of a bit that is set in both bitsets, or -1 if there is no such bit
    int getIndexOfMutualBit(const Bitset & with) const{
        uint64_t tmp = 0;
        int max = size();
        if(max > with.size()){
            max = with.size();
        }
        int max_i = max / BITSET_ELEMENT_SIZE + 1;
        for(int i = 0; i < max_i; i++){
            uint64_t a = buf[i];
            uint64_t b = with.buf[i];
            tmp = (a & b);
            if (tmp != 0){
                return i*BITSET_ELEMENT_SIZE + bit_pos64(tmp);
            }
        }
        return -1;
    }
};
//=================================================================================================
}

#endif /* BITSET_H_ */
