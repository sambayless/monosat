/****************************************************************************************[IntMap.h]
Copyright (c) 2011, Niklas Sorensson
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
//This is the IntMap from Minisat's MTL, adapted to use std::vector.
#ifndef DGL_Minisat_IntMap_h
#define DGL_Minisat_IntMap_h

#include <cassert>
#include <vector>

namespace dgl {
namespace alg {
template<class T>
struct MkIndexDefault {
    int operator()(T t) const{return (int) t;}
};

template<class K, class V, class MkIndex = MkIndexDefault<K>>
class IntMap {
    std::vector<V> map;
    MkIndex index;
public:
    explicit IntMap(MkIndex _index = MkIndex()) : index(_index){}

    bool has(K k) const{return index(k) < map.size();}

    const V& operator[](K k) const{
        assert(has(k));
        return map[index(k)];
    }

    V& operator[](K k){
        assert(has(k));
        return map[index(k)];
    }

    const V* begin() const{return &map[0];}

    const V* end() const{return &map[map.size()];}

    V* begin(){return &map[0];}

    V* end(){return &map[map.size()];}

    void reserve(K key, V pad){
        int ns = index(key) + 1;
        if(map.size() < ns){map.resize(ns, pad);}
    }

    void reserve(K key){
        int ns = index(key) + 1;
        if(map.size() < ns){map.resize(ns);}
    }

    void insert(K key, V val, V pad){
        reserve(key, pad);
        operator[](key) = val;
    }

    void insert(K key, V val){
        reserve(key);
        operator[](key) = val;
    }

    void clear(bool dispose = false){map.clear(dispose);}

    void moveTo(IntMap& to){
        map.moveTo(to.map);
        to.index = index;
    }

    void copyTo(IntMap& to) const{
        map.copyTo(to.map);
        to.index = index;
    }
};


template<class K=int, class MkIndex = MkIndexDefault<K>>
class IntSet {
    IntMap<K, char, MkIndex> in_set;
    std::vector<K> xs;

public:
    // Size operations:
    int size(void) const{return xs.size();}

    void clear(){
        for(int i = 0; i < xs.size(); i++)
            in_set[xs[i]] = 0;
        xs.clear();
    }

    // Stack interface:

    void push(const K& elem){
        insert(elem);
    }

    void pop(void){
        assert(size() > 0);
        assert(in_set[xs.last()]);
        in_set[xs.last()] = 0;
        xs.pop_back();
    }


    const K& last(void) const{
        return xs.last();
    }

    K& last(void){
        return xs.last();
    }

    // Allow inspecting the internal vector:
    const std::vector<K>&
    toVec() const{return xs;}

    // Vector interface:
    K operator[](int index) const{return xs[index];}


    void insert(K k){
        in_set.reserve(k, 0);
        if(!in_set[k]){
            in_set[k] = 1;
            xs.push_back(k);
        }
        assert(this->has(k));
    }

    void insertAll(std::vector<K>& from){
        for(K& l:from){
            insert(l);
        }
    }

    bool has(K k){
        in_set.reserve(k, 0);
        return in_set[k];
    }

    bool contains(K k){return has(k);}

    //stl-style begin and end, to support C++11 range-based for loops
    K* begin(){
        return &xs[0];
    }

    K* end(){
        return &xs[xs.size()];
    }

};

//=================================================================================================
}
}
#endif
