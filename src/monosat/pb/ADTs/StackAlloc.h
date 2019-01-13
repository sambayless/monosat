/************************************************************************************[StackAlloc.h]
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

#ifndef StackAlloc_h
#define StackAlloc_h

//=================================================================================================
#include "Global.h"

namespace Monosat {
namespace PB {
template<class T>
struct Allocator {
    virtual T* alloc(int nwords) = 0;
};


// STACK ALLOCATOR 
//
// 'cap' is capacity, 'lim' is malloc limit (use 'malloc()' for larger sizer 
// than this -- unless enough space left on the current stack).
//
template<class T, int cap = 10000, int lim = cap / 10>
class StackAlloc : public Allocator<T> {
    T* data;
    StackAlloc* prev;
    int index;

    T* alloc_helper(int n);

    StackAlloc(T* d, StackAlloc* p, int i) : data(d), prev(p), index(i){}


    void init(void){
        data = xmalloc<T>(cap);
        index = 0;
        prev = NULL;
    }

public:
    StackAlloc(void){init();}

    virtual ~StackAlloc(){

    }

    T* alloc(int n){
        if(index + n <= cap){
            T* tmp = data + index;
            index += n;
            return tmp;
        }else return alloc_helper(n);
    }

    void freeAll(void);

    void clear(void){
        freeAll();
        init();
    }
};


template<class T, int cap, int lim>
T* StackAlloc<T, cap, lim>::alloc_helper(int n){
    if(n > lim){
        StackAlloc* singleton = new StackAlloc<T, cap, lim>(xmalloc<T>(n), prev, -1);
        prev = singleton;
        return singleton->data;
    }else{
        StackAlloc* copy = new StackAlloc<T, cap, lim>(data, prev, index);
        data = xmalloc<T>(cap);
        index = n;
        prev = copy;
        assert(n <= cap);
        return data;
    }
}


// Call this before allocator is destructed if you do not want to keep the data.
//
template<class T, int cap, int lim>
void StackAlloc<T, cap, lim>::freeAll(void){
    StackAlloc* tmp, * ptr;
    for(ptr = this; ptr != NULL; ptr = ptr->prev)
        xfree(ptr->data);
    for(ptr = prev; ptr != NULL;)
        tmp = ptr->prev, delete ptr, ptr = tmp;
}

}
}
//=================================================================================================
#endif
