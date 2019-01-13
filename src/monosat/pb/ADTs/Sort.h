/******************************************************************************************[Sort.h]
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

/**************************************************************************************************

Template based sorting routines: sort, sortUnique (remove duplicates). Can be applied either on
'vec's or on standard C arrays (pointers).  

**************************************************************************************************/


#ifndef PB_Sort_h
#define PB_Sort_h

//#include <cstdlib>
#include "monosat/mtl/Vec.h"
#include "monosat/mtl/Rnd.h"
//=================================================================================================
namespace Monosat {
namespace PB {

template<class T>
struct LessThan_default {
    bool operator()(T x, T y){return x < y;}
};


//=================================================================================================


template<class T, class LessThan>
void selectionSort(T* array, int size, LessThan lt){
    int i, j, best_i;
    T tmp;

    for(i = 0; i < size - 1; i++){
        best_i = i;
        for(j = i + 1; j < size; j++){
            if(lt(array[j], array[best_i]))
                best_i = j;
        }
        tmp = array[i];
        array[i] = array[best_i];
        array[best_i] = tmp;
    }
}

template<class T>
static inline void selectionSort(T* array, int size){
    PB::selectionSort(array, size, LessThan_default<T>());
}


template<class T, class LessThan>
void sort(T* array, int size, LessThan lt, double& seed){
    if(size <= 15)
        PB::selectionSort(array, size, lt);

    else{
        T pivot = array[Monosat::irand(seed, size)];
        T tmp;
        int i = -1;
        int j = size;

        for(;;){
            do i++;while(lt(array[i], pivot));
            do j--;while(lt(pivot, array[j]));

            if(i >= j) break;

            tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }

        PB::sort(array, i, lt, seed);
        PB::sort(&array[i], size - i, lt, seed);
    }
}

template<class T, class LessThan>
void sort(T* array, int size, LessThan lt){
    double seed = 91648253;
    PB::sort(array, size, lt, seed);
}

template<class T>
static inline void sort(T* array, int size){
    PB::sort(array, size, LessThan_default<T>());
}


template<class T, class LessThan>
void sortUnique(T* array, int& size, LessThan lt){
    int i, j;
    T last;

    if(size == 0) return;

    PB::sort(array, size, lt);

    i = 1;
    last = array[0];
    for(j = 1; j < size; j++){
        if(lt(last, array[j])){
            last = array[i] = array[j];
            i++;
        }
    }

    size = i;
}

template<class T>
static inline void sortUnique(T* array, int& size){
    sortUnique(array, size, LessThan_default<T>());
}


//=================================================================================================
// For 'vec's:


template<class T, class LessThan>
void sort(Monosat::vec<T>& v, LessThan lt){
    PB::sort((T*) v, v.size(), lt);
}

template<class T>
void sort(Monosat::vec<T>& v){
    PB::sort(v, LessThan_default<T>());
}


template<class T, class LessThan>
void sortUnique(Monosat::vec<T>& v, LessThan lt){
    int size = v.size();
    T* data = v.release();
    sortUnique(data, size, lt);
    v.~vec();
    new(&v) Monosat::vec<T>(data, size);
}

template<class T>
void sortUnique(Monosat::vec<T>& v){
    sortUnique(v, LessThan_default<T>());
}

}
}
//=================================================================================================
#endif
