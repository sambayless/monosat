#ifndef VecAlloc_h
#define VecAlloc_h

#include <typeinfo>
#include "Global.h"
//=================================================================================================
namespace Monosat {
namespace PB {
template<class T, int chunk_size = 100>
class VecAlloc {

    union Slot {
        char data[sizeof(T)];     // (would have liked type 'T' here)
        Slot* next;
    };

    Slot* table;
    int index;
    Slot* recycle;
#ifdef DEBUG_PB
    int     nallocs;
#endif

    void newTable(void){
        Slot* t = xmalloc<Slot>(chunk_size);
        t[0].next = table;
        table = t;
        index = 1;
    }

public:
    VecAlloc(void){
        recycle = nullptr;
        table = nullptr;
#ifdef DEBUG_PB
        nallocs = 0;
#endif
        newTable();
    }

    ~VecAlloc(void){
#ifdef DEBUG_PB
        //if (nallocs != 0) fprintf(stderr, "WARNING! VecAlloc detected leak of %d unit(s) of type '%s'.\n", nallocs, typeid(T).name());
        if (nallocs != 0) fprintf(stderr, "WARNING! VecAlloc detected leak of %d unit(s) of size %lu.\n", nallocs, sizeof(T));
#endif
        Slot* curr, * next;
        curr = table;
        while(curr != nullptr)
            next = curr[0].next,
                    xfree(curr),
            curr = next;
    }

    T* alloc(void){
#ifdef DEBUG_PB
        nallocs++;
#endif
        if(recycle == nullptr){
            if(index >= chunk_size)
                newTable();
            return (T*) &table[index++];
        }else{
            T* tmp = (T*) recycle;
            recycle = (*recycle).next;
            return tmp;
        }
    }

    void free(T* ptr){
#ifdef DEBUG_PB
        nallocs--;
#endif
        ((Slot*) ptr)->next = recycle;
        recycle = (Slot*) ptr;
    }
};
}
}
//=================================================================================================

#endif
