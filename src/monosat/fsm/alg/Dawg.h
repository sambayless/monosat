//
// Created by sam on 29/05/16.
//

#ifndef MONOSAT_DAWG_H
#define MONOSAT_DAWG_H
#include <stdexcept>
#include "monosat/mtl/Vec.h"
namespace Monosat {
struct FSMDawg {
    bool removed = false;
    bool seen = false;
    int id = -1;
    /*struct DawgTransition{
        int letter=-1;
        Dawg * dest=nullptr;
    };
    vec<DawgTransition> used_transitions*/;
    Monosat::vec<FSMDawg *> transitions;

    void addTransition(int letter, FSMDawg *d) {
        transitions.growTo(letter + 1, nullptr);
        if (transitions[letter]) {
            throw std::runtime_error("Duplicate dawg transition");
        }
        transitions[letter] = d;
    }

    FSMDawg *getTransition(int letter) {
        if (letter >= transitions.size()) {
            return nullptr;
        } else {
            return transitions[letter];
        }
    }
};
};
#endif //MONOSAT_DAWG_H
