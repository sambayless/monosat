//
// Created by sam on 13/03/16.
//

#ifndef MONOSAT_NFAACCEPTOR_H
#define MONOSAT_NFAACCEPTOR_H

#include "monosat/mtl/Vec.h"
#include "monosat/fsm/alg/NFATypes.h"

namespace Monosat {
class NFAAcceptor {
public:
    NFAAcceptor(){

    }

    virtual ~NFAAcceptor(){

    }

    virtual void update() = 0;

    virtual void setTrackStringAcceptance(int str, int state, bool trackPositiveAcceptance,
                                          bool trackNegativeAcceptance) = 0;

    virtual int numUpdates() const = 0;;

    virtual bool acceptsString(int string, int state) = 0;

    virtual bool getPath(int string, int state, vec<NFATransition>& path) = 0;

    //Get an accepting path in the FSM, but don't include multiple copies of the same FSM transition
    virtual bool getAbstractPath(int string, int state, vec<NFATransition>& path, bool reversed = false) = 0;
};
};
#endif //MONOSAT_NFAACCEPTOR_H
