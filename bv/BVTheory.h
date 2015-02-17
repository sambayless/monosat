/*
 * BVTheory.h
 *
 *  Created on: Feb 17, 2015
 *      Author: sam
 */

#ifndef BVTHEORY_H_
#define BVTHEORY_H_

#include "core/Theory.h"

namespace Monosat {
/**
 * Abstract interface to SMT theory solvers, as accessed by the BV solver
 */
class BVTheory:public Theory {
public:
virtual enqueueBV(int bvID)=0;

}


#endif /* BVTHEORY_H_ */
