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
class BVTheory {
public:
	virtual ~BVTheory() {
	}
	virtual int getTheoryIndexBV()=0;
	virtual void enqueueBV(int bvID)=0;
	virtual void backtrackBV(int bvID)=0;
	//Optional interface, if optimized BV clause learning is supported
	virtual void rewindBV(int bvID){

	}
};

};
#endif /* BVTHEORY_H_ */
