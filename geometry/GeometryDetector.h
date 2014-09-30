/****************************************************************************************[Solver.h]
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

#ifndef GEOMETRY_DETECTOR_H_
#define GEOMETRY_DETECTOR_H_
#include "core/SolverTypes.h"

#include "mtl/Vec.h"
using namespace Monosat;

class GeometryDetector{
public:
	int detectorID;

	int unassigned_positives=0;
	int unassigned_negatives=0;

	int getID(){
		return detectorID;
	}

	virtual bool propagate(vec<Lit> & conflict)=0;
	virtual void buildReason(Lit p, vec<Lit> & reason, CRef marker)=0;
	virtual bool checkSatisfied(){
		return true;
	}
	virtual void backtrackUntil(int level){

	}
	virtual void printStats(){

	}
	virtual void printSolution(){

	}
	virtual void backtrackUntil(Lit p){

	}

	virtual void preprocess(){

	}
	virtual Lit decide(){
		return lit_Undef;
	}
	//virtual vec<int> & getLitMap();
	GeometryDetector(int _detectorID):detectorID(_detectorID){};
	virtual ~GeometryDetector(){

	}
	virtual void addVar(Var v){
		unassigned_negatives++;
		unassigned_positives++;
	}
	virtual void setOccurs(Lit l, bool occurs){
		if(!occurs){
			 if(sign(l))
			   unassigned_negatives--;
			 else
			   unassigned_positives--;
		}else{
			 if(sign(l))
			   unassigned_negatives++;
			 else
			   unassigned_positives++;
		}
		assert(unassigned_positives>=0);
		assert(unassigned_negatives>=0);
	}
	virtual void assign(Lit l){
		 if(sign(l))
		   unassigned_negatives--;
		 else
		   unassigned_positives--;
		 assert(unassigned_positives>=0);
		 assert(unassigned_negatives>=0);

	}
	virtual	 void unassign(Lit l){
		 if(sign(l))
		   unassigned_negatives++;
		 else
		   unassigned_positives++;
	}
};



#endif /* GEOMETRY_DETECTOR_H_ */
