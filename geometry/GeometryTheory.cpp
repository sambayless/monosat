/*
 * GeometryTheory.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "GeometryTheory.h"
#include "ConvexHullDetector.h"
using namespace Minisat;
void GeometryTheorySolver::createConvexHull(int id, const int D){

	GeometryDetector * d;
	int detectorID = detectors.size();

	//should probably swap this out for some fancy template recursion...
	switch (D){
		case 2:
			d = new ConvexHullDetector<2,double>(detectorID,this,drand(rnd_seed));break;
		case 3:
			//d = new ConvexHullDetector<3,double>(detectorID,this,drand(rnd_seed));break;
		default:
			assert(false);
			fprintf(stderr,"Unsupported number of dimensions (%d), exiting\n",D);fflush(stderr);exit(1);
	};
	detectors.push(d);
}
