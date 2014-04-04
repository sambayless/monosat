/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef CONNECTED_COMPONENTS_H_
#define CONNECTED_COMPONENTS_H_

#include "mtl/Vec.h"

class ConnectedComponents{
public:





	virtual ~ConnectedComponents(){};

	virtual void printStats(){

	}

	virtual void update( )=0;


	virtual int numComponents()=0;

	//Get an arbitrary element from the given set
	//virtual int getElement(int set)=0;
	//Get the component this element belongs to
	virtual int getComponent(int node)=0;
	virtual bool connected(int from, int to)=0;
};

#endif /* MINIMUMSPANNINGTREE_H_ */
