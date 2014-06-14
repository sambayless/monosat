/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef CONNECTED_COMPONENTS_H_
#define CONNECTED_COMPONENTS_H_

#include <vector>
#include <cstdio>
class ConnectedComponents{
public:

	struct NullConnectedComponentsStatus{
		void setConnected(int u, int v, bool connected){

		}
		void setComponents(int components){

		}
	};

	static NullConnectedComponentsStatus nullConnectedComponentsStatus;


	virtual ~ConnectedComponents(){};

	virtual void printStats(){

	}

	virtual void update( )=0;

	virtual void addConnectedCheck(int u, int v){
		fprintf(stderr,"Connect checks not supported, aborting\n");
		exit(3);
	}

	virtual int numComponents()=0;

	//Get an arbitrary element from the given set
	//virtual int getElement(int set)=0;
	//Get the component this element belongs to
	virtual int getComponent(int node)=0;
	virtual bool connected(int from, int to)=0;
};

#endif /* MINIMUMSPANNINGTREE_H_ */
