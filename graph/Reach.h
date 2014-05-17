/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef REACH_H_
#define REACH_H_

#include "mtl/Vec.h"

struct NullReachStatus{
	void setReachable(int u, bool reachable){

	}
	bool isReachable(int u) const{
		return false;
	}

	void setMininumDistance(int u, bool reachable, int distance){

	}
};
extern NullReachStatus nullReachStatus;

class Reach{
public:






	virtual ~Reach(){};

	virtual void setSource(int s)=0;
	virtual int getSource()=0;
	//virtual addSource(int s)=0;

	virtual void update( )=0;


	virtual bool connected_unsafe(int t)=0;
	virtual bool connected_unchecked(int t)=0;
	virtual bool connected( int t)=0;
	virtual int distance( int t)=0;
	virtual int distance_unsafe(int t)=0;
	virtual int previous( int node)=0;
	virtual int incomingEdge( int node)=0;
	//The maximum distance to compute up to.
	virtual void setMaxDistance(int maxDistance){

	}
};

#endif /* REACH_H_ */
