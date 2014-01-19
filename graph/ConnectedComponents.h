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

	bool marked;

	int stats_full_updates;
	int stats_fast_updates;
	int stats_fast_failed_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;
	int stats_num_skipable_deletions;
	double mod_percentage;

	double stats_full_update_time;
	double stats_fast_update_time;



	virtual ~ConnectedComponents(){};

	virtual void update( )=0;


	virtual int numComponents()=0;

	//Get an arbitrary element from the given set
	virtual int getElement(int set)=0;
	//Get the component this element belongs to
	virtual int getComponent(int node)=0;
	virtual bool connected(int from, int to)=0;
};

#endif /* MINIMUMSPANNINGTREE_H_ */
