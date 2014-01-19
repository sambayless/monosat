/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef MINIMUMSPANNINGTREE_H_
#define MINIMUMSPANNINGTREE_H_

#include "mtl/Vec.h"

class MinimumSpanningTree{
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



	virtual ~MinimumSpanningTree(){};



	virtual void update( )=0;

	virtual bool dbg_uptodate()=0;
	virtual bool dbg_mst()=0;
	virtual int weight()=0;
	virtual Minisat::vec<int> & getSpanningTree()=0;
	virtual int getParent(int node)=0;
	virtual int getParentEdge(int node)=0;
	virtual bool edgeInTree(int edgeid)=0;
	virtual int numComponents()=0;
	virtual int getComponent(int node)=0;
	virtual int getRoot(int component=0)=0;
};

#endif /* MINIMUMSPANNINGTREE_H_ */
