/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef ALLPAIRS_H_
#define ALLPAIRS_H_

#include "mtl/Vec.h"
namespace Minisat{
class AllPairs{
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

	virtual void drawFull()=0;

	virtual ~AllPairs(){};

	virtual void addSource(int s)=0;

	/*virtual Minisat::vec<int> & getChanged()=0;
	virtual void clearChanged()=0;

*/
	virtual void update( )=0;

	virtual bool dbg_uptodate()=0;
	virtual bool dbg_path(int from,int to)=0;
	virtual bool connected_unsafe(int from,int t)const=0;
	virtual bool connected_unchecked(int from,int t)const=0;
	virtual bool connected(int from,int t)=0;
	virtual int distance(int from,int t)=0;
	virtual int distance_unsafe(int from,int t)=0;
	//virtual int previous(int from,int t)=0;
	virtual void getPath(int source, int to, vec<int> & path_store)=0;
};
};

#endif /* REACH_H_ */
