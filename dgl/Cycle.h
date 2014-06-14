/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef CYCLE_H_
#define CYCLE_H_

#include <vector>
namespace dgl{
class Cycle{
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

	virtual ~Cycle(){};

	virtual bool hasDirectedCycle()=0;
	virtual bool hasUndirectedCycle()=0;
	virtual void update( )=0;
	virtual std::vector<int> & getUndirectedCycle()=0;
	virtual std::vector<int> & getDirectedCycle()=0;
};
};
#endif /* REACH_H_ */
