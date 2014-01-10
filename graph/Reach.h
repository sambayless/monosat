/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef REACH_H_
#define REACH_H_



class Reach{
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

	virtual ~Reach(){};

	virtual void setSource(int s)=0;
	virtual int getSource()=0;

	/*virtual Minisat::vec<int> & getChanged()=0;
	virtual void clearChanged()=0;

*/
	virtual void update( )=0;

	virtual bool dbg_uptodate()=0;
	virtual bool dbg_path(int to)=0;
	virtual bool connected_unsafe(int t)const=0;
	virtual bool connected_unchecked(int t)const=0;
	virtual bool connected(int t)=0;
	virtual int distance(int t)=0;
	virtual int distance_unsafe(int t)=0;
	virtual int previous(int t)=0;
};

#endif /* REACH_H_ */
