/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef ALLPAIRS_H_
#define ALLPAIRS_H_

#include <vector>
namespace dgl{
class AllPairs{
public:
	virtual ~AllPairs(){};

	virtual void addSource(int s)=0;

	/*virtual Minisat::vec<int> & getChanged()=0;
	virtual void clearChanged()=0;

*/
	virtual void update( )=0;


	virtual bool connected_unsafe(int from,int t)=0;
	virtual bool connected_unchecked(int from,int t)=0;
	virtual bool connected(int from,int t)=0;
	virtual int distance(int from,int t)=0;
	virtual int distance_unsafe(int from,int t)=0;
	//virtual int previous(int from,int t)=0;
	virtual void getPath(int source, int to, std::vector<int> & path_store)=0;
};
};

#endif /* REACH_H_ */
