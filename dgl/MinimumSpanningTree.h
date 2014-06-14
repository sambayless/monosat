/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef MINIMUMSPANNINGTREE_H_
#define MINIMUMSPANNINGTREE_H_

#include <vector>


class MinimumSpanningTree{
public:

	struct NullStatus{
		void setMinimumSpanningTree(int min_weight){

		}

		void inMinimumSpanningTree(int edge, bool in_tree){

		}
	};
	static NullStatus nullStatus;



	virtual ~MinimumSpanningTree(){};

	virtual void printStats(){

	}

	virtual void update( )=0;

	virtual bool dbg_uptodate()=0;
	virtual bool dbg_mst()=0;
	//Total weight of the mst (or infinite, if the graph is disconnected)
	virtual int weight()=0;
	//Sum of the weight of the mst of each tree in the forest
	virtual int forestWeight()=0;
	virtual Minisat::vec<int> & getSpanningTree()=0;
	virtual int getParent(int node)=0;
	virtual int getParentEdge(int node)=0;
	virtual bool edgeInTree(int edgeid)=0;
	virtual int numComponents()=0;
	virtual int getComponent(int node)=0;
	virtual int getRoot(int component=0)=0;
};

#endif /* MINIMUMSPANNINGTREE_H_ */
