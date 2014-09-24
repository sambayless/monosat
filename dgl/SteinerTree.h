/*
 * Reach.h
 *
 *  Created on: 2013-07-26
 *      Author: sam
 */

#ifndef STEINERTREE_H_
#define STEINERTREE_H_

#include <vector>

namespace dgl{


template<typename Weight=int>
class SteinerTree{
public:

	struct NullStatus{

		void setMinimumSteinerTree(Weight min_weight){

			}

	};
	static NullStatus nullStatus;



	virtual ~SteinerTree(){};

	virtual void printStats(){

	}

	virtual void update( )=0;

	//Total weight of the steiner tree (or infinite, if the graph is disconnected)
	virtual Weight& weight()=0;

	virtual bool disconnected()=0;

	virtual void getSteinerTree(std::vector<int> & edges)=0;
};
template<typename Weight>
typename SteinerTree<Weight>::NullStatus  SteinerTree<Weight>::nullStatus;
};
#endif /* MINIMUMSPANNINGTREE_H_ */
