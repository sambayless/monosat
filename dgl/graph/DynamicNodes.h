/*
 * DynamicNodes.h
 *
 *  Created on: Jun 17, 2014
 *      Author: sam
 */
#ifndef DYNAMIC_NODES
#define DYNAMIC_NODES
#include <vector>
class DynamicNodes{
	std::vector<bool> nodeStatus;
public:
	void addNode(int nodeID){
		assert(nodeID>=0);
		if(nodeStatus.size()<nodeID)
			nodeStatus.resize(nodeID+1);
		nodeStatus[nodeID]=true;
	}

	int nodes()const{
		return nodeStatus.size();
	}

	bool nodeEnabled(int n)const{
		return nodeStatus[n];
	}
	void setNodeEnabled(int n, bool enabled){
		nodeStatus[n]=enabled;
	}
};

#endif
