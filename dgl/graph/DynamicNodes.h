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
	int n_enabled=0;
	int n_actualNodes=0;
public:
	void addNode(int nodeID){
		assert(nodeID>=0);
		if(nodeStatus.size()<=nodeID)
			nodeStatus.resize(nodeID+1);
		nodeStatus[nodeID]=true;
		n_enabled++;
		n_actualNodes++;
	}

	int numEnabled()const{
		return n_enabled;
	}

	int nodes()const{
		return nodeStatus.size();
	}

	bool nodeEnabled(int n)const{
		return nodeStatus[n];
	}
	void setNodeEnabled(int n, bool enabled){
		if(nodeStatus[n]!=enabled){
			nodeStatus[n]=enabled;
			if(enabled){
				n_enabled++;
			}else{
				n_enabled--;
			}
		}
		assert(n_enabled>=0);
		assert(n_enabled<=n_actualNodes);
	}
};

#endif
