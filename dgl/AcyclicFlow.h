#include <vector>
#include <cassert>
#include "alg/LinkCutCost.h"
namespace dgl {

//Implementation of the acyclic flow algorithm from Sleator and Tarjan (1983)
template<typename Weight, class FlowFunction>
class AcyclicFlow{
	DynamicGraph<Weight>& g;
	LinkCutCost forest;

	//This stores, for each node in the tree, the edge of g connecting it to its parent (or -1 if there is no such edge)
	std::vector<int> parent_edge;
	AcyclicFlow(DynamicGraph<Weight>& g):g(g){
		for(int n = 0;n<g.nodes();n++){
			forest.addNode();
			parent_edge.push_back(-1);
		}
	}
	void getAcyclicFlow(int s, int t,FlowFunction & f,	std::vector<int> & acyclic_flow){

		acyclic_flow.resize(g.edges());
		for (int edgeID =0;edgeID<g.edges();edgeID++){
			if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)){
				int flow =f(edgeID);//restricted to integer weight flows for, because the link cust cost implementation is restricted to ints.
				acyclic_flow[edgeID]=flow;
				assert(flow>=0);

			}
		}
		while(true){
			int v = forest.findRoot(s);//can we optimize this to findRoot(v) after the first iteration?
			//check if any edges with non-zero flow leave v
			for (int i = 0;i<g.nIncident(v);i++){
				int edgeID = g.incident(v,i).id;
				int to = g.incident(v,i).node;
				if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID) && acyclic_flow[edgeID]>0){


					assert( g.getEdge(edgeID).from==v);
					assert(acyclic_flow[edgeID]>0);

					if (forest.findRoot(to)==v){
						//step 3: cycle found
						int flow = forest.minCost(to);
						if (acyclic_flow[edgeID]<flow){
							flow = acyclic_flow[edgeID];
						}
						assert(acyclic_flow[edgeID]>=flow);
						acyclic_flow[edgeID]-=flow;
						forest.updateCostOfPathToRoot(to,-flow);
						assert(forest.minCost(to)==0 || acyclic_flow[edgeID]==0);//the flow cycle must have been broken

						//delete tree edges with no remaining flow
						int mincost_ancestor = forest.ancecstorFindWeight(to,0);
						//int u = forest.minCost(to);
						while(mincost_ancestor!=-1){
							assert( forest.getCost(mincost_ancestor)==0);
							int parent_edgeID = parent_edge[mincost_ancestor];
							if(parent_edgeID>=0){
								assert(g.hasEdge(parent_edgeID));
								assert(g.getEdge(parent_edgeID).to==mincost_ancestor);
								assert(acyclic_flow[parent_edgeID]==flow);
								acyclic_flow[parent_edgeID]=0;
								//record a current flow of zero?
							}
							forest.cut(mincost_ancestor);
							parent_edgeID[mincost_ancestor]=-1;
							//should also reduce the edge cost of (mincost_ancestor,parent(mincost_ancestor)) here...
							mincost_ancestor = forest.ancecstorFindWeight(to,0);
						}
						assert(forest.minCost(to)>0);
					}

					assert(acyclic_flow[edgeID]>0);
					forest.link(v,to,acyclic_flow[edgeID]);
					parent_edge[to]=edgeID;
					break;//go back to step 1
				}
			}
			//step 2: all paths from v to t are acyclic.
			if(v==s){
				//compute the current flow of every tree edge and stop
				for(int n = 0;n<g.nodes();n++){
					if (!forest.isRoot(n)){
						int edgeID = parent_edge[n];
						assert(edgeID>=0);
						int flow = forest.getCost(n);
						acyclic_flow[edgeID]=flow;
						forest.cut(n);//can this flow extraction be done more efficiently?
					}
				}
				return;
			}else{

				for (int i = 0;i<g.nIncoming(v);i++){
					int edgeID = g.incoming(v,i).id;
					int from = g.incident(v,i).node;
					if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID) && acyclic_flow[edgeID]>0){
						//delete from the graph every edge entering v.
						acyclic_flow[edgeID]=0;
						if (!forest.isRoot(from) && parent_edge[from]==edgeID){
							//For each such edge (u,v) that is a tree edge, cut (u).
							//record the current flow
							int flow = forest.getCost(from);
							forest.cut(from);
							parent_edge[from]=-1;
							assert(acyclic_flow[edgeID]>=from);
							acyclic_flow[edgeID]=from;//Is this correct?
						}
					}
				}
			}
		}
	}
};


