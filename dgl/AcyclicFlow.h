#include <vector>
#include <cassert>
//#include "MaxFlow.h"
#include "alg/LinkCutCost.h"
namespace dgl {

//Implementation of the acyclic flow algorithm from Sleator and Tarjan (1983)
template<typename Weight> //, class Flow=MaxFlow<Weight> >
class AcyclicFlow{
	DynamicGraph<Weight>& g;
	LinkCutCost<Weight> forest;

	//This stores, for each node in the tree, the edge of g connecting it to its parent (or -1 if there is no such edge)
	std::vector<int> parent_edge;
	std::vector<int> parents;
	std::vector<int> disabled_edges;//disabled internally to the acyclic flow algorithm, to keep them from being revisited
public:
	AcyclicFlow(DynamicGraph<Weight>& g):g(g){

	}

	void getAcyclicFlow(int s, int t,	std::vector<Weight> & flow_store){
		assert(flow_store.size()==g.edges());
		while(parent_edge.size()<g.nodes()){
			forest.addNode();
			parent_edge.push_back(-1);
			parents.push_back(-1);
		}
#ifndef NDEBUG
		for(int i = 0;i<g.nodes();i++){
			assert(forest.isRoot(i));
			assert(parent_edge[i]==-1);
			assert(parents[i]==-1);
		}
		for(int i = 0;i<g.edges();i++){
			assert(flow_store[i]>=0);
		}
#endif

		disabled_edges.resize(g.edges());
		for (int edgeID =0;edgeID<g.edges();edgeID++){
			disabled_edges[edgeID]=false;
		}
		while(true){
			int v = forest.findRoot(s);//can we optimize this to findRoot(v) after the first iteration?
			bool done=true;
			//check if any edges with non-zero flow leave v
			for (int i = 0;i<g.nIncident(v);i++){
				int edgeID = g.incident(v,i).id;
				int to = g.incident(v,i).node;
				if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID) && !disabled_edges[edgeID] && flow_store[edgeID]>0){


					assert( g.getEdge(edgeID).from==v);
					assert(flow_store[edgeID]>0);

					if (forest.findRoot(to)==v){
						//step 3: cycle found
						//Weight flow = forest.minCost(to);
						int mincost_ancestor = forest.ancecstorFindMin(to);
						Weight flow = forest.getCost(mincost_ancestor);
						if (flow_store[edgeID]<flow){
							flow = flow_store[edgeID];
						}
						assert(flow_store[edgeID]>=flow);
						flow_store[edgeID]-=flow;
						forest.updateCostOfPathToRoot(to,-flow);
						assert(flow_store[edgeID]==0 || forest.getCost(forest.ancecstorFindMin(to))==0);//the flow cycle must have been broken

						//delete tree edges with no remaining flow

						//int u = forest.minCost(to);
						while(mincost_ancestor!=-1 &&  forest.getCost(mincost_ancestor)==0){
							assert( forest.getCost(mincost_ancestor)==0);
							int parent_edgeID = parent_edge[mincost_ancestor];
							if(parent_edgeID>=0){
								assert(g.hasEdge(parent_edgeID));
								assert(g.getEdge(parent_edgeID).to==mincost_ancestor);
								assert(flow_store[parent_edgeID]==flow);
								flow_store[parent_edgeID]=0;
								//record a current flow of zero?
							}
							forest.cut(mincost_ancestor);
							parent_edge[mincost_ancestor]=-1;
							parents[mincost_ancestor]=-1;
							//should also reduce the edge cost of (mincost_ancestor,parent(mincost_ancestor)) here...
							mincost_ancestor = forest.ancecstorFindWeight(to,0);
						}

					}else{

						assert(flow_store[edgeID]>0);
						forest.link(v,to,flow_store[edgeID]);
						parent_edge[v]=edgeID;
						parents[v]=to;
						done=false;
						break;//go back to step 1
					}
				}
			}
			if(done){
				//step 2: all paths from v to t are acyclic.
				if(v==s){
					//compute the current flow of every tree edge and stop
					for(int n = 0;n<g.nodes();n++){
						if (!forest.isRoot(n)){
							int edgeID = parent_edge[n];
							assert(edgeID>=0);
							Weight flow = forest.getCost(n);
							flow_store[edgeID]=flow;
							forest.cut(n);//can this flow extraction be done more efficiently?
							parent_edge[n]=-1;
							parents[n]=-1;
						}
					}
					return;
				}else{

					for (int i = 0;i<g.nIncoming(v);i++){
						int edgeID = g.incoming(v,i).id;
						int from = g.incoming(v,i).node;
						if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID) && !disabled_edges[edgeID] && flow_store[edgeID]>0){
							//delete from the graph every edge entering v.
							disabled_edges[edgeID]=true;
							if (!forest.isRoot(from) && parent_edge[from]==edgeID){
								//For each such edge (u,v) that is a tree edge, cut (u).
								//record the current flow
								Weight flow = forest.getCost(from);
								forest.cut(from);
								parent_edge[from]=-1;
								parents[from]=-1;

								flow_store[edgeID]=flow;//Is this correct?

							}
						}
					}
				}
			}
		}
	}
};

};
