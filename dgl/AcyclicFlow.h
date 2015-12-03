#ifndef ACYCLIC_FLOW_H
#define ACYCLIC_FLOW_H

#include <vector>
#include <cassert>
#include "DFSCycle.h"
#include "EdmondsKarpAdj.h"
//#include "MaxFlow.h"
#include "alg/LinkCutCost.h"
#include <sstream>
#include <stdexcept>
namespace dgl {

//Implementation of the acyclic flow algorithm from Sleator and Tarjan (1983)
template<typename Weight> //, class Flow=MaxFlow<Weight> >
class AcyclicFlow{
	DynamicGraph<Weight>& g;
	LinkCutCost<Weight> forest;
	std::vector<Weight> new_flow;
	//This stores, for each node in the tree, the edge of g connecting it to its parent (or -1 if there is no such edge)
	std::vector<int> parent_edge;
	std::vector<int> parents;
	//std::vector<int> disabled_edges;//disabled internally to the acyclic flow algorithm, to keep them from being revisited
public:
	AcyclicFlow(DynamicGraph<Weight>& g):g(g){

	}

	void getAcyclicFlow(int s, int t,	std::vector<Weight> & flow_store){
		assert(flow_store.size()==g.edges());
		if (s<0 || s>=g.nodes() || t<0 ||t>=g.nodes()){
			return;
		}

		if (g.outfile) {
			fprintf(g.outfile, "acyclic_flow %d %d:  ",s,t);
			std::stringstream ss;
			for (int i = 0;i<g.edges();i++){
				ss<< " " << flow_store[i] << " ";
			}
			fprintf(g.outfile, "%s\n",ss.str().c_str());
		}

#ifndef NDEBUG

		Weight expected_flow=-1;
		{
			DynamicGraph<Weight> test;
			//test.outfile = fopen("/tmp/TEST_ACYCLIC_FLOW", "w");
			test.addNodes(g.nodes());
			for (int i = 0;i<flow_store.size();i++){
				if(i==4025){
					int a =1;
				}
				Weight f = flow_store[i];
				bool h = g.hasEdge(i);
				bool e = g.edgeEnabled(i);
				if(flow_store[i]>0 && g.hasEdge(i) && g.edgeEnabled(i)){

					test.addEdge(g.getEdge(i).from,g.getEdge(i).to,i,flow_store[i]);
				}
			}
			//test.drawFull(true);
			EdmondsKarpAdj<Weight> mf(test,s,t);
			expected_flow = mf.maxFlow(s,t);
		}

		for(int v = 0;v<g.nodes();v++){
			Weight sum_in=0;
			Weight sum_out=0;

			for (int i = 0;i<g.nIncident(v);i++){
					int edgeID = g.incident(v,i).id;
					int to = g.incident(v,i).node;
					if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)  && flow_store[edgeID]>0){
						sum_out+=flow_store[edgeID];
					}
			}
			for (int i = 0;i<g.nIncoming(v);i++){
					int edgeID = g.incoming(v,i).id;
					int to = g.incoming(v,i).node;
					if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)  && flow_store[edgeID]>0){
						sum_in+=flow_store[edgeID];
					}
			}
			if(sum_in!=sum_out){
				if (v==s && sum_out-sum_in==expected_flow){
					continue;
				}
				if (v==t && sum_in-sum_out==expected_flow){
					continue;
				}
				throw std::invalid_argument("Bad flow for acyclic removal");
			}
		}

#endif
		new_flow.clear();
		new_flow.resize(flow_store.size(),0);
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

		/*disabled_edges.resize(g.edges());
		for (int edgeID =0;edgeID<g.edges();edgeID++){
			disabled_edges[edgeID]=false;
		}*/
		//&& !disabled_edges[edgeID]
		while(true){
			//step 1
			int v = forest.findRoot(s);//can we optimize this to findRoot(v) after the first iteration?
			bool done=true;
			if(v==586 || v==537 ){
				int a=1;
			}
			//check if any edges with non-zero flow leave v
			for (int i = 0;i<g.nIncident(v);i++){
				int edgeID = g.incident(v,i).id;
				int to = g.incident(v,i).node;
				if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)  && flow_store[edgeID]>0){
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
		/*				printf("digraph{\n");
						for(int j =0;j<g.nodes();j++){
							if(parent_edge[j]!=-1){
								printf("%d -> %d label=\"v%d\"\n",j,parents[j],parent_edge[j]);
							}
						}
						printf("\n}\n");*/

						//forest.print_forest();
						forest.updateCostOfPathToRoot(to,-flow);
						//forest.print_forest();

						assert(flow_store[edgeID]==0 || forest.getCost(forest.ancecstorFindMin(to))==0);//the flow cycle must have been broken

						//step 4:
						//delete tree edges with no remaining flow

						//int u = forest.minCost(to);
						while(mincost_ancestor!=-1 &&  forest.getCost(mincost_ancestor)==0){
							if(mincost_ancestor==1280){
								int a=1;

							}
							assert( forest.getCost(mincost_ancestor)==0);
							int parent_edgeID = parent_edge[mincost_ancestor];
							if(parent_edgeID>=0){
								assert(g.hasEdge(parent_edgeID));
								assert(g.getEdge(parent_edgeID).from==mincost_ancestor);
								assert(flow_store[parent_edgeID]>=flow);
								flow_store[parent_edgeID]=0;
								//record a current flow of zero?
							}
							forest.cut(mincost_ancestor);
							parent_edge[mincost_ancestor]=-1;
							parents[mincost_ancestor]=-1;
							//should also reduce the edge cost of (mincost_ancestor,parent(mincost_ancestor)) here...
							mincost_ancestor = forest.ancecstorFindWeight(to,0);
						}
						//forest.print_forest();
						done=false;
						break;//go back to step 1
					}else{
						if(edgeID==4025){
							int a =1;
						}
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
					//step 2a
					for(int n = 0;n<g.nodes();n++){
						if (!forest.isRoot(n)){
							int edgeID = parent_edge[n];
							if(edgeID==4025){
								int a =1;
							}
							assert(edgeID>=0);
							Weight flow = forest.getCost(n);
							//new_flow[edgeID]+=flow; //the paper suggests we should include these disconnected flows,
							//but it seems to me that any such flows are actually part of disconnected cycles that no longer are reachable from
							//source, and that we should as such leave them out.
							forest.cut(n);//can this flow extraction be done more efficiently?
							parent_edge[n]=-1;
							parents[n]=-1;
						}
					}
					for(int i = 0;i<flow_store.size();i++)
						flow_store[i]=new_flow[i];
#ifndef NDEBUG
					//verify that no cycle remains
					DynamicGraph<Weight> test;
					//test.outfile = fopen("/tmp/TEST_ACYCLIC_FLOW", "w");
					test.addNodes(g.nodes());
					for (int i = 0;i<flow_store.size();i++){
						if(i==4025){
							int a =1;
						}
						Weight f = flow_store[i];
						bool h = g.hasEdge(i);
						bool e = g.edgeEnabled(i);
						if(flow_store[i]>0 && g.hasEdge(i) && g.edgeEnabled(i)){
							Weight w = flow_store[i];
							if(g.getEdge(i).to==1280){
								int a =1;
							}
							test.addEdge(g.getEdge(i).from,g.getEdge(i).to,i,flow_store[i]);
						}
					}
					//test.drawFull(true);
					DFSCycle<Weight> cycle(test);
					if(cycle.hasDirectedCycle()){
						//fflush(g.outfile);
			     		throw std::logic_error( "Cycle remains in acyclic flow" );
					}
					EdmondsKarpAdj<Weight> mf(test,s,t);
					if(expected_flow>=0 && mf.maxFlow(s,t)!=expected_flow){
						//fflush(g.outfile);
						throw std::logic_error( "Flow is altered by cycle removal");
					}

					for(int v = 0;v<g.nodes();v++){
						Weight sum_in=0;
						Weight sum_out=0;

						for (int i = 0;i<g.nIncident(v);i++){
								int edgeID = g.incident(v,i).id;
								int to = g.incident(v,i).node;
								if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)  && flow_store[edgeID]>0){
									sum_out+=flow_store[edgeID];
								}
						}
						for (int i = 0;i<g.nIncoming(v);i++){
								int edgeID = g.incoming(v,i).id;
								int to = g.incoming(v,i).node;
								if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)  && flow_store[edgeID]>0){
									sum_in+=flow_store[edgeID];
								}
						}
						if(sum_in!=sum_out){
							if (v==s && sum_out-sum_in==expected_flow){
								continue;
							}
							if (v==t && sum_in-sum_out==expected_flow){
								continue;
							}
							//fflush(g.outfile);
							throw std::logic_error("Bad flow after acyclic removal");
						}
					}
#endif
					return;
				}else{
					//step 2b
					for (int i = 0;i<g.nIncoming(v);i++){
						int edgeID = g.incoming(v,i).id;
						int from = g.incoming(v,i).node;
						//&& !disabled_edges[edgeID]
						if (g.hasEdge(edgeID) && g.edgeEnabled(edgeID)  && flow_store[edgeID]>0){
							//according to the paper, you should delete from the graph every edge entering v here.
							//but I don't see how that can work, because you can delete edges to the sink that are not yet examined...
							//disabled_edges[edgeID]=true;
							if (!forest.isRoot(from) && parent_edge[from]==edgeID){
								//For each such edge (u,v) that is a tree edge, cut (u).
								//record the current flow
								Weight flow = forest.getCost(from);
								assert(flow_store[edgeID]>=flow);
								flow_store[edgeID]=0;
								forest.cut(from);
								parent_edge[from]=-1;
								parents[from]=-1;
								if(edgeID==4025){
									int a =1;
								}
								new_flow[edgeID]=flow;//Is this correct?

							}
						}
					}
				}
			}
		}
	}
};

};
#endif
