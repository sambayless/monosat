/****************************************************************************************[Solver.h]
The MIT License (MIT)

Copyright (c) 2014, Sam Bayless

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef DINICS_LINKCUT_H
#define DINICS_LINKCUT_H


#include "MaxFlow.h"
#include <vector>
#include "core/Config.h"
#include "EdmondsKarpAdj.h"
#include "alg/LinkCutCost.h"
#include <algorithm>
#include <climits>
//Implementation of the Link/Cut Tree (aka Dynamic Tree) max-flow algorithm from Sleator and Tarjan, 1983).
namespace dgl{
template< class Capacity>
class DinitzLinkCut:public MaxFlow<int>{

public:

	std::vector<int> F;

	struct LocalEdge{
		int from;
		int id;
		bool backward=false;
		LocalEdge(int from=-1, int id=-1, bool backward=false):from(from),id(id),backward(backward){

		}
	};
	int curflow;
    int last_modification;
    int last_deletion;
    int last_addition;
    int source=-1;
    int sink=-1;
    int history_qhead;
    int last_history_clear;
    //std::vector<LocalEdge> prev;
    //std::vector<int> M;
    std::vector<int> dist;
    std::vector<int> pos;//position in the combined forward and backward adjacency list of each node in the DFS.
    DynamicGraph& g;
    Capacity & capacity;
    int INF;
    int src;
    int dst;
    struct ParentEdge{
    	bool backward;
    	int edgeID;
    };
    std::vector<ParentEdge> parentEdge;

   /* struct EdgeInTree{
    	union{
    		char data;
    		struct{
			char in_tree:1;
			char in_tree_backward:1;
			char disabled:1;
			char backward_disabled:1;
    		};
    	};
    	EdgeInTree():data(0){}
    };
    std::vector<EdgeInTree> tree_edges;//for each edge.*/
    LinkCutCost forest;
    std::vector<int> Q;
    struct Link{
    	int u;
    	int v;
    	bool backward:1;
    	int edgeID:31;
    };
    std::vector<Link> toLink;
#ifdef DEBUG_MAXFLOW
	EdmondsKarpAdj<Capacity> ek;
#endif

	double augtime=0;
	double augtime_search=0;
	double augtime_cleanup=0;
	double bfstime=0;
	double totaltime=0;
	long stats_augmenting_rounds =0;
	long stats_rounds=0;
	long stats_backtracks = 0;
	long stats_avoided_backtracks=0;
	void printStats(){
		printf("Dinics Link Cut:\n");
		printf("Total time: %f\n",totaltime);
		printf("BFS Time: %f\n",bfstime);
		printf("Augmenting Path Time: %f (search: %f, cleanup: %f)\n", augtime, augtime_search,augtime_cleanup);
		printf("Rounds: %d, Augmenting Rounds: %d\n",stats_rounds,stats_augmenting_rounds);
		printf("Backtracks %d (%d avoided)\n",stats_backtracks,stats_avoided_backtracks);
	}

public:
    DinitzLinkCut(DynamicGraph& _g,Capacity & cap,int source=-1,int sink=-1):g(_g),capacity(cap),source(source),sink(sink),INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
		,ek(_g,cap,source,sink)
#endif
    {
    	  curflow=0;
      	last_modification=-1;
      	last_deletion=-1;
      	last_addition=-1;

      	history_qhead=-1;
      	last_history_clear=-1;

    }
    int getSource() const{
    	return source;
    }
    int getSink() const{
    	return sink;
    }
    void setCapacity(int u, int w, int c){
    	//C.resize(g.edges());
    	//C[ ]=c;

    }
    void setAllEdgeCapacities(int c){

    }
    void dbg_print_graph(int from, int to){

   #ifndef NDEBUG
    	return;

       			printf("digraph{\n");
       			for(int i = 0;i<g.nodes();i++){
       				if(i==from){
       					printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
       				}else if (i==to){
       					printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
       				}else if (dist[i]==-1){
       					printf("n%d [style=filled, fillcolor=gray]\n", i);
       				}else
       					printf("n%d\n", i);
       			}

       			for(int i = 0;i<g.edges();i++){
       				if(g.edgeEnabled(i)){
   						auto & e = g.all_edges[i];
   						const char * s = "black";

   						bool link=false;
   						bool backward=false;
   						for(auto & e:toLink){
   							if(e.edgeID==i){
   								link=true;
   								backward = e.backward;
   								if(backward){
   									s="orange";
   								}else{
   									s="green";
   								}
   								break;
   							}
   						}
   						int hasParent=false;
   						bool backwardParent=false;
   						if(parentEdge[e.from].edgeID==i){
   							hasParent=true;
   							s="red";
   						}else if (parentEdge[e.to].edgeID==i){
   							s="blue";
   							//hasParent=true;
   							backwardParent=true;
   						}

   						//assert(!(hasParent && link));

   						if(dist[e.to]==dist[e.from]+1){
   							s="blue";
   						}

   						std::cout<<"n" << e.from <<" -> n" << e.to << " [label=\"" << i <<": " <<  F[i]<<"/" << capacity[i]  << "\" color=\"" << s<<"\"]\n";
       				}
       			}

       			printf("}\n");
   #endif

       		}
    bool buildLevelGraph(int src, int dst) {
    	double start_time = rtime(3);
    	dist.clear();
    	dist.resize(g.nodes(),-1);
        dist[src] = 0;
        Q.push_back(src);
        //Build the level graph using a simple BFS
        for (int i = 0; i < Q.size(); i++) {
            int u = Q[i];
            for (int j = 0;j<g.nIncident(u);j++){
            	int edgeID = g.incident(u,j).id;
               	if(!g.edgeEnabled(edgeID))
   						continue;
               	int v =  g.incident(u,j).node;
                if (dist[v] < 0 && F[edgeID] < capacity[edgeID]) {
                    dist[v] = dist[u] + 1;
                    Q.push_back(v);
                }
            }
            for (int j = 0;j<g.nIncoming(u);j++){
            	int edgeID = g.incoming(u,j).id;
               	if(!g.edgeEnabled(edgeID))
   						continue;
               	int v =  g.incoming(u,j).node;
               	//this is a backward edge, so it has capacity exactly if the forward edge has flow
                if (dist[v] < 0 && F[edgeID]) {
                    dist[v] = dist[u] + 1;
                    Q.push_back(v);
                }
            }
        }
        Q.clear();
        bfstime+=rtime(3)-start_time;
        return dist[dst] >= 0;
    }
    bool dbg_hasLink(int u){
#ifndef NDEBUG
    	for(auto e:toLink){
    		if(e.u==u)
    			return true;
    	}
#endif
    	return false;
    }
    bool dbg_isLinkRoot(int v){
#ifndef NDEBUG
    	int u =  src;
    	int i = 0;
    	while(true){
    		u = forest.findRoot(u);
    		if(i<toLink.size()){
    			auto & e = toLink[i++];
    			assert(e.u==u);
    			u = e.v;
    		}else{
    			break;
    		}
    	}
    	/*if(u!=v){
    		dbg_print_graph(src,dst);
    	}*/
    	assert(u==v);
    	return u==v;
#endif
    	return true;
    }

    int findAugmentingPath(int src, int dst){
    	dbg_print_graph(src,dst);
    	//static int outerit = 0;
    	//static int innerit = 0;
    	toLink.clear();
        	int f = 0;
        	int u =forest.findRoot(src);
        	while(true){
        		/*if(++outerit==10){
        			int a=1;
        		}*/
    			bool found=true;
    			bool foundPath = false;
    			//double starttime=rtime(4);
    			dbg_print_graph(src,dst);

    			 while(found){
    				 //++innerit;
    					 found=false;
    					 assert(dbg_isLinkRoot(u));
    					 //assert(parentEdge[u].edgeID==-1);
    					if (u==dst){
    						foundPath=true;
    						break;
    					}else{
    						//extend the path
    						 for (;pos[u]<g.nIncident(u);pos[u]++){
    							//auto & edge = g.adjacency[u][pos[u]];
    							auto & edge = g.incident(u,pos[u]);
    							int edgeID = edge.id;
    							int v =  edge.node;
    							if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
    									continue;
    							if (F[edgeID] < capacity[edgeID]) {
    								// assert(parentEdge[u].edgeID==-1);
    								 assert(!dbg_hasLink(u));
    								toLink.push_back({u,v,false,edgeID});
    								u = forest.findRoot(v);
    								found=true;
    								//pos[u]++;
    								break;
    							}
    						}
    						 if(!found){
								for (;pos[u]-g.nIncident(u) <g.nIncoming(u);pos[u]++){
									//auto & edge = g.inverted_adjacency[u][pos[u]-g.nIncident(u)];
									auto & edge = g.incoming(u,pos[u]-g.nIncident(u));
									int edgeID = edge.id;
									int v = edge.node;
									if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
											continue;

									//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
									if (F[edgeID]) {
										// assert(parentEdge[u].edgeID==-1);
										 assert(!dbg_hasLink(u));
										toLink.push_back({u,v,true,edgeID});
										u = forest.findRoot(v);
										found=true;
										//pos[u]++;
										break;
									}
    							}
    						}
    					}

    			}
    			 dbg_isLinkRoot(u);
    			 dbg_print_graph(src,dst);

    			// double ctime = rtime(4);
    			// augtime_search+=ctime-starttime;
    			assert(!found);
    			 if(foundPath){
    				 stats_augmenting_rounds++;
    				 dbg_print_graph(src,dst);
    				 //ok, link the path together now that we know there is a path.
    				 for(int i = 0;i<toLink.size();i++){
    					 auto & e = toLink[i];
    					 int u = e.u;
    					 int v = e.v;
    					 int edgeID = e.edgeID;
    					 bool backward = e.backward;
    					 //assert(parentEdge[u].edgeID==-1);
    					 assert(edgeID>=0);
    					 assert(forest.findRoot(src)==u);
    					 assert(capacity[edgeID]>=F[edgeID]);
    					 assert(F[edgeID]>=0);
    					 if(backward){
    						 forest.link(u,v,F[edgeID]);
    						// parent_edge_backward[u]=true;
							 parentEdge[u].edgeID=edgeID;
							 parentEdge[u].backward=true;
    					 }else{
    						 forest.link(u,v,capacity[edgeID]-F[edgeID]);
    						 //parent_edge_backward[u]=false;
    						 parentEdge[u].backward=false;
							 parentEdge[u].edgeID=edgeID;
    					 }
    					 assert(forest.findRoot(src)==forest.findRoot(v));
    				 }
    				 static int iter = 0;
    				 ++iter;
    				 toLink.clear();
    				 dbg_print_graph(src,dst);
    				 assert(forest.findRoot(src)==dst);

    				 //found s-t augmenting path.
    				 //Find the minimum capacity of any edge from src to dst on that path
    				 int c = forest.minCost(src);
    				 //subtract that capacity from all edges on the path
    				 f+=c;
    				 forest.updateCostOfPathToRoot(src,-c);
    				 //delete edges with no remaining capacity
    				 int minC=0;
    				 while(minC==0){
    					 int u = forest.ancecstorFindMin(src);
    					 assert(forest.getCost(u)==minC);
    					 int edgeID = parentEdge[u].edgeID;
    					 assert(edgeID>=0);
    					 if(!parentEdge[u].backward){
    						 assert(F[edgeID] + c<=capacity[edgeID]);
    						 F[edgeID] +=c;
    					 }else{
    						 assert(c<=F[edgeID]);
    						 F[edgeID] -=c;
    					 }
    					 forest.cut(u);
    					 parentEdge[u].edgeID=-1;
    					 minC = forest.minCost(src);
    				 }
    				 toLink.clear();
    				 u= forest.findRoot(src);
    				 dbg_print_graph(src,dst);

    			 }else{
    				 //couldn't find any s-t augmenting path.
    				 //remove the vertex from the graph
    				 if(u==src){
    					 //done - no paths remain.
    					 //Clean up the state of the tree:
    					 for(int u = 0;u<g.nodes();u++){
    						 if(parentEdge[u].edgeID>=0){
    							 assert(!forest.isRoot(u));
    							 int c = forest.getCost(u);

    							 int edgeID = parentEdge[u].edgeID;assert(c<=capacity[edgeID]);
    							 if(!parentEdge[u].backward){
    								 F[edgeID]=capacity[edgeID]-c;
    							 }else{
    								 F[edgeID]=c;
    							 }

    							 forest.cut(u);
    							 parentEdge[u].edgeID=-1;
    						 }
    						 forest.undeleteNode(u);
    					 }

    					 toLink.clear();

    					 break;
    				 }else{
    					 dist[u]=-1;//prevents u from being explored again.
    					 if(toLink.size() && u==toLink.back().v){
    						 stats_avoided_backtracks++;
    						 u = toLink.back().u;
    						 toLink.pop_back();
    					 }else{
    						stats_backtracks++;
    						forest.removeNode(u);
							for(int i = 0;i<g.nIncident(u);i++){
								auto & edge = g.incident(u,i,false);
								int edgeID = edge.id;
								int v = edge.node;
								if(!g.edgeEnabled(edgeID) )
										continue;
								if(parentEdge[v].edgeID==edgeID){
									//need to remember the remaining flow on this edge...
									assert(F[edgeID]>0);//else this edge wouldn't be in the tree
									int residual_capacity = forest.getCost(v);
									F[edgeID]=residual_capacity;
									assert(F[edgeID]>=0);
									forest.cut(v);//this is a backward edge into u
									parentEdge[v].edgeID=-1;
								}
							}
							for(int i = 0;i<g.nIncoming(u);i++){
								auto & edge = g.incoming(u,i,false);
								int edgeID = edge.id;
								int v = edge.node;
								if(!g.edgeEnabled(edgeID)  )
										continue;
								if(parentEdge[v].edgeID==edgeID){
									int residual_capacity = forest.getCost(v);
									F[edgeID]=capacity[edgeID] - residual_capacity;
									assert(F[edgeID]>=0);
									forest.cut(v);
									parentEdge[v].edgeID=-1;
								}
							}
							if(toLink.size()){
								u = forest.findRoot(toLink.back().v);//should this be u or v?
							}else{
								u = forest.findRoot(src);
							}
    					 }
    				 }
    			 }
    			// augtime_cleanup+=rtime(4)- ctime;
        	}

    		return f;
        }

    int findAugmentingPath_pointers(int src, int dst){
    	int f = 0;
    	while(true){

			bool found=true;
			bool foundPath = false;
			double starttime=rtime(4);
			//dbg_print_graph(src,dst);
			int u=forest.findRoot(src);
			 while(found){
					 found=false;
					 //u = forest.findRoot(src);

					if (u==dst){
						foundPath=true;
						break;
					}else{
						//extend the path

						 for (;pos[u]<g.nIncident(u);pos[u]++){
							//auto & edge = g.adjacency[u][pos[u]];
							auto & edge = g.incident(u,pos[u]);
							int edgeID = edge.id;
							int v =  edge.node;
							if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
									continue;

							if (F[edgeID] < capacity[edgeID]) {
								forest.link(u,v,capacity[edgeID]-F[edgeID]);
								//tree_edges[edgeID].in_tree=true;
								parentEdge[u].backward=false;
								parentEdge[u].edgeID=edgeID;
								assert(forest.findRoot(u)==v);
								u = forest.findRoot(v);
								found=true;
								break;
							}
						}
						 if(!found){
								for (;pos[u]-g.nIncident(u) <g.nIncoming(u);pos[u]++){
									//auto & edge = g.inverted_adjacency[u][pos[u]-g.nIncident(u)];
									auto & edge = g.incoming(u,pos[u]-g.nIncident(u));
									int edgeID = edge.id;
									int v = edge.node;
									if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
											continue;

									//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
									if (F[edgeID]) {
										forest.link(u,v,F[edgeID]);
										//tree_edges[edgeID].in_tree_backward=true;
										parentEdge[u].backward=true;
										parentEdge[u].edgeID=edgeID;
										assert(forest.findRoot(u)==v);
										u = forest.findRoot(v);
										found=true;
										break;
									}
							}
						}
					}

			}
			//forest.dbg_print_forest(true);
			 double ctime = rtime(4);
			 augtime_search+=ctime-starttime;
			assert(!found);
			 if(foundPath){
				 stats_augmenting_rounds++;
				 //found s-t augmenting path.
				 //Find the minimum capacity of any edge from src to dst on that path
				 int c = forest.minCost(src);
				 //subtract that capacity from all edges on the path
				 f+=c;
				 forest.updateCostOfPathToRoot(src,-c);
				 //delete edges with no remaining capacity
				 //forest.dbg_print_forest(true);
				 int minC=0;
				 while(minC==0){

					 int u = forest.ancecstorFindMin(src);
					//forest.dbg_print_forest(true);
					 assert(forest.getCost(u)==minC);
					 int edgeID = parentEdge[u].edgeID;
					 if(!parentEdge[u].backward){
						 assert(F[edgeID] + c<=capacity[edgeID]);
						 F[edgeID] +=c;
					 }else{

						 assert(c<=F[edgeID]);
						 F[edgeID] -=c;
					 }
					 forest.cut(u);
					parentEdge[u].edgeID=-1;
					 minC = forest.minCost(src);
					//forest.dbg_print_forest(true);
				 }

			 }else{
				 //couldn't find any s-t augmenting path.
				 //remove the vertex from the graph
				 if(u==src){
					 //done - no paths remain.
					 //Clean up the state of the tree:
					 for(int u = 0;u<g.nodes();u++){
						 if(parentEdge[u].edgeID>=0){
							 int c = forest.getCost(u);
							 int edgeID = parentEdge[u].edgeID;
							 if(!parentEdge[u].backward){
								 F[edgeID]=capacity[edgeID]-c;
							 }else{

								 F[edgeID]=c;
							 }

							 forest.cut(u);
							 parentEdge[u].edgeID=-1;
						 }
					 }
				/*	 for(int i = 0;i<disabled.size();i++){
						 disabled[i]=false;
					 }
*/
					 break;
				 }else{
					// disabled[u]=true;
					 dist[u]=-1;
					 for(int i = 0;i<g.nIncident(u);i++){
						auto & edge = g.incident(u,i);
						int edgeID = edge.id;
						int v = edge.node;
						if(!g.edgeEnabled(edgeID) )
								continue;
						if(parentEdge[v].edgeID==edgeID){
							//need to remember the remaining flow on this edge...
							assert(F[edgeID]>0);//else this edge wouldn't be in the tree
							int residual_capacity = forest.getCost(v);
							F[edgeID]=residual_capacity;
							assert(F[edgeID]>=0);
							forest.cut(v);//this is a backward edge into u
							parentEdge[v].edgeID=-1;
						}
					}
					for(int i = 0;i<g.nIncoming(u);i++){
						auto & edge = g.incoming(u,i);
						int edgeID = edge.id;
						int v = edge.node;
						if(!g.edgeEnabled(edgeID)  )
								continue;
						if(parentEdge[v].edgeID==edgeID){
							int residual_capacity = forest.getCost(v);
							F[edgeID]=capacity[edgeID] - residual_capacity;
							assert(F[edgeID]>=0);
							forest.cut(v);
							parentEdge[v].edgeID=-1;
						}
					}
				 }
			 }
			 augtime_cleanup+=rtime(4)- ctime;
    	}

		return f;
    }

    int findAugmentingPath_original(int src, int dst){

        	int f = 0;
        	while(true){
        		stats_augmenting_rounds++;
    			bool found=true;
    			bool foundPath = false;

    			//dbg_print_graph(src,dst);
    			int u=forest.findRoot(src);
    			 while(found){
    					 found=false;

    					if (u==dst){
    						foundPath=true;
    						break;
    					}else{
    						//extend the path

    						for(int i = 0;i<g.nIncident(u);i++){
								auto & edge = g.incident(u,i);
    							int edgeID = edge.id;
    							int v =  edge.node;
    							if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
    									continue;

    							if (F[edgeID] < capacity[edgeID]) {
    								forest.link(u,v,capacity[edgeID]-F[edgeID]);
    								//tree_edges[edgeID].in_tree=true;
    								parentEdge[u].backward=false;
    								parentEdge[u].edgeID=edgeID;
    								assert(forest.findRoot(u)==v);
    								found=true;
    								break;
    							}
    						}
    						if(!found){
    							for(int i = 0;i<g.nIncoming(u);i++){
									auto & edge = g.incoming(u,i);
    								int edgeID = edge.id;
    								int v = edge.node;
    								if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
    										continue;

    								//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
    								if (F[edgeID]) {
    									forest.link(u,v,F[edgeID]);
    									parentEdge[u].backward=true;
    									//tree_edges[edgeID].in_tree_backward=true;
    									parentEdge[u].edgeID=edgeID;
    									assert(forest.findRoot(u)==v);
    									found=true;
    									break;
    								}
    							}
    						}
    					}
    					u = forest.findRoot(src);
    			}
    			//forest.dbg_print_forest(true);

    			assert(!found);
    			 if(foundPath){
    				 //found s-t augmenting path.
    				 //Find the minimum capacity of any edge from src to dst on that path
    				 int c = forest.minCost(src);
    				 //subtract that capacity from all edges on the path
    				 f+=c;
    				 forest.updateCostOfPathToRoot(src,-c);
    				 //delete edges with no remaining capacity
    				 //forest.dbg_print_forest(true);
    				 int minC=0;
    				 while(minC==0){

    					 int u = forest.ancecstorFindMin(src);
    					//forest.dbg_print_forest(true);
    					 assert(forest.getCost(u)==minC);
    					 int edgeID = parentEdge[u].edgeID;
    					 if(!parentEdge[u].backward){
    					 //if(tree_edges[edgeID].in_tree){
    						 assert(F[edgeID] + c<=capacity[edgeID]);
    						 F[edgeID] +=c;
    					 }else{
    						 //assert(tree_edges[edgeID].in_tree_backward);
    						 assert(c<=F[edgeID]);
    						 F[edgeID] -=c;
    					 }
    					 forest.cut(u);
    					parentEdge[u].edgeID=-1;
    					 minC = forest.minCost(src);
    					//forest.dbg_print_forest(true);
    				 }

    			 }else{
    				 //couldn't find any s-t augmenting path.
    				 //remove the vertex from the graph
    				 if(u==src){
    					 //done - no paths remain.
    					 //Clean up the state of the tree:
    					 for(int u = 0;u<g.nodes();u++){
    						 if(parentEdge[u].edgeID>=0){
    							 int c = forest.getCost(u);
    							 int edgeID = parentEdge[u].edgeID;
    							 if(!parentEdge[u].edgeID){
    							// if(tree_edges[edgeID].in_tree){
    								 F[edgeID]=capacity[edgeID]-c;
    							 }else{
    								 //assert(tree_edges[edgeID].in_tree_backward);
    								 F[edgeID]=c;
    							 }

    							 forest.cut(u);
    							 parentEdge[u].edgeID=-1;
    						 }
    					 }
    				/*	 for(int i = 0;i<disabled.size();i++){
    						 disabled[i]=false;
    					 }*/
    		/*			 for(int i = 0;i<tree_edges.size();i++){
    						 tree_edges[i] = EdgeInTree();
    					 }*/
    					 break;
    				 }else{
    					 dist[u]=-1;
    					 //disabled[u]=true;
    					//forest.cut(u);
    					//for(auto edge:g.adjacency[u]){
    					 for(int i = 0;i<g.nIncident(u);i++){
    						auto & edge = g.incident(u,i);
    						int edgeID = edge.id;
    						int v = edge.node;
    						if(!g.edgeEnabled(edgeID))
    								continue;
    						if(parentEdge[v].edgeID==edgeID){
    							//need to remember the remaining flow on this edge...
    							assert(F[edgeID]>0);//else this edge wouldn't be in the tree
    							int residual_capacity = forest.getCost(v);
    							F[edgeID]=residual_capacity;
    							assert(F[edgeID]>=0);
    							forest.cut(v);//this is a backward edge into u
    							parentEdge[v].edgeID=-1;
    						}
    					}
    					 for(int i = 0;i<g.nIncoming(u);i++){
							auto & edge = g.incoming(u,i);
    						int edgeID = edge.id;
    						int v = edge.node;
    						if(!g.edgeEnabled(edgeID)  )
    								continue;
    						if(parentEdge[v].edgeID==edgeID){
    							int residual_capacity = forest.getCost(v);
    							F[edgeID]=capacity[edgeID] - residual_capacity;
    							assert(F[edgeID]>=0);
    							forest.cut(v);
    							parentEdge[v].edgeID=-1;
    						}
    					}
    				 }
    			 }

        	}

    		return f;
        }
	long num_updates=0;
	int numUpdates()const{
		return num_updates;
	}
   const int update(){
		return maxFlow(source,sink);
	}
   std::vector<int> changed_edges;
   std::vector<int> &  getChangedEdges(){
    	return changed_edges;
    }
    void clearChangedEdges(){

    }
    void setSource(int s){
      	if(source==s){
      		return;
      	}
      	source=s;
      	last_modification= g.modifications-1;
      }
      void setSink(int t){
      	if(sink==t){
      		return;
      	}
      	sink=t;
      	last_modification=g.modifications-1;
      }
    const  int maxFlow(int s, int t){
    	int f = 0;
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"f %d %d\n", s,t);
			fflush(g.outfile);
		}
#endif

      	if(last_modification>0 && g.modifications==last_modification){
			return curflow;
		}
      	src=s;
      	dst=t;
    	F.clear();
    	F.resize(g.all_edges.size());
    	dist.clear();
    	dist.resize(g.nodes());
    	f=0;
    	forest.reset();
    	while(forest.nNodes()<g.nodes())
    		forest.addNode();
      	pos.clear();pos.resize(g.nodes());

		parentEdge.clear();
		parentEdge.resize(g.nodes(),{false,-1});

	    while (buildLevelGraph(s,t)) {
	    	stats_rounds++;
	    	double start_time = rtime(3);
	    	dbg_print_graph(s,t);
	    	for(int i = 0;i<pos.size();i++)
	    		pos[i]=0;

	    	for(int i = 0;i<parentEdge.size();i++){
	    		parentEdge[i]={false,-1};
	    	}

			if(opt_dinics_recursive){
				int delta = findAugmentingPath_original(src,dst);
				f += delta;
			}else{
				int delta = findAugmentingPath(src,dst);
				f += delta;
			}
			dbg_print_graph(s,t);
			augtime+=rtime(3)- start_time;
	    }

#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
    	assert(f==expected_flow);
#endif

    	curflow=f;
    	num_updates++;
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
        return f;
    }


    std::vector<bool> seen;
    std::vector<bool> visited;
    const int minCut( std::vector<MaxFlowEdge> & cut){
    	return minCut(source,sink,cut);
    }
    const  int minCut(int s, int t, std::vector<dgl::MaxFlowEdge> & cut){
    	const int f = maxFlow(s,t);
    	//ok, now find the cut
    	Q.clear();
		Q.push_back(s);
		seen.clear();
		seen.resize(g.nodes());
		seen[s]=true;

		//explore the residual graph
		for(int j = 0;j<Q.size();j++){
		   int u = Q[j];

			for(int i = 0;i<g.nIncident(u);i++){
				if(!g.edgeEnabled(g.incident(u,i).id))
					continue;
				int v = g.incident(u,i).node;
				int id = g.incident(u,i).id;
				if(capacity[id] - F[id] == 0){
					cut.push_back(MaxFlowEdge{u,v,id});//potential element of the cut
				}else if(!seen[v]){
					Q.push_back(v);
					seen[v]=true;
				}
			}
			for(int i = 0;i<g.nIncoming(u);i++){
				if(!g.edgeEnabled(g.incoming(u,i).id))
					continue;
				int v = g.incoming(u,i).node;
				int id = g.incoming(u,i).id;
				if(F[id] == 0){

				}else if(!seen[v]){
					Q.push_back(v);
					seen[v]=true;
				}
			}
		}
		//Now keep only the edges from a seen vertex to an unseen vertex
		int i, j = 0;
		for( i = 0;i<cut.size();i++){
			if(!seen[cut[i].v] && seen[cut[i].u]){
				cut[j++]=cut[i];
			}
		}
		cut.resize(j);
#ifndef NDEBUG
		int dbg_sum = 0;
		for(int i = 0;i<cut.size();i++){
			int id = cut[i].id;
			assert(F[id]==capacity[id]);
			dbg_sum+=F[id];
		}
		assert(dbg_sum==f);
#endif
    	return f;
    }
    const  int getEdgeCapacity(int id){
     	assert(g.edgeEnabled(id));
     	return capacity[id];
     }
    const  int getEdgeFlow(int id){
    	assert(g.edgeEnabled(id));
    	return F[id];// reserve(id);
    }
    const  int getEdgeResidualCapacity(int id){
    	assert(g.edgeEnabled(id));
    	return  capacity[id]-F[id];// reserve(id);
    }
};
};
#endif

