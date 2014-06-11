#ifndef DINICS_LINKCUT_H
#define DINICS_LINKCUT_H


#include "MaxFlow.h"
#include "mtl/Vec.h"
#include "core/Config.h"
#include "graph/EdmondsKarpAdj.h"
#include "graph/LinkCutCost.h"
#include <algorithm>
#include <climits>
using namespace std;
using namespace Minisat;

template< class Capacity ,class EdgeStatus=vec<bool> >
class DinitzLinkCut:public MaxFlow{

public:

	vec<int> F;

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

    int history_qhead;
    int last_history_clear;
    //vec<LocalEdge> prev;
    //vec<int> M;
    vec<int> dist;
    vec<int> pos;//position in the combined forward and backward adjacency list of each node in the DFS.
    DynamicGraph<EdgeStatus>& g;
    Capacity & capacity;
    int INF;
    int src;
    int dst;
    struct ParentEdge{
    	bool backward;
    	int edgeID;
    };
    vec<ParentEdge> parentEdge;

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
    vec<EdgeInTree> tree_edges;//for each edge.*/
    LinkCutCost forest;
    vec<int> Q;
    struct Link{
    	int u;
    	int v;
    	bool backward:1;
    	int edgeID:31;
    };
    vec<Link> toLink;
#ifdef DEBUG_MAXFLOW
	EdmondsKarpAdj<Capacity, EdgeStatus> ek;
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
    DinitzLinkCut(DynamicGraph<EdgeStatus>& _g,Capacity & cap):g(_g),capacity(cap),INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
		,ek(_g,cap)
#endif
    {
    	  curflow=0;
      	last_modification=-1;
      	last_deletion=-1;
      	last_addition=-1;

      	history_qhead=-1;
      	last_history_clear=-1;

    }
    void setCapacity(int u, int w, int c){
    	//C.growTo(g.edges);
    	//C[ ]=c;

    }
    void setAllEdgeCapacities(int c){

    }
    void dbg_print_graph(int from, int to){
   #ifndef NDEBUG


       			printf("digraph{\n");
       			for(int i = 0;i<g.nodes;i++){
       				if(i==from){
       					printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
       				}else if (i==to){
       					printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
       				}else if (dist[i]==-1){
       					printf("n%d [style=filled, fillcolor=gray]\n", i);
       				}else
       					printf("n%d\n", i);
       			}

       			for(int i = 0;i<g.edges;i++){
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

   						assert(!(hasParent && link));

   						/*if(dist[e.to]==dist[e.from]+1){
   							s="blue";
   						}*/
   						/*if(value(e.v)==l_True)
   							s="blue";
   						else if (value(e.v)==l_False)
   							s="red";*/
   						printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"%s\"]\n", e.from,e.to, i, F[i],g.weights[i] , s);
       				}
       			}

       			printf("}\n");
   #endif
       		}
    bool buildLevelGraph(int src, int dst) {
    	double start_time = rtime(3);
    	dist.clear();
    	dist.growTo(g.nodes,-1);
        dist[src] = 0;
        Q.push(src);
        //Build the level graph using a simple BFS
        for (int i = 0; i < Q.size(); i++) {
            int u = Q[i];
            for (int j = 0;j<g.adjacency[u].size();j++){
            	int edgeID = g.adjacency[u][j].id;
               	if(!g.edgeEnabled(edgeID))
   						continue;
               	int v =  g.adjacency[u][j].node;
                if (dist[v] < 0 && F[edgeID] < capacity[edgeID]) {
                    dist[v] = dist[u] + 1;
                    Q.push(v);
                }
            }
            for (int j = 0;j<g.inverted_adjacency[u].size();j++){
            	int edgeID = g.inverted_adjacency[u][j].id;
               	if(!g.edgeEnabled(edgeID))
   						continue;
               	int v =  g.inverted_adjacency[u][j].node;
               	//this is a backward edge, so it has capacity exactly if the forward edge has flow
                if (dist[v] < 0 && F[edgeID]) {
                    dist[v] = dist[u] + 1;
                    Q.push(v);
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
    	if(u!=v){
    		dbg_print_graph(src,dst);
    	}
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
    						 for (;pos[u]<g.adjacency[u].size();pos[u]++){
    							auto & edge = g.adjacency[u][pos[u]];
    							int edgeID = edge.id;
    							int v =  edge.node;
    							if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
    									continue;
    							if (F[edgeID] < capacity[edgeID]) {
    								// assert(parentEdge[u].edgeID==-1);
    								 assert(!dbg_hasLink(u));
    								toLink.push({u,v,false,edgeID});
    								u = forest.findRoot(v);
    								found=true;
    								//pos[u]++;
    								break;
    							}
    						}
    						 if(!found){
								for (;pos[u]-g.adjacency[u].size() <g.inverted_adjacency[u].size();pos[u]++){
									auto & edge = g.inverted_adjacency[u][pos[u]-g.adjacency[u].size()];
									int edgeID = edge.id;
									int v = edge.node;
									if((dist[v] != dist[u] + 1) ||  !g.edgeEnabled(edgeID) )
											continue;

									//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
									if (F[edgeID]) {
										// assert(parentEdge[u].edgeID==-1);
										 assert(!dbg_hasLink(u));
										toLink.push({u,v,true,edgeID});
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
    						 F[edgeID] =c;
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
    					 for(int u = 0;u<g.nodes;u++){
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
    						 forest.undeleteNode(u);
    					 }

    					 toLink.clear();

    					 break;
    				 }else{
    					 dist[u]=-1;//prevents u from being explored again.
    					 if(toLink.size() && u==toLink.last().v){
    						 stats_avoided_backtracks++;
    						 u = toLink.last().u;
    						 toLink.pop();
    					 }else{
    						stats_backtracks++;
    						forest.removeNode(u);
							/*for(auto edge:g.adjacency[u]){
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
							for(auto edge:g.inverted_adjacency[u]){
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
							}*/
							if(toLink.size()){
								u = forest.findRoot(toLink.last().u);
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

						 for (;pos[u]<g.adjacency[u].size();pos[u]++){
							auto & edge = g.adjacency[u][pos[u]];
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
								for (;pos[u]-g.adjacency[u].size() <g.inverted_adjacency[u].size();pos[u]++){
									auto & edge = g.inverted_adjacency[u][pos[u]-g.adjacency[u].size()];
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
						 F[edgeID] =c;
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
					 for(int u = 0;u<g.nodes;u++){
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
					for(auto edge:g.adjacency[u]){
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
					for(auto edge:g.inverted_adjacency[u]){
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

    						for(auto edge:g.adjacency[u]){
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
    							for(auto edge:g.inverted_adjacency[u]){
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
    					 if(!parentEdge[u].edgeID){
    					 //if(tree_edges[edgeID].in_tree){
    						 assert(F[edgeID] + c<=capacity[edgeID]);
    						 F[edgeID] +=c;
    					 }else{
    						 //assert(tree_edges[edgeID].in_tree_backward);
    						 assert(c<=F[edgeID]);
    						 F[edgeID] =c;
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
    					 for(int u = 0;u<g.nodes;u++){
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
    					for(auto edge:g.adjacency[u]){
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
    					for(auto edge:g.inverted_adjacency[u]){
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

    int maxFlow(int s, int t){
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
    	F.growTo(g.all_edges.size());
    	dist.clear();
    	dist.growTo(g.nodes);
    	f=0;
    	forest.reset();
    	while(forest.nNodes()<g.nodes)
    		forest.addNode();
      	pos.clear();pos.growTo(g.nodes);

		parentEdge.clear();
		parentEdge.growTo(g.nodes,{false,-1});

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
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
        return f;
    }


    vec<bool> seen;
    vec<bool> visited;

    int minCut(int s, int t, vec<Edge> & cut){
    	int f = maxFlow(s,t);
    	//ok, now find the cut
    	Q.clear();
    	Q.push(s);
    	seen.clear();
    	seen.growTo(g.nodes);
    	seen[s]=true;

    	for(int j = 0;j<Q.size();j++){
		   int u = Q[j];

    		for(int i = 0;i<g.adjacency[u].size();i++){
    			if(!g.edgeEnabled(g.adjacency[u][i].id))
    				continue;
    			int v = g.adjacency[u][i].node;
    			int id = g.adjacency[u][i].id;
    			if(capacity[id] - F[id] == 0){
    				cut.push(Edge{u,v,id});
    			}else if(!seen[v]){
    				Q.push(v);
    				seen[v]=true;
    			}
    		}
    	}
    	//Now remove any edges that lead to vertices that we ended up visiting
    	int i, j = 0;
    	for( i = 0;i<cut.size();i++){
    		if(!seen[cut[i].v]){
    			cut[j++]=cut[i];
    		}
    	}
    	cut.shrink(i-j);
    	return f;
    }
    int getEdgeCapacity(int id){
     	assert(g.edgeEnabled(id));
     	return capacity[id];
     }
    int getEdgeFlow(int id){
    	assert(g.edgeEnabled(id));
    	return F[id];// capacity(id);
    }
    int getEdgeResidualCapacity(int id){
    	assert(g.edgeEnabled(id));
    	return  capacity[id]-F[id];// capacity(id);
    }
};
#endif

