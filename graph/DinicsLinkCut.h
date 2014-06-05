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
class DinicsLinkCut:public MaxFlow{

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
    //vec<int> dist;
    //vec<int> pos;//position in the combined forward and backward adjacency list of each node in the DFS.
    DynamicGraph<EdgeStatus>& g;
    Capacity & capacity;
    int INF;
    int src;
    int dst;
    vec<int> parentEdge;
    vec<bool> disabled;
    struct EdgeInTree{
    	char in_tree:1;
    	char in_tree_backward:1;
    	char disabled:1;
    	char backward_disabled:1;
    	EdgeInTree():in_tree(0),in_tree_backward(0),disabled(0),backward_disabled(0){}
    };
    vec<EdgeInTree> tree_edges;//for each edge.
    LinkCutCost forest;

#ifdef DEBUG_MAXFLOW
	EdmondsKarpAdj<Capacity, EdgeStatus> ek;
#endif
    vec<int> Q;

public:
    DinicsLinkCut(DynamicGraph<EdgeStatus>& _g,Capacity & cap):g(_g),capacity(cap),INF(0xF0F0F0)
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

       		static int it = 0;
       		if(++it==6){
       			int a =1;
       		}
       		printf("Graph %d\n", it);
       			printf("digraph{\n");
       			for(int i = 0;i<g.nodes;i++){
       				if(i==from){
       					printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
       				}else if (i==to){
       					printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
       				}else
       					printf("n%d\n", i);
       			}

       			for(int i = 0;i<g.edges;i++){
       				if(g.edgeEnabled(i)){
   						auto & e = g.all_edges[i];
   						const char * s = "black";
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
   /* bool buildLevelGraph(int src, int dst) {
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
        return dist[dst] >= 0;
    }*/


    int maxFlow(int s, int t){
    	int f = 0;
#ifdef RECORD
		if(g.outfile && mincutalg==MinCutAlg::ALG_EDKARP_ADJ){
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
/*    	dist.clear();
    	dist.growTo(g.nodes);
    	M.growTo(g.nodes);
    	prev.growTo(g.nodes);*/
    	tree_edges.growTo(g.nodes);
    	f=0;
    	forest.reset();
    	while(forest.nNodes()<g.nodes)
    		forest.addNode();
    	disabled.clear();
    	disabled.growTo(g.nodes);

    	parentEdge.growTo(g.nodes,-1);
    	int u ;
    	while(true){
    		bool found=true;
    		bool foundPath = false;

    		dbg_print_graph(src,dst);
    		u=forest.findRoot(src);
   		 	 while(found){
   		 		 	 found=false;
   		 		 	 //u = forest.findRoot(src);
					if (u==dst){
						foundPath=true;
						break;
					}else{
						//extend the path

						for(auto edge:g.adjacency[u]){
							int edgeID = edge.id;
							int v =  edge.node;
							if(!g.edgeEnabled(edgeID) || disabled[v] || tree_edges[edgeID].disabled)
									continue;

							if (F[edgeID] < capacity[edgeID]) {
								forest.link(u,v,capacity[edgeID]-F[edgeID]);
								tree_edges[edgeID].in_tree=true;
								parentEdge[u]=edgeID;
								assert(forest.findRoot(u)==v);
								found=true;
								break;
							}
						}
						if(!found){
							for(auto edge:g.inverted_adjacency[u]){
								int edgeID = edge.id;
								int v = edge.node;
								if(!g.edgeEnabled(edgeID)|| disabled[v] || tree_edges[edgeID].backward_disabled)
										continue;

								//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
								if (F[edgeID]) {
									forest.link(u,v,F[edgeID]);
									tree_edges[edgeID].in_tree_backward=true;
									parentEdge[u]=edgeID;
									assert(forest.findRoot(u)==v);
									found=true;
									break;
								}
							}
						}
					}
					u = forest.findRoot(src);
			}
   		 	forest.dbg_print_forest(true);

   		 	assert(!found);
   		 	 if(foundPath){
   		 		 //found s-t augmenting path.
   		 		 //Find the minimum capacity of any edge from src to dst on that path
   		 		 int c = forest.minCost(src);
   		 		 //subtract that capacity from all edges on the path
   		 		 f+=c;
   		 		 forest.updateCostOfPathToRoot(src,-c);
   		 		 //delete edges with no remaining capacity
   		 		 forest.dbg_print_forest(true);
   		 		 int minC=0;
   		 		 while(minC==0){

   		 			 int u = forest.ancecstorFindMin(src);
   		 			forest.dbg_print_forest(true);
   		 			 assert(forest.getCost(u)==minC);
   		 			 int edgeID = parentEdge[u];
   		 			 assert(F[edgeID] + c<=capacity[edgeID]);
   		 			 F[edgeID] +=c;
   		 			 forest.cut(u);
   		 			parentEdge[u]=-1;
   		 			 minC = forest.minCost(src);
   		 			forest.dbg_print_forest(true);
   		 		 }

   		 	 }else{
   		 		 //couldn't find any s-t augmenting path.
   		 		 //remove the vertex from the graph
   		 		 if(u==src){
   		 			 //done - no paths remain.
   		 			 //Clean up the state of the tree:
   		 			 for(int u = 0;u<g.nodes;u++){
   		 				 if(parentEdge[u]>=0){
							 int c = forest.getCost(u);
							 int edgeID = parentEdge[u];
							 if(tree_edges[edgeID].in_tree){
								 F[edgeID]=capacity[edgeID]-c;
							 }else{
								 assert(tree_edges[edgeID].in_tree_backward);
								 F[edgeID]=c;
							 }

							 forest.cut(u);
							 parentEdge[u]=-1;
   		 				 }
   		 			 }
   		 			 for(int i = 0;i<disabled.size();i++){
   		 				 disabled[i]=false;
   		 			 }
   		 			 for(int i = 0;i<tree_edges.size();i++){
   		 				 tree_edges[i] = EdgeInTree();
   		 			 }
   		 			 break;
   		 		 }else{
   		 			 disabled[u]=true;
   		 			//forest.cut(u);
   		 			for(auto edge:g.adjacency[u]){
						int edgeID = edge.id;
						int v = edge.node;
						if(!g.edgeEnabled(edgeID)|| disabled[v]  || tree_edges[edgeID].backward_disabled)
								continue;
						if(tree_edges[edgeID].in_tree_backward){
							//need to remember the remaining flow on this edge...
							assert(F[edgeID]>0);//else this edge wouldn't be in the tree
							int residual_capacity = forest.getCost(v);
							F[edgeID]=residual_capacity;
							assert(F[edgeID]>=0);
							forest.cut(v);//this is a backward edge into u
							parentEdge[v]=-1;
						}
					}
   					for(auto edge:g.inverted_adjacency[u]){
						int edgeID = edge.id;
						int v = edge.node;
						if(!g.edgeEnabled(edgeID)|| disabled[v]  || tree_edges[edgeID].disabled)
								continue;
						if(tree_edges[edgeID].in_tree){
							int residual_capacity = forest.getCost(v);
							F[edgeID]=capacity[edgeID] - residual_capacity;
							assert(F[edgeID]>=0);
							forest.cut(v);
							parentEdge[v]=-1;
						}
					}
   		 		 }
   		 	 }


    	}


#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
    	assert(f==expected_flow);
#endif

        //dbg_print_graph(s,t);
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
    //	visited.clear();
    	//visited.growTo(g.nodes);
    //	visited[s]=true;
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

