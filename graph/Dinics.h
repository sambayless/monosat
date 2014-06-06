#ifndef DINICS_H
#define DINICS_H


#include "MaxFlow.h"
#include "mtl/Vec.h"
#include "core/Config.h"
#include "graph/EdmondsKarpAdj.h"
#include <algorithm>
#include <climits>
using namespace std;
using namespace Minisat;

template< class Capacity ,class EdgeStatus=vec<bool> >
class Dinitz:public MaxFlow{

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
    vec<LocalEdge> prev;
    vec<int> M;
    vec<int> dist;
    vec<int> pos;//position in the combined forward and backward adjacency list of each node in the DFS.
    DynamicGraph<EdgeStatus>& g;
    Capacity & capacity;
    int INF;
    int src;
    int dst;
#ifdef DEBUG_MAXFLOW
	EdmondsKarpAdj<Capacity, EdgeStatus> ek;
#endif
    vec<int> Q;

public:
    Dinitz(DynamicGraph<EdgeStatus>& _g,Capacity & cap):g(_g),capacity(cap),INF(0xF0F0F0)
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
    	//setAllEdgeCapacities(1);
    }
    void setCapacity(int u, int w, int c){
    	//C.growTo(g.edges);
    	//C[ ]=c;

    }
    void setAllEdgeCapacities(int c){

    }
    void dbg_print_graph(int from, int to){
   #ifndef NDEBUG
    	return;
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
   						if(dist[e.to]==dist[e.from]+1){
   							s="blue";
   						}
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
    }

    int findAugmentingPath(int u) {
    	int m =0;
    	assert(Q.size()==0);

        Q.push(src);
        M[src]=INT_MAX;
        while(Q.size()){
        	int u = Q.last();
		    if (u == dst)
		    	return M[u];
		    bool found=false;
			for (;pos[u]<g.adjacency[u].size();pos[u]++){
				int edgeID = g.adjacency[u][pos[u]].id;
				if(!g.edgeEnabled(edgeID))
						continue;
				int v =  g.adjacency[u][pos[u]].node;
				if (dist[v] == dist[u] + 1 && F[edgeID] < capacity[edgeID]) {

					M[v] = min(M[u], capacity[edgeID] - F[edgeID]);
					prev[v]=LocalEdge(u,edgeID,false);
					if(v==dst){
						//M[v] = min(M[u], capacity[edgeID] - F[id]);
						m=M[dst];
						pos[u]++;
						break;
					}else{
						Q.push(v);
						found=true;
						pos[u]++;
						break;
					}
				}
			}
			if(!found){
				for (;pos[u]-g.adjacency[u].size() <g.inverted_adjacency[u].size();pos[u]++){
					int edgeID = g.inverted_adjacency[u][pos[u]-g.adjacency[u].size()].id;
					if(!g.edgeEnabled(edgeID))
							continue;
					int v =  g.inverted_adjacency[u][pos[u]-g.adjacency[u].size()].node;
					//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
					if (dist[v] == dist[u] + 1 && F[edgeID]) {
						M[v] = min(M[u], F[edgeID]);
						prev[v]=LocalEdge(u,edgeID,false);
						if(v==dst){
							m=M[dst];
							pos[u]++;
							break;
						}else{
							Q.push(v);
							found=true;
							pos[u]++;
							break;
						}
					}
				}
			}
			if(!found){
				Q.pop();
			}
        }
        Q.clear();
        if(m>0){
        	//we found an augmenting flow, so update all the edge flows correspondingly.
			int v = dst;
		   while (v!=  src){
			   int u = prev[v].from;
			   int id = prev[v].id;
			   if(prev[v].backward){
				   F[id] = F[id] - m;
			   }else
				   F[id] = F[id] + m;
			   v = u;
		   }
        }
        return m;
    }

    int findAugmentingPath_recursive(int u, int f) {
            if (u == dst)
                return f;

            for (;pos[u]<g.adjacency[u].size();pos[u]++){
    			int edgeID = g.adjacency[u][pos[u]].id;
    			if(!g.edgeEnabled(edgeID))
    					continue;
    			int v =  g.adjacency[u][pos[u]].node;
    			if (dist[v] == dist[u] + 1 && F[edgeID] < capacity[edgeID]) {
    				int df = findAugmentingPath_recursive(v, min(f, capacity[edgeID] - F[edgeID]));
    				if (df > 0) {
    					F[edgeID] += df;
    					return df;
    				}
    			}
            }
            for (;pos[u]-g.adjacency[u].size() <g.inverted_adjacency[u].size();pos[u]++){
    			int edgeID = g.inverted_adjacency[u][pos[u]-g.adjacency[u].size()].id;
    			if(!g.edgeEnabled(edgeID))
    					continue;
    			int v =  g.inverted_adjacency[u][pos[u]-g.adjacency[u].size()].node;
    			//these are backwards edges, which have capacity exactly if the forward edge has non-zero flow
    			if (dist[v] == dist[u] + 1 && F[edgeID]) {
    				int df = findAugmentingPath_recursive(v, min(f, F[edgeID]));
    				if (df > 0) {
    					F[edgeID] -= df;
    					return df;
    				}
    			}
            }
            return 0;
        }

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
    	dist.clear();
    	dist.growTo(g.nodes);
    	M.growTo(g.nodes);
    	prev.growTo(g.nodes);
    	f=0;

		while (buildLevelGraph(s,t)) {
			dbg_print_graph(s,t);
			pos.clear();pos.growTo(g.nodes);//does this have to be done inside the inner well loop?
			if(opt_dinics_recursive){
				while (int delta = findAugmentingPath_recursive(s, INT_MAX)){
					f += delta;
					dbg_print_graph(s,t);
				}
			}else{
				while (int delta = findAugmentingPath(s)){
					f += delta;
					dbg_print_graph(s,t);
				}
			}
		}


#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
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

