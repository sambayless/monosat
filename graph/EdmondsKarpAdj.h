#ifndef EDMONDS_KARP_ADJ_H
#define EDMONDS_KARP_ADJ_H

// Maximum Flow- Edmonds-Karp Algorithm
// by Iraklis Karagkiozoglou <i.kar@windowslive.com>

#include "MaxFlow.h"
#include "mtl/Vec.h"
using namespace std;
using namespace Minisat;

template< class Capacity ,class EdgeStatus=vec<bool>>
class EdmondsKarpAdj:public MaxFlow{


/*
    algorithm EdmondsKarp
        input:
            C[1..n, 1..n] (Capacity matrix)
            E[1..n, 1..?] (Neighbour lists)
            s             (Source)
            t             (Sink)
        output:
            f             (Value of maximum flow)
            F             (A matrix giving a legal flow with the maximum value)*/
    //vec<vec<int> > F;//(Residual capacity from u to v is C[u,v] - F[u,v])
    //vec<vec<int> > C;

	vec<int> F;
	vec<int> rev;
	//vec<int> C;

	struct LocalEdge{
		int from;
		int id;
	};
    vec<LocalEdge> prev;
    vec<int> M;
    DynamicGraph<EdgeStatus>& g;
    Capacity capacity;
    int INF;
#ifdef DEBUG_MAXFLOW
    	EdmondsKarp<EdgeStatus> ek;
#endif
    /*
     *            input:
               C, E, s, t, F
           output:
               M[t]          (Capacity of path found)
               P             (Parent table)
     */
    vec<int> Q;




    int  BreadthFirstSearch(int s, int t){
     	for(int i = 0;i<g.nodes;i++)
     		prev[i].from=-1;
     	prev[s].from = -2;
    	Q.clear();
           Q.push(s);

       	for(int j = 0;j<Q.size();j++){
   			   int u = Q[j];

               for (int i = 0;i<g.adjacency[u].size();i++){
            	   if(!g.edgeEnabled(g.adjacency[u][i].id))
						continue;
            	   int id = g.adjacency[u][i].id;
            	   int v = g.adjacency[u][i].node;
                   ///(If there is available capacity, and v is not seen before in search)

            	   int f = F[id];
            	   int c = capacity(id);

            	 //  int fr = F[id];
                   if (((c - F[id]) > 0) && (prev[v].from == -1)){
                       prev[v] = {u,id};
                       M[v] = min(M[u], c - F[id]);
                       if (v != t)
                           Q.push(v);
                       else
                           return M[t];
                   }
               }
           }
           return 0;


	   }
public:
    EdmondsKarpAdj(DynamicGraph<EdgeStatus>& _g,Capacity & cap):g(_g),capacity(cap),INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
    	,ek(_g)
#endif
    {
    	setAllEdgeCapacities(1);
    }
    void setCapacity(int u, int w, int c){
    	//C.growTo(g.edges);
    	//C[ ]=c;

    }
    void setAllEdgeCapacities(int c){

    }
    int maxFlow(int s, int t){
    	int f = 0;

    	//C.growTo(g.nodes);
#ifdef DEBUG_MAXFLOW
    	for(int i = 0;i<g.all_edges.size();i++){
    		int id = g.all_edges[i].id;
    		int cap = capacity(id);
    		int from =  g.all_edges[i].from;
    		int to =  g.all_edges[i].to;

    		ek.setCapacity(from,to,cap);
    	}
#endif
    	if(rev.size()<g.all_edges.size()){
    		rev.clear();

    		rev.growTo(g.all_edges.size());
    		for(int i = 0;i<g.all_edges.size();i++){
    			rev[i]=-1;
    			int from = g.all_edges[i].from;
    			int to = g.all_edges[i].to;
    			for(int j = 0;j<g.adjacency[to].size();j++){
    				if(g.adjacency[to][j].node == from){
    					rev[i]=g.adjacency[to][j].id;
    					break;
    				}
    			}
    		}
    	}
    	F.clear();
    	F.growTo(g.all_edges.size());
    	prev.growTo(g.nodes);
    	M.growTo(g.nodes);

    	for(int i = 0;i<g.nodes;i++){
    		prev[i].from =-1;
    		M[i]=0;
    	}
    	prev[s].from = -2;
    	 M[s] = INF;
        while(true){
        	int m= BreadthFirstSearch(s,t);

            if (m == 0)
                break;

            f = f + m;

            int v = t;
            while (v!=  s){
                int u = prev[v].from;
                int id = prev[v].id;
                F[id] = F[id] + m;
                if(rev[id]>-1){
                	F[rev[id]]-=m;
                }
               // F[v][u] = F[v][u] - m;
                v = u;
            }
        }
#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
    	assert(f==expected_flow);
#endif

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
    			if(capacity(id) - F[id] == 0){
    				cut.push(Edge{u,v});
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
};
#endif

