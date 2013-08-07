#ifndef EDMONDS_KARP_H
#define EDMONDS_KARP_H

// Maximum Flow- Edmonds-Karp Algorithm
// by Iraklis Karagkiozoglou <i.kar@windowslive.com>

#include "MaxFlow.h"
#include "mtl/Vec.h"
using namespace std;
using namespace Minisat;

template<class EdgeStatus=vec<bool> >
class EdmondsKarp:public MaxFlow{


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
    vec<vec<int> > F;//(Residual capacity from u to v is C[u,v] - F[u,v])
    vec<vec<int> > C;
    vec<int> P;
    vec<int> M;
    DynamicGraph<EdgeStatus>& g;
    int INF;
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
			P[i]=-1;
     	P[s] = -2;
    	Q.clear();
           Q.push(s);

           //while (Q.size() > 0){
        	   //int u =Q.last(); Q.pop();
           for(int j = 0;j<Q.size();j++){
        	   int u = Q[j];

               for (int i = 0;i<g.adjacency[u].size();i++){
            	   if(!g.edgeEnabled(g.adjacency[u][i].id))
						continue;
            	   int v = g.adjacency[u][i].node;
                   ///(If there is available capacity, and v is not seen before in search)
            	   int c = C[u][v];
            	   int f = F[u][v];

                   if (((C[u][v] - F[u][v]) > 0) && (P[v] == -1)){
                       P[v] = u;
                       M[v] = min(M[u], C[u][v] - F[u][v]);
                       if (v != t)
                           Q.push(v);
                       else
                           return M[t];
                   }
               }
               //Must also try reverse edges
               for (int i = 0;i<g.inverted_adjacency[u].size();i++){
				   if(!g.edgeEnabled(g.inverted_adjacency[u][i].id))
						continue;
				   int v = g.inverted_adjacency[u][i].node;
					///(If there is available capacity, and v is not seen before in search)
				   int c = C[u][v];
				   int f = F[u][v];

					if (((C[u][v] - F[u][v]) > 0) && (P[v] == -1)){
						P[v] = u;
						M[v] = min(M[u], C[u][v] - F[u][v]);
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
    EdmondsKarp(DynamicGraph<EdgeStatus>& _g):g(_g),INF(0xF0F0F0)
    {
    	setAllEdgeCapacities(1);
    }
    void setCapacity(int u, int w, int c){
    	if(C.size()<g.nodes){
    		C.growTo(g.nodes);
    		for(int i = 0;i<g.nodes;i++){
    			C[i].growTo(g.nodes);
    		}
    	}
    	C[u][w]=c;
    }
    void setAllEdgeCapacities(int c){
    	for(int i = 0;i<g.nodes;i++){
    		for(int j = 0;j<g.adjacency[i].size();j++){
    			if(!g.edgeEnabled(g.adjacency[i][j].id))
						continue;
    			setCapacity(i,g.adjacency[i][j].node,c);
    		}
    	}
    }
    int maxFlow(int s, int t){
    	int f = 0;
    	C.growTo(g.nodes);
    	F.growTo(g.nodes);
    	P.growTo(g.nodes);
    	M.growTo(g.nodes);

    	for(int i = 0;i<g.nodes;i++){
    		P[i]=-1;
    		M[i]=0;
    		F[i].growTo(g.nodes);
    		C[i].growTo(g.nodes);
    		for(int j = 0;j<g.nodes;j++){
    			F[i][j]=0;
    		}
    	}
    	P[s] = -2;
    	 M[s] = INF;
        while(true){
        	int m= BreadthFirstSearch(s,t);

            if (m == 0)
                break;

            f = f + m;

            int v = t;
            while (v!=  s){
                int u = P[v];
                F[u][v] = F[u][v] + m;
                F[v][u] = F[v][u] - m;
                v = u;
            }
        }
        return f;
    }
    int maxFlow(int s, int t, int max_length){
       	int f = 0;
       	C.growTo(g.nodes);
       	F.growTo(g.nodes);
       	P.growTo(g.nodes);
       	M.growTo(g.nodes);

       	for(int i = 0;i<g.nodes;i++){
       		P[i]=-1;
       		M[i]=0;
       		F[i].growTo(g.nodes);
       		C[i].growTo(g.nodes);
       		for(int j = 0;j<g.nodes;j++){
       			F[i][j]=0;
       		}
       	}
       	P[s] = -2;
       	 M[s] = INF;
           while(true){
           	int m= BreadthFirstSearch(s,t);

               if (m == 0)
                   break;

               f = f + m;
               if(f>max_length){
            	   return f;
               }
               int v = t;
               while (v!=  s){
                   int u = P[v];
                   F[u][v] = F[u][v] + m;
                   F[v][u] = F[v][u] - m;
                   v = u;
               }
           }
           return f;
       }

    vec<bool> seen;
    vec<bool> visited;

    bool minCut(int s, int t, int max_length, vec<Edge> & cut){
    	int f = maxFlow(s,t,max_length);
    	if(f>max_length){
    		return false;
    	}
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

				if(C[u][v] - F[u][v] == 0){
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
		return true;
    }

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

    			if(C[u][v] - F[u][v] == 0){
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

