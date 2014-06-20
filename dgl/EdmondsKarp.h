#ifndef EDMONDS_KARP_H
#define EDMONDS_KARP_H

// Maximum Flow- Edmonds-Karp Algorithm
// by Iraklis Karagkiozoglou <i.kar@windowslive.com>

#include "MaxFlow.h"
#include <vector>
#include <algorithm>
#include "graph/DynamicGraph.h"
namespace dgl{

class EdmondsKarp:public MaxFlow{
public:

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
    std::vector<std::vector<int> > F;//(Residual capacity from u to v is C[u,v] - F[u,v])
    std::vector<std::vector<int> > C;
    std::vector<int> P;
    std::vector<int> M;

    int curflow;

    int last_modification;
    int last_deletion;
    int last_addition;

    int history_qhead;
    int last_history_clear;

    DynamicGraph& g;
    int INF;
    /*
     *            input:
               C, E, s, t, F
           output:
               M[t]          (Capacity of path found)
               P             (Parent table)
     */
    std::vector<int> Q;

    int  BreadthFirstSearch(int s, int t){
     	for(int i = 0;i<g.nodes();i++)
			P[i]=-1;
     	P[s] = -2;
    	Q.clear();
           Q.push_back(s);

           //while (Q.size() > 0){
        	   //int u =Q.back(); Q.pop_back();
           for(int j = 0;j<Q.size();j++){
        	   int u = Q[j];

               for (int i = 0;i<g.nIncident(u);i++){
            	   if(!g.edgeEnabled(g.incident(u,i).id))
						continue;
            	   int v = g.incident(u,i).node;
                   ///(If there is available capacity, and v is not seen before in search)
            	   int c = C[u][v];
            	   int f = F[u][v];

                   if (((C[u][v] - F[u][v]) > 0) && (P[v] == -1)){
                       P[v] = u;
                       M[v] = std::min(M[u], C[u][v] - F[u][v]);
                       if (v != t)
                           Q.push_back(v);
                       else
                           return M[t];
                   }
               }
               //Must also try reverse edges
               for (int i = 0;i<g.nIncoming(u);i++){
				   if(!g.edgeEnabled(g.incoming(u,i).id))
						continue;
				   int v = g.incoming(u,i).node;
					///(If there is available capacity, and v is not seen before in search)
				   int c = C[u][v];
				   int f = F[u][v];

					if (((C[u][v] - F[u][v]) > 0) && (P[v] == -1)){
						P[v] = u;
						M[v] = std::min(M[u], C[u][v] - F[u][v]);
						if (v != t)
							Q.push_back(v);
						else
							return M[t];
					}
				}
           }
           return 0;


	   }
public:
    EdmondsKarp(DynamicGraph& _g):g(_g),INF(0xF0F0F0)
    {
        curflow=-1;

    	last_modification=-1;
    	last_deletion=-1;
    	last_addition=-1;

    	history_qhead=-1;
    	last_history_clear=-1;
    	setAllEdgeCapacities(1);
    }
    void setCapacity(int u, int w, int c){
    	if(C.size()<g.nodes()){
    		C.resize(g.nodes());
    		for(int i = 0;i<g.nodes();i++){
    			C[i].resize(g.nodes());
    		}
    	}
    	C[u][w]=c;
    }
    void setAllEdgeCapacities(int c){
    	for(int i = 0;i<g.nodes();i++){
    		for(int j = 0;j<g.nIncident(i);j++){
    			if(!g.edgeEnabled(g.incident(i,j).id))
						continue;
    			setCapacity(i,g.incident(i,j).node,c);
    		}
    	}
    }
    int maxFlow(int s, int t){
    	if(last_modification>0 && g.modifications==last_modification){

    			return curflow;
    		}
    	int f = 0;
    	C.resize(g.nodes());
    	F.resize(g.nodes());
    	P.resize(g.nodes());
    	M.resize(g.nodes());

    	for(int i = 0;i<g.nodes();i++){
    		P[i]=-1;
    		M[i]=0;
    		F[i].resize(g.nodes());
    		C[i].resize(g.nodes());
    		for(int j = 0;j<g.nodes();j++){
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
        curflow=f;

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
        return f;
    }
    int maxFlow(int s, int t, int max_length){
       	int f = 0;
       	C.resize(g.nodes());
       	F.resize(g.nodes());
       	P.resize(g.nodes());
       	M.resize(g.nodes());

       	for(int i = 0;i<g.nodes();i++){
       		P[i]=-1;
       		M[i]=0;
       		F[i].resize(g.nodes());
       		C[i].resize(g.nodes());
       		for(int j = 0;j<g.nodes();j++){
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

    std::vector<bool> seen;
    std::vector<bool> visited;

    bool minCut(int s, int t, int max_length, std::vector<Edge> & cut){
    	int f = maxFlow(s,t,max_length);
    	if(f>max_length){
    		return false;
    	}
		//ok, now find the cut
		Q.clear();
		Q.push_back(s);
		seen.clear();
		seen.resize(g.nodes());
		seen[s]=true;
	//	visited.clear();
		//visited.resize(g.nodes());
	//	visited[s]=true;
		for(int j = 0;j<Q.size();j++){
		   int u = Q[j];

			for(int i = 0;i<g.nIncident(u);i++){
				if(!g.edgeEnabled(g.incident(u,i).id))
					continue;
				int v = g.incident(u,i).node;
				int id = g.incident(u,i).id;
				if(C[u][v] - F[u][v] == 0){
					cut.push_back(Edge{u,v,id});
				}else if(!seen[v]){
					Q.push_back(v);
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
		cut.resize(j);
		return true;
    }

    int getEdgeFlow(int edgeid){
    	assert(g.edgeEnabled(edgeid));
    	int u = g.all_edges[edgeid].from;
    	int v = g.all_edges[edgeid].to;
    	return F[u][v];
    }
    int getEdgeCapacity(int id){
    	assert(g.edgeEnabled(id));
    	int u = g.all_edges[id].from;
    	int v = g.all_edges[id].to;
    	return C[u][v];
       }

      int getEdgeResidualCapacity(int id){
    	  assert(g.edgeEnabled(id));
		int u = g.all_edges[id].from;
		int v = g.all_edges[id].to;
		return C[u][v]-F[u][v];

      }
    int minCut(int s, int t, std::vector<Edge> & cut){
    	int f = maxFlow(s,t);
    	//ok, now find the cut
    	Q.clear();
    	Q.push_back(s);
    	seen.clear();
    	seen.resize(g.nodes());
    	seen[s]=true;
    //	visited.clear();
    	//visited.resize(g.nodes());
    //	visited[s]=true;
    	for(int j = 0;j<Q.size();j++){
		   int u = Q[j];

    		for(int i = 0;i<g.nIncident(u);i++){
    			if(!g.edgeEnabled(g.incident(u,i).id))
    				continue;
    			int v = g.incident(u,i).node;

    			if(C[u][v] - F[u][v] == 0){
    				cut.push_back(Edge{u,v});
    			}else if(!seen[v]){
    				Q.push_back(v);
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
    	cut.resize(j);
    	return f;
    }


};

};
#endif

