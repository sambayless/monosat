#ifndef EDMONDS_KARP_ADJ_H
#define EDMONDS_KARP_ADJ_H

// Maximum Flow- Edmonds-Karp Algorithm
// by Iraklis Karagkiozoglou <i.kar@windowslive.com>

#include "MaxFlow.h"
#include <vector>
#include "core/Config.h"
#include "dgl/EdmondsKarp.h"


template< class Capacity ,class EdgeStatus=vec<bool> >
class EdmondsKarpAdj:public MaxFlow{

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
    DynamicGraph<EdgeStatus>& g;
    Capacity & capacity;
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
            	   int c = capacity[id];

            	 //  int fr = F[id];
                   if (((c - F[id]) > 0) && (prev[v].from == -1)){
                       prev[v] = LocalEdge(u,id,false);
                       M[v] = min(M[u], c - F[id]);
                       if (v != t)
                           Q.push(v);
                       else
                           return M[t];
                   }
               }

               for (int i = 0;i<g.inverted_adjacency[u].size();i++){
            	   int id = g.inverted_adjacency[u][i].id;
            	   if(!g.edgeEnabled(id))
						continue;

				   int v = g.inverted_adjacency[u][i].node;

				   int f = 0;
				   int c = F[id];

					  if (((c - f) > 0) && (prev[v].from == -1)){
						  prev[v] = LocalEdge(u,id,true);
						  M[v] = min(M[u], c - f);
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

    int maxFlow(int s, int t){
    	int f = 0;
#ifdef RECORD
		if(g.outfile && mincutalg==MinCutAlg::ALG_EDKARP_ADJ){
			fprintf(g.outfile,"f %d %d\n", s,t);
			fflush(g.outfile);
		}
#endif

    	//C.growTo(g.nodes);
/*
#ifdef DEBUG_MAXFLOW
    	for(int i = 0;i<g.all_edges.size();i++){
    		int id = g.all_edges[i].id;
    		int cap = capacity[id];
    		int from =  g.all_edges[i].from;
    		int to =  g.all_edges[i].to;

    		ek.setCapacity(from,to,cap);
    	}
#endif
*/
      	if(last_modification>0 && g.modifications==last_modification){
/*
#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
    	assert(curflow==expected_flow);
#endif
*/
        			return curflow;
        		}
    	/*if(rev.size()<g.all_edges.size()){
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
    	}*/
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
    			if(prev[v].backward){
					F[id] = F[id] - m;
				}else
					F[id] = F[id] + m;
                v = u;
            }

        }
/*
#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
    	assert(f==expected_flow);
#endif
*/
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

