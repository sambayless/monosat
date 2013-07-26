#ifndef MAX_FLOW_IBFS_H
#define MAX_FLOW_IBFS_H

    // Maximum Flow- Edmonds-Karp Algorithm
    // by Iraklis Karagkiozoglou <i.kar@windowslive.com>
    #include <iostream>
    #include <climits>
#include "ibfs.h"
#include "mtl/Vec.h"
#ifdef DEBUG_MAXFLOW
#include "EdwardsKarp.h"
#endif

    using namespace std;
    using namespace Minisat;

class IBFS{

    DynamicGraph& g;
#ifdef DEBUG_MAXFLOW
    	EdmondsKarp ek;
#endif
    int INF;
    vec<vec<int> > C;
    IBFSGraph<int,int,int>* ibfs ;
public:
    IBFS(DynamicGraph& _g):g(_g),INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
    	,ek(_g)
#endif
    {
    	ibfs=NULL;
    }
    void setCapacity(int u, int w, int c){
      	if(C.size()<g.nodes){
        		C.growTo(g.nodes);
        		for(int i = 0;i<g.nodes;i++){
        			C[i].growTo(g.nodes);
        		}
        	}
        	C[u][w]=c;
#ifdef DEBUG_MAXFLOW
    	ek.setCapacity(u,w,c);
#endif
    }
    void setAllEdgeCapacities(int c){
    	for(int i = 0;i<g.nodes;i++){
    		for(int j = 0;j<g.adjacency[i].size();j++){
    			setCapacity(i,g.adjacency[i][j],c);
    		}
    	}
    }
    int maxFlow(int s, int t){
    	if(ibfs)
    		delete(ibfs);
    	int edges = 0;
    	for(int i = 0;i<C.size();i++){
        		for(int j = 0;j<C[i].size();j++){
        			if(C[i][j]>0){
        				edges++;
        			}
        		}
        	}

    	ibfs =  new IBFSGraph<int,int,int>(g.nodes,edges) ;
    	while(ibfs->nNodes < g.nodes){
			ibfs->add_node(1);
		}

    	for(int i = 0;i<C.size();i++){
    		for(int j = 0;j<C[i].size();j++){
    			if(C[i][j]>0){
    				ibfs->add_edge(i,j,C[i][j],0);
    			}
    		}
    	}


    	int f = 0;

    	//add edges to the source and sink from s and t (as opposed to just designating s and t the source and sink)
    	ibfs->add_tweights(s,INF,0);
    	ibfs->add_tweights(t,0,INF);
#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif
    	f=ibfs->maxflow();
#ifdef DEBUG_MAXFLOW
    	assert(f==expected_flow);
#endif
        return f;
    }

    vec<int> q;
    vec<bool> seen;
    vec<bool> visited;
    struct Edge {
    	int u;
    	int v;
    };
    int minCut(int s, int t, vec<Edge> & cut){
    	int f = maxFlow(s,t);
    	//ok, now find the cut
    	for(int u = 0;u<g.nodes;u++){
    		IBFSGraph<int,int,int>::termtype utype = ibfs->what_segment(u,IBFSGraph<int,int,int>::SOURCE);
    		if(utype ==IBFSGraph<int,int,int>::SOURCE ){
				for(int j = 0;j<g.adjacency[u].size();j++){
					int v = g.adjacency[u][j];
					if( ibfs->what_segment(v,IBFSGraph<int,int,int>::SOURCE) ==IBFSGraph<int,int,int>::SINK ){
						//then this is on the cut
						cut.push(Edge{u,v});
					}
				}
    		}
    	}

    	return f;
    }





/*




    typedef struct node_t node_t;
    typedef struct edge_t edge_t;
    struct edge_t {
        int cap, ni,i;
    };

    DynamicGraph g;
    vec<vec<int> > cap;
    vec<vec<int> > f;
    int edmondskarp(int source, int sink, int n){
        int max = 0;
        static vec<int> q;
        static vec<int> mins;
        q.growTo(g.nodes);
        mins.growTo(g.nodes);
        vec<int> pre;
        vec<int> ni;
        while(true){
            int   h=0,t=0,c,i,j,min=g.nodes;

            int u; int w;
            q[t++]=source;
            while(t>h && pre[sink]==-1){
                c = q[h++];
                for(int i = 0;i<g.adjacency[c];i++){
                	int d = g.adjacency[c][i];
                	int index = ni[c][d];
                	if(cap[c][d]>0 && pre[j] == -1){
                		q[t++] = index;
                		pre[index] = i;
                		if(mins[c]>cap[c][d])
                			mins[index]=cap[c][d];
                		else
                			mins[index]=mins[c];
                	}
                }

            }
            if(pre[sink]==NULL) break;
            u=pre[sink];
            while(u!=NULL) {
            	cap[u]
                    (*u).cap-=mins[sink];
                    u=pre[(*u).i];
            }
            max+=mins[sink];
        }
        return max;
    }
*/
};
#endif

