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

#ifndef MAX_FLOW_IBFS_H
#define MAX_FLOW_IBFS_H


#include "alg/ibfs.h"
#include <vector>
#ifdef DEBUG_MAXFLOW
#include "EdmondsKarp.h"
#endif

namespace dgl{

class IBFS:public MaxFlow<int>{

    DynamicGraph& g;
    int INF;
#ifdef DEBUG_MAXFLOW
    	EdmondsKarp<int> ek;
#endif

    std::vector<std::vector<int> > C;
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
      	if(C.size()<g.nodes()){
        		C.resize(g.nodes());
        		for(int i = 0;i<g.nodes();i++){
        			C[i].resize(g.nodes());
        		}
        	}
        	C[u][w]=c;
#ifdef DEBUG_MAXFLOW
    	ek.setCapacity(u,w,c);
#endif
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

    	ibfs =  new IBFSGraph<int,int,int>(g.nodes(),edges) ;
    	while(ibfs->nNodes < g.nodes()){
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


    int minCut(int s, int t, std::vector<MaxFlowEdge> & cut){
    	int f = maxFlow(s,t);
    	//ok, now find the cut
    	for(int u = 0;u<g.nodes();u++){
    		IBFSGraph<int,int,int>::termtype utype = ibfs->what_segment(u,IBFSGraph<int,int,int>::SOURCE);
    		if(utype ==IBFSGraph<int,int,int>::SOURCE ){
				for(int j = 0;j<g.nIncident(u);j++){
					if(!g.edgeEnabled(g.incident(u,j).id))
						continue;
					int v = g.incident(u,j).node;
					int id =  g.incident(u,j).id;
					if( ibfs->what_segment(v,IBFSGraph<int,int,int>::SOURCE) ==IBFSGraph<int,int,int>::SINK ){
						//then this is on the cut
						cut.push_back(MaxFlowEdge{u,v,id});
					}
				}
    		}
    	}

    	return f;
    }


    int getEdgeFlow(int edgeid){
      	assert(g.edgeEnabled(edgeid));
      	int u = g.all_edges[edgeid].from;
      	int v = g.all_edges[edgeid].to;
      	return ibfs->edgeflow(u,v);
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
  		return C[u][v]-ibfs->edgeflow(u,v);

        }


/*




    typedef struct node_t node_t;
    typedef struct edge_t edge_t;
    struct edge_t {
        int cap, ni,i;
    };

    DynamicGraph g;
    std::vector<std::vector<int> > cap;
    std::vector<std::vector<int> > f;
    int edmondskarp(int source, int sink, int n){
        int max = 0;
        static std::vector<int> q;
        static std::vector<int> mins;
        q.resize(g.nodes());
        mins.resize(g.nodes());
        std::vector<int> pre;
        std::vector<int> ni;
        while(true){
            int   h=0,t=0,c,i,j,min=g.nodes();

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
};
#endif

