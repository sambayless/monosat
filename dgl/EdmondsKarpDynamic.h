#ifndef EDMONDS_KARP_DYNAMIC_H
#define EDMONDS_KARP_DYNAMIC_H

//dynamic edmonds_karp, following http://cstheory.stackexchange.com/questions/9938/incremental-maximum-flow-in-dynamic-graphs

#include "MaxFlow.h"
#include <vector>
#include "EdmondsKarp.h"
#include "EdmondsKarpAdj.h"

namespace dgl{
template< class Capacity  >
class EdmondsKarpDynamic:public MaxFlow{
	int f = 0;
	std::vector<int> F;

	int shortCircuitFlow=0;

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
    std::vector<LocalEdge> prev;
    std::vector<int> M;
    DynamicGraph& g;
    Capacity & capacity;
    int INF;
#ifdef DEBUG_MAXFLOW
    	EdmondsKarpAdj<Capacity, EdgeStatus> ek;
#endif
    /*
     *            input:
               C, E, s, t, F
           output:
               M[t]          (Capacity of path found)
               P             (Parent table)
     */
    std::vector<int> Q;


    std::vector<bool> edge_enabled;

    int  BreadthFirstSearch(int s, int t,int bound=-1){
     	for(int i = 0;i<g.nodes;i++){//this has to go...
     		prev[i].from=-1;
     	}
     	prev[s].from = -2;
    	Q.clear();
           Q.push_back(s);
           bool found = false;
         int old_m = M[s];
       	for(int j = 0;j<Q.size();j++){
   			   int u = Q[j];

               for (int i = 0;i<g.adjacency[u].size();i++){
            	   if(!edge_enabled[g.adjacency[u][i].id])
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
                           Q.push_back(v);
                       else{
                    	   found=true;
                    	   break;
                       }

                           //return M[t];
                   }
               }
               for (int i = 0;i<g.inverted_adjacency[u].size();i++){
				   int id = g.inverted_adjacency[u][i].id;
				   if(!edge_enabled[(g.inverted_adjacency[u][i].id)])
						continue;

				   int v = g.inverted_adjacency[u][i].node;
					  ///(If there is available capacity, and v is not seen before in search)

				   int f = 0;
				   int c = F[id];

				 //  int fr = F[id];
					  if (((c - f) > 0) && (prev[v].from == -1)){
						  prev[v] = LocalEdge(u,id,true);
						  M[v] = min(M[u], c - f);
				/*		  */
						  if (v != t)
							  Q.push_back(v);
						  else{
	                    	   found=true;
	                    	   break;
	                       }
					  }
				  }
           }
       	M[s]=old_m;
       		if(found){
       			return M[t];
       		}
           return 0;


	   }


public:
    EdmondsKarpDynamic(DynamicGraph& _g,Capacity & cap):g(_g),reserve(cap),INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
    	,ek(_g,cap)
#endif
    {
    	  curflow=0;
      	last_modification=-1;
      	last_deletion=-1;
      	last_addition=-1;

      	history_qhead=0;
      	last_history_clear=-1;
    	//setAllEdgeCapacities(1);
    }
    void setCapacity(int u, int w, int c){
    	//C.resize(g.edges);
    	//C[ ]=c;

    }
    void setAllEdgeCapacities(int c){

    }
    int maxFlow(int s, int t){
    	//see http://cstheory.stackexchange.com/a/10186
    	static int it = 0;
    	if(++it==32){
    		int a=1;
    	}
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"f %d %d\n", s,t);
			fflush(g.outfile);
		}
#endif

    	//C.resize(g.nodes);
#ifdef DEBUG_MAXFLOW
    	for(int i = 0;i<g.all_edges.size();i++){
    		int id = g.all_edges[i].id;
    		int cap = capacity[id];
    		int from =  g.all_edges[i].from;
    		int to =  g.all_edges[i].to;

    		ek.setCapacity(from,to,cap);
    	}
#endif
      	if(last_modification>0 && g.modifications==last_modification){
#ifdef DEBUG_MAXFLOW
      		int expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
      		assert(curflow==expected_flow);
#endif
        	return curflow;
        }else if (last_modification<=0 || g.historyclears!=last_history_clear  || g.changed()){
        	F.clear();
        	F.resize(g.all_edges.size());
			prev.resize(g.nodes);
			M.resize(g.nodes);
			f=0;
			for(int i = 0;i<g.nodes;i++){
				prev[i].from =-1;
				M[i]=0;
			}
			prev[s].from = -2;
			 M[s] = INF;
			 edge_enabled.resize(g.edges);
			 for(int i = 0;i<g.edges;i++)
				 edge_enabled[i]= g.isEdge(i) && g.edgeEnabled(i);

			 dbg_print_graph(s,t,-1, -1);
			 f = maxFlow_p(s,t);
			 dbg_print_graph(s,t,-1, -1);
        }

#ifdef DEBUG_MAXFLOW
    	for(int i = 0;i<g.edges;i++){
			if(edge_enabled[i]){
				int fi = F[i];
				assert(F[i]<= capacity[i]);
			}else{
				assert(F[i]==0);
			}
		}
#endif
    	bool added_Edges=false;
    	bool needsReflow = false;

    	for (int i = history_qhead;i<g.history.size();i++){
    			int edgeid = g.history[i].id;
    			if(g.history[i].addition && g.edgeEnabled(edgeid)){
    				added_Edges=true;
    				edge_enabled[edgeid]=true;
    			}else if (!g.history[i].addition &&  !g.edgeEnabled(edgeid)){
    				//assert(edge_enabled[edgeid]);
    				edge_enabled[edgeid]=false;
    				int fv = F[edgeid]; //g.all_edges[edgeid].from;
    				if(fv==0){
    					//do nothing.
    				}else{
    					//ok, check if the maxflow from u to v has not lowered now that we've removed this edge.
    					//if it hasn't, then we are still safe
    					int u = g.all_edges[edgeid].from;
    					int v = g.all_edges[edgeid].to;

    					if (fv<0){
    						std::swap(u,v);
    						fv = -fv;
    					}
    					assert(fv>0);
    					int flow = maxFlow_residual(u,v, fv);//fix this!
    					assert(flow<=fv);
    					if(flow==fv){
    						//then we are ok.
    					}else{
    						//the total flow in the network has to be decreased by delta.
    						int delta = fv- flow;
    						assert(delta>0);
    						needsReflow=true;
    						//temporarily connect s and t by an arc of infinite capacity and run maxflow algorithm again from vin to vout
    						flow = maxFlow_p(u,v,s,t,delta);
    						//f-=flow;//this doesn't work
    					}
    					F[edgeid]=0;
    				}
    			}
    		}

    		//recompute the s-t flow
    		if (needsReflow){
    			f = 0;
    			for(auto edge:g.adjacency[s]){
    				if (edge_enabled[edge.id]){
    					f+=F[edge.id];
    				}else{
    					assert(F[edge.id]==0);
    				}
    			}
#ifndef NDEBUG
    			for(auto edge:g.inverted_adjacency[s]){
    				//There shouldn't be any backwards flow to s in a maximum flow
    				assert(F[edge.id]==0);
				}
#endif
    		}
    		dbg_check_flow(s,t);
    		if(added_Edges)
    			f = maxFlow_p(s,t);

#ifdef DEBUG_MAXFLOW
    	int expected_flow =ek.maxFlow(s,t);
#endif
   		dbg_print_graph(s,t,-1, -1);
#ifdef DEBUG_MAXFLOW
    	assert(f==expected_flow);
    	for(int i = 0;i<g.edges;i++){
			if(g.edgeEnabled(i)){
				int fi = F[i];
				assert(F[i]<= capacity[i]);
			}else{
				assert(F[i]==0);
			}
		}
#endif

#ifndef NDEBUG
    	assert(edge_enabled.size()==g.edges);
    	for(int i = 0;i<g.edges;i++)
    		assert(edge_enabled[i]==g.edgeEnabled(i));
    	dbg_check_flow(s,t);
#endif

    	curflow=f;
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
        return f;
    }


private:

    int maxFlow_residual(int s, int t, int bound){
#ifndef NDEBUG
    	DynamicGraph<> d;
		d.addNodes(g.nodes);
		//std::vector<int> R;
		for(int i = 0;i<g.edges;i++){
			if(g.isEdge(i)){
				if(edge_enabled[i]){
					int r =capacity[i]-F[i];
					d.addEdge(g.all_edges[i].from,g.all_edges[i].to,g.all_edges[i].id,r);
				}else  {
					d.addEdge(g.all_edges[i].from,g.all_edges[i].to,g.all_edges[i].id,0);
					d.disableEdge(g.all_edges[i].id);
				}
			}
		}
		for(int i = 0;i<g.edges;i++){
			if(edge_enabled[i]){
				if(F[i]>0){
					d.addEdge(g.all_edges[i].to,g.all_edges[i].from,-1,F[i]);
				}
			}
		}
#endif
    	int new_flow = 0;
        	 while(true){
        		 dbg_print_graph(s,t,-1, -1);
    			int m= BreadthFirstSearch(s,t,bound);
    			if(bound >=0 && new_flow+m>bound){
					m=bound-new_flow;
				}

    			if (m <= 0)
    				break;

    			new_flow = new_flow + m;

    			int v = t;
    			while (v!=  s){
    				int u = prev[v].from;
    				int id = prev[v].id;
    				if(prev[v].backward){
    					F[id] = F[id] - m;
    				}else
    					F[id] = F[id] + m;
    			/*	if(rev[id]>-1){
    					F[rev[id]]-=m;
    				}*/
    			   // F[v][u] = F[v][u] - m;
    				assert(F[id]<= capacity[id]);
    				v = u;
    			}
    			dbg_print_graph(s,t,-1, -1);
    		}
#ifndef NDEBUG
			EdmondsKarpAdj<std::vector<int>, EdgeStatus> ek_check(d,d.weights);

        		int expect =  ek_check.maxFlow(s,t);
        		assert(new_flow<=expect);
    #endif

        	 return new_flow;
        }

    int maxFlow_p(int s, int t){
    	dbg_print_graph(s,t,-1, -1);
    	 while(true){
			int m= BreadthFirstSearch(s,t);

			if (m == 0)
				break;

			f = f + m;

			int v = t;
			while (v!=  s){
				int u = prev[v].from;
				int id = prev[v].id;
				if(id==29){
					int a=1;
				}
				if(prev[v].backward){
					F[id] = F[id] - m;
			/*		if(rev[id]>-1){
						if(rev[id]==29){
							int a=1;
						}
						F[rev[id]]+=m;
					}*/
				}else{
					F[id] = F[id] + m;

				/*	if(rev[id]>-1){
						F[rev[id]]-=m;
					}*/
				}
				assert(F[id]<= capacity[id]);
			   // F[v][u] = F[v][u] - m;
				v = u;
			}

		}

#ifndef NDEBUG
    		//EdmondsKarp<EdgeStatus> ek_check(g);
    	 	 EdmondsKarpAdj<std::vector<int>, EdgeStatus> ek_check(g,g.weights);
    		int expect =  ek_check.maxFlow(s,t);
    		assert(f==expect);
#endif
    	 return f;
    }

    int  BreadthFirstSearch(int s, int t, int shortCircuitFrom, int shortCircuitTo, int shortCircuitCapacity, int & shortCircuitFlow){
         	for(int i = 0;i<g.nodes;i++)//this has to go...
         		prev[i].from=-1;
         	prev[s].from = -2;
        	Q.clear();
            Q.push_back(s);

           	for(int j = 0;j<Q.size();j++){
       			   int u = Q[j];

       			   if(u==shortCircuitFrom){
       				   int v = shortCircuitTo;
       				   int f = shortCircuitFlow;
       				   int c = shortCircuitCapacity;
       				   if (((c - f) > 0)){
       					   prev[v] = {u,-1};
						   M[v] = min(M[u], c -f);
						   if (v != t)
							   Q.push_back(v);
						   else
							   return M[t];
					   }
       			   }

                   for (int i = 0;i<g.adjacency[u].size();i++){
                	   if(!edge_enabled[g.adjacency[u][i].id])
    						continue;
                	   int id = g.adjacency[u][i].id;
                	   int v = g.adjacency[u][i].node;
                       ///(If there is available capacity, and v is not seen before in search)
                	   if(id==27 || id==29){
                							   int a=1;
                						   }
                	   int f = F[id];
                	   int c = capacity[id];

                	 //  int fr = F[id];
                       if (((c - f) > 0) && (prev[v].from == -1)){
                           prev[v] = LocalEdge(u,id,false);
                           M[v] = min(M[u], c - f);
                           if (v != t)
                               Q.push_back(v);
                           else
                               return M[t];
                       }
                   }

                   for (int i = 0;i<g.inverted_adjacency[u].size();i++){
                	   int id = g.inverted_adjacency[u][i].id;
					   if(!edge_enabled[(g.inverted_adjacency[u][i].id)])
							continue;
					   if(id==27 || id==29){
						   int a=1;
					   }
					   int v = g.inverted_adjacency[u][i].node;
						  ///(If there is available capacity, and v is not seen before in search)

					   int f = 0;
					   int c = F[id];

					 //  int fr = F[id];
						  if (((c - f) > 0) && (prev[v].from == -1)){
							  prev[v] = LocalEdge(u,id,true);
							  M[v] = min(M[u], c - f);
							  if (v != t)
								  Q.push_back(v);
							  else
								  return M[t];
						  }
					  }

               }
               return 0;


    	   }
    void dbg_print_graph(int from, int to, int shortCircuitFrom=-1, int shortCircuitTo=-1){
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
    				if(edge_enabled[i]){
						auto & e = g.all_edges[i];
						const char * s = "black";
						/*if(value(e.v)==l_True)
							s="blue";
						else if (value(e.v)==l_False)
							s="red";*/
						printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"%s\"]\n", e.from,e.to, i, F[i],g.weights[i] , s);
    				}
    			}

    			if(shortCircuitFrom>=0){
    				printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"black\"]\n", shortCircuitFrom,shortCircuitTo, -1, shortCircuitFlow,-1 );
    			}

    			printf("}\n");
#endif
    		}

    int maxFlow_p(int s, int t, int shortCircuitFrom, int shortCircuitTo, int bound){
    	//
    	int newFlow = 0;
#ifndef NDEBUG
    	 	DynamicGraph<> d;
    	 	d.addNodes(g.nodes);
    	 	//std::vector<int> R;
    	 	for(int i = 0;i<g.edges;i++){
    	 		if(edge_enabled[i]){
    	 			int r =capacity[i]-F[i];
    	 			d.addEdge(g.all_edges[i].from,g.all_edges[i].to,g.all_edges[i].id,r);
    	 		}else if(g.isEdge(i)){
    	 			d.addEdge(g.all_edges[i].from,g.all_edges[i].to,g.all_edges[i].id,0);
    	 			d.disableEdge(g.all_edges[i].id);
    	 		}
    	 	}
    		for(int i = 0;i<g.edges;i++){
				if(edge_enabled[i]){
					if(F[i]>0){
						d.addEdge(g.all_edges[i].to,g.all_edges[i].from,-1,F[i]);
					}
				}
			}

    		int stflow = 0;
			for(auto edge:g.inverted_adjacency[s]){
				if (edge_enabled[edge.id]){
					stflow+=F[edge.id];
				}
			}

#endif

    	shortCircuitFlow=0;
        	 while(true){
        		 dbg_print_graph(s,t,shortCircuitFrom, shortCircuitTo);
    			int m= BreadthFirstSearch(s,t,shortCircuitFrom,shortCircuitTo,bound,shortCircuitFlow);
    			if(bound >=0 && newFlow+m>bound){
						m=bound-newFlow;
					}
    			if (m == 0)
    				break;

    			newFlow = newFlow + m;

    			int v = t;
    			while (v!=  s){
    				int u = prev[v].from;
    				if(u==s){
    					int a=1;
    				}
    				int id = prev[v].id;
    				if(id>=0){
    					if(prev[v].backward){
    						F[id] = F[id] - m;
						/*	if(rev[id]>-1){
								F[rev[id]]+=m;
							}*/
    					}else{
							F[id] = F[id] + m;
					/*		if(rev[id]>-1){
								F[rev[id]]-=m;
							}*/
    					}
        				assert(id>=0);
        				assert(id<F.size());
        				assert(id<capacity.size());
        				assert(F[id]<= capacity[id]);
    				}else{
    					u = shortCircuitFrom;
    				}

    			   // F[v][u] = F[v][u] - m;
    				v = u;
    			}
    		}
#ifndef NDEBUG
        	 dbg_print_graph(s,t,shortCircuitFrom, shortCircuitTo);
    	 	if(shortCircuitFrom>=0){
    	 		d.addEdge(shortCircuitFrom,shortCircuitTo,d.edges,100);
    	 	}
    		EdmondsKarpAdj<std::vector<int>, EdgeStatus> ek_check(d,d.weights);

    		int expect =  ek_check.maxFlow(s,t);


    		assert(newFlow<=expect);

    		for(int i = 0;i<g.edges;i++){

					int fi = F[i];
					assert(F[i]<= capacity[i]);

    		}

#endif
        	 return newFlow;
        }

    void dbg_check_flow(int s, int t){
#ifndef NDEBUG
    	for(int u = 0;u<g.nodes;u++){
    		int inflow = 0;
    		int outflow = 0;
            for (int i = 0;i<g.inverted_adjacency[u].size();i++){
				   int id = g.inverted_adjacency[u][i].id;
				   assert(id<edge_enabled.size());
				   if(!edge_enabled[(g.inverted_adjacency[u][i].id)])
						continue;

				   inflow+=F[id];
            }
            for (int i = 0;i<g.adjacency[u].size();i++){
				   int id = g.adjacency[u][i].id;
				   assert(id<edge_enabled.size());
				   if(!edge_enabled[(g.adjacency[u][i].id)])
						continue;

				   outflow+=F[id];
            }
            if(u!=s && u != t){
            	assert(inflow==outflow);
            }else if (u==s){
            	assert(outflow==f);
            }else if (u==t){
            	assert(inflow==f);
            }
    	}

#endif
    }

    std::vector<bool> seen;
    std::vector<bool> visited;
public:
    int minCut(int s, int t, std::vector<Edge> & cut){
    	int f = maxFlow(s,t);
    	//ok, now find the cut
    	Q.clear();
    	Q.push_back(s);
    	seen.clear();
    	seen.resize(g.nodes);
    	seen[s]=true;

    	for(int j = 0;j<Q.size();j++){
		   int u = Q[j];

    		for(int i = 0;i<g.adjacency[u].size();i++){
    			if(!g.edgeEnabled(g.adjacency[u][i].id))
    				continue;
    			int v = g.adjacency[u][i].node;
    			int id = g.adjacency[u][i].id;
    			if(capacity[id] - F[id] == 0){
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
    	return f;
    }
    int getEdgeCapacity(int id){
     	assert(g.edgeEnabled(id));
     	return capacity[id];
     }
    int getEdgeFlow(int id){
    	assert(g.edgeEnabled(id));
    	return F[id];// reserve(id);
    }
    int getEdgeResidualCapacity(int id){
    	assert(g.edgeEnabled(id));
    	return  capacity[id]-F[id];// reserve(id);
    }
};
};
#endif

