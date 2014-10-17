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

#ifndef MAXFLOW_KOHLI_TORR_H
#define MAXFLOW_KOHLI_TORR_H

//Wrapper around Kohli and Torr's  Dynamic Graph Cuts algorithm (version 2)
//Note that the Kohli Torr implementation itself is under the GPL (Version 2); only this wrapper code is MIT licensed.

#include "MaxFlow.h"
#include <vector>
#include "alg/dyncut/graph.h"
#include "EdmondsKarpAdj.h"
#include <algorithm>
namespace dgl{
template< class Capacity,typename Weight  >
class KohliTorr:public MaxFlow<Weight>{
	Weight f = 0;
	/**
	 * Note: The Kohli Torr implementation does _not_ support multiple edges between the same nodes.
	 */
	kohli_torr::Graph<Weight,Weight,Weight> * kt=nullptr;
	struct LocalEdge{
		int from;
		int id;
		bool backward=false;
		LocalEdge(int from=-1, int id=-1, bool backward=false):from(from),id(id),backward(backward){

			}
	};
	bool dynamic=true;
	Weight curflow;
    int last_modification;
    int last_deletion;
    int last_addition;
    std::vector<int> tmp_edges;
    int history_qhead;
    int last_history_clear;

    typedef  typename kohli_torr::Graph<Weight,Weight,Weight>::arc_id arc;

    std::vector<int> arc_map; //map from edgeids to arcs (note: this may be many to one)

    DynamicGraph& g;
    Capacity & capacity;
    Weight INF;
#ifdef DEBUG_MAXFLOW
    	EdmondsKarpAdj<Capacity,Weight> ek;
#endif

    Weight max_capacity =0;

    std::vector<bool> edge_enabled;


public:
    KohliTorr(DynamicGraph& _g,Capacity & cap):g(_g),capacity(cap),INF(0xF0F0F0)
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

    }
    void setCapacity(int u, int w, Weight c){

    }
    void setAllEdgeCapacities(Weight c){

    }
	long num_updates=0;
	int numUpdates()const{
		return num_updates;
	}

    const  Weight maxFlow(int s, int t){
    	//see http://cstheory.stackexchange.com/a/10186
    	static int it = 0;
    	if(++it==56){
    		int a=1;
    	}
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"f %d %d\n", s,t);
			fflush(g.outfile);
		}
#endif
		dbg_print_graph(s,t);
    	//C.resize(g.nodes());
#ifdef DEBUG_MAXFLOW
    	for(int i = 0;i<g.all_edges.size();i++){
    		int id = g.all_edges[i].id;
    		Weight cap = capacity[id];
    		int from =  g.all_edges[i].from;
    		int to =  g.all_edges[i].to;

    		ek.setCapacity(from,to,cap);
    	}
#endif
      	if(last_modification>0 && g.modifications==last_modification){
#ifdef DEBUG_MAXFLOW
      		Weight expected_flow =ek.maxFlow(s,t);
#endif

#ifdef DEBUG_MAXFLOW
      		assert(curflow==expected_flow);
#endif
        	return curflow;
        }else if (last_modification<=0 || g.historyclears!=last_history_clear  || g.changed()){
        	edge_enabled.clear();

        	if(!kt){
        		kt = new kohli_torr::Graph<Weight,Weight,Weight> (g.nodes(),g.edges());
        	  	kt->maxflow(false);//just to initialize things.
        	}else{
        		kt->reset();
        	}

        	while(kt->get_node_num()<g.nodes()){
        		int node_id= kt->add_node();
        		assert(node_id==kt->get_node_num()-1);
        	}
        	edge_enabled.resize(g.edges(),false);
        	arc_map.clear();
        	arc_map.resize(g.edges(),-1);
        	for(int edgeID = 0;edgeID<g.edges();edgeID++){

        		if(!g.hasEdge(edgeID))
        			continue;

        		max_capacity+=capacity[edgeID];
        		int from = g.getEdge(edgeID).from;
        		int to = g.getEdge(edgeID).to;
        		edge_enabled[edgeID]=false;
        		if(!kt->has_edge(from,to)){
        			assert(arc_map[edgeID]==-1);
        			int arc_id = kt->add_edge(from,to,0,0);
        			arc_map[edgeID]=arc_id;
        			//set the corresponding arc for each other from-to edge
					for(int i = 0;i<g.nIncident(from,true);i++){
						int edgeid = g.incident(from,i,true).id;
						if(g.getEdge(edgeid).from==from && g.getEdge(edgeid).to==to){
							arc_map[edgeid]=arc_id;
						}else if  (g.getEdge(edgeid).from==to && g.getEdge(edgeid).to==from){
							arc_map[edgeid]=arc_id+1;//the reverse arc is always stored right after the forward arc
						}
					}
        		}
        		if(g.edgeEnabled(edgeID)){
        			edge_enabled[edgeID]=true;
        			kt->edit_edge_inc(from,to,capacity[edgeID],0);
        		}
        	}
        	kt->edit_tweights(s,max_capacity,0);
        	kt->edit_tweights(t,0,max_capacity);


        }

      	assert(kt);

    	bool added_Edges=false;
    	bool needsReflow = false;

    	for (int i = history_qhead;i<g.history.size();i++){
    			int edgeid = g.history[i].id;
    			if(g.history[i].addition && g.edgeEnabled(edgeid) && !edge_enabled[edgeid]){
    				added_Edges=true;
    				edge_enabled[edgeid]=true;
    				kt->edit_edge_inc(g.getEdge(edgeid).from,g.getEdge(edgeid).to,capacity[edgeid],0);

    			}else if (!g.history[i].addition &&  !g.edgeEnabled(edgeid) && edge_enabled[edgeid]){
    				assert(edge_enabled[edgeid]);
    				edge_enabled[edgeid]=false;
    				kt->edit_edge_inc(g.getEdge(edgeid).from,g.getEdge(edgeid).to,-capacity[edgeid],0);
    			}
    		}

    	f= kt->maxflow(dynamic);


#ifdef DEBUG_MAXFLOW
    		Weight expected_flow =ek.maxFlow(s,t);
#endif
   		dbg_print_graph(s,t);


#ifndef NDEBUG

    	for(int i = 0;i<g.edges();i++)
    		assert(edge_enabled[i]==g.edgeEnabled(i));
    	dbg_check_flow(s,t);
#endif

    	curflow=f;
    	num_updates++;
		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
        return f;
    }


private:



    void dbg_print_graph(int from, int to){
#ifndef NDEBUG
    	if(edge_enabled.size()<g.edges())
    		return;
    		static int it = 0;
    		if(++it==6){
    			int a =1;
    		}
    		printf("Graph %d\n", it);
    			printf("digraph{\n");
    			for(int i = 0;i<g.nodes();i++){
    				if(i==from){
    					printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
    				}else if (i==to){
    					printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
    				}else
    					printf("n%d\n", i);
    			}

    			for(int i = 0;i<g.edges();i++){
    				if(edge_enabled[i]){
						auto & e = g.all_edges[i];
						const char * s = "black";
						std::cout<<"n" << e.from <<" -> n" << e.to << " [label=\"" << i <<": " << getEdgeFlow(e.id) <<"/" << capacity[i]  << "\" color=\"" << s<<"\"]\n";
						//printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"%s\"]\n", e.from,e.to, i, F[i],capacity[i] , s);
    				}
    			}

    			printf("}\n");
#endif
    		}

    void dbg_check_flow(int s, int t){

    }

    std::vector<bool> seen;
    std::vector<bool> visited;
public:
    const Weight minCut(int s, int t, std::vector<MaxFlowEdge> & cut){
    	Weight f = maxFlow(s,t);
    	assert(false);
    	//ok, now find the cut
/*
    	Q.clear();
    	Q.push_back(s);
    	seen.clear();
    	seen.resize(g.nodes());
    	seen[s]=true;
    	dbg_print_graph(s,t);
    	//explore the residual graph
    	for(int j = 0;j<Q.size();j++){
		   int u = Q[j];

    		for(int i = 0;i<g.nIncident(u);i++){
    			if(!g.edgeEnabled(g.incident(u,i).id))
    				continue;
    			int v = g.incident(u,i).node;
    			int id = g.incident(u,i).id;
    			if(capacity[id] - F[id] == 0){
    				cut.push_back(MaxFlowEdge{u,v,id});//potential element of the cut
    			}else if(!seen[v]){
    				Q.push_back(v);
    				seen[v]=true;
    			}
    		}
    		for(int i = 0;i<g.nIncoming(u);i++){
				if(!g.edgeEnabled(g.incoming(u,i).id))
					continue;
				int v = g.incoming(u,i).node;
				int id = g.incoming(u,i).id;
				if(F[id] == 0){

				}else if(!seen[v]){
					Q.push_back(v);
					seen[v]=true;
				}
			}
    	}
    	//Now keep only the edges from a seen vertex to an unseen vertex
    	int i, j = 0;
    	for( i = 0;i<cut.size();i++){
    		if(!seen[cut[i].v] && seen[cut[i].u]){
    			cut[j++]=cut[i];
    		}
    	}
    	cut.resize(j);
#ifndef NDEBUG
		Weight dbg_sum = 0;
		for(int i = 0;i<cut.size();i++){
			int id = cut[i].id;
			assert(F[id]==capacity[id]);
			dbg_sum+=F[id];
		}
		assert(dbg_sum==f);
#endif
*/
    	return f;
    }
    const Weight getEdgeCapacity(int id){
     	assert(g.edgeEnabled(id));
     	return capacity[id];
     }
    const Weight getEdgeFlow(int flow_edge){

    	//we need to pick, possibly arbitrarily (but deterministically), which of the edges have flow
    	int arc_id = arc_map[flow_edge];
    	assert(arc_id>=0);
    	arc a = kt->get_arc(arc_id);
    	Weight start_cap = kt->get_ecap(a) ;
    	Weight end_cap = kt->get_rcap(a);
    	Weight remaining_flow = kt->get_ecap(a) - kt->get_rcap(a);
		int from = g.getEdge(flow_edge).from;
		int to = g.getEdge(flow_edge).to;

		for(int i = 0;i<g.nIncident(from,true);i++){
			int edgeid = g.incident(from,i,true).id;
			if(g.edgeEnabled(edgeid) &&  ((g.getEdge(edgeid).from==from && g.getEdge(edgeid).to==to))){
				assert(arc_map[edgeid] == arc_id);
				Weight & edge_cap = capacity[edgeid];
				if(edgeid==flow_edge){
					if(remaining_flow>=edge_cap)
						return edge_cap;
					else
						return remaining_flow;
				}
				if(remaining_flow>=edge_cap){
					remaining_flow-=edge_cap;
				}
				if(remaining_flow<=0){
					return 0;
				}
			}
		}

		//there may be multiple edges. Assign to the low-id's first?
		/*tmp_edges.clear();
		//we can probably arrange to avoid this, but for now, sorting by edge id to ensure deterministic results.
		for(int i = 0;i<g.nIncident(from,true);i++){
			int edgeid = g.incident(from,i,true).id;
			if(g.edgeEnabled(edgeid) &&  ((g.getEdge(edgeid).from==from && g.getEdge(edgeid).to==to) || (g.getEdge(edgeid).from==to && g.getEdge(edgeid).to==from))){
				assert(arc_map[edgeid]==arc_id);
				tmp_edges.push_back(edgeid);
			}
		}
		if(tmp_edges.size()<=1){
			if(tmp_edges.size()==0){
				return 0;
			}else{
				return remaining_flow;
			}
		}
		if(tmp_edges.size()>1){
			 std::sort (tmp_edges.begin(), tmp_edges.end());
		}

		for(int i = 0;i<tmp_edges.size();i++){
			int e = tmp_edges[i];
			Weight & edge_cap = capacity[e];
			if(e==edgeID){
				if(remaining_flow>=edge_cap)
					return edge_cap;
				else
					return remaining_flow;
			}
			if(remaining_flow>=edge_cap){
				remaining_flow-=edge_cap;
			}
			if(remaining_flow<=0){
				return 0;
			}
		}*/
    	return 0;
    }
    const Weight getEdgeResidualCapacity(int id){
    	return getEdgeCapacity(id)-getEdgeFlow(id);
    }
};
};
#endif

