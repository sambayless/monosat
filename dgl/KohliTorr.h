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
#include "DynamicGraph.h"
#include "MaxFlow.h"
#include <vector>
#include "alg/dyncut/graph.h"
#include "EdmondsKarpDynamic.h"
#include <algorithm>

namespace dgl {
template<class Capacity, typename Weight>
class KohliTorr: public MaxFlow<Weight> {
	Weight f = 0;
	/**
	 * Note: The Kohli Torr implementation does _not_ support multiple edges between the same nodes.
	 */
	kohli_torr::Graph<Weight, Weight, Weight> * kt = nullptr;
	struct LocalEdge {
		int from;
		int id;
		bool backward = false;
		LocalEdge(int from = -1, int id = -1, bool backward = false) :
				from(from), id(id), backward(backward) {
			
		}
	};
	int source = -1;
	int sink = -1;
	bool dynamic = true;
	Weight curflow;
	int last_modification;
	int last_deletion;
	int last_addition;
	std::vector<int> tmp_edges;
	int history_qhead;
	int last_history_clear;
	//bool backward_maxflow=false;
	bool kt_preserve_order = false;

	std::vector<std::vector<int>> multi_edges;

	typedef typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id arc;

	std::vector<int> arc_map; //map from edgeids to arcs (note: this may be many to one)
	std::vector<std::vector<int>> edge_map; //map from arc ids to edgeids
	DynamicGraph& g;
	Capacity & capacity;
	Weight INF;

	static const auto KT_SOURCE = kohli_torr::Graph<Weight, Weight, Weight>::SOURCE;
	static const auto KT_SINK = kohli_torr::Graph<Weight, Weight, Weight>::SINK;
#ifdef DEBUG_MAXFLOW
	EdmondsKarpDynamic<Capacity,Weight> ek;
#endif
	
	Weight max_capacity = 1; //start this at 1, not 0
	
	std::vector<bool> edge_enabled;
	std::vector<int> changed_edges;
	bool flow_needs_recalc = true;

public:
	double stats_calc_time = 0;
	double stats_flow_time = 0;
	long stats_flow_calcs = 0;
	KohliTorr(DynamicGraph& _g, Capacity & cap, int source, int sink, bool kt_preserve_order = false) :
			g(_g), capacity(cap), source(source), sink(sink), kt_preserve_order(kt_preserve_order), INF(0xF0F0F0)
#ifdef DEBUG_MAXFLOW
	,ek(_g,cap,source,sink)
#endif
	{
		curflow = 0;
		last_modification = -1;
		last_deletion = -1;
		last_addition = -1;
		
		history_qhead = 0;
		last_history_clear = -1;
	}
	
	int getSource() const {
		return source;
	}
	int getSink() const {
		return sink;
	}
	
	void setSource(int s) {
		if (source == s) {
			return;
		}
		if (kt) {
			
			if (dynamic) {
				kt->edit_tweights(source, 0, 0);
				//if(!backward_maxflow){
				kt->edit_tweights(s, max_capacity, 0);
				/*}else{
				 kt->edit_tweights(s,0,max_capacity);
				 }*/
			} else {
				kt->edit_tweights_wt(source, 0, 0);
				//if(!backward_maxflow){
				kt->edit_tweights_wt(s, max_capacity, 0);
				/*}else{
				 kt->edit_tweights_wt(s,0,max_capacity);
				 }*/
			}
		}
		source = s;
		last_modification = g.modifications - 1;
		flow_needs_recalc = true;
	}
	
	void setSink(int t) {
		if (sink == t) {
			return;
		}
		if (kt) {
			
			if (dynamic) {
				kt->edit_tweights(sink, 0, 0);
				//if(!backward_maxflow){
				kt->edit_tweights(t, 0, max_capacity);
				/*}else{
				 kt->edit_tweights(t,max_capacity,0);
				 }*/
			} else {
				kt->edit_tweights_wt(sink, 0, 0);
				//if(!backward_maxflow){
				kt->edit_tweights(t, 0, max_capacity);
				/*}else{
				 kt->edit_tweights(t,max_capacity,0);
				 }*/
			}
		}
		sink = t;
		last_modification = std::min(last_modification, g.modifications - 1);
		flow_needs_recalc = true;
	}
	
	void setCapacity(int u, int w, Weight c) {
		
	}
	void setAllEdgeCapacities(Weight c) {
		
	}
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	
	const Weight update() {
		int s = source;
		int t = sink;
		
		//see http://cstheory.stackexchange.com/a/10186
		static int it = 0;
		if (++it == 56) {
			int a = 1;
		}
#ifdef RECORD
		if (g.outfile) {
			fprintf(g.outfile, "f %d %d\n", s, t);
			fflush(g.outfile);
		}
#endif
		
		//C.resize(g.nodes());
#ifdef DEBUG_MAXFLOW
		for(int i = 0;i<g.all_edges.size();i++) {
			if(!g.hasEdge(i))
			continue;
			int id = g.all_edges[i].id;
			Weight cap = capacity[id];
			int from = g.all_edges[i].from;
			int to = g.all_edges[i].to;

			ek.setCapacity(from,to,cap);
		}
#endif
		if (last_modification > 0 && g.modifications == last_modification) {
#ifdef DEBUG_MAXFLOW
			Weight expected_flow =ek.maxFlow(source,sink);
#endif
			
#ifdef DEBUG_MAXFLOW
			assert(curflow==expected_flow);
			bassert(curflow == expected_flow);

#endif
			return curflow;
		} else if (!kt || last_modification <= 0 || kt->get_node_num() != g.nodes()
				|| edge_enabled.size() != g.edges()) {
			clearChangedEdges();
			edge_enabled.clear();
			flow_needs_recalc = true;
			if (kt) {
				delete (kt);
				kt = nullptr;
			}
			if (!kt) {
				kt = new kohli_torr::Graph<Weight, Weight, Weight>(g.nodes(), g.edges());
				kt->preserve_backward_order = kt_preserve_order;
				kt->maxflow(false); //just to initialize things.
			}
			
			while (kt->get_node_num() < g.nodes()) {
				int node_id = kt->add_node();
				assert(node_id == kt->get_node_num() - 1);
			}
			
			edge_enabled.resize(g.edges(), false);
			multi_edges.clear();
			multi_edges.resize(g.edges());
			arc_map.clear();
			arc_map.resize(g.edges(), -1);
			edge_map.clear();
			
			for (int edgeID = 0; edgeID < g.edges(); edgeID++) {
				
				if (!g.hasEdge(edgeID) || g.selfLoop(edgeID))
					continue;
				
				max_capacity += capacity[edgeID];
				int from = g.getEdge(edgeID).from;
				int to = g.getEdge(edgeID).to;
				edge_enabled[edgeID] = false;
				if (from == to)
					continue; //skip self edges.
				if (!kt->has_edge(from, to)) {
					assert(arc_map[edgeID] == -1);
					int arc_id = kt->add_edge(from, to, 0, 0);
					
					if (edge_map.size() <= arc_id + 1)
						edge_map.resize(arc_id + 2);
					
					arc_map[edgeID] = arc_id;
					//set the corresponding arc for each other from-to edge
					for (int i = 0; i < g.nIncident(from, true); i++) {
						int edgeid = g.incident(from, i, true).id;
						if (g.getEdge(edgeid).from == from && g.getEdge(edgeid).to == to) {
							edge_map[arc_id].push_back(edgeid);
							arc_map[edgeid] = arc_id;
						} else if (g.getEdge(edgeid).from == to && g.getEdge(edgeid).to == from) {
							edge_map[arc_id + 1].push_back(edgeid);
							arc_map[edgeid] = arc_id + 1; //the reverse arc is always stored right after the forward arc
						}
					}
				}
				if (g.edgeEnabled(edgeID)) {
					edge_enabled[edgeID] = true;
					//if(!backward_maxflow){
					kt->edit_edge_inc(from, to, capacity[edgeID], 0);
					/*}else{
					 kt->edit_edge_inc(to,from,capacity[edgeID],0);
					 }*/
				}
			}
			if (dynamic) {
				//if(!backward_maxflow){
				kt->edit_tweights(s, max_capacity, 0);
				kt->edit_tweights(t, 0, max_capacity);
				/*	}else{
				 kt->edit_tweights(t,max_capacity,0);
				 kt->edit_tweights(s,0,max_capacity);
				 }*/
			} else {
				//if(!backward_maxflow){
				kt->edit_tweights_wt(s, max_capacity, 0);
				kt->edit_tweights_wt(t, 0, max_capacity);
				/*	}else{
				 kt->edit_tweights_wt(t,max_capacity,0);
				 kt->edit_tweights_wt(s,0,max_capacity);
				 }*/
			}
		} else if (g.historyclears != last_history_clear || g.changed()) {
			flow_needs_recalc = true;
			for (int edgeid = 0; edgeid < g.edges(); edgeid++) {
				if (!g.hasEdge(edgeid) || g.selfLoop(edgeid))
					continue;
				if (g.edgeEnabled(edgeid) && !edge_enabled[edgeid]) {
					
					edge_enabled[edgeid] = true;
					//if(!backward_maxflow){
					kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, capacity[edgeid], 0);
					/*}else{
					 kt->edit_edge_inc(g.getEdge(edgeid).to,g.getEdge(edgeid).from,capacity[edgeid],0);
					 }*/
				} else if (!g.edgeEnabled(edgeid) && edge_enabled[edgeid]) {
					assert(edge_enabled[edgeid]);
					edge_enabled[edgeid] = false;
					//if(!backward_maxflow){
					kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, -capacity[edgeid], 0);
					/*	}else{
					 kt->edit_edge_inc(g.getEdge(edgeid).to,g.getEdge(edgeid).from,-capacity[edgeid],0);
					 }*/
				}
			}
			history_qhead = g.history.size();
		}
		flow_needs_recalc = true;
		assert(kt);
		
		for (int i = history_qhead; i < g.history.size(); i++) {
			int edgeid = g.history[i].id;
			if (g.selfLoop(edgeid))
				continue; //skip self loops
			if (g.history[i].addition && g.edgeEnabled(edgeid) && !edge_enabled[edgeid]) {
				
				edge_enabled[edgeid] = true;
				//if(!backward_maxflow){
				kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, capacity[edgeid], 0);
				/*}else{
				 kt->edit_edge_inc(g.getEdge(edgeid).to,g.getEdge(edgeid).from,capacity[edgeid],0);
				 }*/
			} else if (!g.history[i].addition && !g.edgeEnabled(edgeid) && edge_enabled[edgeid]) {
				assert(edge_enabled[edgeid]);
				edge_enabled[edgeid] = false;
				//if(!backward_maxflow){
				kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, -capacity[edgeid], 0);
				/*}else{
				 kt->edit_edge_inc(g.getEdge(edgeid).to,g.getEdge(edgeid).from,-capacity[edgeid],0);
				 }*/
			}
		}
		
		f = kt->maxflow(dynamic);
		
#ifdef DEBUG_MAXFLOW
		Weight expected_flow =ek.maxFlow(source,sink);
		bassert(f == expected_flow);

#endif
		//dbg_print_graph(s,t,true);
		
#ifndef NDEBUG
		
		for (int i = 0; i < g.edges(); i++) {
			if (g.getEdge(i).from == g.getEdge(i).to)
				continue; //skip self edges.
			assert(edge_enabled[i] == g.edgeEnabled(i));
		}
		
#endif
		dbg_check_flow(s, t);
		curflow = f;
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.history.size();
		last_history_clear = g.historyclears;
		return f;
	}
	
	std::vector<int> & getChangedEdges() {
		calc_flow();
		while (kt->changed_edges.size()) {
			int arc = kt->changed_edges.back();
			kt->changed_edges.pop_back();
			kt->unmarkFlowEdge(kt->get_arc(arc));
			assert(arc < edge_map.size());
			for (int edgeID : edge_map[arc]) {
				changed_edges.push_back(edgeID);
			}
		}
		return changed_edges;
	}
	
	void clearChangedEdges() {
		changed_edges.clear();
	}
private:
	
	inline void collect_multi_edges(int for_edge) {
		if (multi_edges[for_edge].size() == 0) {
			int from = g.getEdge(for_edge).from;
			int to = g.getEdge(for_edge).to;
			//for(int i = g.nIncident(from,true)-1;i>=0;i--){
			for (int i = 0; i < g.nIncident(from, true); i++) {
				int edgeid = g.incident(from, i, true).id;
				if ((g.getEdge(edgeid).from == from && g.getEdge(edgeid).to == to)) {
					multi_edges[for_edge].push_back(edgeid);
				}
			}
		}
	}
	
	void dbg_print_graph(int from, int to, bool only_flow = false) {
#ifndef NDEBUG
		
		if (edge_enabled.size() < g.edges())
			return;
		static int it = 0;
		if (++it == 6) {
			int a = 1;
		}
		printf("Graph %d\n", it);
		printf("digraph{\n");
		for (int i = 0; i < g.nodes(); i++) {
			if (i == from) {
				printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
			} else if (i == to) {
				printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
			}
			if (kt->what_segment(i, kohli_torr::Graph<Weight, Weight, Weight>::SOURCE)
					== kohli_torr::Graph<Weight, Weight, Weight>::SOURCE) {
				printf("n%d [label=\"%d\", style=filled, fillcolor=blue]\n", i, i);
			} else if (kt->what_segment(i, kohli_torr::Graph<Weight, Weight, Weight>::SOURCE)
					== kohli_torr::Graph<Weight, Weight, Weight>::SINK) {
				printf("n%d [label=\"%d\", style=filled, fillcolor=red]\n", i, i);
			} else
				printf("n%d\n", i);
		}
		printf("outer_source\n");
		printf("outer_sink\n");
		for (int i = 0; i < g.edges(); i++) {
			if (edge_enabled[i]) {
				auto & e = g.all_edges[i];
				const char * s = "black";
				if (only_flow && getEdgeFlow(e.id) == 0)
					continue;
				std::cout << "n" << e.from << " -> n" << e.to << " [label=\"" << i << ": " << getEdgeFlow(e.id) << "/"
						<< capacity[i] << "\" color=\"" << s << "\"]\n";
				//printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"%s\"]\n", e.from,e.to, i, F[i],capacity[i] , s);
			}
		}
		for (int i = 0; i < g.nodes(); i++) {
			Weight w = kt->get_trcap(i);
			//Weight c = kt->t_cap(i);
			/**
			 * 	// if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
			 // otherwise         -tr_cap is residual capacity of the arc node->SINK
			 */
			/*if(w<0){
			 std::cout<<"n" << from <<" -> n" << i << " [label=\"" << w <<  "\" color=\"blue" <<"\"]\n";
			 }else if(w>0){
			 std::cout<<"n" << i <<" -> n" << to << " [label=\""<< w  <<  "\" color=\"red" << "\"]\n";
			 }*/
			if (w < 0) {
				std::cout << "outer_source" << " -> n" << i << " [label=\"" << w << "\" color=\"blue" << "\"]\n";
			} else if (w > 0) {
				std::cout << "n" << i << " -> outer_sink " << " [label=\"" << w << "\" color=\"red" << "\"]\n";
			}
		}
		printf("}\n");
#endif
	}
	void bassert(bool c) {
#ifndef NDEBUG
		assert(c);
		if (!c) {
			exit(4);
		}
#endif
	}
	
	/*  void dbg_print_matrix(){
	 int width = 11;
	 int height = 11;
	 int node = 0;
	 for(int y = 0;y<height;y++){
	 for(int x = 0;x<width;x++){
	 for (int i = 0;i<g.nIncident(node);i++){
	 auto & e = g.incident(node,i);
	 int flow = this->getEdgeFlow(e.id);
	 if(flow!=0){
	 printf("%d ", flow);
	 }else{
	 printf(" ")
	 }

	 }

	 node = node+1;
	 }
	 }


	 }*/

	void dbg_check_flow(int s, int t) {
#ifndef NDEBUG
		//check that the flow is legal
		for (int i = 0; i < g.edges(); i++) {
			Weight flow = getEdgeFlow(i);
			bassert(flow >= 0);
			bassert(flow <= capacity[i]);
			if (flow != 0) {
				bassert(g.edgeEnabled(i));
			}
		}
		
		for (int n = 0; n < g.nodes(); n++) {
			Weight flow_in = 0;
			Weight flow_out = 0;
			for (int i = 0; i < g.nIncident(n); i++) {
				int edge = g.incident(n, i).id;
				if (g.edgeEnabled(edge)) {
					Weight flow = getEdgeFlow(edge);
					bassert(flow >= 0);
					flow_out += flow;
				}
			}
			for (int i = 0; i < g.nIncoming(n); i++) {
				int edge = g.incoming(n, i).id;
				if (g.edgeEnabled(edge)) {
					Weight flow = getEdgeFlow(edge);
					bassert(flow >= 0);
					flow_in += flow;
				}
			}
			if (n == s) {
				//if(!backward_maxflow){
				bassert(flow_in == 0);
				bassert(flow_out == f);
				/*}else{
				 bassert(flow_in==0);
				 bassert(flow_out==f);
				 }*/
			} else if (n == t) {
				//if(!backward_maxflow){
				bassert(flow_out == 0);
				bassert(flow_in == f);
				/*}else{
				 bassert(flow_out==0);
				 bassert(flow_in==f);
				 }*/
			} else {
				bassert(flow_in == flow_out);
			}
			
		}
		
#endif
	}
	
	std::vector<bool> seen;
	std::vector<bool> visited;

	inline void calc_flow() {
		if (!flow_needs_recalc)
			return;
		//double startflowtime = Monosat::rtime(0);
		flow_needs_recalc = false;
		stats_flow_calcs++;
		//apply edmonds karp to the current flow.
		Weight maxflow = kt->maxflow(true, nullptr);
		double startcalctime = Monosat::rtime(0);
		//if(backward_maxflow){
		
		/*kt->clear_t_edges(sink,source);

		 dbg_check_flow(source,sink);*/
		/*}else{*/

		kt->clear_t_edges(source, sink);
		
		dbg_check_flow(source, sink);
		//}
		
		assert(kt->maxflow(true, nullptr) == maxflow);
		//stats_calc_time+=  Monosat::rtime(0)-startcalctime;
		//stats_flow_time+=  startcalctime-startflowtime;
		//printf("flow calc time %f %f\n", stats_flow_time,stats_calc_time);
	}
	
	std::vector<int> Q;

public:
	
	const Weight minCut(std::vector<MaxFlowEdge> & cut) {
		Weight f = this->maxFlow();
		
		int s = source;
		int t = sink;
#ifdef RECORD
		if (g.outfile) {
			fprintf(g.outfile, "m %d %d\n", s, t);
			fflush(g.outfile);
		}
#endif
		cut.clear();
		dbg_print_graph(s, t);
		/*   	if(f==0)
		 return 0;*/

		auto SOURCE = kohli_torr::Graph<Weight, Weight, Weight>::SOURCE; // backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SINK : kohli_torr::Graph<Weight,Weight,Weight>::SOURCE;
		auto SINK = kohli_torr::Graph<Weight, Weight, Weight>::SINK; //backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SOURCE :  kohli_torr::Graph<Weight,Weight,Weight>::SINK;
		
		//KT allows for all the nodes to be source or sink, if the cut is placed at the final source->'outer source' or 'outer sink' -> sink edge.
		//Ideally, this shouldn't happen here, because those should have infinite weight edges.
		if (kt->what_segment(s, SOURCE) == SINK) {
			exit(1);
			return -1;
		} else if (kt->what_segment(t, SOURCE) == SOURCE) {
			exit(1);
			return -1;
		}
		
		for (int n = 0; n < g.nodes(); n++) {
			auto t = kt->what_segment(n, SOURCE);
			if (t == SINK) {
				//check to see if any neighbouring edges are in the source segment.
				for (int i = 0; i < g.nIncoming(n); i++) {
					//any edge that crosses from the source side into the sink side is on the cut.
					auto & e = g.incoming(n, i);
					if (g.edgeEnabled(e.id)) {
						auto segment = kt->what_segment(e.node, SOURCE);
						if (segment == SOURCE) {
							//then this edge is on the cut
							cut.push_back(MaxFlowEdge { n, e.node, e.id });
						}
					}
				}
			}
		}
		
		//kohli-torr is unusual in that it maintains the _mincut_, explicitly.
		
		/*calc_flow();
		 cut.clear();
		 Q.clear();
		 Q.push_back(s);
		 seen.clear();
		 seen.resize(g.nodes());
		 seen[s]=true;
		 //
		 //explore the residual graph
		 for(int j = 0;j<Q.size();j++){
		 int u = Q[j];

		 for(int i = 0;i<g.nIncident(u);i++){
		 if(!g.edgeEnabled(g.incident(u,i).id))
		 continue;
		 int v = g.incident(u,i).node;
		 int id = g.incident(u,i).id;
		 if(getEdgeResidualCapacity(id) == 0){
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
		 if(getEdgeFlow(id) == 0){

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
		 */
#ifndef NDEBUG
		Weight dbg_sum = 0;
		for (int i = 0; i < cut.size(); i++) {
			int id = cut[i].id;
			bassert(getEdgeFlow(id) == capacity[id]);
			dbg_sum += getEdgeFlow(id);
		}
		bassert(dbg_sum == f);
#endif
		
		return f;
	}
	const Weight getEdgeCapacity(int id) {
		assert(g.edgeEnabled(id));
		return capacity[id];
	}
	
public:
	const Weight getEdgeFlow(int flow_edge) {
		if (g.getEdge(flow_edge).from == g.getEdge(flow_edge).to)
			return 0;    	//self edges have no flow.
			
		calc_flow();
		collect_multi_edges(flow_edge);
		//we need to pick, possibly arbitrarily (but deterministically), which of the edges have flow
		int arc_id = arc_map[flow_edge];
		assert(arc_id >= 0);
		arc a;
		//if(!backward_maxflow){
		a = kt->get_arc(arc_id);
		/*	}else{
		 a= kt->get_reverse( kt->get_arc(arc_id));
		 }*/
		Weight start_cap = kt->get_ecap(a);
		Weight end_cap = kt->get_rcap(a);
		Weight remaining_flow = kt->get_ecap(a) - kt->get_rcap(a);
		if (remaining_flow <= 0)
			return 0;
		for (int edgeid : multi_edges[flow_edge]) {
			
			assert(g.getEdge(edgeid).from == g.getEdge(flow_edge).from);
			assert(g.getEdge(edgeid).to == g.getEdge(flow_edge).to);
			
			if (g.edgeEnabled(edgeid)) {
				assert(arc_map[edgeid] == arc_id);
				Weight edge_cap = capacity[edgeid];
				if (edgeid == flow_edge) {
					if (remaining_flow >= edge_cap)
						return edge_cap;
					else
						return remaining_flow;
				}
				if (remaining_flow >= edge_cap) {
					remaining_flow -= edge_cap;
				}
				if (remaining_flow <= 0) {
					return 0;
				}
			}
		}
		
		return 0;
	}
	const Weight getEdgeResidualCapacity(int id) {
		return getEdgeCapacity(id) - getEdgeFlow(id);
	}
};
}
;
#endif

