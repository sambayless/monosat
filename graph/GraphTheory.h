/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef DGRAPH_H_
#define DGRAPH_H_

#include "core/Theory.h"
#include "Graph.h"
#include "Reach.h"
#include "Dijkstra.h"
#include "Connectivity.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "MaxFlow.h"
#include "IBFS.h"
#include "EdmondsKarp.h"
#include "utils/System.h"
#ifdef DEBUG_GRAPH
#include "TestGraph.h"
#endif

#ifdef DEBUG_SOLVER
#include "TestGraph.h"
#endif
namespace Minisat{



class GraphTheorySolver:public GraphTheory{
private:


	Lit False;
	Lit True;
	int local_q;

	Solver * S;
#ifdef DEBUG_GRAPH
	Solver * dbg;
	TestGraph *dbg_graph;
#endif
#ifdef DEBUG_SOLVER
	TestGraph * shadow_dbg;
#endif


	vec<lbool> edge_assignments;

	typedef DefaultEdgeStatus PositiveEdgeStatus;
	typedef DefaultEdgeStatus NegativeEdgeStatus;
	typedef DefaultEdgeStatus CutEdgeStatus;

	/*struct PositiveEdgeStatus{

		const GraphTheorySolver & outer;

	    bool operator [] (int index) const {
	    	return outer.edge_assignments[index]==l_True;
	    }
	    void setStatus(int id,bool value){
	    	//Do nothing
	    }

	    int size(){
	    	return outer.edge_assignments.size();
	    }
	    void growTo(int size) const{
	    	assert(outer.edge_assignments.size()>=size);
	    }
	    PositiveEdgeStatus(GraphTheorySolver & _outer):outer(_outer){

	    }
	};*/

/*	struct NegativeEdgeStatus{

		const GraphTheorySolver & outer;
		bool ignore;
	    const bool operator [] (int index) const { return outer.edge_assignments[index]!=l_False; }
	    bool&       operator [] (int index)       { return ignore; }
	    int size(){
	    	return outer.edge_assignments.size();
	    }
	    void growTo(int size) const{
	    	assert(outer.edge_assignments.size()>=size);
	    }
	    NegativeEdgeStatus(GraphTheorySolver & _outer):outer(_outer), ignore(false){

	    }
	};*/

	PositiveEdgeStatus g_status;
	NegativeEdgeStatus antig_status;
	CutEdgeStatus cutGraph_status;

	DynamicGraph<PositiveEdgeStatus> g;
	DynamicGraph<NegativeEdgeStatus> antig;
	DynamicGraph<CutEdgeStatus> cutGraph;



	struct ReachDetector{
		int within;
		int source;
		Reach * positive_reach_detector;
		Reach * negative_reach_detector;
		Dijkstra<PositiveEdgeStatus>* positive_dist_detector;
		Dijkstra<NegativeEdgeStatus>* negative_dist_detector;
		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;


		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

		Lit getLit(int node){

			return reach_lits[node];
		}

		ReachDetector():within(-1),source(-1),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_dist_detector(NULL),negative_dist_detector(NULL){}
	};

	struct ReachInfo{
		int source;
		ReachDetector * detector;
		ReachInfo():source(-1),detector(NULL){}
	};

	vec<ReachInfo> reach_info;
public:
	vec<ReachDetector*> reach_detectors;

	vec<CRef> reach_markers;
	vec<CRef> non_reach_markers;

	vec<int> marker_map;

	struct Edge{
		Var v;
		int from;
		int to;
	};
	vec<MaxFlow::Edge> cut;

	//Full matrix
	vec<vec<Edge> > edges;

	//Just a list of the edges
	vec<Edge> edge_list;

	vec<vec<Edge> > inv_adj;

	Var min_edge_var;
	int num_edges;
	struct Assignment{
		bool isEdge:1;
		bool assign:1;
		int from:30;
		int to;
		Var var;
	};

	vec<Assignment> trail;
	vec<int> trail_lim;
	MaxFlow * mc;

    vec<char> seen;
	vec<int> to_visit;


public:

	double mctime;
	double reachtime;
	double unreachtime;
	double pathtime;
	double propagationtime;
	double reachupdatetime;
	double unreachupdatetime;
	double stats_initial_propagation_time;
	int stats_mc_calls;
	GraphTheorySolver(Solver * S_):S(S_),g(g_status),antig(antig_status) ,cutGraph(cutGraph_status){
		True = mkLit(S->newVar(),false);
			False=~True;
			S->addClause(True);
			num_edges=0;
			local_q=0;
			mctime=0;
			stats_mc_calls=0;
			reachtime=0;
			unreachtime=0;
			pathtime=0;
			propagationtime=0;
			reachupdatetime=0;
			unreachupdatetime=0;
			stats_initial_propagation_time=0;
			if(mincutalg==ALG_IBFS){
				mc = new IBFS<CutEdgeStatus>(cutGraph);
			}else{
				mc = new EdmondsKarp<CutEdgeStatus>(cutGraph);
			}

#ifdef DEBUG_GRAPH
		dbg=new Solver();
		dbg_graph = new TestGraph(dbg);
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg = new TestGraph(S->dbg_solver);
#endif
	}

	void printStats(){

		printf("Graph stats:\n");
		printf("Prop Time: %f (initial: %f)\n", propagationtime,stats_initial_propagation_time);
		printf("Reach Time: %f (update Time: %f)\n", reachtime,reachupdatetime);
		printf("Unreach Time: %f (Update Time: %f)\n", unreachtime,unreachupdatetime);
		printf("Path Time: %f\n", pathtime);
		printf("Min-cut Time: %f (%d calls, %f average)\n", mctime, stats_mc_calls,(mctime/(stats_mc_calls ? stats_mc_calls:1)));

		int stats_full_updates=0;
		int stats_fast_updates=0;
		int stats_skip_deletes=0;
		int stats_skipped_updates=0;
		int skipable_deletions = 0;
		double stats_full_update_time=0;
		double stats_fast_update_time=0;

		for(int i = 0;i<reach_detectors.size();i++){
			skipable_deletions+=reach_detectors[i]->positive_reach_detector->stats_num_skipable_deletions;
			skipable_deletions+=reach_detectors[i]->negative_reach_detector->stats_num_skipable_deletions;

			stats_full_updates+=reach_detectors[i]->positive_reach_detector->stats_full_updates;
			stats_fast_updates+=reach_detectors[i]->positive_reach_detector->stats_fast_updates;
			stats_skip_deletes+=reach_detectors[i]->positive_reach_detector->stats_skip_deletes;
			stats_skipped_updates+=reach_detectors[i]->positive_reach_detector->stats_skipped_updates;

			stats_full_update_time+=reach_detectors[i]->positive_reach_detector->stats_full_update_time;
			stats_fast_update_time+=reach_detectors[i]->positive_reach_detector->stats_fast_update_time;

			stats_full_updates+=reach_detectors[i]->negative_reach_detector->stats_full_updates;
			stats_fast_updates+=reach_detectors[i]->negative_reach_detector->stats_fast_updates;
			stats_skip_deletes+=reach_detectors[i]->negative_reach_detector->stats_skip_deletes;
			stats_skipped_updates+=reach_detectors[i]->negative_reach_detector->stats_skipped_updates;

			stats_full_update_time+=reach_detectors[i]->negative_reach_detector->stats_full_update_time;
			stats_fast_update_time+=reach_detectors[i]->negative_reach_detector->stats_fast_update_time;
		}

		printf("Dijkstra Full Updates: %d (time: %f, average: %f)\n",stats_full_updates, stats_full_update_time,(stats_full_update_time/(stats_full_updates ? stats_full_updates:1)));
		printf("Dijkstra Fast Updates: %d (time: %f, average: %f)\n",stats_fast_updates, stats_fast_update_time,(stats_fast_update_time/(stats_fast_updates ? stats_fast_updates:1)));
		printf("Dijkstra Skipped Updates: %d (deletionSkips: %d, skipable deletes %d)\n",stats_skipped_updates, stats_skip_deletes,skipable_deletions);

	}

     ~GraphTheorySolver(){};
	 int newNode(){
#ifdef DEBUG_GRAPH
		 dbg_graph->newNode();
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg->newNode();
#endif
		 edges.push();
		 for(int i = 0;i<edges.size();i++)
			 edges[i].growTo(edges.size());
		 inv_adj.push();
		 reach_info.push();
		 antig.addNode();
		 cutGraph.addNode();


		 seen.growTo(nNodes());

		return g.addNode();
	}
	 void newNodes(int n){
		 for(int i = 0;i<n;i++)
			 newNode();
	 }
	int nNodes(){
		return g.nodes;
	}
	bool isNode(int n){
		return n>=0 && n<nNodes();
	}

#ifdef DEBUG_GRAPH
		bool dbg_clause(const vec<Lit> &conflict){

			static vec<Lit> c;
			c.clear();
			for(int i = 0;i<conflict.size();i++)
				c.push(~conflict[i]);
		//	assert(dbg->solve());
			bool res = dbg->solve(c);
			assert(~res);

			return true;
		}
		bool dbg_propgation(Lit l){

			static vec<Lit> c;
			c.clear();
			for(int i = 0;i<S->trail.size() ;i++){
				Lit l = S->trail[i];
				Var v = var(l);
				if(v>=min_edge_var && v<min_edge_var+num_edges){
					if(edge_list[v-min_edge_var].v<0)
								continue;
					Edge e = edge_list[v-min_edge_var];
					c.push(l);
				}
			}
			c.push(~l);
			bool res = dbg->solve(c);
			assert(~res);

			return true;
		}
#endif
	void	dbg_sync_reachability(){
#ifdef DEBUG_GRAPH

		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			for(int j = 0;j< d->reach_lits.size();j++){
				Lit l = d->reach_lits[j];
				if(l!=lit_Undef){
					int node = d->getNode(var(l));

					if(d->positive_reach_detector->connected(node)){
						assert(S->value(l)==l_True);
					}else if(! d->negative_reach_detector->connected(node)){
						assert(S->value(l)==l_False);
					}
				}

			}
		}
#endif
	}
	void dbg_sync(){
#ifdef DEBUG_GRAPH
		static vec<lbool> assigned;
		assigned.clear();
		for(int i = 0;i<edge_list.size();i++)
			assigned.push(l_Undef);
		int lev = 0;
		int j =0;
		for(int i = 0;i<trail.size();i++){
			while(j<trail_lim.size() && i>=trail_lim[j]){
				lev++;
				j++;
			}

			Assignment e = trail[i];
			if(e.isEdge){

				Lit l = mkLit( e.var,!e.assign);
				assert(S->value(l)==l_True);
				int expected_level = S->level(var(l));
				assert(S->level(var(l))==lev);
				int edge_num = e.var-min_edge_var;
				if(edge_list[edge_num].v<0)
							continue;
				assert(assigned[edge_num]==l_Undef);

				assigned[edge_num] = sign(l)?l_False:l_True;
			}else{
				//Reachability assignment


			}
		}

		for(int i = 0;i<assigned.size();i++){
			assert(edge_assignments[i]== assigned[i]);
		}

		for(int i = 0;i<S->trail.size();i++){

			Lit l = S->trail[i];
			Var v = var(l);

			int lev = S->level(v);

			if(v>= min_edge_var && v<min_edge_var+num_edges){
				int edge_num = v-min_edge_var;
				if(edge_list[edge_num].v<0)
					continue;
				lbool assigned_val=assigned[edge_num];
				assert(assigned_val== (sign(l)?l_False:l_True));
			}

		}



#endif
	}

	void backtrackUntil(int level){
		static int it = 0;
		if(++it==28){
			int a =1;
		}
		//need to remove and add edges in the two graphs accordingly.
		if(trail_lim.size()>level){
			int stop = trail_lim[level];
			for(int i = trail.size()-1;i>=trail_lim[level];i--){
				Assignment e = trail[i];
				if(e.isEdge){
					assert(S->value(e.var)==l_Undef);
					int edge_num = e.var-min_edge_var;
					assert(edge_assignments[edge_num]!=l_Undef);
					edge_assignments[edge_num]=l_Undef;
					if(e.assign){
						g.disableEdge(e.from,e.to, edge_num);
					}else{
						antig.enableEdge(e.from,e.to,edge_num);
						assert(antig.hasEdge(e.from,e.to));
					}
				}
			}
			trail.shrink(trail.size()-stop);
			trail_lim.shrink(trail_lim.size()-level);
			assert(trail_lim.size()==level);


		}

		if(local_q>S->qhead)
			local_q=S->qhead;
		assert(dbg_graphsUpToDate());
		/*for(int i = 0;i<reach_detectors.size();i++){
			if(reach_detectors[i]->positive_reach_detector)
				reach_detectors[i]->positive_reach_detector->update();
			if(reach_detectors[i]->negative_reach_detector)
				reach_detectors[i]->negative_reach_detector->update();
		}*/

		dbg_sync();


	};

	void backtrackUntil(Lit p){
			//need to remove and add edges in the two graphs accordingly.
			int i = trail.size()-1;
			for(;i>=0;i--){
				Assignment e = trail[i];
				if(e.isEdge){
					int edge_num = e.var-min_edge_var;
					assert(edge_assignments[edge_num]!=l_Undef);
					edge_assignments[edge_num]=l_Undef;
					if(e.assign){
						g.disableEdge(e.from,e.to, edge_num);
					}else{
						antig.enableEdge(e.from,e.to,edge_num);
						assert(antig.hasEdge(e.from,e.to));
					}
				}else{
					if(var(p)==e.var){
						assert(sign(p)!=e.assign);
						break;
					}
				}
			}

			trail.shrink(trail.size()-(i+1));
			//while(trail_lim.size() && trail_lim.last()>=trail.size())
			//	trail_lim.pop();

	/*		for(int i = 0;i<reach_detectors.size();i++){
					if(reach_detectors[i]->positive_reach_detector)
						reach_detectors[i]->positive_reach_detector->update();
					if(reach_detectors[i]->negative_reach_detector)
						reach_detectors[i]->negative_reach_detector->update();
				}*/
		};

	void newDecisionLevel(){
		trail_lim.push(trail.size());
	};

	void buildReason(Lit p, vec<Lit> & reason){
		CRef marker = S->reason(var(p));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef- marker;
		int d = marker_map[pos];

		backtrackUntil(p);

		reason.push(p);
		assert(d!=0);
		if(d>0){
		//	double startpathtime = cpuTime();
			d--;
			/*Dijkstra & detector = *reach_detectors[d]->positive_dist_detector;
			//the reason is a path from s to p(provided by d)
			//p is the var for a reachability detector in dijkstra, and corresponds to a node
			detector.update();
			Var v = var(p);
			int u =reach_detectors[d]->getNode(v); //reach_detectors[d]->reach_lit_map[v];
			assert(detector.connected(u));
			int w;
			while(( w= detector.previous(u)) > -1){
				Lit l = mkLit( edges[w][u].v,true );
				assert(S->value(l)==l_False);
				reason.push(l);
				u=w;
			}*/
			Var v = var(p);
			int u =reach_detectors[d]->getNode(v);
			buildReachReason(v,opt_conflict_shortest_path ? *reach_detectors[d]->positive_dist_detector: *reach_detectors[d]->positive_reach_detector,reason);

#ifdef DEBUG_GRAPH
		 assert(dbg_clause(reason));

#endif
			//double elapsed = cpuTime()-startpathtime;
		//	pathtime+=elapsed;
		}else{
			d=-d-1;

			//the reason is a cut separating p from s;
			//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

			//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
			//assign the mincut edge weights if they aren't already assigned.

			//set weights

			//compute the mincut
			/*Var v = var(p);
				int t =  v- var(reach_lits[d][0]);
				assert(!detector.connected(t));

			cut.clear();
			mc.minCut (detector.getSource(), t,cut);
			for(int i = 0;i<cut.size();i++){
				reason.push(mkLit( edges[cut[i].u][cut[i].v].v,false ));
			}*/
			Var v = var(p);
			int t = reach_detectors[d]->getNode(v); // v- var(reach_lits[d][0]);
			buildNonReachReason(t,d,reason);
			/*for(int i = 0;i<g.nodes;i++){
				for(int j = 0;j<g.adjacency[i];j++){
					int v = g.adjacency[i][j];
					if(mc.parity[i]!=mc.parity[v]){
						//this edge is in the min-cut
						reason.push(mkLit( edges[i][v].v,true ));
					}
				}
			}*/
		}

	}
	void buildReachReason(int node,Reach & d,vec<Lit> & conflict){
		//drawFull();
		assert(dbg_reachable(d.getSource(),node));

		double starttime = cpuTime();
		d.update();
		int u = node;
		int p;
		while(( p = d.previous(u)) != -1){
			Edge edg = edges[p][u];
			Var e =edges[p][u].v;
			lbool val = S->value(e);
			assert(S->value(e)==l_True);
			conflict.push(mkLit(e, true));
			u = p;
		}
		double elapsed = cpuTime()-starttime;
		pathtime+=elapsed;
#ifdef DEBUG_GRAPH
		 assert(dbg_clause(conflict));

#endif
	}


	bool dbg_reachable(int from, int to){
#ifdef DEBUG_GRAPH
		DefaultEdgeStatus tmp;
		DynamicGraph<> g(tmp);
		for(int i = 0;i<nNodes();i++){
			g.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
				continue;
			Edge e  = edge_list[i];
			if(S->assigns[e.v]==l_True){
				g.addEdge(e.from,e.to);
			}
		}

		Dijkstra<> d(from,g);

		return d.connected(to);
#else
		return true;
#endif
	}

	bool dbg_notreachable(int from, int to){

#ifdef DEBUG_GRAPH
		//drawFull(from,to);
		DefaultEdgeStatus tmp;
		DynamicGraph<> g(tmp);
		for(int i = 0;i<nNodes();i++){
			g.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge e  = edge_list[i];
			if(S->assigns[e.v]!=l_False){
				g.addEdge(e.from,e.to);
			}
		}

		Dijkstra<> d(from,g);

		return !d.connected(to);
#else
		return true;
#endif
	}

	bool dbg_graphsUpToDate(){
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
				continue;
			Edge e = edge_list[i];
			lbool val = S->value(e.v);
			assert(edge_assignments[i]==val);
			if(val==l_True || val==l_Undef){
				assert(antig.edgeEnabled(i));
				assert(antig.hasEdge(e.from,e.to));
			}else{
				assert(!antig.edgeEnabled(i));
				assert(!antig.hasEdge(e.from,e.to));
			}
			if(val==l_True){
				assert(g.edgeEnabled(i));
				assert(g.hasEdge(e.from,e.to));
			}else{
				assert(!g.edgeEnabled(i));
				assert(!g.hasEdge(e.from,e.to));
			}
		}

#endif
		return true;
	}

	void buildNonReachReason(int node, int detector ,vec<Lit> & conflict){
		static int it = 0;
		++it;
		int u = node;
		//drawFull( non_reach_detectors[detector]->getSource(),u);
		assert(dbg_notreachable( reach_detectors[detector]->source,u));
		double starttime = cpuTime();
		cutGraph.clearHistory();
		stats_mc_calls++;
		if(opt_conflict_min_cut){
			//ok, set the weights for each edge in the cut graph.
			//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
			for(int u = 0;u<cutGraph.adjacency.size();u++){
				for(int j = 0;j<cutGraph.adjacency[u].size();j++){
					int v = cutGraph.adjacency[u][j].to;
					Var var = edges[u][v].v;
					/*if(S->value(var)==l_False){
						mc.setCapacity(u,v,1);
					}else{*/
						mc->setCapacity(u,v,0xF0F0F0);
					//}
				}
			}

			//find any edges assigned to false, and set their capacity to 1
			for(int i =0;i<trail.size();i++){
				if(trail[i].isEdge && !trail[i].assign){
					mc->setCapacity(trail[i].from, trail[i].to,1);
				}
			}

			cut.clear();

			int f =mc->minCut(reach_detectors[detector]->source,node,cut);
			assert(f<0xF0F0F0); assert(f==cut.size());//because edges are only ever infinity or 1
			for(int i = 0;i<cut.size();i++){
				MaxFlow::Edge e = cut[i];

				Lit l = mkLit( edges[e.u][e.v].v,false);
				assert(S->value(l)==l_False);
				conflict.push(l);
			}
		}else{
			//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
			//or perhaps we can learn the actual 1-uip cut?

			    to_visit.clear();
			    to_visit.push(node);
			    seen.clear();
				seen.growTo(nNodes());
			    seen[node]=true;

			    do{

			    	assert(to_visit.size());
			    	int u = to_visit.last();
			    	assert(u!=reach_detectors[detector]->source);
			    	to_visit.pop();
			    	assert(seen[u]);
			    	assert(!reach_detectors[detector]->negative_reach_detector->connected_unsafe(u));
			    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
			    	for(int i = 0;i<inv_adj[u].size();i++){
			    		int v = inv_adj[u][i].v;
			    		int from = inv_adj[u][i].from;
			    		assert(from!=u);
			    		assert(inv_adj[u][i].to==u);
			    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
			    		int edge_num = v-min_edge_var;

			    		if(edge_assignments[edge_num]==l_False){
			    			//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
			    			//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
			    			//if(!seen[from])
			    				conflict.push(mkLit(v,false));
			    		}else{
			    			assert(from!=reach_detectors[detector]->source);
			    			//even if it is undef? probably...
			    			if(!seen[from]){
			    				seen[from]=true;
			    				to_visit.push(from);
			    			}
			    		}
			    	}
			    }  while (to_visit.size());


		}
		double elapsed = cpuTime()-starttime;
			mctime+=elapsed;

#ifdef DEBUG_GRAPH
		 assert(dbg_clause(conflict));
#endif
	}
	bool propagateTheory(vec<Lit> & conflict){
		static int itp = 0;
		if(	++itp==13){
			int a =1;
		}

		bool any_change = false;
		double startproptime = cpuTime();
		static vec<int> detectors_to_check;

		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.

		//At level 0, need to propagate constant reaches/source nodes/edges...

		while(local_q<S->qhead){
			Lit l = S->trail[local_q++];
			Var v = var(l);

			int lev = S->level(v);
			while(lev>trail_lim.size()){
				newDecisionLevel();
			}
			if(v>= min_edge_var && v<min_edge_var+edge_list.size()){

				//this is an edge assignment
				int edge_num = v-min_edge_var;
				if(edge_list[edge_num].v<0)
					continue;
				int from = edge_list[edge_num].from;
				int to = edge_list[edge_num].to;
				trail.push({true,!sign(l), from,to,v});
				assert(edge_assignments[edge_num]==l_Undef);
				edge_assignments[edge_num]=sign(l) ? l_False:l_True;
				Assignment e = trail.last();
				assert(e.from==from);
				assert(e.to==to);

				if (!sign(l)){

					g.enableEdge(from,to,edge_num);
				}else{
					antig.disableEdge(from,to,edge_num);

				}
			}
		}
		stats_initial_propagation_time += cpuTime() - startproptime;
		dbg_sync();
			assert(dbg_graphsUpToDate());

			for(int d = 0;d<reach_detectors.size();d++){

				{
					double startdreachtime = cpuTime();

					reach_detectors[d]->positive_reach_detector->update();
					double reachUpdateElapsed = cpuTime()-startdreachtime;
					reachupdatetime+=reachUpdateElapsed;
				/*	for(int j =0;j<reach_lits[d].size();j++){
						Lit l = reach_lits[d][j];
						if( reach_detectors[d]->connected(j)){*/
					for(int u = 0;u<reach_detectors[d]->reach_lits.size();u++){
							//int u = reach_detectors[d]->positive_reach_detector->getChanged()[j];
							if(u>= reach_detectors[d]->reach_lits.size() || reach_detectors[d]->reach_lits[u] == lit_Undef)
								continue;
							//assert(reach_detectors[d]->reach_lits[u]!=lit_Undef);
							Lit l = (reach_detectors[d]->reach_lits[u]); // reach_lits[d][u];

							if( reach_detectors[d]->positive_reach_detector->connected_unsafe(u)){
								assert(dbg_reachable( reach_detectors[d]->source,u));
							if(S->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
								assert(dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
								if(S->dbg_solver)
									S->dbg_check_propagation(l);
#endif
								trail.push({false,true,d,0,var(l)});
								S->uncheckedEnqueue(l,reach_markers[d]) ;
							}else if (S->value(l)==l_False){
								double elapsed = cpuTime()-startdreachtime;
								reachtime+=elapsed;
								//conflict
								//The reason is all the literals in the shortest path in to s in d
								conflict.push(l);

								buildReachReason(u,opt_conflict_shortest_path ? *reach_detectors[d]->positive_dist_detector: *reach_detectors[d]->positive_reach_detector,conflict);
								//add it to s
								//return it as a conflict
#ifdef DEBUG_GRAPH
								for(int i = 0;i<conflict.size();i++)
											 assert(S->value(conflict[i])==l_False);
#endif
#ifdef DEBUG_SOLVER
								if(S->dbg_solver)
									S->dbg_check(conflict);
#endif


							//	reach_detectors[d]->positive_reach_detector->clearChanged();
								return false;
							}else{
								int  a=1;
							}
						}else{
#ifdef DEBUG_GRAPH
							assert(!dbg_reachable( reach_detectors[d]->source,u));
#endif
						}

					}
					//reach_detectors[d]->positive_reach_detector->clearChanged();
					double elapsed = cpuTime()-startdreachtime;
								reachtime+=elapsed;
				}

				{

					double startunreachtime = cpuTime();

					reach_detectors[d]->negative_reach_detector->update();
					double unreachUpdateElapsed = cpuTime()-startunreachtime;
					unreachupdatetime+=unreachUpdateElapsed;
					assert(dbg_graphsUpToDate());
					for(int u = 0;u<reach_detectors[d]->reach_lits.size();u++){
							//int u = reach_detectors[d]->negative_reach_detector->getChanged()[j];
							if(u>= reach_detectors[d]->reach_lits.size()||  reach_detectors[d]->reach_lits[u] == lit_Undef)
									continue;
							assert(reach_detectors[d]->reach_lits[u]!=lit_Undef);
							Lit l = ~ (reach_detectors[d]->reach_lits[u]); // ~reach_lits[d][u];
							if(! reach_detectors[d]->negative_reach_detector->connected_unsafe(u)){
								assert(dbg_notreachable( reach_detectors[d]->source,u));

							if(S->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
								assert(dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
								if(S->dbg_solver)
									S->dbg_check_propagation(l);
#endif
								trail.push({false,false,d,0,var(l)});
								S->uncheckedEnqueue(l,non_reach_markers[d]) ;
							}else if (S->value(l)==l_False){
								double elapsed = cpuTime()-startunreachtime;
												unreachtime+=elapsed;
								//conflict
								//The reason is a cut separating s from t
								conflict.push(l);
								buildNonReachReason(u,d,conflict);
#ifdef DEBUG_SOLVER
								if(S->dbg_solver)
									S->dbg_check(conflict);
#endif
								//add it to s
								//return it as a conflict
#ifdef DEBUG_GRAPH
								for(int i = 0;i<conflict.size();i++)
											 assert(S->value(conflict[i])==l_False);
#endif
							//	reach_detectors[d]->negative_reach_detector->clearChanged();
								return false;
							}
						}

					}
				//	reach_detectors[d]->negative_reach_detector->clearChanged();
					double elapsed = cpuTime()-startunreachtime;
					unreachtime+=elapsed;
				}
			}




		g.clearHistory();
		antig.clearHistory();

		detectors_to_check.clear();

		double elapsed = cpuTime()-startproptime;
					propagationtime+=elapsed;
					dbg_sync();
					dbg_sync_reachability();
		return true;
	};
	bool solveTheory(vec<Lit> & conflict){return true;};
	void drawFull(int from, int to){
			printf("digraph{\n");
			for(int i = 0;i<nNodes();i++){
				if(i==from){
					printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
				}else if (i==to){
					printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
				}else
					printf("n%d\n", i);
			}

			for(int i = 0;i<edge_list.size();i++){
				if(edge_list[i].v<0)
							continue;
				Edge & e = edge_list[i];
				const char * s = "black";
				if(S->value(e.v)==l_True)
					s="blue";
				else if (S->value(e.v)==l_False)
					s="red";
				printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from,e.to, e.v, s);
			}

			printf("}\n");
		}
	void drawFull(){
		printf("digraph{\n");
		for(int i = 0;i<nNodes();i++){
			printf("n%d\n", i);
		}

		for(int i = 0;i<edge_list.size();i++){
			Edge & e = edge_list[i];
			const char * s = "black";
			if(S->value(e.v)==l_True)
				s="blue";
			else if (S->value(e.v)==l_False)
				s="red";
			else{
				int  a=1;
			}
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from,e.to, e.v, s);
		}

		printf("}\n");
	}

	bool check_solved(){
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge & e = edge_list[i];
			lbool val = S->value(e.v);
			if(val==l_Undef){
				return false;
			}

			if(val==l_True){
				if(!g.hasEdge(e.from,e.to)){
					return false;
				}
				if(!antig.hasEdge(e.from,e.to)){
					return false;
				}
			}else{
				if(g.hasEdge(e.from,e.to)){
					return false;
				}
				if(antig.hasEdge(e.from,e.to)){
					return false;
				}

			}


		}
		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			for(int j = 0;j< d->reach_lits.size();j++){
				Lit l = d->reach_lits[j];
				if(l!=lit_Undef){
					int node = d->getNode(var(l));

					if(S->value(l)==l_True){
						if(!d->positive_reach_detector->connected(node)){
							return false;
						}
					}else if (S->value(l)==l_False){
						if( d->negative_reach_detector->connected(node)){
							return false;
						}
					}else{
						if(d->positive_reach_detector->connected(node)){
							return false;
						}
						if(!d->negative_reach_detector->connected(node)){
							return false;
						}
					}
				}
			}
		}
		return true;
	}

	bool dbg_solved(){
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge & e = edge_list[i];
			lbool val = S->value(e.v);
			assert(val!=l_Undef);

			if(val==l_True){
				assert(g.hasEdge(e.from,e.to));
				assert(antig.hasEdge(e.from,e.to));
			}else{
				assert(!g.hasEdge(e.from,e.to));
				assert(!antig.hasEdge(e.from,e.to));
			}


		}
		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			for(int j = 0;j< d->reach_lits.size();j++){
				Lit l = d->reach_lits[j];
				if(l!=lit_Undef){
					int node = d->getNode(var(l));

					if(S->value(l)==l_True){
						assert(d->positive_reach_detector->connected(node));
					}else if (S->value(l)==l_False){
						assert(! d->negative_reach_detector->connected(node));
					}else{
						assert(!d->positive_reach_detector->connected(node));
						assert(d->negative_reach_detector->connected(node));
					}
				}
			}
		}
#endif
		return true;
	}

	void drawCurrent(){

	}

	Lit newEdge(int from,int to, Var v = var_Undef)
    {

		if(v==var_Undef)
			v = S->newVar();
#ifdef DEBUG_GRAPH
		 dbg_graph->newEdge(from,to,v);
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg->newEdge(from,to,v);
#endif
		if(num_edges>0){
		}else
			min_edge_var=v;

		int index = v-min_edge_var;

		while(edge_list.size()<=index){
			edge_list.push({-1,-1,-1});
			edge_assignments.push(l_Undef);
		}

		inv_adj[to].push({v,from,to});

		num_edges++;
		edge_list[index].v =v;
		edge_list[index].from=from;
		edge_list[index].to =to;

		edges[from][to]= {v,from,to};
		g.addEdge(from,to,index);
		g.disableEdge(from,to, index);
		antig.addEdge(from,to,index);
		cutGraph.addEdge(from,to,index);
    	return mkLit(v,false);
    }

	void reaches(int from, int to, Var reach_var,int within_steps=-1){

#ifdef DEBUG_GRAPH
		 dbg_graph->reaches(from,  to,reach_var,within_steps);
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg->reaches(from,  to,reach_var,within_steps);
#endif
			assert(from<g.nodes);

			if (reach_info[from].source<0){

				reach_markers.push(S->newReasonMarker(this));
				int mnum = CRef_Undef- reach_markers.last();
				marker_map.growTo(mnum+1);
				marker_map[mnum] = reach_markers.size();
				//marker_map.insert(reach_markers.last(),reach_markers.size());

				non_reach_markers.push(S->newReasonMarker(this));
				//marker_map[non_reach_markers.last()]=-non_reach_markers.size();
				//marker_map.insert(non_reach_markers.last(),non_reach_markers.size());

				mnum = CRef_Undef- non_reach_markers.last();
				marker_map.growTo(mnum+1);
				marker_map[mnum] = -non_reach_markers.size();

				reach_detectors.push(new ReachDetector());
				if(reachalg==ALG_CONNECTIVITY){
					reach_detectors.last()->positive_reach_detector = new Connectivity<PositiveEdgeStatus>(from,g);
					reach_detectors.last()->negative_reach_detector = new Connectivity<NegativeEdgeStatus>(from,antig);
					if(opt_conflict_shortest_path)
						reach_detectors.last()->positive_dist_detector = new Dijkstra<PositiveEdgeStatus>(from,g);
				}else{

					reach_detectors.last()->positive_reach_detector = new Dijkstra<PositiveEdgeStatus>(from,g);
					reach_detectors.last()->negative_reach_detector = new Dijkstra<NegativeEdgeStatus>(from,antig);
					reach_detectors.last()->positive_dist_detector = (Dijkstra<PositiveEdgeStatus>*)reach_detectors.last()->positive_reach_detector;
					//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
				}
				//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
				reach_detectors.last()->source=from;

				reach_info[from].source=from;
				reach_info[from].detector=reach_detectors.last();

				reach_detectors.last()->within=within_steps;

				reach_detectors.last()->positive_reach_detector->update();
				reach_detectors.last()->negative_reach_detector->update();

				reach_detectors.last()->first_reach_var = reach_var;

			}

			ReachDetector * d = reach_info[from].detector;
			assert(d);



			while( d->reach_lits.size()<=to)
				d->reach_lits.push(lit_Undef);

			while(S->nVars()<=reach_var)
				S->newVar();

			Lit reachLit=mkLit(reach_var,false);

			if(d->reach_lits[to]==lit_Undef){
				d->reach_lits[to] = reachLit;

				while(d->reach_lit_map.size()<= reach_var- d->first_reach_var ){
					d->reach_lit_map.push(-1);
				}

				d->reach_lit_map[reach_var-d->first_reach_var]=to;
			}else{
				Lit r = d->reach_lits[to];
				//force equality between the new lit and the old reach lit, in the SAT solver
				S->addClause(~r, reachLit);
				S->addClause(r, ~reachLit);
			}




	    }

	void reachesAny(int from, Var firstVar,int within_steps=-1){

		for(int i = 0;i<g.nodes;i++){
			reaches(from,i,firstVar+i,within_steps);
		}

	  }
	void reachesAny(int from, vec<Lit> & reachlits_out,int within_steps=-1){
		for(int i = 0;i<g.nodes;i++){
			Var reachVar = S->newVar();
			reaches(from,i,reachVar,within_steps);
			reachlits_out.push(mkLit(reachVar,false));
		}
    }

};

};

#endif /* DGRAPH_H_ */
