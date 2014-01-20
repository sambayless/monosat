/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef DGRAPH_H_
#define DGRAPH_H_

#include "utils/System.h"
#include "core/Theory.h"
#include "Graph.h"
#include "Reach.h"
#include "Dijkstra.h"
#include "Connectivity.h"
#include "Distance.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "MaxFlow.h"
#include "IBFS.h"
#include "EdmondsKarp.h"
#include "EdmondsKarpAdj.h"
#include "Chokepoint.h"
#include "WeightedDijkstra.h"
#include "GraphTheoryTypes.h"
#include "utils/System.h"
#include "core/Solver.h"

#ifdef DEBUG_GRAPH
#include "TestGraph.h"
#endif
#include "AllPairsDetector.h"
#include "ReachDetector.h"
#include "DistanceDetector.h"
#include "MSTDetector.h"
namespace Minisat{

class GraphTheorySolver;


#ifdef DEBUG_SOLVER
#include "TestGraph.h"
#endif







class GraphTheorySolver:public GraphTheory{
public:

	double rnd_seed;
	Lit False;
	Lit True;
	int local_q;

	Solver * S;
	int id;
#ifdef DEBUG_GRAPH
	Solver * dbg;
	TestGraph *dbg_graph;
#endif
#ifdef DEBUG_SOLVER
	TestGraph * shadow_dbg;
#endif


	vec<lbool> edge_assignments;

	MSTDetector * mstDetector;
	vec<ReachabilityConstraint> unimplemented_graph_constraints;


	PositiveEdgeStatus g_status;
	NegativeEdgeStatus antig_status;
	CutEdgeStatus cutGraph_status;

	DynamicGraph<PositiveEdgeStatus> g;
	DynamicGraph<NegativeEdgeStatus> antig;
	DynamicGraph<CutEdgeStatus> cutGraph;

	Var min_edge_var;
	int num_edges;

	vec<Assignment> trail;
	vec<int> trail_lim;






	struct ReachInfo{
		int source;
		bool distance;
		Detector * detector;

		ReachInfo():source(-1),detector(NULL){}
	};


	vec<ReachInfo> dist_info;
	vec<ReachInfo> reach_info;
public:
	vec<Detector*> detectors;
	vec<ReachDetector*> reach_detectors;
	vec<DistanceDetector*> distance_detectors;


	vec<int> marker_map;


	vec<MaxFlow::Edge> cut;

	//Full matrix
	vec<vec<Edge> > edges;

	//Just a list of the edges
	vec<Edge> edge_list;

	vec<vec<Edge> > inv_adj;

	//vector of the weights for each edge
	vec<int> edge_weights;

	MaxFlow * mc;
	//MaxFlow * reachprop;

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
	double stats_decision_time;
	double stats_reason_initial_time;
	double stats_reason_time;
	int num_learnt_paths;
	int learnt_path_clause_length;
	int num_learnt_cuts;
	int learnt_cut_clause_length;

	int stats_mc_calls;
	vec<Lit> reach_cut;

	struct CutStatus{
		GraphTheorySolver & outer;
		int operator () (int id) const {

			if(outer.edge_assignments[id]==l_False){
				return 1;
			}else{
				return 0xF0F0F0;
			}
		}
		CutStatus(GraphTheorySolver & _outer):outer(_outer){}

	} cutStatus;

	struct PropCutStatus{
		GraphTheorySolver & outer;
		int operator () (int id) const {

			if(outer.edge_assignments[id]==l_Undef){
				return 1;
			}else{
				assert(outer.edge_assignments[id]==l_True);
				return 0xF0F0F0;
			}
		}
		PropCutStatus(GraphTheorySolver & _outer):outer(_outer){}

	}propCutStatus;

	GraphTheorySolver(Solver * S_, int _id=-1):S(S_),id(_id),g(g_status),antig(antig_status) ,cutGraph(cutGraph_status),cutStatus(*this),propCutStatus(*this){
		mstDetector = NULL;
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
			stats_decision_time=0;
			stats_reason_initial_time=0;
			stats_reason_time=0;
			 num_learnt_paths=0;
			 learnt_path_clause_length=0;
			 num_learnt_cuts=0;
			 learnt_cut_clause_length=0;

			 rnd_seed=opt_random_seed;

			if(mincutalg==ALG_IBFS){
				mc = new IBFS<CutEdgeStatus>(cutGraph);

			}else if (mincutalg == ALG_EDKARP_ADJ){

				mc = new EdmondsKarpAdj<CutStatus, CutEdgeStatus>(cutGraph,cutStatus);
				//reachprop = new EdmondsKarpAdj<PropCutStatus, NegativeEdgeStatus>(antig,propCutStatus);
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

	int getGraphID(){
		return id;
	}
	void printStats(){

		printf("Graph stats:\n");
		printf("Decision Time: %f\n", stats_decision_time);
		printf("Prop Time: %f (initial: %f)\n", propagationtime,stats_initial_propagation_time);
		printf("Conflict Time: %f (initial: %f)\n", stats_reason_time,stats_reason_initial_time);

		printf("Reach Time: %f (update Time: %f)\n", reachtime,reachupdatetime);
		printf("Unreach Time: %f (Update Time: %f)\n", unreachtime,unreachupdatetime);
		printf("Path Time: %f (#Paths: %d, AvgLength %f, total: %d)\n", pathtime, num_learnt_paths, (learnt_path_clause_length /  ((float) num_learnt_paths+1)),learnt_path_clause_length);
		printf("Min-cut Time: %f (%d calls, %f average, #Cuts: %d, AvgLength %f, total: %d)\n", mctime, stats_mc_calls,(mctime/(stats_mc_calls ? stats_mc_calls:1)),  num_learnt_cuts, (learnt_cut_clause_length /  ((float) num_learnt_cuts+1)),learnt_cut_clause_length);

		int stats_full_updates=0;
		int stats_fast_updates=0;
		int stats_failed_fast_updates=0;
		int stats_skip_deletes=0;
		int stats_skipped_updates=0;
		int skipable_deletions = 0;
		double stats_full_update_time=0;
		double stats_fast_update_time=0;

		for(int i = 0;i<reach_detectors.size();i++){
			skipable_deletions+=reach_detectors[i]->positive_reach_detector->stats_num_skipable_deletions;
			skipable_deletions+=reach_detectors[i]->negative_reach_detector->stats_num_skipable_deletions;

			stats_failed_fast_updates+=reach_detectors[i]->positive_reach_detector->stats_fast_failed_updates;
			stats_failed_fast_updates+=reach_detectors[i]->negative_reach_detector->stats_fast_failed_updates;

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
		printf("Dijkstra Fast Updates: %d (time: %f, average: %f, failed: %d)\n",stats_fast_updates, stats_fast_update_time,(stats_fast_update_time/(stats_fast_updates ? stats_fast_updates:1)),stats_failed_fast_updates);
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
		 dist_info.push();
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
			d->dbg_sync_reachability();
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

	Lit decideTheory(){
		if(!opt_decide_graph)
			return lit_Undef;
		double start = cpuTime();
		for(int i = 0;i<detectors.size();i++){
			Detector * r = detectors[i];
			Lit l =r->decide();
			if(l!=lit_Undef)
				return l;

		}
		stats_decision_time += cpuTime() - start;
		return lit_Undef;
	}

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
		double initial_start = cpuTime();
		backtrackUntil(p);

		double start = cpuTime();

		assert(d<detectors.size());
		detectors[d]->buildReason(p,reason,marker);
		double finish = cpuTime();
		stats_reason_time+=finish-start;
		stats_reason_initial_time+=start-initial_start;

	}



	bool dbg_reachable(int from, int to){
#ifdef DEBUG_GRAPH
		DefaultEdgeStatus tmp;
		/*DynamicGraph<> gtest(tmp);
		for(int i = 0;i<nNodes();i++){
			gtest.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
				continue;
			Edge e  = edge_list[i];
			if(S->assigns[e.v]==l_True){
				gtest.addEdge(e.from,e.to);
				assert(g.edgeEnabled(e.v));
			}else{
				assert(!g.edgeEnabled(e.v));
			}
		}*/

		Dijkstra<> d(from,g);
		d.update();
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

	int getEdgeID(Var v){
		assert(v>= min_edge_var && v<min_edge_var+edge_list.size());

						//this is an edge assignment
		int edge_num = v-min_edge_var;
		return edge_num;
	}

	void preprocess(){
		for (int i = 0;i<detectors.size();i++){
			detectors[i]->preprocess();
		}
	}

	bool propagateTheory(vec<Lit> & conflict){
		static int itp = 0;
		if(	++itp==2){
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

			for(int d = 0;d<detectors.size();d++){
				assert(conflict.size()==0);
				bool r =detectors[d]->propagate(trail,conflict);
				if(!r){
					propagationtime+= cpuTime()-startproptime;
					return false;
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
		for(int i = 0;i<detectors.size();i++){
			if(!detectors[i]->checkSatisfied()){
				return false;
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
		assert(check_solved());
		/*for(int i = 0;i<reach_detectors.size();i++){
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
		}*/
#endif
		return true;
	}

	void drawCurrent(){

	}

	CRef newReasonMarker(int detectorID){
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef- reasonMarker;
		marker_map.growTo(mnum+1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}

	Lit newEdge(int from,int to, Var v = var_Undef, int weight=1)
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
			edge_list.push({-1,-1,-1,-1,1});
			edge_assignments.push(l_Undef);
		}

		inv_adj[to].push({v,from,to,index,weight});

		num_edges++;
		edge_list[index].v =v;
		edge_list[index].from=from;
		edge_list[index].to =to;
		edge_list[index].edgeID=index;
		edge_list[index].weight=weight;

		edge_weights.push(weight);

		edges[from][to]= {v,from,to,index,weight};
		g.addEdge(from,to,index,weight);
		g.disableEdge(from,to, index);
		antig.addEdge(from,to,index,weight);
		cutGraph.addEdge(from,to,index,weight);
    	return mkLit(v,false);
    }
	int getEdgeID(int from, int to){
		assert(edges[from][to].edgeID>=0);
		return edges[from][to].edgeID;
	}
	int getWeight(int edgeID){
		return edge_list[edgeID].weight;
	}
	void reachesWithinSteps(int from, int to, Var reach_var, int within_steps){

	#ifdef DEBUG_GRAPH
			 dbg_graph->reaches(from,  to,reach_var,within_steps);
	#endif
	#ifdef DEBUG_SOLVER
			if(S->dbg_solver)
				shadow_dbg->reaches(from,  to,reach_var,within_steps);
	#endif
				assert(from<g.nodes);
				if(within_steps<=-1)
					within_steps = g.nodes;

				if (dist_info[from].source<0){
					DistanceDetector * d =new DistanceDetector(detectors.size(), this,g,antig,from,within_steps,drand(rnd_seed));
					detectors.push(d);
					//reach_detectors.push(reach_detectors.last());
					distance_detectors.push(d);
					//detectors.push(new DistanceDetector(detectors.size(), this,g,antig,from,within_steps,drand(rnd_seed)));
						//reach_detectors.push(reach_detectors.last());

					assert(detectors.last()->getID()==detectors.size()-1);





					//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
					//reach_detectors.last()->source=from;

					dist_info[from].source=from;
					dist_info[from].detector=detectors.last();

					//reach_detectors.last()->within=within_steps;

				}

				DistanceDetector * d = (DistanceDetector*)dist_info[from].detector;
				assert(d);

				d->addLit(from,to,reach_var,within_steps);


		    }

	void implementConstraints(){
		if(opt_allpairs_percentage>=1){
			for(int i = 0;i<unimplemented_graph_constraints.size();i++){
				ReachabilityConstraint c = unimplemented_graph_constraints[i];
				reaches_private(c.from,c.to,c.reach_var,c.distance);
			}
		}else if (opt_allpairs_percentage==0){
			for(int i = 0;i<unimplemented_graph_constraints.size();i++){
				ReachabilityConstraint c = unimplemented_graph_constraints[i];
				allpairs(c.from,c.to,c.reach_var,c.distance);
			}
		}else{
			vec<bool> seen;
			int count=0;
			seen.growTo(nNodes());
			for(int i = 0;i<unimplemented_graph_constraints.size();i++){
						ReachabilityConstraint c = unimplemented_graph_constraints[i];
						if(!seen[c.from]){
							seen[c.from]=true;
							count++;
						}
					}
			double frac = ((double)count)/((double)nNodes());

			if (opt_verb>0 && frac>=opt_allpairs_percentage){
				printf("Allpairs solver triggered for graph %d by percentage of source nodes: %d/%d=%f>%f\n",getGraphID() ,count,nNodes(),frac,(double)opt_allpairs_percentage);
			}

			for(int i = 0;i<unimplemented_graph_constraints.size();i++){
									ReachabilityConstraint c = unimplemented_graph_constraints[i];
									if(frac>=opt_allpairs_percentage)
										allpairs(c.from,c.to,c.reach_var,c.distance);
									else
										reaches_private(c.from,c.to,c.reach_var,c.distance);
								}
		}
		unimplemented_graph_constraints.clear();
	}
	void allpairs(int from, int to, Var reach_var,int within_steps=-1){
				//for now, reachesWithinSteps to be called instead

		#ifdef DEBUG_GRAPH
				 dbg_graph->reaches(from,  to,reach_var,within_steps);
		#endif
		#ifdef DEBUG_SOLVER
				if(S->dbg_solver)
					shadow_dbg->reaches(from,  to,reach_var,within_steps);
		#endif
					assert(from<g.nodes);
					if(within_steps>g.nodes)
						within_steps=-1;

					if (reach_info[from].source<0){


						detectors.push(new AllPairsDetector(detectors.size(), this,g,antig,drand(rnd_seed)));
							//reach_detectors.push(reach_detectors.last());






						assert(detectors.last()->getID()==detectors.size()-1);





						//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
						//reach_detectors.last()->source=from;

						reach_info[from].source=from;
						reach_info[from].detector=detectors.last();

						//reach_detectors.last()->within=within_steps;

					}

					AllPairsDetector * d =(AllPairsDetector*) reach_info[from].detector;
					assert(d);

					d->addLit(from,to,reach_var,within_steps);


			    }
	void reaches_private(int from, int to, Var reach_var,int within_steps=-1){
			//for now, reachesWithinSteps to be called instead
			if(within_steps>=0 || opt_force_distance_solver){
				reachesWithinSteps(from,to,reach_var,within_steps);
				return;
			}

	#ifdef DEBUG_GRAPH
			 dbg_graph->reaches(from,  to,reach_var,within_steps);
	#endif
	#ifdef DEBUG_SOLVER
			if(S->dbg_solver)
				shadow_dbg->reaches(from,  to,reach_var,within_steps);
	#endif
				assert(from<g.nodes);
				if(within_steps>g.nodes)
					within_steps=-1;

				if (reach_info[from].source<0){



						ReachDetector*rd = new ReachDetector(detectors.size(), this,g,antig,from,drand(rnd_seed));
						detectors.push(rd);
						reach_detectors.push(rd);



					assert(detectors.last()->getID()==detectors.size()-1);


					//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
					//reach_detectors.last()->source=from;

					reach_info[from].source=from;
					reach_info[from].detector=detectors.last();

					//reach_detectors.last()->within=within_steps;

				}

				ReachDetector * d = (ReachDetector*) reach_info[from].detector;
				assert(d);

				d->addLit(from,to,reach_var,within_steps);


		    }

	void reaches(int from, int to, Var reach_var,int within_steps=-1){
			unimplemented_graph_constraints.push({from,to,within_steps,reach_var});
			//to allow us to alter the solving algorithm based on the number and type of constraints, we aren't implementing them here directly any more - instead,
			//we just store the constraints in this vector, then implement them later when 'implementConstraints' is called.
	    }

	void reachesAny(int from, Var firstVar,int within_steps=-1){
		for(int i = 0;i<g.nodes;i++){
			reaches(from,i,firstVar+i,within_steps);
		}
	}

	void reachesAny(int from, vec<Lit> & reachlits_out,int within_steps=-1){
		for(int i = 0;i<g.nodes;i++){
			Var reachVar = S->newVar();
			//reaches(from,i,reachVar,within_steps);
			reaches(from,i,reachVar,within_steps);
			reachlits_out.push(mkLit(reachVar,false));
		}
    }
	//v will be true if the minimum weight is <= the specified value
	void minimumSpanningTree(Var v, int minimum_weight){
		if(!mstDetector){
			mstDetector = new MSTDetector(detectors.size(),this, g, antig, this->edge_weights,drand(rnd_seed));
			detectors.push(mstDetector);
		}
		mstDetector->addWeightLit(v,minimum_weight);
	}
	void edgeInMinimumSpanningTree(Var edgeVar, Var var){
		if(!mstDetector){
			mstDetector = new MSTDetector(detectors.size(),this, g, antig, this->edge_weights,drand(rnd_seed));
			detectors.push(mstDetector);
		}
		int edgeid =getEdgeID(edgeVar);
		assert(edgeid>=0);
		mstDetector->addTreeEdgeLit(edgeid,var);
	}
};

};

#endif /* DGRAPH_H_ */
