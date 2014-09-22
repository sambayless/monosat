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

#include "dgl/Reach.h"
#include "dgl/Dijkstra.h"
#include "dgl/BFS.h"
#include "dgl/UnweightedDistance.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "dgl/MaxFlow.h"
#include "dgl/IBFS.h"
#include "dgl/EdmondsKarp.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/Chokepoint.h"
#include "dgl/WeightedDijkstra.h"
#include "GraphTheoryTypes.h"
#include "utils/System.h"
#include "core/Solver.h"

#include "AllPairsDetector.h"
#include "ReachDetector.h"
#include "ConnectDetector.h"
#include "DistanceDetector.h"
#include "MSTDetector.h"
#include "MaxflowDetector.h"
#include "ConnectedComponentsDetector.h"
#include "CycleDetector.h"
#include "SteinerDetector.h"
#include <vector>
#include <gmpxx.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>


using namespace dgl;
namespace Minisat{
template<typename Weight>
class GraphTheorySolver;



template<typename Weight>
class GraphTheorySolver:public Theory{
public:

	double rnd_seed;
	Lit False;
	Lit True;
	int local_q;
private:
	Solver * S;
public:
	int id;


	bool all_edges_unit=true;
	vec<lbool> assigns;

	MSTDetector<Weight> * mstDetector;
	vec<ReachabilityConstraint> unimplemented_reachability_constraints;
	vec<ConnectivityConstraint> unimplemented_connectivity_constraints;




	DynamicGraph g;
	DynamicGraph antig;
	DynamicGraph cutGraph;

	//Var min_edge_var;
	//int num_edges;

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
	vec<ReachInfo> connect_info;
public:
	vec<Detector*> detectors;
	vec<ReachDetector<Weight>*> reach_detectors;
	vec<ConnectDetector<Weight>*> connect_detectors;
	vec<DistanceDetector<Weight>*> distance_detectors;
	vec<MaxflowDetector<Weight>*> flow_detectors;
	ConnectedComponentsDetector<Weight>* component_detector;
	CycleDetector<Weight> * cycle_detector;
	vec<SteinerDetector<Weight>*>  steiner_detectors;

	vec<int> marker_map;


	std::vector<MaxFlowEdge> cut;

	//Full matrix
	vec<vec<Edge> > edges;



	//Just a list of the edges
	vec<Edge> edge_list;
	vec<vec<Edge> > undirected_adj;
	vec<vec<Edge> > inv_adj;

	//vector of the weights for each edge
	std::vector<Weight> edge_weights;
	//Some algorithms support computations using rational weights (or capacities).
	std::vector<mpq_class> rational_weights;
	bool requiresPropagation;
	MaxFlow<int> * mc;
	//MaxFlow * reachprop;

    vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData{
		int isEdge:1;
		int occursPositive:1;
		int occursNegative:1;
		int detector_edge:29;//the detector this variable belongs to, or its edge number, if it is an edge variable
		Var solverVar;
	};

	vec<VarData> vars;
	int theory_index;
public:

	double mctime;
	double reachtime;
	double unreachtime;
	double pathtime;
	double propagationtime;
	long stats_propagations;
	long stats_num_conflicts;
	long stats_decisions;
	long stats_num_reasons;

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
	int stats_pure_skipped;
	int stats_mc_calls;
	long stats_propagations_skipped;
	vec<Lit> reach_cut;

	struct CutStatus{
		GraphTheorySolver & outer;
		int operator [] (int id) const {

			if(outer.value(outer.edge_list[id].v) ==l_False){
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

			if(outer.value(outer.edge_list[id].v) ==l_Undef){
				return 1;
			}else{
				assert(outer.value(outer.edge_list[id].v)==l_True);
				return 0xF0F0F0;
			}
		}
		PropCutStatus(GraphTheorySolver & _outer):outer(_outer){}

	}propCutStatus;

	GraphTheorySolver(Solver * S_, int _id=-1):S(S_),id(_id),cutStatus(*this),propCutStatus(*this){
		mstDetector = NULL;
		//True = mkLit(S->newVar(),false);
			//False=~True;
			//S->addClause(True);
#ifdef RECORD
			{
				char t[30];
				sprintf(t,"TEST_GRAPH%d",id);
				g.outfile=fopen(t,"w");
			}
			{
						char t[30];
						sprintf(t,"TEST_ANTI_GRAPH%d",id);
						antig.outfile=fopen(t,"w");
					}
#endif
			local_q=0;
			theory_index=0;
			mctime=0;
			stats_mc_calls=0;
			reachtime=0;
			unreachtime=0;
			pathtime=0;
			propagationtime=0;
			stats_propagations=0;
			stats_num_conflicts=0;
			stats_num_reasons=0;
			stats_decisions = 0;
			stats_propagations_skipped=0;
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
			 component_detector=NULL;
			 requiresPropagation=true;
			 rnd_seed=opt_random_seed;

			if(mincutalg==MinCutAlg::ALG_IBFS){
				mc = new IBFS(cutGraph);

			}else if (mincutalg == MinCutAlg::ALG_EDKARP_ADJ){

				mc = new EdmondsKarpAdj<CutStatus,int>(cutGraph, cutStatus);
				//reachprop = new EdmondsKarpAdj<PropCutStatus, NegativeEdgeStatus>(antig,propCutStatus);
			}else{
				mc = new EdmondsKarp<int>(cutGraph);
			}


	}

	void printStats(int detailLevel){
		if(detailLevel>0){
			for(Detector * d:detectors)
				d->printStats();
		}

		printf("Graph stats:\n");
/*		printf("Decision Time: %f\n", stats_decision_time);
		printf("Prop Time: %f (initial: %f)\n", propagationtime,stats_initial_propagation_time);
		printf("Conflict Time: %f (initial: %f)\n", stats_reason_time,stats_reason_initial_time);*/

	/*	printf("Reach Time: %f (update Time: %f)\n", reachtime,reachupdatetime);
		printf("Unreach Time: %f (Update Time: %f)\n", unreachtime,unreachupdatetime);
		printf("Path Time: %f (#Paths: %d, AvgLength %f, total: %d)\n", pathtime, num_learnt_paths, (learnt_path_clause_length /  ((float) num_learnt_paths+1)),learnt_path_clause_length);
		printf("Min-cut Time: %f (%d calls, %f average, #Cuts: %d, AvgLength %f, total: %d)\n", mctime, stats_mc_calls,(mctime/(stats_mc_calls ? stats_mc_calls:1)),  num_learnt_cuts, (learnt_cut_clause_length /  ((float) num_learnt_cuts+1)),learnt_cut_clause_length);
*/

		printf("Propagations: %ld (%f s, avg: %f s, %ld skipped)\n",stats_propagations ,propagationtime, (propagationtime)/((double)stats_propagations+1),stats_propagations_skipped);
		printf("Decisions: %ld (%f s, avg: %f s)\n",stats_decisions,stats_decision_time, (stats_decision_time)/((double)stats_decisions+1));
		printf("Conflicts: %ld\n",stats_num_conflicts);
		printf("Reasons: %ld (%f s, avg: %f s)\n",stats_num_reasons,stats_reason_time, (stats_reason_time)/((double)stats_num_reasons+1));
		fflush(stdout);
	}

	inline int getTheoryIndex(){
	    	return theory_index;
	    }
	  inline  void setTheoryIndex(int id){
	    	theory_index=id;
	    }
	inline int getGraphID(){
		return id;
	}
	inline bool isEdgeVar(Var v){
		assert(v<vars.size());
		return vars[v].isEdge;
	}
	inline int getEdgeID(Var v){
		assert(isEdgeVar(v));
		return vars[v].detector_edge;
	}
	inline int getDetector(Var v){
		assert(!isEdgeVar(v));
		return vars[v].detector_edge;
	}

	inline Var getEdgeVar(int edgeID){
		Var v = edge_list[edgeID].v;
		assert(v<vars.size());
		assert(vars[v].isEdge);
		return v;
	}

	void makeEqual(Lit l1, Lit l2){
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1,o2);
		S->addClause(o1, ~o2);
	}
	void makeEqualInSolver(Lit l1, Lit l2){
		S->addClause(~l1,l2);
		S->addClause(l1, ~l2);
	}
	void addClause(Lit l1){
		Lit o1 = toSolver(l1);
		S->addClause(o1);
	}
	void addClause(Lit l1, Lit l2){
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(o1,o2);
	}
	void addClause(Lit l1, Lit l2, Lit l3){
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		Lit o3 = toSolver(l3);
		S->addClause(o1,o2,o3);
	}
	void addClause(vec<Lit> & c){
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		S->addClause(tmp_clause);
	}
	void addClauseSafely(vec<Lit> & c){
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);

		S->addClauseSafely(tmp_clause);
	}
	/*void addConflictClause(vec<Lit> & c){
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		CRef ignore;
		S->addConflictClause(tmp_clause,ignore);
	}*/
	Var newVar(int forDetector=-1, bool connectToTheory=false){
		Var s= S->newVar();
	/*	Var v = vars.size();
		vars.push();
		vars[v].isEdge=false;
		vars[v].occursPositive=true;
		vars[v].occursNegative=true;
		vars[v].detector_edge=-1;
		vars[v].solverVar=s;
		if(connectToTheory){
			S->setTheoryVar(s,getTheoryIndex(),v);
		}*/
		return newVar(s,forDetector,false,connectToTheory);
	}
	Var newVar(Var solverVar, int detector, bool isEdge=false, bool connectToTheory=true){
		while(S->nVars()<=solverVar)
				S->newVar();
		Var v = vars.size();
		vars.push();
		vars[v].isEdge=isEdge;
		vars[v].detector_edge=detector;
		vars[v].solverVar=solverVar;
		assigns.push(l_Undef);
		if(connectToTheory){
			S->setTheoryVar(solverVar,getTheoryIndex(),v);
			assert(toSolver(v)==solverVar);
		}
		if(!isEdge && detector>=0)
			detectors[detector]->addVar(v);
		return v;
	}
	inline int level(Var v){
		return S->level(toSolver(v));
	}
	inline int decisionLevel(){
		return trail_lim.size(); //S->decisionLevel();
	}
	inline int nVars()const{
		return vars.size();//S->nVars();
	}
	inline Var toSolver(Var v){
		//return v;
		assert(v<vars.size());
		//assert(S->hasTheory(vars[v].solverVar));
		//assert(S->getTheoryVar(vars[v].solverVar)==v);
		return vars[v].solverVar;
	}

	inline Lit toSolver(Lit l){
		//assert(S->hasTheory(vars[var(l)].solverVar));
		//assert(S->getTheoryVar(vars[var(l)].solverVar)==var(l));
		return mkLit(vars[var(l)].solverVar,sign(l));
	}

	void toSolver(vec<Lit> & c){
		for(int i = 0;i<c.size();i++){
			c[i]=toSolver(c[i]);
		}
	}

	inline lbool value(Var v){
		if(assigns[v]!=l_Undef)
			assert(S->value(toSolver(v))==assigns[v]);

		return assigns[v]; //S->value(toSolver(v));
	}
	inline lbool value(Lit l){
		if(assigns[var(l)]!=l_Undef){
			assert(S->value(toSolver(l))==  (assigns[var(l)]^ sign(l)));
		}
		return assigns[var(l)]^ sign(l);;//S->value(toSolver(l));
	}
	inline lbool dbg_value(Var v){
		return S->value(toSolver(v));
	}
	inline lbool dbg_value(Lit l){
		return S->value(toSolver(l));
	}
	inline bool enqueue(Lit l, CRef reason){
		assert(assigns[var(l)]==l_Undef);

		Lit sl = toSolver(l);
		if( S->enqueue(sl,reason)){
			enqueueTheory(l);
			return true;
		}else{
			return false;
		}
	}



     ~GraphTheorySolver(){};
	 int newNode(){

		 edges.push();
		 for(int i = 0;i<edges.size();i++)
			 edges[i].growTo(edges.size());
		 inv_adj.push();
		 undirected_adj.push();
		 reach_info.push();
		 connect_info.push();
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
		return g.nodes();
	}
	bool isNode(int n){
		return n>=0 && n<nNodes();
	}


		bool dbg_propgation(Lit l){
#ifndef NDEBUG
			static vec<Lit> c;
			c.clear();
			for(int i = 0;i<S->trail.size() ;i++){
				if(!S->hasTheory(S->trail[i]) || S->getTheoryID(S->trail[i])!= getTheoryIndex())
					continue;
				Lit l = S->getTheoryLit(S->trail[i]);
				Var v =  var(l);
				if(isEdgeVar(v)){
					Edge & e = edge_list[getEdgeID(v)];
					c.push(l);
				}

		/*		if(v>=min_edge_var && v<min_edge_var+num_edges){
					if(edge_list[v-min_edge_var].v<0)
								continue;
					Edge e = edge_list[v-min_edge_var];
					c.push(l);
				}*/
			}
			c.push(~l);
			//bool res = dbg->solve(c);
			//assert(~res);
#endif
			return true;
		}

	void	dbg_sync_reachability(){
/*
#ifdef DEBUG_GRAPH

		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			d->dbg_sync_reachability();
		}
#endif
*/
	}
	void dbg_sync(){
#ifdef DEBUG_DIJKSTRA

		for(int i = 0;i<assigns.size();i++){
			lbool val = assigns[i];

			if(val!=l_Undef){
				bool found=false;
				for(int j = 0;j<trail.size();j++){
					if(trail[j].var==i){
						assert(!found);
						assert(trail[j].assign== (val==l_True));
						found=true;
					}
				}
				assert(found);
			}else{
				for(int j = 0;j<trail.size();j++){
					if(trail[j].var==i){
						assert(false);
					}
				}
			}
			if(val!=l_Undef)
					assert(val==S->value(toSolver(i)));
		}
		/*static vec<lbool> assigned;
		assigned.clear();*.


		/*for(int i = 0;i<edge_list.size();i++)
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
				assert(value(l)==l_True);
				int expected_level = level(var(l));
				assert(level(var(l))==lev);
				int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
				if(edge_list[edge_num].v<0)
							continue;
				assert(assigned[edge_num]==l_Undef);

				assigned[edge_num] = sign(l)?l_False:l_True;
			}else{
				//Reachability assignment


			}
		}

		for(int i = 0;i<assigned.size();i++){
			assert(ass[i]== assigned[i]);
		}

		for(int i = 0;i<S->trail.size();i++){

			Lit l = S->trail[i];
			if(S->getTheoryID(l)==getTheoryIndex()){
				l = S->getTheoryLit(l);
			Var v = var(l);

			int lev = level(v);
			if(isEdgeVar(v)){
				int edge_num = getEdgeID(v);
				if(edge_list[edge_num].v<0)
					continue;
				lbool assigned_val=assigned[edge_num];
				assert(assigned_val== (sign(l)?l_False:l_True));
			}
			if(v>= min_edge_var && v<min_edge_var+num_edges){
				int edge_num = v-min_edge_var;
				if(edge_list[edge_num].v<0)
					continue;
				lbool assigned_val=assigned[edge_num];
				assert(assigned_val== (sign(l)?l_False:l_True));
			}

			}
		}*/


#endif
	}
	void dbg_full_sync(){
	#ifdef DEBUG_GRAPH
			dbg_sync();
			for(int i =0;i<edge_list.size();i++){

				if(edge_list[i].edgeID>=0 && g.edgeEnabled(i)){
					assert(value(edge_list[i].v)==l_True);
				}else if (edge_list[i].edgeID>=0 &&  !antig.edgeEnabled(i)){
					assert(value(edge_list[i].v)==l_False);
				}
			}

	#endif
		}
	
	

	void backtrackUntil(int level){
		static int it = 0;

		bool changed=false;
		//need to remove and add edges in the two graphs accordingly.
		if(trail_lim.size()>level){

			int stop = trail_lim[level];
			for(int i = trail.size()-1;i>=trail_lim[level];i--){

				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if(e.isEdge){
					assert(dbg_value(e.var)==l_Undef);
					int edge_num = getEdgeID(e.var); //e.var-min_edge_var;

					if(e.assign){
						g.disableEdge(e.from,e.to, edge_num);
					}else{
						antig.enableEdge(e.from,e.to,edge_num);
						assert(antig.hasEdge(e.from,e.to));
					}
				}else{
				  //This is a reachability literal				  
				  detectors[getDetector(e.var)]->unassign(mkLit(e.var,!e.assign));
				}
				assigns[e.var]=l_Undef;
				changed=true;
			}
			trail.shrink(trail.size()-stop);
			trail_lim.shrink(trail_lim.size()-level);
			assert(trail_lim.size()==level);


		}
		if(changed){
			requiresPropagation=true;
			g.markChanged();
			antig.markChanged();
			cutGraph.markChanged();
		}
/*		if(local_q>S->qhead)
			local_q=S->qhead;*/
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
		double start = rtime(1);

		dbg_full_sync();
		for(int i = 0;i<detectors.size();i++){
			Detector * r = detectors[i];
			Lit l =r->decide();
			if(l!=lit_Undef){
				stats_decisions++;
				return toSolver(l);
			}
		}
		stats_decision_time += rtime(1) - start;
		return lit_Undef;
	}

	void backtrackUntil(Lit p){
			//need to remove and add edges in the two graphs accordingly.
			int i = trail.size()-1;
			for(;i>=0;i--){
				Assignment e = trail[i];
				if(e.isEdge){
					int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
					assert(assigns[e.var]!=l_Undef);
					assigns[e.var]=l_Undef;
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
					assigns[e.var]=l_Undef;
					detectors[getDetector(e.var)]->unassign(mkLit(e.var,!e.assign));
				}
			}

			trail.shrink(trail.size()-(i+1));
			if(i>0){
				requiresPropagation=true;
				g.markChanged();
				antig.markChanged();
				cutGraph.markChanged();
			}

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
		CRef marker = S->reason(var(toSolver(p)));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef- marker;
		int d = marker_map[pos];
		//double initial_start = rtime(1);
		double start = rtime(1);
		backtrackUntil(p);



		assert(d<detectors.size());
		detectors[d]->buildReason(p,reason,marker);
		toSolver(reason);
		double finish = rtime(1);
		stats_reason_time+=finish-start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;

	}



	bool dbg_reachable(int from, int to, bool undirected=false){
#ifdef DEBUG_DIJKSTRA

		if(undirected){
			UnweightedDijkstra<Reach::NullStatus, true> d(from,g);
			d.update();
			return d.connected(to);
		}else{
			UnweightedDijkstra<> d(from,g);
			d.update();
			return d.connected(to);
		}

#else
		return true;
#endif
	}

	bool dbg_notreachable(int from, int to, bool undirected=false){

#ifndef NDEBUG
		//drawFull(from,to);

		DynamicGraph g;
		for(int i = 0;i<nNodes();i++){
			g.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge e  = edge_list[i];
			if(value(e.v)!=l_False){
				g.addEdge(e.from,e.to);
			}
		}
		if(undirected){
			UnweightedDijkstra<Reach::NullStatus, true> d(from,g);

			return !d.connected(to);
		}else{
			UnweightedDijkstra<Reach::NullStatus, false> d(from,g);

					return !d.connected(to);
		}
#endif
		return true;

	}

	bool dbg_graphsUpToDate(){
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
				continue;
			Edge e = edge_list[i];
			lbool val = value(e.v);

			if(val==l_True || val==l_Undef){
				assert(antig.edgeEnabled(i));
				//assert(antig.hasEdge(e.from,e.to));
			}else{
				assert(!antig.edgeEnabled(i));
				//assert(!antig.hasEdge(e.from,e.to));
			}
			if(val==l_True){
				assert(g.edgeEnabled(i));
				//assert(g.hasEdge(e.from,e.to));
			}else{
				assert(!g.edgeEnabled(i));
				//assert(!g.hasEdge(e.from,e.to));
			}
		}

#endif
		return true;
	}

/*	int getEdgeID(Var v){
		assert(v>= min_edge_var && v<min_edge_var+edge_list.size());

						//this is an edge assignment
		int edge_num = v-min_edge_var;
		return edge_num;
	}*/

	void preprocess(){
		for (int i = 0;i<detectors.size();i++){
			detectors[i]->preprocess();
		}
	}
	void setLiteralOccurs(Lit l, bool occurs){
		if(isEdgeVar(var(l))){
			//don't do anything
		}else{
			//this is a graph property detector var
			if(!sign(l) && vars[var(l)].occursPositive!=occurs)
				detectors[getDetector(var(l))]->setOccurs(l,occurs);
			else if(sign(l) && vars[var(l)].occursNegative!=occurs)
				detectors[getDetector(var(l))]->setOccurs(l,occurs);
		}

	}

	void enqueueTheory(Lit l){
		Var v = var(l);

		int lev = level(v);

		assert(decisionLevel()<=lev);

		while(lev>trail_lim.size()){
			newDecisionLevel();
		}

		if(assigns[var(l)]!=l_Undef){
			return;//this is already enqueued.
		}
		assert(assigns[var(l)]==l_Undef);
		assigns[var(l)]=sign(l) ? l_False:l_True;
		requiresPropagation=true;
		//printf("enqueue %d\n", dimacs(l));
#ifndef NDEBUG

		{
			for(int i = 0;i<trail.size();i++){
				assert(trail[i].var !=v);
			}
		}
#endif
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"enqueue %d\n", dimacs(l));

			fprintf(g.outfile,"\n");
			fflush(g.outfile);
		}
		if(antig.outfile){
			fprintf(antig.outfile,"enqueue %d\n", dimacs(l));
			fprintf(antig.outfile,"\n");
				fflush(antig.outfile);
			}
#endif

		//if(v>= min_edge_var && v<min_edge_var+edge_list.size())
		if(isEdgeVar(var(l))){

			//this is an edge assignment
			int edge_num = getEdgeID(var(l)); //v-min_edge_var;
			assert(edge_list[edge_num].v==var(l));

			int from = edge_list[edge_num].from;
			int to = edge_list[edge_num].to;
			trail.push({true,!sign(l), from,to,v});

			Assignment e = trail.last();
			assert(e.from==from);
			assert(e.to==to);

			if (!sign(l)){
				g.enableEdge(from,to,edge_num);
			}else{
				antig.disableEdge(from,to,edge_num);
			}

		}else{

			trail.push({false,!sign(l), 0,0,v});
			//this is an assignment to a non-edge atom. (eg, a reachability assertion)
			detectors[getDetector(var(l))]->assign(l);
		}

	};
	bool propagateTheory(vec<Lit> & conflict){
		static int itp = 0;
		if(	++itp==62279){
			int a =1;
		}
		stats_propagations++;
		dbg_sync();
		if(!requiresPropagation){
			stats_propagations_skipped++;
			assert(dbg_graphsUpToDate());
			return true;
		}

		bool any_change = false;
		double startproptime = rtime(1);
		static vec<int> detectors_to_check;

		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.

		//At level 0, need to propagate constant reaches/source nodes/edges...


		//stats_initial_propagation_time += rtime(1) - startproptime;
		dbg_sync();
		assert(dbg_graphsUpToDate());

		for(int d = 0;d<detectors.size();d++){
			assert(conflict.size()==0);
			bool r =detectors[d]->propagate(conflict);
			if(!r){
				stats_num_conflicts++;
				toSolver(conflict);
				propagationtime+= rtime(1)-startproptime;
				return false;
			}
		}



		dbg_full_sync();

		requiresPropagation=false;
		g.clearChanged();
		antig.clearChanged();
		cutGraph.clearChanged();

		g.clearHistory();
		antig.clearHistory();
		cutGraph.clearHistory();
		detectors_to_check.clear();

		double elapsed = rtime(1)-startproptime;
		propagationtime+=elapsed;
		dbg_sync();
		dbg_sync_reachability();
		return true;
	};

	bool solveTheory(vec<Lit> & conflict){
		requiresPropagation=true;//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	};

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
				if(value(e.v)==l_True)
					s="blue";
				else if (value(e.v)==l_False)
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
			if(value(e.v)==l_True)
				s="blue";
			else if (value(e.v)==l_False)
				s="red";
			else{
				int  a=1;
			}
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from,e.to, e.v, s);
		}

		printf("}\n");
	}

	bool check_solved(){
		if(opt_print_graph){
			drawFull();
		}
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge & e = edge_list[i];
			lbool val = value(e.v);
			if(val==l_Undef){
				return false;
			}

			if(val==l_True){
			/*	if(!g.hasEdge(e.from,e.to)){
					return false;
				}
				if(!antig.hasEdge(e.from,e.to)){
					return false;
				}*/
				if(!g.edgeEnabled(e.edgeID)){
					return false;
				}
				if(!antig.edgeEnabled(e.edgeID)){
					return false;
				}
			}else{
				/*if(g.hasEdge(e.from,e.to)){
					return false;
				}*/
				if(g.edgeEnabled(e.edgeID)){
					return false;
				}
				if(antig.edgeEnabled(e.edgeID)){
					return false;
				}
				/*if(antig.hasEdge(e.from,e.to)){
					return false;
				}*/

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
			lbool val = value(e.v);
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

					if(value(l)==l_True){
						assert(d->positive_reach_detector->connected(node));
					}else if (value(l)==l_False){
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
	int nEdges(){
		return edge_list.size();
	}
	CRef newReasonMarker(int detectorID){
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef- reasonMarker;
		marker_map.growTo(mnum+1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}

	Lit newEdge(int from,int to, Var outerVar = var_Undef, Weight weight=1)
    {
		assert(outerVar!=var_Undef);
	/*	if(outerVar==var_Undef)
			outerVar = S->newVar();*/


		all_edges_unit&=(weight==1);
		int index = edge_list.size();
		edge_list.push();
		Var v = newVar(outerVar,index,true);

/*
		if(num_edges>0){
		}else
			min_edge_var=v;

		int index = v-min_edge_var;*/
/*
		while(edge_list.size()<=index){
			edge_list.push({-1,-1,-1,-1,-1,1});
			assigns.push(l_Undef);
		}*/
		undirected_adj[to].push({v,outerVar,from,to,index});
		undirected_adj[from].push({v,outerVar,to,from,index});
		inv_adj[to].push({v,outerVar,from,to,index});

		//num_edges++;
		edge_list[index].v =v;
		edge_list[index].outerVar =outerVar;
		edge_list[index].from=from;
		edge_list[index].to =to;
		edge_list[index].edgeID=index;
		//edge_list[index].weight=weight;
		if(edge_weights.size()<=index){
			edge_weights.resize(index+1);
		}
		edge_weights[index]=weight;

		edges[from][to]= {v,outerVar,from,to,index};
		g.addEdge(from,to,index);
		g.disableEdge(from,to, index);
		antig.addEdge(from,to,index);
		cutGraph.addEdge(from,to,index);

    	return mkLit(v,false);
    }
	int getEdgeID(int from, int to){
		assert(edges[from][to].edgeID>=0);
		return edges[from][to].edgeID;
	}
	Weight getWeight(int edgeID){
		return edge_weights[edgeID];
	}
	void reachesWithinSteps(int from, int to, Var reach_var, int within_steps){

				assert(from<g.nodes());
				if(within_steps<=-1)
					within_steps = g.nodes();

				if (dist_info[from].source<0){
					DistanceDetector<Weight> * d =new DistanceDetector<Weight>(detectors.size(), this,edge_weights,g,antig,from,drand(rnd_seed));
					detectors.push(d);
					distance_detectors.push(d);
					assert(detectors.last()->getID()==detectors.size()-1);
					dist_info[from].source=from;
					dist_info[from].detector=detectors.last();
				}

				DistanceDetector<Weight>  * d = (DistanceDetector<Weight>*)dist_info[from].detector;
				assert(d);

				d->addUnweightedShortestPathLit(from,to,reach_var,within_steps);


		    }
	void reachesWithinDistance(int from, int to, Var reach_var, Weight distance){

					assert(from<g.nodes());


					if (dist_info[from].source<0){
						DistanceDetector<Weight> * d =new DistanceDetector<Weight>(detectors.size(), this,edge_weights,g,antig,from,drand(rnd_seed));
						detectors.push(d);
						distance_detectors.push(d);
						assert(detectors.last()->getID()==detectors.size()-1);
						dist_info[from].source=from;
						dist_info[from].detector=detectors.last();
					}

					DistanceDetector<Weight>  * d = (DistanceDetector<Weight>*)dist_info[from].detector;
					assert(d);

					d->addWeightedShortestPathLit(from,to,reach_var,distance);


			    }
	void implementConstraints(){
		if(!S->okay())
			return;
		if(opt_allpairs_percentage>=1){
			for(int i = 0;i<unimplemented_reachability_constraints.size();i++){
				ReachabilityConstraint c = unimplemented_reachability_constraints[i];
				reaches_private(c.from,c.to,c.reach_var,c.distance);
			}
			for(int i = 0;i<unimplemented_connectivity_constraints.size();i++){
				ConnectivityConstraint c = unimplemented_connectivity_constraints[i];
				connects_private(c.from,c.to,c.connect_var,c.distance);
			}
		}else if (opt_allpairs_percentage==0){
			for(int i = 0;i<unimplemented_reachability_constraints.size();i++){
				ReachabilityConstraint c = unimplemented_reachability_constraints[i];
				allpairs(c.from,c.to,c.reach_var,c.distance);
			}
			for(int i = 0;i<unimplemented_connectivity_constraints.size();i++){
				ConnectivityConstraint c = unimplemented_connectivity_constraints[i];
				allpairs_undirected(c.from,c.to,c.connect_var,c.distance);
			}
		}else{
			{
				vec<bool> seen;
				int count=0;
				seen.growTo(nNodes());
				for(int i = 0;i<unimplemented_reachability_constraints.size();i++){
							ReachabilityConstraint c = unimplemented_reachability_constraints[i];
							if(!seen[c.from]){
								seen[c.from]=true;
								count++;
							}
						}
				double frac = ((double)count)/((double)nNodes());

				if (opt_verb>0 && frac>=opt_allpairs_percentage){
					printf("Allpairs solver triggered for graph %d by percentage of source nodes: %d/%d=%f>%f\n",getGraphID() ,count,nNodes(),frac,(double)opt_allpairs_percentage);
				}

				for(int i = 0;i<unimplemented_reachability_constraints.size();i++){
					ReachabilityConstraint c = unimplemented_reachability_constraints[i];
					if(frac>=opt_allpairs_percentage)
						allpairs(c.from,c.to,c.reach_var,c.distance);
					else
						reaches_private(c.from,c.to,c.reach_var,c.distance);
				}
			}

			{
				vec<bool> seen;
				int count=0;
				seen.growTo(nNodes());
				for(int i = 0;i<unimplemented_connectivity_constraints.size();i++){
					ConnectivityConstraint c = unimplemented_connectivity_constraints[i];
							if(!seen[c.from]){
								seen[c.from]=true;
								count++;
							}
						}
				double frac = ((double)count)/((double)nNodes());

				if (opt_verb>0 && frac>=opt_allpairs_percentage){
					printf("Allpairs-undirected solver triggered for graph %d by percentage of source nodes: %d/%d=%f>%f\n",getGraphID() ,count,nNodes(),frac,(double)opt_allpairs_percentage);
				}

				for(int i = 0;i<unimplemented_connectivity_constraints.size();i++){
					ConnectivityConstraint c = unimplemented_connectivity_constraints[i];
					if(frac>=opt_allpairs_percentage)
						allpairs_undirected(c.from,c.to,c.connect_var,c.distance);
					else
						connects_private(c.from,c.to,c.connect_var,c.distance);
				}
			}
		}
		unimplemented_reachability_constraints.clear();


	}
	void allpairs_undirected(int from, int to, Var reach_var,int within_steps=-1){

	}

	void allpairs(int from, int to, Var reach_var,int within_steps=-1){
				//for now, reachesWithinSteps to be called instead

					assert(from<g.nodes());
					if(within_steps>g.nodes())
						within_steps=-1;

					if (reach_info[from].source<0){


						detectors.push(new AllPairsDetector<Weight>(detectors.size(), this,g,antig,drand(rnd_seed)));
							//reach_detectors.push(reach_detectors.last());






						assert(detectors.last()->getID()==detectors.size()-1);





						//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
						//reach_detectors.last()->source=from;

						reach_info[from].source=from;
						reach_info[from].detector=detectors.last();

						//reach_detectors.last()->within=within_steps;

					}

					AllPairsDetector<Weight> * d =(AllPairsDetector<Weight>*) reach_info[from].detector;
					assert(d);

					d->addLit(from,to,reach_var,within_steps);


			    }
	void connects_private(int from, int to, Var reach_var,int within_steps=-1){
		//for now, reachesWithinSteps to be called instead
		if(within_steps>=0 || opt_force_distance_solver){
			//reachesWithinSteps(from,to,reach_var,within_steps);
			printf("Not supported yet\n");
			exit(1);
			return;
		}


			assert(from<g.nodes());
			if(within_steps>g.nodes())
				within_steps=-1;

			if (connect_info[from].source<0){

					ConnectDetector<Weight>*rd = new ConnectDetector<Weight>(detectors.size(), this,g,antig,from,drand(rnd_seed));
					detectors.push(rd);
					connect_detectors.push(rd);

				assert(detectors.last()->getID()==detectors.size()-1);


				connect_info[from].source=from;
				connect_info[from].detector=detectors.last();

				//reach_detectors.last()->within=within_steps;

			}

			ConnectDetector<Weight> * d = (ConnectDetector<Weight>*) connect_info[from].detector;
			assert(d);
			assert(within_steps==-1);
			d->addLit(from,to,reach_var);


	}

	void reaches_private(int from, int to, Var reach_var,int within_steps=-1){
			//for now, reachesWithinSteps to be called instead
			if(within_steps>=0 || opt_force_distance_solver){
				reachesWithinSteps(from,to,reach_var,within_steps);
				return;
			}

				assert(from<g.nodes());
				if(within_steps>g.nodes())
					within_steps=-1;

				if (reach_info[from].source<0){



						ReachDetector<Weight>*rd = new ReachDetector<Weight>(detectors.size(), this,g,antig,from,drand(rnd_seed));
						detectors.push(rd);
						reach_detectors.push(rd);



					assert(detectors.last()->getID()==detectors.size()-1);


					//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
					//reach_detectors.last()->source=from;

					reach_info[from].source=from;
					reach_info[from].detector=detectors.last();

					//reach_detectors.last()->within=within_steps;

				}

				ReachDetector<Weight> * d = (ReachDetector<Weight>*) reach_info[from].detector;
				assert(d);
				assert(within_steps==-1);
				d->addLit(from,to,reach_var);


		    }

	//Undirected reachability query
	void connects(int from, int to, Var connect_var, int within_steps=-1){
		unimplemented_connectivity_constraints.push({from,to,within_steps,connect_var});
	}

	void reaches(int from, int to, Var reach_var,int within_steps=-1){
			unimplemented_reachability_constraints.push({from,to,within_steps,reach_var});
			//to allow us to alter the solving algorithm based on the number and type of constraints, we aren't implementing them here directly any more - instead,
			//we just store the constraints in this vector, then implement them later when 'implementConstraints' is called.
	    }

	void reachesAny(int from, Var firstVar,int within_steps=-1){
		for(int i = 0;i<g.nodes();i++){
			reaches(from,i,firstVar+i,within_steps);
		}
	}

	void reachesAny(int from, vec<Lit> & reachlits_out,int within_steps=-1){
		for(int i = 0;i<g.nodes();i++){
			Var reachVar = S->newVar();
			//reaches(from,i,reachVar,within_steps);
			reaches(from,i,reachVar,within_steps);
			reachlits_out.push(mkLit(reachVar,false));
		}
    }
	//v will be true if the minimum weight is <= the specified value
	void minimumSpanningTree(Var v, Weight minimum_weight){
		if(!mstDetector){
			mstDetector = new MSTDetector<Weight>(detectors.size(),this, g, antig, edge_weights,drand(rnd_seed));
			detectors.push(mstDetector);
		}
		mstDetector->addWeightLit(v,minimum_weight);
	}
	void edgeInMinimumSpanningTree(Var edgeVar, Var var){
		if(!mstDetector){
			mstDetector = new MSTDetector<Weight>(detectors.size(),this, g, antig, edge_weights,drand(rnd_seed));
			detectors.push(mstDetector);
		}
		if(!S->hasTheory(edgeVar) || (S->getTheoryID(edgeVar)!= getTheoryIndex()) || ! isEdgeVar(S->getTheoryVar(edgeVar)) ){
			fprintf(stderr,"%d is not an edge variable for theory %d! Aborting\n",edgeVar+1, getTheoryIndex());
			exit(1);
		}
		edgeVar = S->getTheoryVar(edgeVar);
		int edgeid =getEdgeID(edgeVar);
		assert(edgeid>=0);
		if(edge_list[edgeid].v==var_Undef){
			printf("MST edge constraint for undefined edge %d with variable %d, aborting.\n", edgeid,edgeVar+1);
			exit(1);
		}
		mstDetector->addTreeEdgeLit(edgeid,var);
	}
	void maxFlow(int from, int to, int max_flow, Var v){

		for (int i =0;i< flow_detectors.size();i++){
			if(flow_detectors[i]->source==from && flow_detectors[i]->target == to){
				flow_detectors[i]->addFlowLit(max_flow,v);
				return;
			}
		}
		MaxflowDetector<Weight> *f = new MaxflowDetector<Weight>(detectors.size(),this,edge_weights, g, antig,from,to,drand(rnd_seed)) ;
		flow_detectors.push(f);
		detectors.push(f);
		f->addFlowLit(max_flow,v);

	}
	void minConnectedComponents(int min_components, Var v){
		if(!component_detector){
			component_detector = new  ConnectedComponentsDetector<Weight>(detectors.size(),this, g, antig,drand(rnd_seed));
			detectors.push(component_detector);
		}
		component_detector->addConnectedComponentsLit(v,min_components);
	}
	void detectCycle(bool directed, Var v){
		if(!cycle_detector){
			cycle_detector = new  CycleDetector<Weight>(detectors.size(),this, g, antig,true,drand(rnd_seed));
			detectors.push(cycle_detector);
		}
		cycle_detector->addCycleDetectorLit(directed,v);
	}


	void addSteinerTree(const vec<std::pair<int, Var> > & terminals, int steinerTreeID){
		steiner_detectors.growTo(steinerTreeID+1);
		assert(!steiner_detectors[steinerTreeID]);
		steiner_detectors[steinerTreeID]= new SteinerDetector<Weight>(detectors.size(),this,edge_weights, g, antig,drand(rnd_seed));
		detectors.push(steiner_detectors[steinerTreeID]);
		for(int i =0;i<terminals.size();i++){
			steiner_detectors[steinerTreeID]->addTerminalNode(terminals[i].first,terminals[i].second);
		}
	}

	void addSteinerWeightConstraint(int steinerTreeID, Weight weight, Var outerVar){
		if(steinerTreeID >= steiner_detectors.size()){
			fprintf(stderr,"invalid steinerTreeID %d\n", steinerTreeID);
			exit(1);
		}
		steiner_detectors[steinerTreeID]->addWeightLit(weight,outerVar);
	}

/*	void inTerminalSet(int node, int terminalSet, Var outerVar){
		terminalSets.growTo(terminalSet+1);
		while(terminalSets[terminalSet].nodes()<node){
			terminalSets[terminalSet].addNode();
		}
		Var v= newVar(outerVar,node,true);
		Lit l = mkLit(v,false);
		if(terminalSets[terminalSet].getNodeVar(node)<0){
			terminalSets[terminalSet].setNodeVar(v);
		}
	}*/

	void printSolution(){
		if(S->model.size()==0)
			return;
						int width = sqrt(nNodes());
						if(opt_width>0){
							width=opt_width;
						}
						int height =width;
						if(opt_height>0){
							height = opt_height;
						}
						int bits = 1;
						if(opt_bits>0)
								bits=opt_bits;
						int v = 0;
						//for (int i = 0;i<w;i++){
						//	for(int j = 0;j<w;j++){
						int lasty= 0;

						int maxwidth = log10(pow(2, bits))+1; //highestbit(bits);
						if(getTheoryIndex()==0){
						for(int n = 0;n<height*width*bits;n+=bits){
							int x = n%(width*bits)/bits;
							int y = n/(width*bits);
							if(y > lasty)
								printf("\n");
		#if not defined(__MINGW32__)
								if (!opt_csv && isatty(fileno(stdout))){
		#else
								if(false){
		#endif
									unsigned long val = 0;
									for(int j = 0;j<bits;j++){
										if(S->model[n+j]==l_True){
											val = val + (1<<j);
										}
									}

									//if(val>0){
										int backcolor = 0;
										if(val>0){
											backcolor=log2(val)+1;
										}
										if(backcolor<0){
											int a=1;
										}
										int forecolor = 7;
										if(backcolor>7){
											backcolor=7;
										}
										if(backcolor==3 || backcolor==7){
											forecolor=0;
										}
										printf("\033[1;4%dm\033[1;3%dm%*lu \033[0m",backcolor,forecolor,maxwidth,val);
									//}else{
										//printf("\033[1;44m\033[1;37m%*lu \033[0m",maxwidth,val);
										//printf("\033[1;40m\033[1;30m%*lu \033[0m",maxwidth,val);
									//}
								}else if (opt_csv){
									unsigned long val = 0;
									for(int j = 0;j<bits;j++){
										if(S->model[n+j]==l_True){
											val = val + (1<<j);
										}
									}
									printf("%*lu",maxwidth,val);
									if (x<width-1){
										printf(",");
									}
								}else{
									unsigned long val = 0;
									for(int j = 0;j<bits;j++){
										if(S->model[n+j]==l_True){
											val = val + (1<<j);
										}
									}
									printf(" %*lu ",maxwidth,val);

						/*			if(S.model[n]==l_True)
										printf(" 1");
									else
										printf(" 0");*/
								}

							lasty=y;
						}
						printf("\n\n");
						}
						if(opt_check_solution){
									if(!check_solved()){
										fprintf(stderr,"Error! Solution doesn't satisfy graph properties!\n");
										exit(3);
									}
								}

						if(opt_print_reach){
						 v = 0;
						//for (int i = 0;i<w;i++){
						//	for(int j = 0;j<w;j++){
					/*	 lasty= 0;
						for(int n = 0;n<nNodes();n++){
							int x = n%width;
							int y = n/width;
							if(y > lasty)
								printf("\n");
		#if not defined(__MINGW32__)
								if (isatty(fileno(stdout))){
		#else
								if(false){
		#endif

									if(S.model[n]==l_True)
										printf("\033[1;42m\033[1;37m%3d\033[0m",n);
									else
										printf("\033[1;44m\033[1;37m%3d\033[0m",n);
								}else{

									if(S.model[n]==l_True)
										printf(" 1");
									else
										printf(" 0");
								}

							lasty=y;
						}

						printf("\n");printf("\n");
						*/



							printf("Theory %d\n", getTheoryIndex());

							int nnodes = nNodes();

							int maxw = log10(nNodes() )+1; //highestbit(bits);

							{

								for(int r = 0;r<reach_detectors.size();r++){

									int width = sqrt(nNodes());
									if(opt_width>0){
											width=opt_width;
										}
										int height =width;
										if(opt_height>0){
											height = opt_height;
										}
									int lasty= 0;
									int extra =  nNodes() % width ? (width- nNodes() % width ):0;
									for(int n = 0;n<nNodes();n++){
										int x = n%width;

										int y = (n + extra )/width;
										if(y > lasty)
											printf("\n");

										int v =var( reach_detectors[r]->reach_lits[n]);
		#if not defined(__MINGW32__)
										if (isatty(fileno(stdout)))
		#else
										if(false)
		#endif
										{
												if(value(v)==l_True)
													printf("\033[1;42m\033[1;37m%4d\033[0m", v+1);
												else
													printf("\033[1;44m\033[1;37m%4d\033[0m",v+1);
											}else{

												if(value(v)==l_True)
													printf(" 1");
												else
													printf(" 0");
											}

											lasty=y;
										}
										printf("\n");
									}



									//drawFull();

									assert(dbg_solved());
								}

							{
										for(int r = 0;r<distance_detectors.size();r++){

													int width = sqrt(nNodes());
													if(opt_width>0){
															width=opt_width;
														}
														int height =width;
														if(opt_height>0){
															height = opt_height;
														}
													int lasty= 0;
													int extra =  nNodes() % width ? (width- nNodes() % width ):0;
													for(int n = 0;n<nNodes();n++){
														int x = n%width;

														int y = (n + extra )/width;
														if(y > lasty)
															printf("\n");

														int d = distance_detectors[r]->positive_reach_detector->distance(n);
														printf("%*d ",maxw,d);


															lasty=y;
														}
														printf("\n");


														distance_detectors[r]->printSolution();

													}


									}

							{
								for(int r = 0;r<flow_detectors.size();r++){

								int width = sqrt(nNodes());
								if(opt_width>0){
										width=opt_width;
									}
									int height =width;
									if(opt_height>0){
										height = opt_height;
									}
									int lasty= 0;
									int extra =  nNodes() % width ? (width- nNodes() % width ):0;
									for(int n = 0;n<nNodes();n++){
										int x = n%width;

										int y = (n )/width;
										if(y > lasty)
											printf("\n");
										Weight total_flow = 0;
										for(int e = 0;e<g.edges();e++){
											if(g.getEdge(e).to==n){
												total_flow+=flow_detectors[r]->positive_detector->getEdgeFlow(e);

											}
										}

										//printf("%*d ",maxw,total_flow);
										std::cout<<total_flow<<" ";
											lasty=y;
										}
										printf("\n");

										for(int n = 0;n<nNodes();n++){
												int x = n%width;

												int y = (n )/width;
												if(y > lasty)
													printf("\n");
												int total_flow = 0;
												for(int e = 0;e<g.edges();e++){
													if(g.getEdge(e).to==n){
														Weight flow = flow_detectors[r]->positive_detector->getEdgeFlow(e);
														if(flow>0){
															printf("flow (%d,%d) %d to %d\n",x,y, g.getEdge(e).from,g.getEdge(e).to);
														}
													}
												}


												lasty=y;
											}
									printf("\n");
								}
							}

								if(mstDetector){
										Weight min_weight = mstDetector->positive_reach_detector->weight();
										//printf("Min Spanning Tree Weight: %d\n",min_weight);
										std::cout<<"Min Spanning Tree Weight: " << min_weight <<"\n";
										int width = sqrt(nNodes());
										if(opt_width>0){
												width=opt_width;
											}
											int height =width;
											if(opt_height>0){
												height = opt_height;
											}
										int lasty= 0;
										vec<bool> down_edge;
										int extra =  nNodes() % width ? (width- nNodes() % width ):0;
										for(int n = 0;n<nNodes();n++){
											int x = n%width;

											int y = (n + extra )/width;
											if(y > lasty){
												printf("\n");

												for(int i = 0;i<down_edge.size();i++){
													if(down_edge[i]){
														printf("|");
													}else{
														printf(" ");
													}
													printf(" ");
												}
												down_edge.clear();
												printf("\n");
											}
											printf("*");
											if(x<width-1){
												int edge_left = getEdgeID(n,n+1);
												Var edge_var = edge_list[edge_left].v;
												if(value(edge_var)==l_True &&  mstDetector->positive_reach_detector->edgeInTree(edge_left)){
													printf("-");
												}else{
													printf(" ");
												}
											}

											if(y<height-1){
													int edge_down = getEdgeID(n,n+width);
													Var edge_var = edge_list[edge_down].v;
													bool in_tree = mstDetector->positive_reach_detector->edgeInTree(edge_down);
													if(value(edge_var)==l_True &&  in_tree){
														down_edge.push(true);
													}else{
														down_edge.push(false);
													}
												}

												lasty=y;
											}
											printf("\n");
										}

								if(component_detector){
									int numComponents = component_detector->positive_component_detector->numComponents();
									printf("Number of connected components is: %d\n",numComponents);

								}





						}
		/*        		for(int r = 0;r<reach_detectors.size();r++){

							int width = sqrt(nNodes());
							int lasty= 0;
							int extra =  nNodes() % width ? (width- nNodes() % width ):0;
							for(int n = 0;n<nNodes();n++){
								int x = n%width;

								int y = (n + extra )/width;

								int v =var( reach_detectors[r]->reach_lits[n]);
								if(v==306){
									int a =1;
								}
								if(S.model[v]==l_True){
									assert(S.value(v)==l_True);
									int node = reach_detectors[r]->getNode(v);
									reach_detectors[r]->positive_reach_detector->dbg_path(node);
									int  b=1;

								}


								lasty=y;
							}
							printf("\n");
						}*/
	}

};

};

#endif /* DGRAPH_H_ */
