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

#ifndef FSM_THEORY_H_
#define FSM_THEORY_H_

#include "utils/System.h"
#include "core/Theory.h"
#include "dgl/DynamicGraph.h"
#include "DynamicFSM.h"
#include "FSMReachDetector.h"

#include "core/SolverTypes.h"
#include "mtl/Map.h"


#include "utils/System.h"
#include "core/Solver.h"

#include <vector>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sstream>

using namespace dgl;
namespace Monosat {

struct EdgeLabel{


};


class FSMTheorySolver: public Theory {
public:
	
	double rnd_seed;

private:
	Solver * S;
	int local_q = 0;
public:
	int id;
	bool all_edges_unit = true;
	vec<lbool> assigns;

	int n_labels=1;//number of transition labels. Transition labels start from 0 (which is the non-consuming epsilon transition) and go to n_labels-1.


	vec<EdgeLabel> edge_labels;

	DynamicFSM<8> g_under;
	DynamicFSM<8> g_over;

	/**
	 * The cutgraph is (optionally) used for conflict analysis by some graph theories.
	 * It has two edges for every edge in the real graph (with indices edgeID*2 and edgeID*2+1).
	 * If edge ID is assigned to FALSE, then edge ID*2 is assigned to enabled in the cutgraph, and
	 * edge ID*2+1 is disabled.
	 * Otherwise, if edge ID is unassigned or true, then edge ID*2 is disabled in the cutgraph, and
	 * edge ID*2+1 is enabled.
	 */
	//DynamicGraph cutGraph;

	
	vec<Assignment> trail;
	vec<int> trail_lim;


public:
	vec<FSMDetector*> detectors;
	vec<FSMReachDetector*> reach_detectors;


	vec<int> marker_map;



	bool requiresPropagation = true;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;
	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {
		int isEdge :1;
		int occursPositive :1;
		int occursNegative :1;
		int detector_edge :29;	//the detector this variable belongs to, or its edge number, if it is an edge variable
		Var solverVar;
	};

	vec<VarData> vars;
	int theory_index = 0;
public:
	
	double mctime = 0;
	double reachtime = 0;
	double unreachtime = 0;
	double pathtime = 0;
	double propagationtime = 0;
	long stats_propagations = 0;
	long stats_num_conflicts = 0;
	long stats_decisions = 0;
	long stats_num_reasons = 0;

	double reachupdatetime = 0;
	double unreachupdatetime = 0;
	double stats_initial_propagation_time = 0;
	double stats_decision_time = 0;
	double stats_reason_initial_time = 0;
	double stats_reason_time = 0;
	long num_learnt_paths = 0;
	long learnt_path_clause_length = 0;
	long num_learnt_cuts = 0;
	long learnt_cut_clause_length = 0;
	long stats_pure_skipped = 0;
	long stats_mc_calls = 0;
	long stats_propagations_skipped = 0;
	vec<Lit> reach_cut;



	FSMTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id){


		rnd_seed = opt_random_seed;
	}
	
	void printStats(int detailLevel) {
		if (detailLevel > 0) {
			for (FSMDetector * d : detectors)
				d->printStats();
		}
		
		printf("FSM %d stats:\n", getGraphID());

		fflush(stdout);
	}
	
	void writeTheoryWitness(std::ostream& write_to) {
		
		for (FSMDetector * d : detectors) {
			write_to << "Graph " << this->getGraphID() << ", detector " << d->getID() << ":\n";
			d->printSolution(write_to);
		}
	}
	
	inline int getTheoryIndex() {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}
	inline int getGraphID() {
		return id;
	}
	inline bool isEdgeVar(Var v) {
		assert(v < vars.size());
		return vars[v].isEdge;
	}
	inline int getEdgeID(Var v) {
		assert(isEdgeVar(v));
		return vars[v].detector_edge;
	}
	inline int getDetector(Var v) {
		assert(!isEdgeVar(v));
		return vars[v].detector_edge;
	}
	
	inline Var getEdgeVar(int edgeID) {
		Var v = edge_list[edgeID].v;
		assert(v < vars.size());
		assert(vars[v].isEdge);
		return v;
	}
	
	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1, o2);
		S->addClause(o1, ~o2);
	}
	void makeEqualInSolver(Lit l1, Lit l2) {
		S->addClause(~l1, l2);
		S->addClause(l1, ~l2);
	}
	void addClause(Lit l1) {
		Lit o1 = toSolver(l1);
		S->addClause(o1);
	}
	void addClause(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(o1, o2);
	}
	void addClause(Lit l1, Lit l2, Lit l3) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		Lit o3 = toSolver(l3);
		S->addClause(o1, o2, o3);
	}
	void addClause(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		S->addClause(tmp_clause);
	}
	void addClauseSafely(vec<Lit> & c) {
		tmp_clause.clear();
		c.copyTo(tmp_clause);
		toSolver(tmp_clause);
		
		S->addClauseSafely(tmp_clause);
	}
	
	Var newVar(int forDetector = -1, bool connectToTheory = false) {
		Var s = S->newVar();
		return newVar(s, forDetector, false, connectToTheory);
	}
	Var newVar(Var solverVar, int detector, bool isEdge = false, bool connectToTheory = true) {
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();
		vars.push();
		vars[v].isEdge = isEdge;
		vars[v].detector_edge = detector;
		vars[v].solverVar = solverVar;
		assigns.push(l_Undef);
		if (connectToTheory) {
			S->setTheoryVar(solverVar, getTheoryIndex(), v);
			assert(toSolver(v) == solverVar);
		}
		if (!isEdge && detector >= 0)
			detectors[detector]->addVar(v);
		return v;
	}
	inline int level(Var v) {
		return S->level(toSolver(v));
	}
	inline int decisionLevel() {
		return trail_lim.size(); //S->decisionLevel();
	}
	inline int nVars() const {
		return vars.size(); //S->nVars();
	}
	inline Var toSolver(Var v) {
		//return v;
		assert(v < vars.size());
		//assert(S->hasTheory(vars[v].solverVar));
		//assert(S->getTheoryVar(vars[v].solverVar)==v);
		return vars[v].solverVar;
	}
	
	inline Lit toSolver(Lit l) {
		//assert(S->hasTheory(vars[var(l)].solverVar));
		//assert(S->getTheoryVar(vars[var(l)].solverVar)==var(l));
		return mkLit(vars[var(l)].solverVar, sign(l));
	}
	
	void toSolver(vec<Lit> & c) {
		for (int i = 0; i < c.size(); i++) {
			c[i] = toSolver(c[i]);
		}
	}
	
	inline lbool value(Var v) {
		if (assigns[v] != l_Undef)
			assert(S->value(toSolver(v)) == assigns[v]);
		
		return assigns[v]; //S->value(toSolver(v));
	}
	inline lbool value(Lit l) {
		if (assigns[var(l)] != l_Undef) {
			assert(S->value(toSolver(l)) == (assigns[var(l)] ^ sign(l)));
		}
		return assigns[var(l)] ^ sign(l);
	}
	inline lbool dbg_value(Var v) {
		return S->value(toSolver(v));
	}
	inline lbool dbg_value(Lit l) {
		return S->value(toSolver(l));
	}
	inline bool enqueue(Lit l, CRef reason) {
		assert(assigns[var(l)]==l_Undef);
		
		Lit sl = toSolver(l);
		if (S->enqueue(sl, reason)) {
			enqueueTheory(l);
			return true;
		} else {
			return false;
		}
	}
	
	~FSMTheorySolver() {
	}
	;
	int newNode() {
		
		inv_adj.push();
		undirected_adj.push();

		g_over.addNode();
		cutGraph.addNode();
		seen.growTo(nNodes());
		
		return g_under.addNode();
	}
	void newNodes(int n) {
		for (int i = 0; i < n; i++)
			newNode();
	}
	int nNodes() {
		return g_under.nodes();
	}
	bool isNode(int n) {
		return n >= 0 && n < nNodes();
	}
	
	bool dbg_propgation(Lit l) {
#ifndef NDEBUG
		static vec<Lit> c;
		c.clear();
		for (int i = 0; i < S->trail.size(); i++) {
			if (!S->hasTheory(S->trail[i]) || S->getTheoryID(S->trail[i]) != getTheoryIndex())
				continue;
			Lit l = S->getTheoryLit(S->trail[i]);
			Var v = var(l);
			if (isEdgeVar(v)) {
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
	
	void dbg_sync_reachability() {
		

		
	}
	void dbg_sync() {
#ifndef NDEBUG
		if (opt_conflict_min_cut) {
			for (int i = 0; i < edge_list.size(); i++) {
				if (value(edge_list[i].v) == l_False) {
					assert(cutGraph.edgeEnabled(i * 2));
					assert(!cutGraph.edgeEnabled(i * 2 + 1));
				} else {
					assert(!cutGraph.edgeEnabled(i * 2));
					assert(cutGraph.edgeEnabled(i * 2 + 1));
				}
			}
		}
#endif

	}
	void dbg_full_sync() {
#ifdef DEBUG_GRAPH
		dbg_sync();
		for(int i =0;i<edge_list.size();i++) {

			if(edge_list[i].edgeID>=0 && g_under.edgeEnabled(i)) {
				assert(value(edge_list[i].v)==l_True);
			} else if (edge_list[i].edgeID>=0 && !g_over.edgeEnabled(i)) {
				assert(value(edge_list[i].v)==l_False);
			}
		}

#endif
	}
	
	void backtrackUntil(int level) {
		static int it = 0;
		
		bool changed = false;
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if (e.isEdge) {
					assert(dbg_value(e.var)==l_Undef);
					int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
							
					if (e.assign) {
						g_under.disableEdge(e.from, e.to, edge_num);
						assert(!cutGraph.edgeEnabled(edge_num * 2));
					} else {
						g_over.enableEdge(e.from, e.to, edge_num);
						assert(g_over.hasEdge(e.from, e.to));
						if (opt_conflict_min_cut) {
							assert(cutGraph.edgeEnabled(edge_num * 2));
							cutGraph.disableEdge(e.from, e.to, edge_num * 2);
							assert(!cutGraph.edgeEnabled(edge_num * 2 + 1));
							cutGraph.enableEdge(e.from, e.to, edge_num * 2 + 1);
						}
					}
				} else {
					//This is a reachability literal				  
					detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
				}
				assigns[e.var] = l_Undef;
				changed = true;
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			
			if (changed) {
				requiresPropagation = true;
				/*				g.markChanged();
				 antig.markChanged();
				 cutGraph.markChanged();*/
			}
			
			for (FSMDetector * d : detectors) {
				d->backtrack(level);
			}
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
		
	}
	;

	Lit decideTheory() {
		if (!opt_decide_theories)
			return lit_Undef;
		double start = rtime(1);
		
		dbg_full_sync();
		for (int i = 0; i < detectors.size(); i++) {
			FSMDetector * r = detectors[i];
			Lit l = r->decide(decisionLevel());
			if (l != lit_Undef) {
				stats_decisions++;
				r->stats_decisions++;
				stats_decision_time += rtime(1) - start;
				return toSolver(l);
			}
		}
		stats_decision_time += rtime(1) - start;
		return lit_Undef;
	}
	
	void backtrackUntil(Lit p) {
		//need to remove and add edges in the two graphs accordingly.
		
		int i = trail.size() - 1;
		for (; i >= 0; i--) {
			Assignment e = trail[i];
			if (e.isEdge) {
				int edge_num = getEdgeID(e.var); //e.var-min_edge_var;
				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
					g_under.disableEdge(e.from, e.to, edge_num);
					
					assert(!cutGraph.edgeEnabled(edge_num * 2));
				} else {
					g_over.enableEdge(e.from, e.to, edge_num);
					assert(g_over.hasEdge(e.from, e.to));
					if (opt_conflict_min_cut) {
						assert(cutGraph.edgeEnabled(edge_num * 2));
						cutGraph.disableEdge(e.from, e.to, edge_num * 2);
						assert(!cutGraph.edgeEnabled(edge_num * 2 + 1));
						cutGraph.enableEdge(e.from, e.to, edge_num * 2 + 1);
					}
				}
			} else {
				if (var(p) == e.var) {
					assert(sign(p) != e.assign);
					break;
				}
				assigns[e.var] = l_Undef;
				detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
			}
		}
		
		trail.shrink(trail.size() - (i + 1));
		//if(i>0){
		requiresPropagation = true;
		/*			g.markChanged();
		 antig.markChanged();
		 cutGraph.markChanged();*/
		//}
		for (FSMDetector * d : detectors) {
			d->backtrack(this->decisionLevel());
		}
		//while(trail_lim.size() && trail_lim.last()>=trail.size())
		//	trail_lim.pop();
		dbg_sync();
		/*		for(int i = 0;i<reach_detectors.size();i++){
		 if(reach_detectors[i]->positive_reach_detector)
		 reach_detectors[i]->positive_reach_detector->update();
		 if(reach_detectors[i]->negative_reach_detector)
		 reach_detectors[i]->negative_reach_detector->update();
		 }*/
	}
	;

	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;

	void buildReason(Lit p, vec<Lit> & reason) {
		CRef marker = S->reason(var(toSolver(p)));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef - marker;
		int d = marker_map[pos];
		//double initial_start = rtime(1);
		double start = rtime(1);
		backtrackUntil(p);
		
		assert(d < detectors.size());
		detectors[d]->buildReason(p, reason, marker);
		toSolver(reason);
		double finish = rtime(1);
		stats_reason_time += finish - start;
		stats_num_reasons++;
		//stats_reason_initial_time+=start-initial_start;
		
	}
	
	bool dbg_reachable(int from, int to, bool undirected = false) {
#ifdef DEBUG_DIJKSTRA
		
		if(undirected) {
			UnweightedDijkstra<Reach::NullStatus, true> d(from,g_under);
			d.update();
			return d.connected(to);
		} else {
			UnweightedDijkstra<> d(from,g_under);
			d.update();
			return d.connected(to);
		}

#else
		return true;
#endif
	}
	
	bool dbg_notreachable(int from, int to, bool undirected = false) {
		
#ifndef NDEBUG
		//drawFull(from,to);
		
		DynamicGraph g;
		for (int i = 0; i < nNodes(); i++) {
			g.addNode();
		}
		
		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge e = edge_list[i];
			if (value(e.v) != l_False) {
				g.addEdge(e.from, e.to);
			}
		}
		if (undirected) {
			UnweightedDijkstra<Reach::NullStatus, true> d(from, g);
			
			return !d.connected(to);
		} else {
			UnweightedDijkstra<Reach::NullStatus, false> d(from, g);
			
			return !d.connected(to);
		}
#endif
		return true;
		
	}
	
	bool dbg_graphsUpToDate() {
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++) {
			if(edge_list[i].v<0)
			continue;
			Edge e = edge_list[i];
			lbool val = value(e.v);

			if(val==l_True || val==l_Undef) {
				assert(g_over.edgeEnabled(i));
				//assert(antig.hasEdge(e.from,e.to));
			} else {
				assert(!g_over.edgeEnabled(i));
				//assert(!antig.hasEdge(e.from,e.to));
			}
			if(val==l_True) {
				assert(g_under.edgeEnabled(i));
				//assert(g.hasEdge(e.from,e.to));
			} else {
				assert(!g_under.edgeEnabled(i));
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

	void preprocess() {
		for (int i = 0; i < detectors.size(); i++) {
			detectors[i]->preprocess();
		}
	}
	void setLiteralOccurs(Lit l, bool occurs) {
		if (isEdgeVar(var(l))) {
			//don't do anything
		} else {
			//this is a graph property detector var
			if (!sign(l) && vars[var(l)].occursPositive != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
			else if (sign(l) && vars[var(l)].occursNegative != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
		}
		
	}
	
	void enqueueTheory(Lit l) {
		Var v = var(l);
		
		int lev = level(v);
		
		assert(decisionLevel() <= lev);
		
		while (lev > trail_lim.size()) {
			newDecisionLevel();
		}
		
		if (assigns[var(l)] != l_Undef) {
			return;			//this is already enqueued.
		}
		assert(assigns[var(l)]==l_Undef);
		assigns[var(l)] = sign(l) ? l_False : l_True;
		requiresPropagation = true;
		//printf("enqueue %d\n", dimacs(l));
		
#ifndef NDEBUG
		{
			for (int i = 0; i < trail.size(); i++) {
				assert(trail[i].var != v);
			}
		}
#endif
		
#ifdef RECORD
		if (g_under.outfile) {
			fprintf(g_under.outfile, "enqueue %d\n", dimacs(l));
			
			fprintf(g_under.outfile, "\n");
			fflush(g_under.outfile);
		}
		if (g_over.outfile) {
			fprintf(g_over.outfile, "enqueue %d\n", dimacs(l));
			fprintf(g_over.outfile, "\n");
			fflush(g_over.outfile);
		}
#endif
		
		if (isEdgeVar(var(l))) {
			
			//this is an edge assignment
			int edge_num = getEdgeID(var(l)); //v-min_edge_var;
			assert(edge_list[edge_num].v == var(l));
			
			int from = edge_list[edge_num].from;
			int to = edge_list[edge_num].to;
			trail.push( { true, !sign(l), from, to, v });
			
			Assignment e = trail.last();
			assert(e.from == from);
			assert(e.to == to);
			
			if (!sign(l)) {
				g_under.enableEdge(from, to, edge_num);
			} else {
				g_over.disableEdge(from, to, edge_num);
				if (opt_conflict_min_cut) {
					assert(cutGraph.edgeEnabled(edge_num * 2 + 1));
					assert(!cutGraph.edgeEnabled(edge_num * 2));
					cutGraph.enableEdge(e.from, e.to, edge_num * 2);
					cutGraph.disableEdge(e.from, e.to, edge_num * 2 + 1);
				}
			}
			
		} else {
			
			trail.push( { false, !sign(l), 0, 0, v });
			//this is an assignment to a non-edge atom. (eg, a reachability assertion)
			detectors[getDetector(var(l))]->assign(l);
		}
		
	}
	;
	bool propagateTheory(vec<Lit> & conflict) {
		static int itp = 0;
		if (++itp == 62279) {
			int a = 1;
		}
		stats_propagations++;
		dbg_sync();
		if (!requiresPropagation) {
			stats_propagations_skipped++;
			assert(dbg_graphsUpToDate());
			return true;
		}
		
		bool any_change = false;
		double startproptime = rtime(1);
		//static vec<int> detectors_to_check;
		
		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.
		
		//At level 0, need to propagate constant reaches/source nodes/edges...
		
		//stats_initial_propagation_time += rtime(1) - startproptime;
		dbg_sync();
		assert(dbg_graphsUpToDate());
		
		for (int d = 0; d < detectors.size(); d++) {
			assert(conflict.size() == 0);
			bool r = detectors[d]->propagate(conflict);
			if (!r) {
				stats_num_conflicts++;
				toSolver(conflict);
				propagationtime += rtime(1) - startproptime;
				return false;
			}
		}
		
		dbg_full_sync();
		
		requiresPropagation = false;
		g_under.clearChanged();
		g_over.clearChanged();
		cutGraph.clearChanged();
		
		g_under.clearHistory();
		g_over.clearHistory();
		cutGraph.clearHistory();
		//detectors_to_check.clear();
		
		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;
		dbg_sync();
		dbg_sync_reachability();
		return true;
	}
	;

	bool solveTheory(vec<Lit> & conflict) {
		requiresPropagation = true;		//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}
	;

	void drawFull(int from, int to) {
		printf("digraph{\n");
		for (int i = 0; i < nNodes(); i++) {
			if (i == from) {
				printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
			} else if (i == to) {
				printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
			} else
				printf("n%d\n", i);
		}
		
		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge & e = edge_list[i];
			const char * s = "black";
			if (value(e.v) == l_True)
				s = "blue";
			else if (value(e.v) == l_False)
				s = "red";
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from, e.to, e.v, s);
		}
		
		printf("}\n");
	}
	void drawFull() {
		printf("digraph{\n");
		for (int i = 0; i < nNodes(); i++) {
			printf("n%d\n", i);
		}
		
		for (int i = 0; i < edge_list.size(); i++) {
			Edge & e = edge_list[i];
			const char * s = "black";
			if (value(e.v) == l_True)
				s = "blue";
			else if (value(e.v) == l_False)
				s = "red";
			else {
				int a = 1;
			}
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from, e.to, e.v, s);
		}
		
		printf("}\n");
	}
	
	bool check_solved() {
		if (opt_print_graph) {
			drawFull();
		}
		for (int i = 0; i < edge_list.size(); i++) {
			if (edge_list[i].v < 0)
				continue;
			Edge & e = edge_list[i];
			lbool val = value(e.v);
			if (val == l_Undef) {
				return false;
			}
			
			if (val == l_True) {
				/*	if(!g.hasEdge(e.from,e.to)){
				 return false;
				 }
				 if(!antig.hasEdge(e.from,e.to)){
				 return false;
				 }*/
				if (!g_under.edgeEnabled(e.edgeID)) {
					return false;
				}
				if (!g_over.edgeEnabled(e.edgeID)) {
					return false;
				}
			} else {
				/*if(g.hasEdge(e.from,e.to)){
				 return false;
				 }*/
				if (g_under.edgeEnabled(e.edgeID)) {
					return false;
				}
				if (g_over.edgeEnabled(e.edgeID)) {
					return false;
				}
				/*if(antig.hasEdge(e.from,e.to)){
				 return false;
				 }*/

			}
			
		}
		for (int i = 0; i < detectors.size(); i++) {
			if (!detectors[i]->checkSatisfied()) {
				return false;
			}
		}
		return true;
	}
	
	bool dbg_solved() {
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++) {
			if(edge_list[i].v<0)
			continue;
			Edge & e = edge_list[i];
			lbool val = value(e.v);
			assert(val!=l_Undef);

			if(val==l_True) {
				assert(g_under.hasEdge(e.from,e.to));
				assert(g_over.hasEdge(e.from,e.to));
			} else {
				assert(!g_under.hasEdge(e.from,e.to));
				assert(!g_over.hasEdge(e.from,e.to));
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
	
	void drawCurrent() {
		
	}
	int nEdges() {
		return edge_list.size();
	}
	CRef newReasonMarker(int detectorID) {
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef - reasonMarker;
		marker_map.growTo(mnum + 1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}
	
	Lit newEdge(int from, int to, Var outerVar = var_Undef) {
		assert(outerVar!=var_Undef);
		/*	if(outerVar==var_Undef)
		 outerVar = S->newVar();*/

		all_edges_unit &= (weight == 1);
		int index = edge_list.size();
		edge_list.push();
		Var v = newVar(outerVar, index, true);
		
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
		undirected_adj[to].push( { v, outerVar, from, to, index });
		undirected_adj[from].push( { v, outerVar, to, from, index });
		inv_adj[to].push( { v, outerVar, from, to, index });
		
		//num_edges++;
		edge_list[index].v = v;
		edge_list[index].outerVar = outerVar;
		edge_list[index].from = from;
		edge_list[index].to = to;
		edge_list[index].edgeID = index;
		//edge_list[index].weight=weight;
		if (edge_labels.size() <= index) {
			edge_labels.resize(index + 1);
		}
		edge_labels[index] = weight;
		
		//edges[from][to]= {v,outerVar,from,to,index};
		g_under.addEdge(from, to, index);
		g_under.disableEdge(from, to, index);
		g_over.addEdge(from, to, index);
		cutGraph.addEdge(from, to, index * 2);
		cutGraph.addEdge(from, to, index * 2 + 1);
		cutGraph.disableEdge(from, to, index * 2);
		
#ifdef RECORD
		if (g_under.outfile) {
			std::stringstream wt;
			wt << weight;
			fprintf(g_under.outfile, "edge_weight %d %s\n", index, wt.str().c_str());
			fflush(g_under.outfile);
		}
		if (g_over.outfile) {
			std::stringstream wt;
			wt << weight;
			fprintf(g_over.outfile, "edge_weight %d %s\n", index, wt.str().c_str());
			fflush(g_over.outfile);
		}
		if (cutGraph.outfile) {
			
			fprintf(cutGraph.outfile, "edge_weight %d %d\n", index * 2, 1);
			fprintf(cutGraph.outfile, "edge_weight %d %d\n", index * 2 + 1, 0xFFFF);
			fflush(cutGraph.outfile);
		}
#endif
		return mkLit(v, false);
	}

	void printSolution() {
		
		for (auto * d : detectors) {
			assert(d);
			d->printSolution();
		}
	}
	
};

}
;

#endif /* DGRAPH_H_ */
