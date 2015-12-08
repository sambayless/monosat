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

#ifndef GEOMETRY_THEORY_H_
#define GEOMETRY_THEORY_H_

#include "utils/System.h"
#include "core/Theory.h"
#include "mtl/Map.h"
#include "GeometryTypes.h"
#include "utils/System.h"
#include "core/Solver.h"
#include "core/Config.h"
#include "GeometryDetector.h"
#include "PointSet.h"
#include "ConvexHullDetector.h"
#include "ConvexHullCollisionDetector.h"

#ifndef NDEBUG
#include <cstdio>
#endif

using namespace Monosat;

template<unsigned int D, class T = double>
class GeometryTheorySolver: public Theory {
	
public:
	
	double rnd_seed;

	int local_q=0;
private:
	Solver * S;
public:
	int id;

	vec<lbool> assigns;

	std::vector<PointSet<D, T> > under_sets;
	std::vector<PointSet<D, T> > over_sets;
	struct Assignment {
		bool isPoint :1;
		bool assign :1;
		int pointID;
		Var var;
		Assignment(bool isPoint, bool _assign, int pointID, Var v) :
				isPoint(isPoint), assign(_assign), pointID(pointID), var(v) {
			
		}
		Assignment() :
				isPoint(false), assign(false), pointID(0), var(var_Undef) {
			
		}
	};
	vec<Assignment> trail;
	vec<int> trail_lim;

	std::vector<std::vector<GeometryDetector*>> pointsetDetectors;
	std::vector<GeometryDetector*> detectors;
	std::vector<GeometryDetector*> detectors_rnd;
	std::vector<ConvexHullDetector<D, T>*> convexHullDetectors;
	//std::vector<GeometricSteinerDetector<D,T>*> steinerTreeDetectors;
	ConvexHullCollisionDetector<D, T>* collisionDetector = nullptr;

	vec<int> marker_map;

	//bool anyRequiresPropagation=true;
	vec<bool> requiresPropagation;
	vec<int> toPropagate;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;

	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData {
		int isPoint :1;
		int occursPositive :1;
		int occursNegative :1;
		int detector_point :29; //the detector this variable belongs to, or its point number, if it is an point variable
		Var solverVar;
	};
	struct PointData {
		int id;
		int pointset;
		int pointset_index;
		Var var;
		Point<D, T> point;
	};
	std::vector<VarData> vars;
	std::vector<PointData> points;
	int theory_index=0;
public:
	
	double mctime=0;
	double reachtime=0;
	double unreachtime=0;
	double pathtime=0;
	double propagationtime=0;
	long stats_propagations=0;
	long stats_num_conflicts=0;
	long stats_decisions=0;
	long stats_num_reasons=0;

	double reachupdatetime=0;
	double unreachupdatetime=0;
	double stats_initial_propagation_time=0;
	double stats_decision_time=0;
	double stats_reason_initial_time=0;
	double stats_reason_time=0;
	long num_learnt_paths=0;
	long learnt_path_clause_length=0;
	long num_learnt_cuts=0;
	long learnt_cut_clause_length=0;
	long stats_pure_skipped=0;
	long stats_mc_calls=0;
	long stats_propagations_skipped=0;
	vec<Lit> reach_cut;

	GeometryTheorySolver(Solver * S_, int _id = -1) :
			S(S_), id(_id) {
		

		rnd_seed = opt_random_seed;
		
	}
	
	~GeometryTheorySolver() {
		for (auto * d : detectors) {
			delete (d);
		}
	}
	
	void printStats(int detailLevel) {
		if (detailLevel > 0) {
			for (GeometryDetector * d : detectors)
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

		printf("Propagations: %ld (%f s, avg: %f s, %ld skipped)\n", stats_propagations, propagationtime,
				(propagationtime) / ((double) stats_propagations + 1), stats_propagations_skipped);
		printf("Decisions: %ld (%f s, avg: %f s)\n", stats_decisions, stats_decision_time,
				(stats_decision_time) / ((double) stats_decisions + 1));
		printf("Conflicts: %ld\n", stats_num_conflicts);
		printf("Reasons: %ld (%f s, avg: %f s)\n", stats_num_reasons, stats_reason_time,
				(stats_reason_time) / ((double) stats_num_reasons + 1));
		fflush(stdout);
	}
	inline int nPointSets() const {
		return under_sets.size();
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
	inline bool isPointVar(Var v) {
		assert(v < vars.size());
		return vars[v].isPoint;
	}
	inline int getPointID(Var v) {
		assert(isPointVar(v));
		return vars[v].detector_point;
	}
	
	inline int getPointsetIndex(int pointID) {
		return points[pointID].pointset_index;
	}
	
	inline int getPointset(int pointID) {
		return points[pointID].pointset;
	}
	
	inline int getDetector(Var v) {
		assert(!isPointVar(v));
		return vars[v].detector_point;
	}
	
	inline Var getPointVar(int pointID) {
		Var v = points[pointID].var;
		assert(v < vars.size());
		assert(vars[v].isPoint);
		return v;
	}
	inline Var getPointVar(int pointset, int pointsetIndex) {
		int pointID = under_sets[pointset][pointsetIndex].getID();
		Var v = points[pointID].var;
		assert(v < vars.size());
		assert(vars[v].isPoint);
		return v;
	}
	
	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1, o2);
		S->addClause(o1, ~o2);
	}
	//Both l1 and l2 are SOLVER vars, not theory vars!
	void makeEqualInSolver(Lit l1, Lit l2) {
		S->addClause(~l1, l2);
		S->addClause(l1, ~l2);
	}
	void makeTrueInSolver(Lit l1) {
		S->addClause(l1);
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
	Var newVar(Var solverVar, int detector, bool isPoint = false, bool connectToTheory = true) {
		while (S->nVars() <= solverVar)
			S->newVar();
		Var v = vars.size();
		vars.push_back(VarData());
		vars[v].isPoint = isPoint;
		vars[v].detector_point = detector;
		vars[v].solverVar = solverVar;
		assigns.push(l_Undef);
		if (connectToTheory) {
			S->setTheoryVar(solverVar, getTheoryIndex(), v);
			assert(toSolver(v) == solverVar);
		}
		if (!isPoint && detector >= 0)
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
		return assigns[var(l)] ^ sign(l); //S->value(toSolver(l));
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
	
	bool dbg_propgation(Lit l) {
#ifndef NDEBUG
		static vec<Lit> c;
		c.clear();
		for (int i = 0; i < S->trail.size(); i++) {
			if (!S->hasTheory(S->trail[i]) || S->getTheoryID(S->trail[i]) != getTheoryIndex())
				continue;
			Lit l = S->getTheoryLit(S->trail[i]);
			Var v = var(l);
			if (isPointVar(v)) {
				PointData & e = points[getPointID(v)];
				c.push(l);
			}
			
		}
		c.push(~l);
		//bool res = dbg->solve(c);
		//assert(~res);
#endif
		return true;
	}
	
	void dbg_sync() {
#ifndef NDEBUG
		return;
		int sz = vars.size();
		for (int i = 0; i < vars.size(); i++) {
			Var v = vars[i].solverVar;
			lbool val = value(i);
			lbool exp = dbg_value(i);
			
			if (value(i) != S->value(vars[i].solverVar)) {
				assert(false);
			}
		}
#endif
	}
	void dbg_full_sync() {
#ifndef NDEBUG
		dbg_sync();
		
		for (int i = 0; i < points.size(); i++) {
			int pointset = getPointset(i);
			int pointsetIndex = getPointsetIndex(i);
			if (points[i].id >= 0 && under_sets[pointset].isEnabled(pointsetIndex)) {
				assert(value(points[i].var)==l_True);
			} else if (points[i].id >= 0 && !over_sets[pointset].isEnabled(pointsetIndex)) {
				assert(value(points[i].var)==l_False);
			}
		}
		
#endif
	}
	
	void backtrackUntil(int level) {
		static int it = 0;
		
		bool changed = false;
		//need to remove and add points in the two graphs accordingly.
		if (trail_lim.size() > level) {
			
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				
				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if (e.isPoint) {
					assert(dbg_value(e.var)==value(e.var));
					int point_num = getPointID(e.var); //e.var-min_point_var;
					int pointSet = points[point_num].pointset;
					int pointsetIndex = getPointsetIndex(point_num);
					if (e.assign) {
						under_sets[pointSet].disablePoint(pointsetIndex);
					} else {
						over_sets[pointSet].enablePoint(pointsetIndex);
						//assert(over.hasPoint(e.pointID));
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
			
		}
		//is this correct?
		for (int i : toPropagate) {
			requiresPropagation[i] = false;
		}
		toPropagate.clear();
		
		assert(dbg_graphsUpToDate());
	}
	;
	virtual bool supportsDecisions() {
		return true;
	}
	Lit decideTheory() {
		if (!opt_decide_theories)
			return lit_Undef;
		double start = rtime(1);
		
		dbg_full_sync();
		for (int i = 0; i < detectors.size(); i++) {
			GeometryDetector * r = detectors[i];
			Lit l = r->decide();
			if (l != lit_Undef) {
				stats_decisions++;
				return toSolver(l);
			}
		}
		stats_decision_time += rtime(1) - start;
		return lit_Undef;
	}
	
	void backtrackUntil(Lit p) {
		//need to remove and add points in the two graphs accordingly.
		int i = trail.size() - 1;
		for (; i >= 0; i--) {
			Assignment e = trail[i];
			if (var(p) == e.var) {
				assert(sign(p) != e.assign);
				break;
			}
			if (e.isPoint) {
				int point_num = getPointID(e.var); //e.var-min_point_var;
				int pointSet = points[point_num].pointset;
				int pointsetIndex = getPointsetIndex(point_num);
				assert(assigns[e.var]!=l_Undef);
				assigns[e.var] = l_Undef;
				if (e.assign) {
					under_sets[pointSet].disablePoint(pointsetIndex);
				} else {
					over_sets[pointSet].enablePoint(pointsetIndex);
					
				}
			} else {

				assigns[e.var] = l_Undef;
				detectors[getDetector(e.var)]->unassign(mkLit(e.var, !e.assign));
			}
		}
		
		trail.shrink(trail.size() - (i + 1));
		//is this correct?
		for (int i : toPropagate) {
			requiresPropagation[i] = false;
		}
		toPropagate.clear();
		
		//while(trail_lim.size() && trail_lim.last()>=trail.size())
		//	trail_lim.pop();
		
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
	
	bool dbg_graphsUpToDate() {
#ifdef DEBUG_GRAPH
		return true;
		for(int i = 0;i<points.size();i++) {
			if(points[i].var<0)
			continue;
			PointData e = points[i];
			lbool val = value(e.var);
			int pointSet = e.pointset;
			if(val==l_True || val==l_Undef) {
				assert(over_sets[pointSet].pointEnabled(e.pointset_index));
			} else {
				assert(!over_sets[pointSet].pointEnabled(e.pointset_index));
			}
			if(val==l_True) {
				assert(under_sets[pointSet].pointEnabled(e.pointset_index));

			} else {
				assert(!under_sets[pointSet].pointEnabled(e.pointset_index));

			}
		}

#endif
		return true;
	}
	
	/*	int getPointID(Var v){
	 assert(v>= min_point_var && v<min_point_var+points.size());

	 //this is an point assignment
	 int point_num = v-min_point_var;
	 return point_num;
	 }*/

	void preprocess() {
		for (int i = 0; i < detectors.size(); i++) {
			detectors[i]->preprocess();
		}
	}
	void setLiteralOccurs(Lit l, bool occurs) {
		if (isPointVar(var(l))) {
			//don't do anything
		} else {
			//this is a graph property detector var
			if (!sign(l) && vars[var(l)].occursPositive != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
			else if (sign(l) && vars[var(l)].occursNegative != occurs)
				detectors[getDetector(var(l))]->setOccurs(l, occurs);
		}
		
	}
	void printSolution() {
		for (GeometryDetector* d : detectors) {
			d->printSolution();
		}
	}
	void enqueueTheory(Lit l) {
		Var v = var(l);
		assert(S->value(toSolver(l))==l_True);
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
		
#ifndef NDEBUG
		/*{
		 for(int i = 0;i<trail.size();i++){
		 assert(trail[i].var !=v);
		 }
		 }*/
#endif
		/*
		 #ifdef RECORD
		 if(under.outfile){
		 fprintf(under.outfile,"enqueue %d\n", dimacs(l));

		 fprintf(under.outfile,"\n");
		 fflush(under.outfile);
		 }
		 if(over.outfile){
		 fprintf(over.outfile,"enqueue %d\n", dimacs(l));
		 fprintf(over.outfile,"\n");
		 fflush(over.outfile);
		 }
		 #endif
		 */

		//if(v>= min_point_var && v<min_point_var+points.size())
		if (isPointVar(var(l))) {
			
			//this is an point assignment
			int pointID = getPointID(var(l)); //v-min_point_var;
			assert(points[pointID].var == var(l));
			
			int pointsetIndex = getPointsetIndex(pointID);
			int pointsetID = getPointset(pointID);
			trail.push(Assignment(true, !sign(l), pointID, var(l)));
			
			if (!requiresPropagation[pointsetID]) {
				toPropagate.push(pointsetID);
				//anyRequiresPropagation=true;
				requiresPropagation[pointsetID] = true;
			}
			Assignment e = trail.last();
			int pointSet = points[pointID].pointset;
			if (!sign(l)) {
				under_sets[pointSet].enablePoint(pointsetIndex);
			} else {
				over_sets[pointSet].disablePoint(pointsetIndex);
			}
			
		} else {
			
			trail.push(Assignment(false, !sign(l), 0, v));
			//this is an assignment to a non-point atom. (eg, a reachability assertion)
			detectors[getDetector(var(l))]->assign(l);
		}
		
	}
	;
	bool propagateTheory(vec<Lit> & conflict) {
		static int itp = 0;
		++itp;
		/*if(	++itp>=955880){
		 int a =1;
		 cout<<"pause\n";
		 for(int i = 0;i<vars.size();i++){
		 if(value(i)!= dbg_value(i)){
		 cout << "Error! Theory unsolved or out of sync: theory var " << i;

		 if(isPointVar(i)){
		 cout << " for point " << getPointID(i) << " " << points[ getPointID(i) ].point;
		 }else{
		 cout << " for detector " << getDetector(i);
		 }
		 cout << " has value " << toInt(value(i)) << " but expected value was " << toInt(dbg_value(i));
		 cout<< "!\n";

		 }
		 }
		 }*/
		stats_propagations++;
		dbg_sync();
		if (!toPropagate.size()) {
			stats_propagations_skipped++;
			assert(dbg_graphsUpToDate());
			return true;
		}
		
		bool any_change = false;
		double startproptime = rtime(1);
		static vec<int> detectors_to_check;
		
		conflict.clear();
		
		dbg_sync();
		assert(dbg_graphsUpToDate());
		
		while (toPropagate.size()) {
			int pointset = toPropagate.last();
			
			for (GeometryDetector * d : pointsetDetectors[pointset]) {
				assert(conflict.size() == 0);
				bool r = d->propagate(conflict);
				if (!r) {
					stats_num_conflicts++;
					toSolver(conflict);
					propagationtime += rtime(1) - startproptime;
					return false;
				}
			}
			toPropagate.pop();
			requiresPropagation[pointset] = false;
			under_sets[pointset].clearChanged();
			over_sets[pointset].clearChanged();
			
			under_sets[pointset].clearHistory();
			over_sets[pointset].clearHistory();
		}
		
		if (collisionDetector) {
			bool r = collisionDetector->propagate(conflict);
			if (!r) {
				stats_num_conflicts++;
				toSolver(conflict);
				propagationtime += rtime(1) - startproptime;
				return false;
			}
		}
		
#ifndef NDEBUG
		/*for(int i = 0;i<detectors.size();i++){
		 assert(detectors[i]->checkSatisfied());
		 }*/
#endif
		dbg_full_sync();
		
		detectors_to_check.clear();
		/*	if(	itp>=955879){
		 int a =1;
		 cout<<"pause after\n";
		 for(int i = 0;i<vars.size();i++){
		 if(value(i)!= dbg_value(i)){
		 cout << "Error! Theory unsolved or out of sync: theory var " << i;

		 if(isPointVar(i)){
		 cout << " for point " << getPointID(i) << " " << points[ getPointID(i) ].point;
		 }else{
		 cout << " for detector " << getDetector(i);
		 }
		 cout << " has value " << toInt(value(i)) << " but expected value was " << toInt(dbg_value(i));
		 cout<< "!\n";

		 }
		 }
		 for(int i = 0;i<detectors.size();i++){
		 if(!detectors[i]->checkSatisfied()){
		 cout<< "detector " << i << "unsat\n";
		 exit(3);
		 }
		 }
		 }*/

		double elapsed = rtime(1) - startproptime;
		propagationtime += elapsed;
		dbg_sync();
		
		return true;
	}
	;

	bool solveTheory(vec<Lit> & conflict) {
		
		for (int i = 0; i < nPointSets(); i++) {
			//Just to be on the safe side, force all detectors to propagate... this shouldn't really be required.
			if (!requiresPropagation[i]) {
				toPropagate.push(i);
				requiresPropagation[i] = true;
			}
		}
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	}
	;

	void drawFull() {
		
	}
	
	bool check_solved() {
		dbg_full_sync();
		for (int i = 0; i < vars.size(); i++) {
			if (value(i) != dbg_value(i)) {
				cout << "Error! Theory unsolved or out of sync: theory var " << i;
				
				if (isPointVar(i)) {
					cout << " for point " << getPointID(i) << " " << points[getPointID(i)].point;
				} else {
					cout << " for detector " << getDetector(i);
				}
				cout << " has value " << toInt(value(i)) << " but expected value was " << toInt(dbg_value(i));
				cout << "!\n";
				return false;
			}
		}
		for (int i = 0; i < points.size(); i++) {
			if (points[i].var < 0)
				continue;
			PointData & e = points[i];
			lbool val = value(e.var);
			if (val == l_Undef) {
				cout << "Error! Theory unsolved!\n";
				return false;
			}
			int pointset = e.pointset;
			if (val == l_True) {
				/*	if(!g.haspoint(e.from,e.to)){
				 return false;
				 }
				 if(!over.haspoint(e.from,e.to)){
				 return false;
				 }*/
				if (!under_sets[pointset].pointEnabled(e.pointset_index)) {
					cout << "Error! Theory out of sync!\n";
					return false;
				}
				if (!over_sets[pointset].pointEnabled(e.pointset_index)) {
					cout << "Error! Theory out of sync!\n";
					return false;
				}
			} else {
				/*if(g.haspoint(e.from,e.to)){
				 return false;
				 }*/
				if (under_sets[pointset].pointEnabled(e.pointset_index)) {
					cout << "Error! Theory out of sync!\n";
					return false;
				}
				if (over_sets[pointset].pointEnabled(e.pointset_index)) {
					cout << "Error! Theory out of sync!\n";
					return false;
				}
				/*if(over.haspoint(e.from,e.to)){
				 return false;
				 }*/

			}
			
		}
		for (int i = 0; i < detectors.size(); i++) {
			if (!detectors[i]->checkSatisfied()) {
				cout << "Error! Detector " << i << " unsatisfied\n";
				return false;
			}
		}
		return true;
	}
	
	bool dbg_solved() {
#ifdef DEBUG_GRAPH
		for(int i = 0;i<points.size();i++) {
			if(points[i].var<0)
			continue;
			PointData & e = points[i];
			lbool val = value(e.var);
			assert(val!=l_Undef);
			int pointset = e.pointset;
			if(val==l_True) {
				assert(under_sets[pointset].pointEnabled(e.pointset_index));
				assert(over_sets[pointset].haspoint(e.pointset_index));
			} else {
				assert(!under_sets[pointset].haspoint(e.pointset_index));
				assert(!over_sets[pointset].haspoint(e.pointset_index));
			}

		}
		assert(check_solved());

#endif
		return true;
	}
	
	void drawCurrent() {
		
	}
	int npoints() {
		return points.size();
	}
	CRef newReasonMarker(int detectorID) {
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef - reasonMarker;
		marker_map.growTo(mnum + 1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}
	
	Lit newPoint(int pointSet, Point<D, T> _point, Var outerVar = var_Undef) {
		assert(outerVar!=var_Undef);
		/*	if(outerVar==var_Undef)
		 outerVar = S->newVar();*/

		int index = points.size();
		points.push_back(PointData());
		
		Var v = newVar(outerVar, index, true);
		
		//num_points++;
		points[index].var = v;
		//points[index].outerVar =outerVar;
		points[index].point = _point;
		points[index].point.setID(index);
		points[index].id = index;
		points[index].pointset = pointSet;
		while (under_sets.size() <= pointSet) {
			under_sets.push_back(under_sets.size());
		}
		while (over_sets.size() <= pointSet) {
			over_sets.push_back(over_sets.size());
		}
		
		int pointsetIndex = under_sets[pointSet].addPoint(points[index].point);
		points[index].pointset_index = pointsetIndex;
		under_sets[pointSet].disablePoint(pointsetIndex);
		int pointsetIndex2 = over_sets[pointSet].addPoint(points[index].point);
		assert(pointsetIndex == pointsetIndex2);
		over_sets[pointSet].enablePoint(pointsetIndex);
		if (pointsetDetectors.size() < nPointSets()) {
			pointsetDetectors.resize(nPointSets());
		}
		
		//anyRequiresPropagation=true;
		requiresPropagation.growTo(nPointSets(), false);
		if (!requiresPropagation[pointSet]) {
			requiresPropagation[pointSet] = true;
			toPropagate.push(pointSet);
		}
		return mkLit(v, false);
	}
	
	void implementConstraints() {
		if (!S->okay())
			return;
	}
	
	void convexHullArea(int pointSet, T areaGreaterThan, Var outerVar) {
		if (convexHullDetectors.size() < nPointSets())
			convexHullDetectors.resize(nPointSets());
		if (!convexHullDetectors[pointSet]) {
			int detectorID = detectors.size();
			auto * convexHull = new ConvexHullDetector<D, T>(detectorID, under_sets[pointSet], over_sets[pointSet],
					this, drand(rnd_seed));
			pointsetDetectors[pointSet].push_back(convexHull);
			convexHullDetectors[pointSet] = convexHull;
			detectors.push_back(convexHull);
		}
		convexHullDetectors[pointSet]->addAreaDetectorLit(areaGreaterThan, outerVar);
	}
	
	void convexHullContains(int pointSet, Point<D, T> point, Var outerVar) {
		if (convexHullDetectors.size() < nPointSets())
			convexHullDetectors.resize(nPointSets());
		if (!convexHullDetectors[pointSet]) {
			int detectorID = detectors.size();
			auto * convexHull = new ConvexHullDetector<D, T>(detectorID, under_sets[pointSet], over_sets[pointSet],
					this, drand(rnd_seed));
			pointsetDetectors[pointSet].push_back(convexHull);
			convexHullDetectors[pointSet] = convexHull;
			detectors.push_back(convexHull);
		}
		convexHullDetectors[pointSet]->addPointContainmentLit(point, outerVar);
	}
	
	void pointOnHull(int pointSet, int pointIndex, Var outerVar) {
		if (convexHullDetectors.size() < nPointSets())
			convexHullDetectors.resize(nPointSets());
		if (!convexHullDetectors[pointSet]) {
			int detectorID = detectors.size();
			auto * convexHull = new ConvexHullDetector<D, T>(detectorID, under_sets[pointSet], over_sets[pointSet],
					this, drand(rnd_seed));
			pointsetDetectors[pointSet].push_back(convexHull);
			convexHullDetectors[pointSet] = convexHull;
			detectors.push_back(convexHull);
		}
		
		convexHullDetectors[pointSet]->addPointOnHullLit(pointIndex, outerVar);
	}
	
	/*	void convexHullIntersectsLine(int pointSet,Point<D,T> p,Point<D,T> q,Var outerVar){
	 convexHullDetectors.growTo(nPointSets(),nullptr);
	 if(!convexHullDetectors[pointSet]){
	 int detectorID = detectors.size();
	 auto * convexHull = new ConvexHullDetector<D,T>(detectorID,under_sets[pointSet], over_sets[pointSet],this,drand(rnd_seed));
	 convexHullDetectors[pointSet]=convexHull;
	 detectors.push(convexHull);
	 }

	 convexHullDetectors[pointSet]->addLineIntersection(LineSegment<D,T>(p,q),outerVar);
	 }*/

	void convexHullIntersectsPolygon(int pointSet, std::vector<Point<D, T>> & points, Var outerVar, bool inclusive =
			true) {
		//For now, polygon must be a simple, CONVEX polygon with clockwise winding
		if (convexHullDetectors.size() < nPointSets())
			convexHullDetectors.resize(nPointSets());
		if (!convexHullDetectors[pointSet]) {
			int detectorID = detectors.size();
			auto * convexHull = new ConvexHullDetector<D, T>(detectorID, under_sets[pointSet], over_sets[pointSet],
					this, drand(rnd_seed));
			convexHullDetectors[pointSet] = convexHull;
			pointsetDetectors[pointSet].push_back(convexHull);
			detectors.push_back(convexHull);
		}
		if (computeWinding(points) == Winding::COUNTER_CLOCKWISE) {
			fprintf(stderr, "Winding must be clockwise, aborting\n");
			exit(1);
		}
		if (!isConvex(points)) {
			fprintf(stderr, "Polygon must be convex, aborting\n");
			exit(1);
		}
		NConvexPolygon<D, T> polygon;
		for (auto & p : points) {
			polygon.addVertex(p);
		}
		if (polygon.size() == 1) {
			convexHullDetectors[pointSet]->addPointContainmentLit(polygon[0], outerVar, inclusive);
		} else if (polygon.size() == 2) {
			convexHullDetectors[pointSet]->addLineIntersectionLit(LineSegment<D, T>(polygon[0], polygon[1]), outerVar,
					inclusive);
		} else
			convexHullDetectors[pointSet]->addConvexIntersectionLit(polygon, outerVar, inclusive);
	}
	void convexHullsIntersect(int pointSet1, int pointSet2, Var outerVar, bool inclusive = true) {
		if (!collisionDetector) {
			collisionDetector = new ConvexHullCollisionDetector<D, T>(detectors.size(), this, under_sets, over_sets,
					convexHullDetectors, drand(rnd_seed));
			detectors.push_back(collisionDetector);
		}
		if (!convexHullDetectors[pointSet1]) {
			int detectorID = detectors.size();
			auto * convexHull = new ConvexHullDetector<D, T>(detectorID, under_sets[pointSet1], over_sets[pointSet1],
					this, drand(rnd_seed));
			convexHullDetectors[pointSet1] = convexHull;
			pointsetDetectors[pointSet1].push_back(convexHull);
			detectors.push_back(convexHull);
		}
		if (!convexHullDetectors[pointSet2]) {
			int detectorID = detectors.size();
			auto * convexHull = new ConvexHullDetector<D, T>(detectorID, under_sets[pointSet2], over_sets[pointSet2],
					this, drand(rnd_seed));
			convexHullDetectors[pointSet2] = convexHull;
			pointsetDetectors[pointSet2].push_back(convexHull);
			detectors.push_back(convexHull);
		}
		collisionDetector->addCollisionDetectorLit(pointSet1, pointSet2, outerVar, inclusive);
	}
	
	/*void euclidianSteinerTreeSize(int pointSet, int sizeLessThan, Var outerVar){
	 steinerTreeDetectors.growTo(nPointSets(),nullptr);
	 if(!steinerTreeDetectors[pointSet]){
	 int detectorID = detectors.size();
	 auto * steinerTree = new GeometricSteinerDetector<D,T>(detectorID,this,drand(rnd_seed));
	 pointsetDetectors[pointSet].push_back(steinerTree);
	 steinerTreeDetectors[pointSet]=steinerTree;
	 detectors.push_back(steinerTree);
	 }

	 steinerTreeDetectors[pointSet]->addAreaDetectorLit(sizeLessThan,outerVar);
	 }*/

};

#endif /* DGRAPH_H_ */
