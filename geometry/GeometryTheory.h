/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

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


namespace Minisat{


class GeometryTheorySolver:public Theory{
public:

	double rnd_seed;
	Lit False;
	Lit True;
	int local_q;

	Solver * S;
	int id;

	vec<lbool> edge_assignments;

	Var min_edge_var;
	int num_edges;



public:
	vec<GeometryDetector*> detectors;
	vec<GeometryDetector*> polygonList;

	vec<int> marker_map;


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

	GeometryTheorySolver(Solver * S_, int _id=-1):S(S_),id(_id){

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

	}

	int getGraphID(){
		return id;
	}
	void printStats(){

	}

     ~GeometryTheorySolver(){};

 	CRef newReasonMarker(int detectorID){
 		CRef reasonMarker = S->newReasonMarker(this);
 		int mnum = CRef_Undef- reasonMarker;
 		marker_map.growTo(mnum+1);
 		marker_map[mnum] = detectorID;
 		return reasonMarker;
 	}
	void backtrackUntil(int level){
		for(int i = 0;i<detectors.size();i++){
			detectors[i]->backtrackUntil(level);
		}

	};

	Lit decideTheory(){
		if(!opt_decide_graph)
			return lit_Undef;
		double start = cpuTime();
		for(int i = 0;i<detectors.size();i++){

			Lit l = detectors[i]->decide();
			if(l!=lit_Undef)
				return l;

		}
		stats_decision_time += cpuTime() - start;
		return lit_Undef;
	}

	void backtrackUntil(Lit p){
		for(int i = 0;i<detectors.size();i++){
				detectors[i]->backtrackUntil(p);
			}

	};

	void newDecisionLevel(){

		//trail_lim.push(trail.size());
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

		stats_initial_propagation_time += cpuTime() - startproptime;

		for(int d = 0;d<detectors.size();d++){
			assert(conflict.size()==0);
			bool r =detectors[d]->propagate(S->trail,conflict);
			if(!r){
				propagationtime+= cpuTime()-startproptime;
				return false;
			}
		}

		detectors_to_check.clear();

		double elapsed = cpuTime()-startproptime;
		propagationtime+=elapsed;

		return true;
	};



	bool solveTheory(vec<Lit> & conflict){return true;};

	bool check_solved(){

		for(int i = 0;i<detectors.size();i++){
			if(!detectors[i]->checkSatisfied()){
				return false;
			}
		}
		return true;
	}



	void createConvexHull(int hullID, const int D);
	void addHullPoint(int hullID,Lit l, vec<int> & point){
		static vec<double> tmp;
		tmp.clear();
		for(int i = 0;i<point.size();i++){
			tmp.push(point[i]);
		}
		addHullPoint(hullID,l,tmp);
	}
	void addHullPoint(int hullID, Lit l, vec<double> & point);

	void addPointContainmentLit(int polyID, Lit l, vec<int> & point){
		static vec<double> tmp;
		tmp.clear();
		for(int i = 0;i<point.size();i++){
			tmp.push(point[i]);
		}
		addPointContainmentLit(polyID,l,tmp);
	}
	void addPointContainmentLit(int polyID, Lit l, vec<double> & point){
		polygonList[polyID]->addPointContainmentLit(l,point);
	}
};

};

#endif /* DGRAPH_H_ */
