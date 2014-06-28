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
#include "PointSet.h"
#include "ConvexHullDetector.h"
#include "GeometrySteinerDetector.h"
#ifndef NDEBUG
#include <cstdio>
#endif

using namespace Minisat;

template<unsigned int D, class T=double>
class GeometryTheorySolver:public Theory{

public:

	double rnd_seed;
	Lit False;
	Lit True;
	int local_q;
private:
	Solver * S;
public:
	int id;


	vec<lbool> assigns;

	PointSet<D,T> under;
	PointSet<D,T> over;
	struct Assignment{
		bool isPoint:1;
		bool assign:1;
		int pointID;
		Var var;
		Assignment(bool isPoint,bool _assign,int pointID, Var v):isPoint(isPoint),assign(_assign),pointID(pointID),var(v){

		}
		Assignment():isPoint(false),assign(false),pointID(0),var(var_Undef){

		}
	};
	vec<Assignment> trail;
	vec<int> trail_lim;

	vec<GeometryDetector*> detectors;
	ConvexHullDetector<D,T>* convexHull=nullptr;
 	GeometricSteinerDetector<D,T> * steinerTree=nullptr;
	vec<int> marker_map;

	bool requiresPropagation;

	vec<char> seen;
	vec<int> to_visit;
	vec<Lit> tmp_clause;

	//Data about local theory variables, and how they connect to the sat solver's variables
	struct VarData{
		int isPoint:1;
		int occursPositive:1;
		int occursNegative:1;
		int detector_point:29;//the detector this variable belongs to, or its point number, if it is an point variable
		Var solverVar;
	};
	struct PointData{
		int id;
		Var var;
		Point<D,T> point;
	};
	vec<VarData> vars;
	vec<PointData> points;
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


	GeometryTheorySolver(Solver * S_, int _id=-1):S(S_),id(_id){
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

			 requiresPropagation=true;
			 rnd_seed=opt_random_seed;

	}

	void printStats(int detailLevel){
		if(detailLevel>0){
			for(GeometryDetector * d:detectors)
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
	inline bool isPointVar(Var v){
		assert(v<vars.size());
		return vars[v].isPoint;
	}
	inline int getPointID(Var v){
		assert(isPointVar(v));
		return vars[v].detector_point;
	}
	inline int getDetector(Var v){
		assert(!isPointVar(v));
		return vars[v].detector_point;
	}

	inline Var getPointVar(int pointID){
		Var v = points[pointID].var;
		assert(v<vars.size());
		assert(vars[v].isPoint);
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

	Var newVar(int forDetector=-1, bool connectToTheory=false){
		Var s= S->newVar();

		return newVar(s,forDetector,false,connectToTheory);
	}
	Var newVar(Var solverVar, int detector, bool isPoint=false, bool connectToTheory=true){
		while(S->nVars()<=solverVar)
				S->newVar();
		Var v = vars.size();
		vars.push();
		vars[v].isPoint=isPoint;
		vars[v].detector_point=detector;
		vars[v].solverVar=solverVar;
		assigns.push(l_Undef);
		if(connectToTheory){
			S->setTheoryVar(solverVar,getTheoryIndex(),v);
			assert(toSolver(v)==solverVar);
		}
		if(!isPoint && detector>=0)
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



     ~GeometryTheorySolver(){};

		bool dbg_propgation(Lit l){
#ifndef NDEBUG
			static vec<Lit> c;
			c.clear();
			for(int i = 0;i<S->trail.size() ;i++){
				if(!S->hasTheory(S->trail[i]) || S->getTheoryID(S->trail[i])!= getTheoryIndex())
					continue;
				Lit l = S->getTheoryLit(S->trail[i]);
				Var v =  var(l);
				if(isPointVar(v)){
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


	void dbg_sync(){

	}
	void dbg_full_sync(){
	#ifndef NDEBUG
			dbg_sync();
			for(int i =0;i<points.size();i++){

				if(points[i].id>=0 && under.isEnabled(i)){
					assert(value(points[i].var)==l_True);
				}else if (points[i].id>=0 &&  !over.isEnabled(i)){
					assert(value(points[i].var)==l_False);
				}
			}

	#endif
		}



	void backtrackUntil(int level){
		static int it = 0;

		bool changed=false;
		//need to remove and add points in the two graphs accordingly.
		if(trail_lim.size()>level){

			int stop = trail_lim[level];
			for(int i = trail.size()-1;i>=trail_lim[level];i--){

				Assignment & e = trail[i];
				assert(assigns[e.var]!=l_Undef);
				if(e.isPoint){
					assert(dbg_value(e.var)==l_Undef);
					int point_num = getPointID(e.var); //e.var-min_point_var;

					if(e.assign){
						under.disablePoint( point_num);
					}else{
						over.enablePoint(point_num);
						//assert(over.hasPoint(e.pointID));
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

		}

		assert(dbg_graphsUpToDate());
	};

	Lit decideTheory(){
		if(!opt_decide_graph)
			return lit_Undef;
		double start = rtime(1);

		dbg_full_sync();
		for(int i = 0;i<detectors.size();i++){
			GeometryDetector * r = detectors[i];
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
			//need to remove and add points in the two graphs accordingly.
			int i = trail.size()-1;
			for(;i>=0;i--){
				Assignment e = trail[i];
				if(e.isPoint){
					int point_num = getPointID(e.var); //e.var-min_point_var;
					assert(assigns[e.var]!=l_Undef);
					assigns[e.var]=l_Undef;
					if(e.assign){
						under.disablePoint(point_num);
					}else{
						over.enablePoint(point_num);

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


	bool dbg_graphsUpToDate(){
#ifdef DEBUG_GRAPH
		for(int i = 0;i<points.size();i++){
			if(points[i].var<0)
				continue;
			PointData e = points[i];
			lbool val = value(e.var);

			if(val==l_True || val==l_Undef){
				assert(over.pointEnabled(e.id));

			}else{
				assert(!over.pointEnabled(e.id));

			}
			if(val==l_True){
				assert(under.pointEnabled(e.id));

			}else{
				assert(!under.pointEnabled(e.id));

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

	void preprocess(){
		for (int i = 0;i<detectors.size();i++){
			detectors[i]->preprocess();
		}
	}
	void setLiteralOccurs(Lit l, bool occurs){
		if(isPointVar(var(l))){
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

#ifndef NDEBUG
		{
			for(int i = 0;i<trail.size();i++){
				assert(trail[i].var !=v);
			}
		}
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
		if(isPointVar(var(l))){

			//this is an point assignment
			int pointID = getPointID(var(l)); //v-min_point_var;
			assert(points[pointID].var==var(l));


			trail.push(Assignment(true,!sign(l), pointID,var(l)));

			Assignment e = trail.last();
			if (!sign(l)){
				under.enablePoint(pointID);
			}else{
				over.disablePoint(pointID);
			}

		}else{

			trail.push(Assignment(false,!sign(l), 0,v));
			//this is an assignment to a non-point atom. (eg, a reachability assertion)
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
		under.clearChanged();
		over.clearChanged();


		under.clearHistory();
		over.clearHistory();

		detectors_to_check.clear();

		double elapsed = rtime(1)-startproptime;
		propagationtime+=elapsed;
		dbg_sync();

		return true;
	};

	bool solveTheory(vec<Lit> & conflict){
		requiresPropagation=true;//Just to be on the safe side... but this shouldn't really be required.
		bool ret = propagateTheory(conflict);
		//Under normal conditions, this should _always_ hold (as propagateTheory should have been called and checked by the parent solver before getting to this point).
		assert(ret);
		return ret;
	};


	void drawFull(){

	}

	bool check_solved(){
		for(int i = 0;i<vars.size();i++){
			if(value(i)!= dbg_value(vars[i].solverVar)){
				return false;
			}

		}
		for(int i = 0;i<points.size();i++){
			if(points[i].var<0)
						continue;
			PointData & e = points[i];
			lbool val = value(e.var);
			if(val==l_Undef){
				return false;
			}

			if(val==l_True){
			/*	if(!g.haspoint(e.from,e.to)){
					return false;
				}
				if(!over.haspoint(e.from,e.to)){
					return false;
				}*/
				if(!under.pointEnabled(e.id)){
					return false;
				}
				if(!over.pointEnabled(e.id)){
					return false;
				}
			}else{
				/*if(g.haspoint(e.from,e.to)){
					return false;
				}*/
				if(under.pointEnabled(e.id)){
					return false;
				}
				if(over.pointEnabled(e.id)){
					return false;
				}
				/*if(over.haspoint(e.from,e.to)){
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
		for(int i = 0;i<points.size();i++){
			if(points[i].var<0)
				continue;
			PointData & e = points[i];
			lbool val = value(e.var);
			assert(val!=l_Undef);

			if(val==l_True){
				assert(under.pointEnabled(e.id));
				assert(over.haspoint(e.id));
			}else{
				assert(!under.haspoint(e.id));
				assert(!over.haspoint(e.id));
			}


		}
		assert(check_solved());

#endif
		return true;
	}

	void drawCurrent(){

	}
	int npoints(){
		return points.size();
	}
	CRef newReasonMarker(int detectorID){
		CRef reasonMarker = S->newReasonMarker(this);
		int mnum = CRef_Undef- reasonMarker;
		marker_map.growTo(mnum+1);
		marker_map[mnum] = detectorID;
		return reasonMarker;
	}

	Lit newPoint(Point<D,T> point, Var outerVar = var_Undef)
    {
		assert(outerVar!=var_Undef);
	/*	if(outerVar==var_Undef)
			outerVar = S->newVar();*/


		int index = points.size();
		points.push();
		Var v = newVar(outerVar,index,true);

		//num_points++;
		points[index].var =v;
		//points[index].outerVar =outerVar;
		points[index].point=point;
		points[index].id=index;

		under.addPoint(point,index);
		under.disablePoint(index);
		over.addPoint(point,index);
		over.enablePoint(index);
    	return mkLit(v,false);
    }

	void implementConstraints(){
		if(!S->okay())
			return;
	}

	void printSolution(){

	}

	void convexHullArea(T areaGreaterThan, Var outerVar){
		if(!convexHull){
			int detectorID = detectors.size();
			convexHull = new ConvexHullDetector<D,T>(detectorID,under, over,this,drand(rnd_seed));
			detectors.push(convexHull);
		}
		convexHull->addAreaDetectorLit(areaGreaterThan,outerVar);
	}

	void convexHullContains(Point<D,T> point, Var outerVar){
		if(!convexHull){
			int detectorID = detectors.size();
			convexHull = new ConvexHullDetector<D,T>(detectorID,under, over,this,drand(rnd_seed));
			detectors.push(convexHull);
		}
		convexHull->addPointContainmentLit(point,outerVar);
	}

	void euclidianSteinerTreeSize(int sizeLessThan, Var outerVar){
		if(!steinerTree){
			int detectorID = detectors.size();
			steinerTree = new GeometricSteinerDetector<D,T>(detectorID,this,drand(rnd_seed));
			detectors.push(steinerTree);
		}

		steinerTree->addAreaDetectorLit(sizeLessThan,outerVar);
	}

};


#endif /* DGRAPH_H_ */
