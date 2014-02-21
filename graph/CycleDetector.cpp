/*
 * CycleDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "CycleDetector.h"
#include "GraphTheory.h"
#include <limits>
CycleDetector::CycleDetector(int _detectorID, GraphTheorySolver * _outer,  DynamicGraph<PositiveEdgeStatus> &_g,DynamicGraph<NegativeEdgeStatus> &_antig, bool detect_directed_cycles,double seed):
Detector(_detectorID),outer(_outer),g(_g),antig(_antig),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL){

	undirected_cycle_lit = lit_Undef;
	directed_cycle_lit=lit_Undef;

		//Note: these are _intentionalyl_ swapped
		negative_reach_detector = new DFSCycle<PositiveEdgeStatus>(_g,detect_directed_cycles,1);
		positive_reach_detector = new DFSCycle<NegativeEdgeStatus>(_antig,detect_directed_cycles,1);

		 directed_cycle_marker=outer->newReasonMarker(getID());
			 no_directed_cycle_marker=outer->newReasonMarker(getID());

			 undirected_cycle_marker=outer->newReasonMarker(getID());
			 no_undirected_cycle_marker=outer->newReasonMarker(getID());
		//forced_reach_marker=outer->newReasonMarker(getID());
}
void CycleDetector::addCycleDetectorLit(bool directed, Var v){
	Lit l = mkLit(v,false);
	g.invalidate();
	antig.invalidate();
	if(!directed){
		if(undirected_cycle_lit==lit_Undef){
			undirected_cycle_lit=l;
		}else{
			outer->S->addClause(undirected_cycle_lit, ~l);
			outer->S->addClause(~undirected_cycle_lit, l);
		}
	}else{
		if(directed_cycle_lit==lit_Undef){
			directed_cycle_lit=l;
		}else{
			outer->S->addClause(directed_cycle_lit, ~l);
			outer->S->addClause(~directed_cycle_lit, l);
		}
	}
}

		void CycleDetector::buildNoUndirectedCycleReason(vec<Lit> & conflict){
			//its clear that we can do better than this, but its also not clear how to do so efficiently...
			//for now, learn the trivial clause...
			for(int i = 0;i<outer->edge_list.size();i++){
				Var v = outer->edge_list[i].v;
				if(outer->S->value(v)==l_False){
					conflict.push(mkLit(v,false));
				}
			}
		}
		void CycleDetector::buildNoDirectedCycleReason(vec<Lit> & conflict){
			//its clear that we can do better than this, but its also not clear how to do so efficiently...
			//for now, learn the trivial clause...
			for(int i = 0;i<outer->edge_list.size();i++){
				Var v = outer->edge_list[i].v;
				if(outer->S->value(v)==l_False){
					conflict.push(mkLit(v,false));
				}
			}

		}


		void CycleDetector::buildUndirectedCycleReason(vec<Lit> & conflict){
			assert(positive_reach_detector->hasUndirectedCycle());

			vec<int> & cycle = positive_reach_detector->getUndirectedCycle();
			for(int i = 0;i<cycle.size();i++){
				int e = cycle[i];
				Lit l = mkLit( outer->edge_list[e].v,false);
				assert(outer->S->value(l)==l_True);
				conflict.push(~l);
			}

		}
		void CycleDetector::buildDirectedCycleReason(vec<Lit> & conflict){
			assert(positive_reach_detector->hasDirectedCycle());

			vec<int> & cycle = positive_reach_detector->getDirectedCycle();
			for(int i = 0;i<cycle.size();i++){
				int e = cycle[i];
				Lit l = mkLit( outer->edge_list[e].v,false);
				assert(outer->S->value(l)==l_True);
				conflict.push(~l);
			}

		}

		void CycleDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){

				if(marker==directed_cycle_marker){
					reason.push(p);

					buildDirectedCycleReason(reason);

				}else if(marker==no_directed_cycle_marker){
					reason.push(p);

					buildNoDirectedCycleReason(reason);

				}else if(marker==undirected_cycle_marker){
					reason.push(p);

					buildUndirectedCycleReason(reason);

				}else if(marker==no_undirected_cycle_marker){
					reason.push(p);

					buildNoUndirectedCycleReason(reason);

				}else{
					assert(false);
				}
		}

		bool CycleDetector::propagate(vec<Assignment> & trail,vec<Lit> & conflict){


		double startdreachtime = cpuTime();

		positive_reach_detector->update();
		double reachUpdateElapsed = cpuTime()-startdreachtime;
		outer->reachupdatetime+=reachUpdateElapsed;

		double startunreachtime = cpuTime();
		negative_reach_detector->update();
		double unreachUpdateElapsed = cpuTime()-startunreachtime;
		outer->unreachupdatetime+=unreachUpdateElapsed;


		if(directed_cycle_lit!=lit_Undef){

			if(positive_reach_detector->hasDirectedCycle()){
				Lit l = directed_cycle_lit;

				if(outer->S->value(l)==l_True){
					//do nothing
				}else if(outer->S->value(l)==l_Undef){
					trail.push(Assignment(false,true,detectorID,0,var(l)));
					outer->S->uncheckedEnqueue(l,directed_cycle_marker) ;
				}else if (outer->S->value(l)==l_False){
					conflict.push(l);
					buildDirectedCycleReason(conflict);
					return false;
				}
			}else if(!negative_reach_detector->hasDirectedCycle()){
				Lit l = ~directed_cycle_lit;

				if(outer->S->value(l)==l_True){
					//do nothing
				}else if(outer->S->value(l)==l_Undef){
					trail.push(Assignment(false,false,detectorID,0,var(l)));
					outer->S->uncheckedEnqueue(l,no_directed_cycle_marker) ;
				}else if (outer->S->value(l)==l_False){
					conflict.push(l);
					buildNoDirectedCycleReason(conflict);
					return false;
				}
			}

		}else if(undirected_cycle_lit!=lit_Undef){

			if(positive_reach_detector->hasUndirectedCycle()){
				Lit l = undirected_cycle_lit;

				if(outer->S->value(l)==l_True){
					//do nothing
				}else if(outer->S->value(l)==l_Undef){
					trail.push(Assignment(false,true,detectorID,0,var(l)));
					outer->S->uncheckedEnqueue(l,undirected_cycle_marker) ;
				}else if (outer->S->value(l)==l_False){
					conflict.push(l);
					buildUndirectedCycleReason(conflict);
					return false;
				}
			}else if(!negative_reach_detector->hasUndirectedCycle()){
				Lit l = ~directed_cycle_lit;

				if(outer->S->value(l)==l_True){
					//do nothing
				}else if(outer->S->value(l)==l_Undef){
					trail.push(Assignment(false,false,detectorID,0,var(l)));
					outer->S->uncheckedEnqueue(l,no_undirected_cycle_marker) ;
				}else if (outer->S->value(l)==l_False){
					conflict.push(l);
					buildNoUndirectedCycleReason(conflict);
					return false;
				}
			}

		}
			return true;
		}

bool CycleDetector::checkSatisfied(){

	return true;
}
Lit CycleDetector::decide(){

	return lit_Undef;
};


