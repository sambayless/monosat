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


#include "CycleDetector.h"
#include "GraphTheory.h"
#include <limits>
using namespace Monosat;
template<typename Weight>
CycleDetector<Weight>::CycleDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,  DynamicGraph &_g,DynamicGraph &_antig, bool detect_directed_cycles,double seed):
Detector(_detectorID),outer(_outer),g(_g),antig(_antig),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL){

	undirected_cycle_lit = lit_Undef;
	directed_cycle_lit=lit_Undef;

		//Note: these are _intentionalyl_ swapped
		negative_reach_detector = new DFSCycle(_g,detect_directed_cycles,1);
		positive_reach_detector = new DFSCycle(_antig,detect_directed_cycles,1);

		 directed_cycle_marker=outer->newReasonMarker(getID());
			 no_directed_cycle_marker=outer->newReasonMarker(getID());

			 undirected_cycle_marker=outer->newReasonMarker(getID());
			 no_undirected_cycle_marker=outer->newReasonMarker(getID());
		//forced_reach_marker=outer->newReasonMarker(getID());
}
template<typename Weight>
void CycleDetector<Weight>::addCycleDetectorLit(bool directed, Var outer_reach_var){
	Var v= outer->newVar(outer_reach_var,getID());
	Lit l = mkLit(v,false);
	g.invalidate();
	antig.invalidate();
	if(!directed){
		if(undirected_cycle_lit==lit_Undef){
			undirected_cycle_lit=l;

		}else{
			outer->makeEqual(undirected_cycle_lit,l);
			/*outer->S->addClause(undirected_cycle_lit, ~l);
			outer->S->addClause(~undirected_cycle_lit, l);*/
		}
	}else{
		if(directed_cycle_lit==lit_Undef){
			directed_cycle_lit=l;

		}else{
			outer->makeEqual(directed_cycle_lit,l);
/*			outer->S->addClause(directed_cycle_lit, ~l);
			outer->S->addClause(~directed_cycle_lit, l);*/
		}
	}
}
template<typename Weight>
		void CycleDetector<Weight>::buildNoUndirectedCycleReason(vec<Lit> & conflict){
			//its clear that we can do better than this, but its also not clear how to do so efficiently...
			//for now, learn the trivial clause...
			for(int i = 0;i<outer->edge_list.size();i++){
				Var v = outer->edge_list[i].v;
				if(outer->value(v)==l_False){
					conflict.push(mkLit(v,false));
				}
			}
		}
template<typename Weight>
		void CycleDetector<Weight>::buildNoDirectedCycleReason(vec<Lit> & conflict){
			//its clear that we can do better than this, but its also not clear how to do so efficiently...
			//for now, learn the trivial clause...
			for(int i = 0;i<outer->edge_list.size();i++){
				Var v = outer->edge_list[i].v;
				if(outer->value(v)==l_False){
					conflict.push(mkLit(v,false));
				}
			}

		}

template<typename Weight>
		void CycleDetector<Weight>::buildUndirectedCycleReason(vec<Lit> & conflict){
			assert(positive_reach_detector->hasUndirectedCycle());

			std::vector<int> & cycle = positive_reach_detector->getUndirectedCycle();
			for(int i = 0;i<cycle.size();i++){
				int e = cycle[i];
				Lit l = mkLit( outer->edge_list[e].v,false);
				assert(outer->value(l)==l_True);
				conflict.push(~l);
			}

		}
template<typename Weight>
		void CycleDetector<Weight>::buildDirectedCycleReason(vec<Lit> & conflict){
			assert(positive_reach_detector->hasDirectedCycle());

			std::vector<int> & cycle= positive_reach_detector->getDirectedCycle();
			for(int i = 0;i<cycle.size();i++){
				int e = cycle[i];
				Lit l = mkLit( outer->edge_list[e].v,false);
				assert(outer->value(l)==l_True);
				conflict.push(~l);
			}

		}
template<typename Weight>
		void CycleDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker){

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
template<typename Weight>
		bool CycleDetector<Weight>::propagate(vec<Lit> & conflict){


		double startdreachtime = rtime(2);

		positive_reach_detector->update();
		double reachUpdateElapsed = rtime(2)-startdreachtime;
		outer->reachupdatetime+=reachUpdateElapsed;

		double startunreachtime = rtime(2);
		negative_reach_detector->update();
		double unreachUpdateElapsed = rtime(2)-startunreachtime;
		outer->unreachupdatetime+=unreachUpdateElapsed;


		if(directed_cycle_lit!=lit_Undef){

			if(positive_reach_detector->hasDirectedCycle()){
				Lit l = directed_cycle_lit;

				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					//trail.push(Assignment(false,true,detectorID,0,var(l)));
					outer->enqueue(l,directed_cycle_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildDirectedCycleReason(conflict);
					return false;
				}
			}else if(!negative_reach_detector->hasDirectedCycle()){
				Lit l = ~directed_cycle_lit;

				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					//trail.push(Assignment(false,false,detectorID,0,var(l)));
					outer->enqueue(l,no_directed_cycle_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildNoDirectedCycleReason(conflict);
					return false;
				}
			}

		}else if(undirected_cycle_lit!=lit_Undef){

			if(positive_reach_detector->hasUndirectedCycle()){
				Lit l = undirected_cycle_lit;

				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					//trail.push(Assignment(false,true,detectorID,0,var(l)));
					outer->enqueue(l,undirected_cycle_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildUndirectedCycleReason(conflict);
					return false;
				}
			}else if(!negative_reach_detector->hasUndirectedCycle()){
				Lit l = ~directed_cycle_lit;

				if(outer->value(l)==l_True){
					//do nothing
				}else if(outer->value(l)==l_Undef){
					//trail.push(Assignment(false,false,detectorID,0,var(l)));
					outer->enqueue(l,no_undirected_cycle_marker) ;
				}else if (outer->value(l)==l_False){
					conflict.push(l);
					buildNoUndirectedCycleReason(conflict);
					return false;
				}
			}

		}
			return true;
		}
template<typename Weight>
bool CycleDetector<Weight>::checkSatisfied(){

	return true;
}
template<typename Weight>
Lit CycleDetector<Weight>::decide(){

	return lit_Undef;
};

template class CycleDetector<int>;
template class CycleDetector<double>;
#include <gmpxx.h>
template class CycleDetector<mpq_class>;
