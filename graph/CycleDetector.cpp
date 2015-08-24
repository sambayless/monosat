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
#include "dgl/PKTopologicalSort.h"
#include "dgl/DFSCycle.h"
#include "core/Config.h"
#include <limits>
using namespace Monosat;
template<typename Weight>
CycleDetector<Weight>::CycleDetector(int _detectorID, GraphTheorySolver<Weight> * _outer, DynamicGraph<Weight>  &g_under,
		DynamicGraph<Weight>  &g_over, bool detect_directed_cycles, double seed) :
		Detector(_detectorID), outer(_outer), g_under(g_under), g_over(g_over), rnd_seed(seed), underapprox_directed_cycle_detector(NULL), overapprox_directed_cycle_detector(
				NULL) {
	
	undirected_acyclic_lit = lit_Undef;
	directed_acyclic_lit = lit_Undef;
	

	if(cyclealg==CycleAlg::ALG_DFS_CYCLE){
		underapprox_directed_cycle_detector = new DFSCycle<Weight,true,true>(g_under, detect_directed_cycles, 1);
		overapprox_directed_cycle_detector = new DFSCycle<Weight,true,true>(g_over, detect_directed_cycles, 1);

		overapprox_undirected_cycle_detector=overapprox_directed_cycle_detector;
		underapprox_undirected_cycle_detector=underapprox_directed_cycle_detector;

	}else if(cyclealg==CycleAlg::ALG_PK_CYCLE){
		underapprox_directed_cycle_detector = new PKToplogicalSort<Weight>(g_under,  1);
		overapprox_directed_cycle_detector = new PKToplogicalSort<Weight>(g_over,  1);

		underapprox_undirected_cycle_detector = new DFSCycle<Weight,false,true>(g_under, false, 1);
		overapprox_undirected_cycle_detector = new DFSCycle<Weight,false,true>(g_over, false, 1);
	}
	directed_cycle_marker = outer->newReasonMarker(getID());
	no_directed_cycle_marker = outer->newReasonMarker(getID());
	
	undirected_cycle_marker = outer->newReasonMarker(getID());
	no_undirected_cycle_marker = outer->newReasonMarker(getID());
	//forced_reach_marker=outer->newReasonMarker(getID());
}
template<typename Weight>
void CycleDetector<Weight>::addAcyclicLit(bool directed, Var outer_reach_var) {
	Var v = outer->newVar(outer_reach_var, getID());
	Lit l = mkLit(v, false);
	g_under.invalidate();
	g_over.invalidate();
	if (!directed) {
		if (undirected_acyclic_lit == lit_Undef) {
			undirected_acyclic_lit = l;
			
		} else {
			outer->makeEqual(undirected_acyclic_lit, l);
			/*outer->S->addClause(undirected_cycle_lit, ~l);
			 outer->S->addClause(~undirected_cycle_lit, l);*/
		}
	} else {
		if (directed_acyclic_lit == lit_Undef) {
			directed_acyclic_lit = l;
			
		} else {
			outer->makeEqual(directed_acyclic_lit, l);
			/*			outer->S->addClause(directed_cycle_lit, ~l);
			 outer->S->addClause(~directed_cycle_lit, l);*/
		}
	}
}
template<typename Weight>
void CycleDetector<Weight>::buildNoUndirectedCycleReason(vec<Lit> & conflict) {
	//its clear that we can do better than this, but its also not clear how to do so efficiently...
	//for now, learn the trivial clause...

	assert(!overapprox_undirected_cycle_detector->hasUndirectedCycle());
	for (int edgeID = 0; edgeID < g_over.edges(); edgeID++) {
		if(g_over.hasEdge(edgeID) && ! g_over.edgeEnabled(edgeID)){
			Var v = outer->getEdgeVar(edgeID);
			if (outer->value(v) == l_False) {
				conflict.push(mkLit(v, false));
			}
		}
	}
}
template<typename Weight>
void CycleDetector<Weight>::buildNoDirectedCycleReason(vec<Lit> & conflict) {
	//its clear that we can do better than this, but its also not clear how to do so efficiently...
	//for now, learn the trivial clause...
	//One thing you could do would be to first do an over-approx cycle detection at level 0, and exclude from here any edges that can't possibly be part of any scc.
	//Is that the best one can do?
	assert(!overapprox_directed_cycle_detector->hasDirectedCycle());
	for (int edgeID = 0; edgeID < g_over.edges(); edgeID++) {
		if(g_over.hasEdge(edgeID) && ! g_over.edgeEnabled(edgeID)){
			Var v = outer->getEdgeVar(edgeID);
			if (outer->value(v) == l_False) {
				conflict.push(mkLit(v, false));
			}
		}
	}
	
}

template<typename Weight>
void CycleDetector<Weight>::buildUndirectedCycleReason(vec<Lit> & conflict) {
	assert(underapprox_undirected_cycle_detector->hasUndirectedCycle());
	
	std::vector<int> & cycle = underapprox_undirected_cycle_detector->getUndirectedCycle();
	for (int i = 0; i < cycle.size(); i++) {
		int edgeID = cycle[i];
		Lit l = mkLit(outer->getEdgeVar(edgeID), false);
		assert(outer->value(l)==l_True);
		conflict.push(~l);
	}
	
}
template<typename Weight>
void CycleDetector<Weight>::buildDirectedCycleReason(vec<Lit> & conflict) {
	assert(underapprox_directed_cycle_detector->hasDirectedCycle());
	
	std::vector<int> & cycle = underapprox_directed_cycle_detector->getDirectedCycle();
	for (int i = 0; i < cycle.size(); i++) {
		int edgeID = cycle[i];
		Lit l = mkLit(outer->getEdgeVar(edgeID), false);
		assert(outer->value(l)==l_True);
		conflict.push(~l);
	}
	
}
template<typename Weight>
void CycleDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker) {
	
	if (marker == directed_cycle_marker) {
		reason.push(p);
		
		buildDirectedCycleReason(reason);
		
	} else if (marker == no_directed_cycle_marker) {
		reason.push(p);
		
		buildNoDirectedCycleReason(reason);
		
	} else if (marker == undirected_cycle_marker) {
		reason.push(p);
		
		buildUndirectedCycleReason(reason);
		
	} else if (marker == no_undirected_cycle_marker) {
		reason.push(p);
		
		buildNoUndirectedCycleReason(reason);
		
	} else {
		assert(false);
	}
}
template<typename Weight>
bool CycleDetector<Weight>::propagate(vec<Lit> & conflict) {
	static int it = 0;
/*	if(++it==687){
		g_under.drawFull();
		int a=1;
	}*/

	double startdreachtime = rtime(2);
	if(directed_acyclic_lit != lit_Undef && outer->value(directed_acyclic_lit)==l_False && outer->level(var(directed_acyclic_lit))==0){
		underapprox_directed_cycle_detector->forceDAG();
	}
	

	
	if (directed_acyclic_lit != lit_Undef) {
		
		if (outer->value(directed_acyclic_lit) !=l_False && underapprox_directed_cycle_detector->hasDirectedCycle()) {

			Lit l = ~directed_acyclic_lit;

			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				//trail.push(Assignment(false,false,detectorID,0,var(l)));
				outer->enqueue(l, no_directed_cycle_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildDirectedCycleReason(conflict);
				return false;
			}
		} else if (outer->value(directed_acyclic_lit) !=l_True && !overapprox_directed_cycle_detector->hasDirectedCycle()) {
			Lit l = directed_acyclic_lit;
			
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				//trail.push(Assignment(false,true,detectorID,0,var(l)));
				outer->enqueue(l, directed_cycle_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildNoDirectedCycleReason(conflict);
				return false;
			}
		}
		
	} else if (undirected_acyclic_lit != lit_Undef) {
		
		if (outer->value(undirected_acyclic_lit) !=l_False && underapprox_undirected_cycle_detector->hasUndirectedCycle()) {

			Lit l = ~directed_acyclic_lit;
			
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				//trail.push(Assignment(false,false,detectorID,0,var(l)));
				outer->enqueue(l, no_undirected_cycle_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildUndirectedCycleReason(conflict);
				return false;
			}
		} else if (outer->value(undirected_acyclic_lit) !=l_True &&  !overapprox_undirected_cycle_detector->hasUndirectedCycle()) {
			Lit l = undirected_acyclic_lit;
			
			if (outer->value(l) == l_True) {
				//do nothing
			} else if (outer->value(l) == l_Undef) {
				//trail.push(Assignment(false,true,detectorID,0,var(l)));
				outer->enqueue(l, undirected_cycle_marker);
			} else if (outer->value(l) == l_False) {
				conflict.push(l);
				buildNoUndirectedCycleReason(conflict);
				return false;
			}
		}
		
	}
	return true;
}
template<typename Weight>
bool CycleDetector<Weight>::checkSatisfied() {
	DFSCycle<Weight,true,false> checkDirected (g_under, true, 1);
	DFSCycle<Weight,false,true> checkUndirected(g_under, false, 1);
	if(directed_acyclic_lit != lit_Undef){
		if(outer->value(directed_acyclic_lit)==l_True && checkDirected.hasDirectedCycle()){
			return false;
		}else if (outer->value(directed_acyclic_lit)==l_False && !checkDirected.hasDirectedCycle()){
			return false;
		}

	}
	if(undirected_acyclic_lit != lit_Undef){
		if(outer->value(undirected_acyclic_lit)==l_True && checkUndirected.hasUndirectedCycle()){
			return false;
		}else if (outer->value(undirected_acyclic_lit)==l_False && !checkUndirected.hasUndirectedCycle()){
			return false;
		}

	}
	return true;
}
template<typename Weight>
Lit CycleDetector<Weight>::decide() {
	
	return lit_Undef;
}
;

template class Monosat::CycleDetector<int> ;
template class Monosat::CycleDetector<long> ;
template class Monosat::CycleDetector<double> ;
#include <gmpxx.h>
template class Monosat::CycleDetector<mpq_class> ;
