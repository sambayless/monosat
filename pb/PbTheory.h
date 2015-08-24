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

#ifndef PBTHEORY_H_
#define PBTHEORY_H_

#include "mtl/Vec.h"
#include "mtl/Map.h"
#include "mtl/Sort.h"
#include "core/SolverTypes.h"
#include "core/Theory.h"
#include <cmath>
namespace Monosat {

static BoolOption opt_dbg_prop("DEBUG", "prove-prop", "", false);

//Note: Although we include this (very basic) pseudo-boolean constraint solver
//as an example of a monotonic theory, it doesn't (so far) perform nearly as well as dedicated PB solvers, such as Minisat+.
//For now, it is included only for the sake of completeness.

class PbTheory: public Theory {
	Solver * S;
	int theory_index;

	int qhead = 0;
	//vec<Assignment> trail;
	vec<Lit> trail;
	vec<int> trail_lim;
	vec<int> inq;
	vec<int> t_weights;
public:
	enum class PbType {
		GT, GE, LT, LE, EQ, NE
	};

	enum class ConstraintSide {
		Upper, Lower, Both
	};
private:
	ConstraintSide invert(ConstraintSide c) {
		if (c == ConstraintSide::Upper)
			return ConstraintSide::Lower;
		else if (c == ConstraintSide::Lower) {
			return ConstraintSide::Upper;
		}
		return c;
	}
	
	struct ConstraintToImplement {
		
		bool implemented = false;
		vec<Lit> clause;
		vec<int> weights;
		Lit rhs_lit;
		int total;
		ConstraintSide side;
		PbType op;
		
	};

	vec<ConstraintToImplement> constraintsToImplement;

	struct VarData {
		//bool onesided;//whether this is a one-sided or two-sided pb constraint (if one sided, then if the output var is false, the constraint is simply unenforced)
		Var solverVar;
		int sign :1;
		int rhs :1;
		int clauseID :30;
		unsigned int weight;
		VarData() :
				solverVar(var_Undef), sign(0), rhs(false), clauseID(-1), weight(0) {
		}
	};
	struct PbElement {
		Lit lit;
		unsigned int weight;
	};
	struct PbClause {
		ConstraintSide side;
		//bool isSatisfied;
		//bool isOverEq; //true iff the sum of the lhs is >= the rhs

		//maintain running counts of the under/over approximations of this clause.
		long under;
		long unassigned;
		bool inQueue;
		PbElement rhs;
		CRef reason;
		vec<PbElement> clause;
		PbClause() :
				side(ConstraintSide::Both), under(0), unassigned(0), inQueue(false), reason(CRef_Undef) {
		}
	};

	//Map reasons to clauseIDs
	Map<CRef, int> reasonMap;
	vec<PbClause> clauses;
	vec<VarData> vars;
	vec<lbool> assigns;
	vec<Lit> clause;
	vec<int> weights;
	vec<int> indices;
	vec<Lit> tmp_clause;
	vec<int> tmp_weights;
	double propagationtime;
	long stats_propagations, stats_propagations_skipped;
	long stats_shrink_removed = 0;
	long stats_reasons = 0;
	long stats_conflicts = 0;
	/*	 Var newVar(){
	 Var s= S->newVar();
	 Var v = vars.size();
	 vars.push();
	 vars[v].solverVar=s;

	 return v;
	 }*/
	Var newVar(Var solverVar, int clauseID, int weight, bool sign, bool isRHS = false) {
		while (S->nVars() <= solverVar)
			S->newVar();
		
		if (S->hasTheory(solverVar)) {
			//introduce a new replacement variable
			Var v2 = S->newVar();
			S->addClause(~mkLit(solverVar), mkLit(v2));
			S->addClause(mkLit(solverVar), ~mkLit(v2));
			solverVar = v2;
		}
		
		Var v = vars.size();
		assigns.push(l_Undef);
		vars.push();
		vars[v].solverVar = solverVar;
		vars[v].clauseID = clauseID;
		vars[v].weight = weight;
		vars[v].sign = sign;
		vars[v].rhs = isRHS;
		
		S->setTheoryVar(solverVar, getTheoryIndex(), v);
		assert(toSolver(v) == solverVar);
		
		return v;
	}
	inline Var toSolver(Var v) {
		//return v;
		assert(v < vars.size());
		assert(S->hasTheory(vars[v].solverVar));
		assert(S->getTheoryVar(vars[v].solverVar) == v);
		return vars[v].solverVar;
	}
	
	inline Lit toSolver(Lit l) {
		assert(S->hasTheory(vars[var(l)].solverVar));
		assert(S->getTheoryVar(vars[var(l)].solverVar) == var(l));
		return mkLit(vars[var(l)].solverVar, sign(l));
	}
	
	void toSolver(vec<Lit> & c) {
		for (int i = 0; i < c.size(); i++) {
			c[i] = toSolver(c[i]);
		}
	}
	int level(Var v) {
		return S->level(toSolver(v));
	}
	inline lbool value(Var v) {
		if (assigns[v] != l_Undef) {
			assert(S->value(toSolver(v)) == assigns[v]);
		}
		return assigns[v]; //S->value(toSolver(v));
	}
	inline lbool value(Lit l) {
		if (assigns[var(l)] != l_Undef) {
			assert(S->value(toSolver(var(l))) == assigns[var(l)]);
		}
		return assigns[var(l)] ^ sign(l); //S->value(toSolver(l));
	}
	inline bool enqueue(Lit l, CRef reason) {
		
		Lit sl = toSolver(l);
		if (S->enqueue(sl, reason)) {
			enqueueTheory(l);
			return true;
		} else {
			return false;
		}
	}
	void makeEqual(Lit l1, Lit l2) {
		Lit o1 = toSolver(l1);
		Lit o2 = toSolver(l2);
		S->addClause(~o1, o2);
		S->addClause(o1, ~o2);
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
		vec<Lit> t;
		c.copyTo(t);
		toSolver(t);
		S->addClause(t);
	}
public:
	
	PbTheory(Solver * S) :
			S(S), theory_index(0), propagationtime(0), stats_propagations(0), stats_propagations_skipped(0) {
		
	}
	~PbTheory() {
	}
	;

	inline int getTheoryIndex() {
		return theory_index;
	}
	inline void setTheoryIndex(int id) {
		theory_index = id;
	}
	void newDecisionLevel() {
		trail_lim.push(trail.size());
	}
	;
	inline int decisionLevel() {
		return S->decisionLevel();
	}
	void enqueueTheory(Lit l) {
		Var v = var(l);
		
		int lev = level(v);
		
		assert(decisionLevel() <= lev);
		if (value(l) != l_Undef)
			return;
		while (lev > trail_lim.size()) {
			newDecisionLevel();
		}
		
		assigns[var(l)] = sign(l) ? l_False : l_True;
		trail.push(l);
		
		int clauseID = vars[v].clauseID;
		
		PbClause & c = clauses[clauseID];
		bool q = false;
		if (v == var(c.rhs.lit)) {
			q = true;
		} else {
			bool s = vars[v].sign;
			c.unassigned -= vars[v].weight;
			if (s == sign(l)) {
				//this is a positive assignment
				c.under += vars[v].weight;
				if (c.under >= c.rhs.weight) {
					q = true;
				}
			} else {
				
				if (c.under + c.unassigned < c.rhs.weight) {
					//this is in conflict; enqueue this clause
					q = true;
				}
			}
			
		}
		
		if (q) {
			if (!c.inQueue) {
				c.inQueue = true;
				inq.push(clauseID);
			}
			assert(inq.contains(clauseID));
		}
		//trail.push({true,toInt(l)});
	}
	;
	bool propagateTheory(vec<Lit> & conflict) {
		
		/*	 		if(trail.size()==0){
		 stats_propagations_skipped++;
		 return true;
		 }*/
		if (inq.size() == 0) {
			dbg_fully_propped();
			return true;
		}
		
		stats_propagations++;
		
		static int iter = 0;
		++iter;
		double startproptime = rtime(2);
		
		//This is wrong! Only need to visit each _clause_ that has any involved literals once per propagation round.
		//for(int i = 0;i<trail.size();i++){
		//int clauseID = trail[i];
		//for(int clauseID = 0;clauseID<clauses.size();clauseID++){
		while (inq.size()) {
			int clauseID = inq.last();
			PbClause & pbclause = clauses[clauseID];
			assert(pbclause.inQueue);
			//if(pbclause.isSatisfied)
			//	continue;
			
			Lit rhs = pbclause.rhs.lit;
			lbool rhs_val = value(rhs);
			if (pbclause.side == ConstraintSide::Upper && rhs_val == l_False) {
				//this clause is free.
				//pbclause.isSatisfied=true;
				pbclause.inQueue = false;
				inq.pop();
				continue;
			} else if (pbclause.side == ConstraintSide::Lower && rhs_val == l_True) {
				pbclause.inQueue = false;
				inq.pop();
				continue;
			}
			
			unsigned int total = pbclause.rhs.weight;
			//compute over and under approximations...
			unsigned int underApprox = pbclause.under;
			int unassignedWeight = pbclause.unassigned;
			unsigned int overApprox = pbclause.under + pbclause.unassigned;
			/*	int n_Free=0;
			 unsigned int unassignedWeight = 0;
			 unsigned int smallestUnassignedWeight = INT32_MAX;
			 unsigned int largestUnassignedWeight =0;
			 Lit largestUnassigned = lit_Undef;
			 for(PbElement e:pbclause.clause){
			 Lit l = e.lit;
			 lbool val = value(l);
			 if(val==l_True){
			 underApprox +=e.weight;
			 }else if(val==l_Undef){
			 n_Free++;
			 unassignedWeight+=e.weight;
			 if(e.weight>largestUnassignedWeight){
			 largestUnassignedWeight=e.weight;
			 largestUnassigned = e.lit;
			 }
			 if (e.weight<=smallestUnassignedWeight){
			 smallestUnassignedWeight=e.weight;
			 }
			 }else{
			 int  a=1;
			 }
			 }
			 */
			//overApprox=underApprox+;
			if (rhs_val == l_True && underApprox >= total) {
				//this is a satisfied constraint
				//pbclause.isSatisfied=true;
				
			} else if (rhs_val == l_False && overApprox < total) {
				//this is a satisfied constraint
				//pbclause.isSatisfied=true;
				
			} else if (rhs_val == l_True && overApprox < total) {
				//conflict
				assert(pbclause.side != ConstraintSide::Lower);
				static int iter = 0;
				++iter;
				conflict.push(~rhs);
				buildSumLTReason(clauseID, conflict);
				dbg_prove(pbclause, conflict);
				dbg_min_conflict(pbclause, conflict);
				toSolver(conflict);
				stats_conflicts++;
				return false;
			} else if (rhs_val == l_False && underApprox >= total) {
				//conflict
				assert(pbclause.side != ConstraintSide::Upper);
				conflict.push(rhs);
				
				buildSumGEReason(clauseID, conflict);
				dbg_prove(pbclause, conflict);
				dbg_min_conflict(pbclause, conflict);
				toSolver(conflict);
				stats_conflicts++;
				return false;
			} else if (underApprox >= total && rhs_val == l_Undef) {
				//pbclause.isSatisfied=true;
				//pbclause.isOverEq=true;
				dbg_prop(pbclause, rhs);
				//trail.push({false,clauseID});
				enqueue(rhs, pbclause.reason);
			} else if (overApprox < total && rhs_val == l_Undef) {
				//pbclause.isSatisfied=true;
				//pbclause.isOverEq=false;
				dbg_prop(pbclause, ~rhs);
				//trail.push({false,clauseID});
				enqueue(~rhs, pbclause.reason);
			} else if (rhs_val == l_False) {
				//pbclause.isOverEq=false;
				//then _all_ unassigned lits whose weight would push the under approx over over the limit must be assigned to false
				assert(pbclause.side != ConstraintSide::Upper);
				
				for (PbElement e : pbclause.clause) {
					Lit l = e.lit;
					if (value(l) == l_Undef && underApprox + e.weight >= total) {
						dbg_prop(pbclause, ~l);
						enqueue(~l, pbclause.reason);
					}
				}
				
			} else if (rhs_val == l_True) {
				assert(pbclause.side != ConstraintSide::Lower);
				//assert(n_Free>0);
				assert(overApprox >= total);
				//if(n_Free==1){
				if (overApprox - 2 < total) {
					//then the largest unassigned weight is forced
					//assert(smallestUnassignedWeight==largestUnassignedWeight);
					/*assert(largestUnassigned!=lit_Undef);
					 //assert(underApprox+largestUnassignedWeight>=total);//else the over approx would have triggered a conflict
					 dbg_prop(pbclause,largestUnassigned);
					 enqueue(largestUnassigned,pbclause.reason);

					 //it may also be the case that the second largest weight is forced.
					 if(n_Free>1){*/
					for (PbElement e : pbclause.clause) {
						Lit l = e.lit;
						lbool val = value(l);
						if (val == l_True) {
							
						} else if (val == l_Undef) {
							if (overApprox - e.weight < total) {
								dbg_prop(pbclause, e.lit);
								enqueue(e.lit, pbclause.reason);
							}
						}
					}
					//}
				}
			}
			clauses[clauseID].inQueue = false;
			inq.pop();
		}
		
#ifndef NDEBUG
		for (PbClause & c : clauses) {
			assert(!c.inQueue);
		}
#endif
		
		double elapsed = rtime(2) - startproptime;
		propagationtime += elapsed;
		dbg_fully_propped();
		return true;
	}
	;

private:
	void buildElementForcedFalseReason(int clauseID, Lit element, vec<Lit> & conflict) {
		//the reason why a lit is forced to false is the set of true literals in the clause, and the rhs
		PbClause & pbclause = clauses[clauseID];
		//assert(pbclause.isSatisfied);
		assert(value(pbclause.rhs.lit)==l_False);
		assert(pbclause.side != ConstraintSide::Upper);
		conflict.push(pbclause.rhs.lit);
		
		int underApprox = 0;
		int unassignedWeight = 0;
		int forcedWeight = 0;
		for (PbElement e : pbclause.clause) {
			Lit l = e.lit;
			if (e.lit == ~element) {
				forcedWeight = e.weight;
				assert(value(e.lit)==l_False);
			} else {
				assert(var(element) != var(e.lit));
			}
			
			if (value(l) == l_True) {
				//if(level(var(l))>0)
				conflict.push(~l);
				underApprox += e.weight;
			}
		}
		assert(forcedWeight + underApprox >= pbclause.rhs.weight);
	}
	void buildElementForcedTrueReason(int clauseID, Lit element, vec<Lit> & conflict) {
		//the reason why a lit is forced to false is the set of true literals in the clause, and the rhs
		PbClause & pbclause = clauses[clauseID];
		assert(pbclause.side != ConstraintSide::Lower);
		//assert(pbclause.isSatisfied);
		lbool v = value(pbclause.rhs.lit);
		assert(value(pbclause.rhs.lit)==l_True);
		conflict.push(~pbclause.rhs.lit);
		int startSize = conflict.size();
		int rhs = pbclause.rhs.weight;
		int overApprox = 0;
		int forcedWeight = 0;
		t_weights.clear();
		for (PbElement e : pbclause.clause) {
			Lit l = e.lit;
			if (e.lit == element) {
				forcedWeight = e.weight;
				assert(value(e.lit)==l_True);
			} else {
				assert(var(element) != var(e.lit));
				//assert(value(l)!=l_Undef);
			}
			
			if (value(l) == l_False) {
				//if(level(var(l))>0)
				assert(var(e.lit) != var(element));
				
				conflict.push(l);
				if (opt_shrink_theory_conflicts)
					t_weights.push(e.weight);
			} else {
				overApprox += e.weight;
			}
			
		}
		
		if (opt_shrink_theory_conflicts) {
			assert(t_weights.size() + startSize == conflict.size());
			int i, j = startSize;
			for (i = startSize; i < conflict.size(); i++) {
				Lit l = conflict[i];
				int w = t_weights[i - startSize];
				if (overApprox - forcedWeight + w < rhs) {
					//then we can safely drop this from the conflict
					stats_shrink_removed++;
					overApprox += w;
				} else {
					conflict[j++] = l;
				}
			}
			conflict.shrink(i - j);
		}

		//assert(overApprox>=pbclause.rhs.weight);
		assert(overApprox - forcedWeight < pbclause.rhs.weight);
		
	}
	void buildSumGEReason(int clauseID, vec<Lit> & conflict) {
		PbClause & pbclause = clauses[clauseID];
		//assert(!pbclause.isSatisfied);
		
		int underApprox = 0;
		int unassignedWeight = 0;
		for (PbElement e : pbclause.clause) {
			Lit l = e.lit;
			if (value(l) == l_True) {
				//if(level(var(l))>0)
				conflict.push(~l);
				underApprox += e.weight;
			} else if (value(l) == l_Undef) {
				unassignedWeight += e.weight;
			}
		}
		assert(underApprox >= pbclause.rhs.weight);
	}
	
	void buildSumLTReason(int clauseID, vec<Lit> & conflict) {
		PbClause & pbclause = clauses[clauseID];
		//assert(!pbclause.isSatisfied);
		t_weights.clear();
		int startSize = conflict.size();
		int rhs = pbclause.rhs.weight;
		int underApprox = 0;
		int unassignedWeight = 0;
		for (PbElement e : pbclause.clause) {
			Lit l = e.lit;
			if (value(l) == l_True) {
				underApprox += e.weight;
			} else if (value(l) == l_Undef) {
				unassignedWeight += e.weight;
			} else if (value(l) == l_False) {
				//if(level(var(l))>0)
				conflict.push(l);
				if (opt_shrink_theory_conflicts)
					t_weights.push(e.weight);
				
			}
		}
		int overApprox = underApprox + unassignedWeight;
		
		assert(underApprox + unassignedWeight < pbclause.rhs.weight);
		if (opt_shrink_theory_conflicts) {
			assert(t_weights.size() + startSize == conflict.size());
			int i, j = startSize;
			for (i = startSize; i < conflict.size(); i++) {
				Lit l = conflict[i];
				int w = t_weights[i - startSize];
				if (overApprox + w < rhs) {
					overApprox += w;
					stats_shrink_removed++;
					//then we can safely drop this from the conflict
				} else {
					conflict[j++] = l;
				}
			}
			conflict.shrink(i - j);
		}

		assert(overApprox < pbclause.rhs.weight);
	}
	
	struct sortByClause {
		vec<Lit> & clause;
		sortByClause(vec<Lit> & clause) :
				clause(clause) {
		}
		bool operator ()(int x, int y) {
			return toInt(clause[x]) < toInt(clause[y]);
		}
	};

	struct sortByVec {
		vec<int> & weights;
		sortByVec(vec<int> & weights) :
				weights(weights) {
		}
		bool operator ()(int x, int y) {
			return weights[x] < weights[y];
		}
	};

	struct sortByVecDec {
		vec<int> & weights;
		sortByVecDec(vec<int> & weights) :
				weights(weights) {
		}
		bool operator ()(int x, int y) {
			return weights[x] > weights[y];
		}
	};
	//Follow the normalization guidelines from the minisatp paper.
	bool normalize(vec<Lit> & clause, vec<int> & weights, Lit & rhs_lit, int &rhs, ConstraintSide side) {
		int i, j = 0;
		assert(clause.size() == weights.size());
		assert(S->decisionLevel() == 0);
		
		//sort the clause and weights by the literals in the clause; so we can easily merge variables
		{
			indices.clear();
			for (int i = 0; i < clause.size(); i++) {
				indices.push(i);
			}
			sort(indices, sortByClause(clause));
			
			clause.copyTo(tmp_clause);
			clause.clear();
			weights.copyTo(tmp_weights);
			weights.clear();
			for (int i = 0; i < indices.size(); i++) {
				int index = indices[i];
				clause.push(tmp_clause[index]);
				weights.push(tmp_weights[index]);
			}
		}
		bool replace_rhs = false;
		Lit last = lit_Undef;
		for (i = 0; i < clause.size(); i++) {
			Lit l = clause[i];
			int w = weights[i];
			
			if (w < 0) {
				//then invert this literal and update the right hand side
				w = -w;
				l = ~l;
				assert(w > 0);
				rhs += w;
			}
			
			if (S->value(l) == l_True) {
				//drop this literal, and subtract is total from the rhs
				rhs -= w;
				continue;
			} else if (S->value(l) == l_False) {
				//then drop this literal
				continue;
			} else if (w == 0) {
				//then drop this literal
				continue;
			} else {
				//don't change
			}
			if (var(l) == var(rhs_lit))
				replace_rhs = true;
			if (w == 0) {
				//drop this literal
			} else if (l == last) {
				//merge these lits
				assert(clause[j] == l);
				weights[j] += w;
			} else if (l == ~last) {
				assert(clause[j] == ~l);
				weights[j] -= w;
				//may now need to invert this previous weight... this is ugly, can we avoid this?
				if (weights[j] < 0) {
					weights[j] = -weights[j];
					clause[j] = ~clause[j];
					assert(weights[j] > 0);
					rhs += weights[j];
				} else if (weights[j] == 0) {
					//drop this literal...
					j--;
				}
			} else {
				weights[j] = w;
				clause[j++] = l;
			}
			
			if (j > 0) {
				last = clause[j - 1];
			} else {
				last = lit_Undef;
			}
		}
		
		clause.shrink(i - j);
		weights.shrink(i - j);
		
		if (rhs <= 0) {
			//then this is a trivially true constraint (even if the total weights are 0, because this is a >= constraint)
			clause.clear();
			weights.clear();
			if (side != ConstraintSide::Upper)
				S->addClause(rhs_lit);
			//if we are one sided, then any value for rhs_lit is acceptable
			return false;
		}
		if (clause.size() == 0) {
			assert(weights.size() == 0);
			if (side != ConstraintSide::Lower)
				S->addClause(~rhs_lit);
			return false;
		}
		
		for (int i = 0; i < weights.size(); i++) {
			assert(weights[i] > 0);
			if (weights[i] > rhs) {
				weights[i] = rhs;
			}
		}
		
		int weightTotal = 0;
		bool allWeights1 = true;
		//as per the minisatp paper, divide the left and right hand sides by the gcd of the left hand side.
		int d = gcd(weights);
		for (int i = 0; i < weights.size(); i++) {
			weights[i] /= d;
			weightTotal += weights[i];
			allWeights1 &= (weights[i] == 1);
		}
		
		int rhs2 = (int) ceil((double) rhs / (double) d);
		assert(rhs / d == rhs2 || rhs / d == rhs2 - 1);
		rhs = rhs2;
		
		allWeights1 &= (rhs == 1);
		
		if (weightTotal < rhs) {
			//then this is a trivially unsat constraint
			clause.clear();
			weights.clear();
			if (side != ConstraintSide::Lower)
				S->addClause(~rhs_lit);
			return false;
		}
		assert(clause.size() > 0);							//else weightTotal is 0;
				
		//It is always the case that at most one element of the clause must be true, OR the value must be false.
		//so we can learn this clause right away
		if (side != ConstraintSide::Lower) {
			clause.push(~rhs_lit);
			S->addClause(clause);
			clause.pop();
		}
		
		if (side != ConstraintSide::Upper) {
			//then enforce the other side of the constraint
			for (int i = 0; i < clause.size(); i++) {
				if (weights[i] >= rhs) {
					S->addClause(rhs_lit, ~clause[i]); //If value is false, then all the lits in the clause that have >= weight as rhs must be false.
				}
			}
		}
		
		if (allWeights1) {
			//this as a pure SAT constraint
			clause.clear();
			weights.clear();
			return false;
		} else {
			assert(clause.size() > 1); 	 			//because we divided by gcd...
		}
		
		if (replace_rhs) {
			//the theory solver doesn't support using the rhs variable in the clause itself, so swap this literal out
			Lit new_rhs = mkLit(S->newVar());
			S->addClause(new_rhs, ~rhs_lit);
			S->addClause(~new_rhs, rhs_lit);
			rhs_lit = new_rhs;
			
		}
		
		// the minisatp paper suggests splitting the clause into a pb part and a pure-sat part by introducing a new variable
		// might want to implement that in the future...
		
		return true;
	}
	
	Lit greaterThanEq_implement(const vec<Lit> & _clause, const vec<int> &_weights, Lit rhs_lit, int total,
			ConstraintSide side = ConstraintSide::Both) {
		_clause.copyTo(clause);
		_weights.copyTo(weights);
		if (normalize(clause, weights, rhs_lit, total, side)) {
			
			for (int w : weights) {
				assert(w > 0);
				
			}
			int clauseID = clauses.size();
			assert(total > 0);
			
			assert(weights.size() == clause.size());
			
			clauses.push();
			CRef reason = S->newReasonMarker(this);
			PbClause & pbclause = clauses.last();
			rhs_lit = mkLit(newVar(var(rhs_lit), clauseID, total, sign(rhs_lit), true), sign(rhs_lit));
			pbclause.rhs.lit = rhs_lit;
			pbclause.rhs.weight = (unsigned int) total;
			pbclause.side = side;
			pbclause.reason = reason;
			reasonMap.insert(reason, clauseID);
			
			vars[var(rhs_lit)].clauseID = clauseID;
			vars[var(rhs_lit)].rhs = true;
			vars[var(rhs_lit)].weight = total;
			for (int i = 0; i < clause.size(); i++) {
				Lit l = clause[i];
				l = mkLit(newVar(var(l), clauseID, weights[i], sign(l)), sign(l));
				
				pbclause.clause.push( { l, (unsigned int) weights[i] });
				pbclause.unassigned += weights[i];
			}
		}
		return rhs_lit;
		
	}
	void dbg_fully_propped() {
#ifndef NDEBUG
		
		for (int i = 0; i < clauses.size(); i++) {
			int clauseID = i;
			
			PbClause & c = clauses[i];
			Lit rhsLit = c.rhs.lit;
			int underApprox = 0;
			int unassigned = 0;
			int overApprox = 0;
			int rhs = c.rhs.weight;
			for (PbElement & e : c.clause) {
				if (value(e.lit) == l_True) {
					underApprox += e.weight;
				} else if (value(e.lit) == l_Undef) {
					unassigned += e.weight;
				}
			}
			overApprox = underApprox + unassigned;
			assert(c.unassigned == unassigned);
			assert(c.under == underApprox);
			if (value(rhsLit) == l_True) {
				if (c.side != ConstraintSide::Lower) {
					assert(overApprox >= rhs);
				}
			} else if (value(rhsLit) == l_False) {
				if (c.side != ConstraintSide::Upper) {
					assert(underApprox < rhs);
				}
			}
			if (underApprox >= rhs) {
				assert(c.side== ConstraintSide::Upper || value(rhsLit)==l_True);
			} else if (overApprox < rhs) {
				assert(c.side== ConstraintSide::Lower || value(rhsLit)==l_False);
			}
			
		}
		
		if (opt_dbg_prop) {
			for (PbClause & c : clauses) {
				dbg_fully_propped(c);
				
			}
		}
#endif
	}
	void dbg_min_conflict(const PbClause & c, vec<Lit> &conflict) {
#ifndef NDEBUG
		if (opt_shrink_theory_conflicts) {
			dbg_prove(c, conflict);
			vec<Lit> t;
			for (int i = 0; i < conflict.size(); i++) {
				t.push(conflict[i]);
			}
			//dbg_unsat(c,t);
			for (int i = 0; i < conflict.size(); i++) {
				if (i != 2)
					continue;
				bool found = false;
				for (int j = 0; j < c.clause.size(); j++)
					found |= var(c.clause[j].lit) == var(conflict[i]);
				if (found) {
					t.clear();
					for (int j = 0; j < conflict.size(); j++) {
						if (j != i) {
							t.push(conflict[j]);
						}
					}
					dbg_unprove(c, t);
				}
			}
		}
#endif
	}
	void dbg_min_conflictb(const PbClause & c, vec<Lit> &conflict) {
#ifndef NDEBUG
		
		dbg_prove(c, conflict);
		vec<Lit> t;
		for (int i = 0; i < conflict.size(); i++) {
			t.push(conflict[i]);
		}
		for (int i = 0; i < c.clause.size(); i++) {
			if (value(c.clause[i].lit) != l_False) {
				printf("%d\n", c.clause[i].weight);
			}
		}
		//dbg_unsat(c,t);
		for (int i = 0; i < conflict.size(); i++) {
			if (i != 2)
				continue;
			bool found = false;
			for (int j = 0; j < c.clause.size(); j++)
				found |= var(c.clause[j].lit) == var(conflict[i]);
			if (found) {
				t.clear();
				for (int j = 0; j < conflict.size(); j++) {
					if (j != i) {
						t.push(conflict[j]);
					}
				}
				dbg_unprove(c, t);
			}
		}
#endif
	}
	void dbg_fully_propped(const PbClause & c) {
#ifndef NDEBUG
		if (!opt_dbg_prop)
			return;
		vec<Lit> prove;
		for (PbElement p : c.clause) {
			if (value(p.lit) == l_False) {
				prove.push(~p.lit);
			} else if (value(p.lit) == l_True) {
				prove.push(p.lit);
			}
		}
		
		if (value(c.rhs.lit) == l_False) {
			prove.push(~c.rhs.lit);
		} else if (value(c.rhs.lit) == l_True) {
			prove.push(c.rhs.lit);
		}
		
		dbg_sat(c, prove);
		
		for (PbElement p : c.clause) {
			// if(var(p.lit)==23 || toInt(p.lit)==23){
			
			if (value(p.lit) == l_Undef) {
				prove.push(p.lit);
				dbg_sat(c, prove);
				prove.pop();
				prove.push(~p.lit);
				dbg_sat(c, prove);
				prove.pop();
				
			}
			// }
		}
		
		if (value(c.rhs.lit) == l_Undef) {
			prove.push(c.rhs.lit);
			dbg_sat(c, prove);
			prove.pop();
			prove.push(~c.rhs.lit);
			dbg_sat(c, prove);
			prove.pop();
		}
#endif
	}
	
	void dbg_prop(const PbClause & c, Lit e) {
#ifndef NDEBUG
		vec<Lit> prove;
		for (PbElement p : c.clause) {
			if (value(p.lit) == l_False) {
				prove.push(p.lit);
			} else if (value(p.lit) == l_True) {
				prove.push(~p.lit);
			}
		}
		
		if (value(c.rhs.lit) == l_False) {
			prove.push(c.rhs.lit);
		} else if (value(c.rhs.lit) == l_True) {
			prove.push(~c.rhs.lit);
		}
		
		prove.push(e);
		dbg_prove(c, prove);
		
#endif
	}
	
	void dbg_prove(const PbClause & c, const vec<Lit> & clause) {
#ifndef NDEBUG
	/*	return;
		bool rhs_val = true;
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				if (l == c.rhs.lit) {
					rhs_val = false;
				} else {
					assert(l == ~c.rhs.lit);
					rhs_val = true;
				}
			}
		}
		FILE *f = fopen("testa.opb", "w");
		for (PbElement & e : c.clause) {
			fprintf(f, "%s%d x%d ","+" , e.weight, var(e.lit) + 1);
		}
		
		if (rhs_val) {
			fprintf(f, " >= %d ;\n", c.rhs.weight);
		} else {
			fprintf(f, " < %d ;\n", c.rhs.weight);
		}
		
		fflush(f);
		int r = system("minisat+ testa.opb >/dev/null") >> 8;
		assert(r == 10);
		
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				
			} else {
				bool found = false;
				bool sign = false;
				for (PbElement & e : c.clause) {
					if (e.lit == l) {
						found = true;
					} else if (e.lit == ~l) {
						sign = true;
						found = true;
					}
				}
				assert(found);
				if (!sign) {
					fprintf(f, "1 x%d < 1 ;\n", var(l) + 1);
				} else {
					fprintf(f, "1 x%d >= 1 ;\n", var(l) + 1);
				}
				
			}
		}
		fflush(f);
		fclose(f);
		r = system("minisat+ testa.opb >/dev/null") >> 8;
		assert(r == 20);*/
#endif
	}
	void dbg_unprove(const PbClause & c, const vec<Lit> & clause) {
#ifndef NDEBUG
		bool rhs_val = true;
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				if (l == c.rhs.lit) {
					rhs_val = false;
				} else {
					assert(l == ~c.rhs.lit);
					rhs_val = true;
				}
			}
		}
		FILE *f = fopen("testa.opb", "w");
		for (PbElement & e : c.clause) {
			fprintf(f, "%s%d x%d ", "+", e.weight, var(e.lit) + 1);
		}
		
		if (rhs_val) {
			fprintf(f, " >= %d ;\n", c.rhs.weight);
		} else {
			fprintf(f, " < %d ;\n", c.rhs.weight);
		}
		
		fflush(f);
		int r = system("minisat+ testa.opb >/dev/null") >> 8;
		assert(r == 10);
		
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				
			} else {
				bool found = false;
				bool sign = false;
				for (PbElement & e : c.clause) {
					if (e.lit == l) {
						found = true;
					} else if (e.lit == ~l) {
						sign = true;
						found = true;
					}
				}
				assert(found);
				if (!sign) {
					fprintf(f, "1 x%d < 1 ;\n", var(l) + 1);
				} else {
					fprintf(f, "1 x%d >= 1 ;\n", var(l) + 1);
				}
				
			}
		}
		fflush(f);
		fclose(f);
		r = system("minisat+ testa.opb >/dev/null") >> 8;
		assert(r == 10);
#endif
	}
	void dbg_unsat(const PbClause & c, const vec<Lit> & clause) {
#ifndef NDEBUG
		bool rhs_val = true;
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				if (l == c.rhs.lit) {
					rhs_val = true;
				} else {
					assert(l == ~c.rhs.lit);
					rhs_val = false;
				}
			}
		}
		FILE *f = fopen("testa.opb", "w");
		for (PbElement & e : c.clause) {
			fprintf(f, "%s%d x%d ",  "+", e.weight, var(e.lit) + 1);
		}
		
		if (rhs_val) {
			fprintf(f, " >= %d ;\n", c.rhs.weight);
		} else {
			fprintf(f, " < %d ;\n", c.rhs.weight);
		}
		
		fflush(f);
		// int r =system("minisat+ testa.opb >/dev/null")>>8;
		// assert(r==10);
		
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				
			} else {
				bool found = false;
				bool sign = false;
				for (PbElement & e : c.clause) {
					if (e.lit == l) {
						found = true;
					} else if (e.lit == ~l) {
						sign = true;
						found = true;
					}
				}
				assert(found);
				if (!sign) {
					fprintf(f, "1 x%d >= 1 ;\n", var(l) + 1);
				} else {
					fprintf(f, "1 x%d < 1 ;\n", var(l) + 1);
				}
				
			}
		}
		fflush(f);
		fclose(f);
		int r = system("minisat+ testa.opb >/dev/null") >> 8;
		assert(r == 20);
#endif
	}
	void dbg_sat(const PbClause & c, const vec<Lit> & clause) {
#ifndef NDEBUG
		bool rhs_val = true;
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				if (l == c.rhs.lit) {
					rhs_val = true;
				} else {
					assert(l == ~c.rhs.lit);
					rhs_val = false;
				}
			}
		}
		FILE *f = fopen("testa.opb", "w");
		for (PbElement & e : c.clause) {
			fprintf(f, "%s%d x%d ",  "+", e.weight, var(e.lit) + 1);
		}
		
		if (rhs_val) {
			fprintf(f, " >= %d ;\n", c.rhs.weight);
		} else {
			fprintf(f, " < %d ;\n", c.rhs.weight);
		}
		
		fflush(f);
		// int r =system("minisat+ testa.opb >/dev/null")>>8;
		// assert(r==10);
		
		for (Lit l : clause) {
			if (var(l) == var(c.rhs.lit)) {
				
			} else {
				bool found = false;
				bool sign = false;
				for (PbElement & e : c.clause) {
					if (e.lit == l) {
						found = true;
					} else if (e.lit == ~l) {
						sign = true;
						found = true;
					}
				}
				assert(found);
				if (!sign) {
					fprintf(f, "1 x%d >= 1 ;\n", var(l) + 1);
				} else {
					fprintf(f, "1 x%d < 1 ;\n", var(l) + 1);
				}
				
			}
		}
		fflush(f);
		fclose(f);
		int r = system("minisat+ testa.opb >/dev/null") >> 8;
		assert(r == 10);
#endif
	}
	
public:
	void buildReason(Lit p, vec<Lit> & reason) {
		backtrackUntil(p);
		assert(value(p)==l_True);
		CRef marker = S->reason(var(toSolver(p)));
		assert(marker != CRef_Undef);
		assert(reasonMap.has(marker));
		stats_reasons++;
		
		int clauseID = reasonMap[marker];
		static int iter = 0;
		if (++iter == 13) {
			int a = 1;
		}
		PbClause & pbclause = clauses[clauseID];
		reason.push(p);
		if (var(p) == var(pbclause.rhs.lit)) {
			if (sign(p) == sign(pbclause.rhs.lit))
				buildSumGEReason(clauseID, reason);
			else {
				buildSumLTReason(clauseID, reason);
			}
			
		} else {
			bool foundNeg = false;
			//it would be preferable to avoid this loop...
			for (int i = 0; i < pbclause.clause.size(); i++) {
				if (pbclause.clause[i].lit == p) {
					break;
				} else if (pbclause.clause[i].lit == ~p) {
					foundNeg = true;
					break;
				}
			}
			//this must be an element of the clause that was forced
			if (foundNeg)
				buildElementForcedFalseReason(clauseID, p, reason);
			else {
				
				buildElementForcedTrueReason(clauseID, p, reason);
			}
		}
		
		dbg_prove(pbclause, reason);
		//if( iter==13){
		dbg_min_conflict(pbclause, reason);
		//}
		toSolver(reason);
	}
	
	void printStats(int detailLevel) {
		printf("PbTheory: %d clauses\n", clauses.size());
		printf("%ld propagations (%ld skipped)\n", stats_propagations, stats_propagations_skipped);
		printf("%ld conflicts, %ld reasons\n", stats_conflicts, stats_reasons);
		printf("Shrink removed %ld lits from conflict clauses\n", stats_shrink_removed);
	}
	void preprocess() {
		
	}
	void implementConstraints() {
		for (ConstraintToImplement &c : constraintsToImplement) {
			if (c.implemented)
				continue;
			c.implemented = true;
			//convert to >= constraints
			if (c.rhs_lit == lit_Undef) {
				//this is an unconditional constraint
				c.rhs_lit = mkLit(S->newVar());
				S->addClause(c.rhs_lit); 		//constraint must hold.
						
			}
			if (c.op == PbType::GE) {
				greaterThanEq_implement(c.clause, c.weights, c.rhs_lit, c.total, c.side);
			} else if (c.op == PbType::GT) {
				greaterThanEq_implement(c.clause, c.weights, c.rhs_lit, c.total + 1, c.side);
			} else if (c.op == PbType::LT) {
				greaterThanEq_implement(c.clause, c.weights, ~c.rhs_lit, c.total, invert(c.side));
			} else if (c.op == PbType::LE) {
				greaterThanEq_implement(c.clause, c.weights, ~c.rhs_lit, c.total + 1, invert(c.side));
			} else if (c.op == PbType::EQ) {
				greaterThanEq_implement(c.clause, c.weights, c.rhs_lit, c.total, c.side);
				greaterThanEq_implement(c.clause, c.weights, ~c.rhs_lit, c.total + 1, invert(c.side));
				
			} else if (c.op == PbType::NE) {
				assert(false); 		//not yet supported..
			}
			
		}
		// constraintsToImplement.clear();
	}
	
	//If the rhs_lit is lit_Undef, then this is an unconditional constraint
	void addConstraint(vec<Lit> & clause, vec<int> & weights, int rhs, Lit rhs_lit = lit_Undef, PbType op = PbType::GE,
			ConstraintSide side = ConstraintSide::Both) {
		constraintsToImplement.push();
		assert(clause.size());
		assert(clause.size() == weights.size());
		clause.copyTo(constraintsToImplement.last().clause);
		weights.copyTo(constraintsToImplement.last().weights);
		constraintsToImplement.last().rhs_lit = rhs_lit;
		constraintsToImplement.last().total = rhs;
		constraintsToImplement.last().op = op;
		constraintsToImplement.last().side = side;
	}
	
	void backtrackUntil(int level) {
		//  bool changed=false;
		for (int clauseID : inq) {
			clauses[clauseID].inQueue = false;
		}
#ifndef NDEBUG
		for (PbClause & c : clauses) {
			assert(!c.inQueue);
		}
#endif
		inq.clear();
		//need to remove and add edges in the two graphs accordingly.
		if (trail_lim.size() > level) {
			int stop = trail_lim[level];
			for (int i = trail.size() - 1; i >= trail_lim[level]; i--) {
				Lit l = trail[i];
				
				assigns[var(l)] = l_Undef;
				if (!vars[var(l)].rhs) {
					Var v = var(l);
					bool s = vars[v].sign;
					PbClause & c = clauses[vars[v].clauseID];
					c.unassigned += vars[v].weight;
					if (s == sign(l)) {
						c.under -= vars[v].weight;
					} else {
						
					}
				}
			}
			trail.shrink(trail.size() - stop);
			trail_lim.shrink(trail_lim.size() - level);
			assert(trail_lim.size() == level);
			if (qhead > trail.size())
				qhead = trail.size();
		}
	}
	
	void backtrackUntil(Lit p) {
		for (int clauseID : inq) {
			clauses[clauseID].inQueue = false;
		}
#ifndef NDEBUG
		for (PbClause & c : clauses) {
			assert(!c.inQueue);
		}
#endif
		inq.clear();
		//need to remove and add edges in the two graphs accordingly.
		if (value(p) != l_Undef) {
			int i;
			for (i = trail.size() - 1; i >= 0; i--) {
				Lit l = trail[i];
				
				if (l == p) {
					i++;
					break;
				} else {
					assert(var(l) != var(p));
				}
				
				assigns[var(l)] = l_Undef;
				if (!vars[var(l)].rhs) {
					Var v = var(l);
					bool s = vars[v].sign;
					PbClause & c = clauses[vars[v].clauseID];
					c.unassigned += vars[v].weight;
					if (s == sign(l)) {
						c.under -= vars[v].weight;
					} else {
						
					}
				}
			}
			trail.shrink(trail.size() - i);
			while (trail_lim.size() && trail_lim.last() >= trail.size())
				trail_lim.pop();
			assert(trail_lim.size() == decisionLevel());
			if (qhead > trail.size())
				qhead = trail.size();
		}
	}
	bool check_solved() {
		for (int i = 0; i < clauses.size(); i++) {
			PbClause & c = clauses[i];
			Lit rhs_lit = c.rhs.lit;
			lbool rhs_val = value(rhs_lit);
			if (c.side == ConstraintSide::Upper && rhs_val == l_True) {
				continue; 		//this constraint is unenforced
			} else if (c.side == ConstraintSide::Lower && rhs_val == l_False) {
				continue; 		//this constraint is unenforced
			}
			
			long lhs = 0;
			for (PbElement & e : c.clause) {
				if (value(e.lit) == l_True) {
					lhs += e.weight;
				}
			}
			if (rhs_val == l_True && lhs < c.rhs.weight && c.side != ConstraintSide::Lower) {
				assert(false);
				return false;
			} else if (rhs_val == l_False && lhs >= c.rhs.weight && c.side != ConstraintSide::Upper) {
				assert(false);
				return false;
			}
			
		}
		
		for (ConstraintToImplement & c : constraintsToImplement) {
			
			Lit rhs_lit = c.rhs_lit;
			lbool rhs_val = S->value(rhs_lit);
			if (c.side == ConstraintSide::Upper && rhs_val == l_False) {
				continue; 		//this constraint is unenforced
			} else if (c.side == ConstraintSide::Lower && rhs_val == l_True) {
				continue; 		//this constraint is unenforced
			}
			
			long lhs = 0;
			for (int i = 0; i < c.clause.size(); i++) {
				Lit l = c.clause[i];
				if (S->value(l) == l_True) {
					lhs += c.weights[i];
				}
			}
			if (c.op == PbType::GE) {
				if (rhs_val == l_True && lhs < c.total && c.side != ConstraintSide::Lower) {
					assert(false);
					return false;
				} else if (rhs_val == l_False && lhs >= c.total && c.side != ConstraintSide::Upper) {
					assert(false);
					return false;
				}
			} else if (c.op == PbType::GT) {
				if (rhs_val == l_True && lhs <= c.total && c.side != ConstraintSide::Lower) {
					assert(false);
					return false;
				} else if (rhs_val == l_False && lhs > c.total && c.side != ConstraintSide::Upper) {
					assert(false);
					return false;
				}
			} else if (c.op == PbType::LE) {
				if (rhs_val == l_True && lhs > c.total && c.side != ConstraintSide::Lower) {
					assert(false);
					return false;
				} else if (rhs_val == l_False && lhs <= c.total && c.side != ConstraintSide::Upper) {
					assert(false);
					return false;
				}
			} else if (c.op == PbType::LT) {
				if (rhs_val == l_True && lhs >= c.total && c.side != ConstraintSide::Lower) {
					assert(false);
					return false;
				} else if (rhs_val == l_False && lhs < c.total && c.side != ConstraintSide::Upper) {
					assert(false);
					return false;
				}
			} else if (c.op == PbType::EQ) {
				if (rhs_val == l_True && lhs != c.total && c.side != ConstraintSide::Lower) {
					assert(false);
					return false;
				} else if (rhs_val == l_False && lhs == c.total && c.side != ConstraintSide::Upper) {
					assert(false);
					return false;
				}
			}
		}
		
		return true;
		
	}
	
	bool solveTheory(vec<Lit> & conflict) {
		return propagateTheory(conflict);
	}
	
	void printSolution() {
		
	}
	
	static unsigned int gcd(unsigned int a, unsigned int b) {
		//Adapted from lingeling
		assert(a > 0);
		assert(b > 0);
		if (a < b)
			std::swap(a, b);
		while (b > 0) {
			unsigned int tmp = b;
			b = a % b;
			a = tmp;
		}
		return a;
	}
	
	static unsigned int gcd(vec<int> & of) {
		if (of.size() == 0) {
			return 0;
		} else if (of.size() == 1) {
			return (of[0]);
		}
		
		int curgcd = gcd(of[0], of[1]);
		for (int i = 2; i < of.size(); i++) {
			curgcd = gcd(curgcd, of[i]);
		}
		return curgcd;
	}
	
};

}
;

#endif /* PBTHEORY_H_ */
