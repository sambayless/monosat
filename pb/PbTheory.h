/*
 * PbTheory.h
 *
 *  Created on: May 3, 2014
 *      Author: sam
 */

#ifndef PBTHEORY_H_
#define PBTHEORY_H_



#include "mtl/Vec.h"
#include "mtl/Map.h"
#include "core/SolverTypes.h"
#include "core/Theory.h"
#include <cmath>
namespace Minisat {


class PbTheory: public Theory{
	Solver * S;
	int theory_index;
	int qhead=0;
	struct Assignment{
		int isLit:1;
		int lit_clauseID:31;
	};

	vec<Assignment> trail;
	vec<int> trail_lim;

	enum class PbType{
		GT, GE,LT,LE,EQ
	};

	struct ConstraintToImplement{
		vec<Lit> clause;
		vec<int> weights;
		Lit rhs_lit;
		int total;
		bool oneSided;
		PbType op;
	};

	vec<ConstraintToImplement> constraintsToImplement;

	struct VarData{
		//bool onesided;//whether this is a one-sided or two-sided pb constraint (if one sided, then if the output var is false, the constraint is simply unenforced)
		Var solverVar;
		int rhs:1;
		int clauseID:31;
		unsigned int weight;
		VarData():solverVar(var_Undef),rhs(false),clauseID(-1),weight(0){}
	};
	struct PbElement{
		Lit lit;
		unsigned int weight;
	};
	struct PbClause{
		bool oneSided;
		bool isSatisfied;
		bool isOverEq; //true iff the sum of the lhs is >= the rhs
		PbElement rhs;
		CRef reason;
		vec<PbElement> clause;
		PbClause():oneSided(false),isSatisfied(false){}
	};
	int firstCRef;
	//Map reasons to clauseIDs
	Map<CRef,int> reasonMap;
	vec<PbClause> clauses;
	vec<VarData> vars;
	vec<Lit> clause;
	vec<int> weights;

	 double propagationtime;
	 long stats_propagations,stats_propagations_skipped;

/*	 Var newVar(){
	 		Var s= S->newVar();
	 		Var v = vars.size();
	 		vars.push();
	 		vars[v].solverVar=s;

	 		return v;
	 	}*/
	 	Var newVar(Var solverVar, int clauseID, int weight, bool isRHS =false){
	 		while(S->nVars()<=solverVar)
	 				S->newVar();

	 		if(S->hasTheory(solverVar)){
	 			//introduce a new replacement variable
			  Var v2 = S->newVar();
				S->addClause(~mkLit(solverVar),mkLit(v2));
				S->addClause(mkLit(solverVar),~mkLit(v2));
				solverVar=v2;
	 		}

	 		Var v = vars.size();
	 		vars.push();
	 		vars[v].solverVar=solverVar;
	 		vars[v].clauseID=clauseID;
	 		vars[v].weight=weight;
	 		vars[v].rhs=isRHS;

	 		S->setTheoryVar(solverVar,getTheoryIndex(),v);
	 		assert(toSolver(v)==solverVar);

	 		return v;
	 	}
	 	inline Var toSolver(Var v){
	 			//return v;
	 			assert(v<vars.size());
	 			assert(S->hasTheory(vars[v].solverVar));
	 			assert(S->getTheoryVar(vars[v].solverVar)==v);
	 			return vars[v].solverVar;
	 		}

	 		inline Lit toSolver(Lit l){
	 			assert(S->hasTheory(vars[var(l)].solverVar));
	 			assert(S->getTheoryVar(vars[var(l)].solverVar)==var(l));
	 			return mkLit(vars[var(l)].solverVar,sign(l));
	 		}

	 		void toSolver(vec<Lit> & c){
	 			for(int i = 0;i<c.size();i++){
	 				c[i]=toSolver(c[i]);
	 			}
	 		}
	 		int level(Var v){
	 			return S->level(toSolver(v));
	 		}
	 		inline lbool value(Var v){
	 			return S->value(toSolver(v));
	 		}
	 		inline lbool value(Lit l){
	 			return S->value(toSolver(l));
	 		}
	 		inline bool enqueue(Lit l, CRef reason){


	 			Lit sl = toSolver(l);
	 			if( S->enqueue(sl,reason)){
	 				enqueueTheory(l);
	 				return true;
	 			}else{
	 				return false;
	 			}
	 		}
	 		void makeEqual(Lit l1, Lit l2){
	 			Lit o1 = toSolver(l1);
	 			Lit o2 = toSolver(l2);
	 			S->addClause(~o1,o2);
	 			S->addClause(o1, ~o2);
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
	 			vec<Lit> t;
	 			c.copyTo(t);
	 			toSolver(t);
	 			S->addClause(t);
	 		}
public:

	 PbTheory(Solver * S):S(S),theory_index(0),propagationtime(0),stats_propagations(0),stats_propagations_skipped(0){

	 }
     ~PbTheory(){};

 	inline int getTheoryIndex(){
 	    	return theory_index;
 	    }
 	  inline  void setTheoryIndex(int id){
 	    	theory_index=id;
 	    }
 		void newDecisionLevel(){
 			trail_lim.push(trail.size());
 		};
 		inline int decisionLevel(){
 			return trail_lim.size(); //S->decisionLevel();
 		}
 	 void enqueueTheory(Lit l){
 	 		Var v = var(l);

 	 		int lev = level(v);

 	 		assert(decisionLevel()<=lev);

 	 		while(lev>trail_lim.size()){
 	 			newDecisionLevel();
 	 		}

			trail.push({true,toInt(l)});
 	 	};
 	 	bool propagateTheory(vec<Lit> & conflict){
 	 		assert(qhead<=trail.size());
 	 		if(qhead==trail.size()){
 	 			stats_propagations_skipped++;
 	 			return true;
 	 		}

 	 		stats_propagations++;

 	 		static int iter = 0;
 	 		double startproptime= rtime(2);

 	 		//This is wrong! Only need to visit each _clause_ that has any involved literals once per propagation round.
 	 		for(;qhead<trail.size();qhead++){

 	 			if(!trail[qhead].isLit)
 	 				continue;
 	 			Lit l2 = toLit(trail[qhead].lit_clauseID);

 	 			if(++iter==6){
 	 				int a=1;
 	 			}
 	 			Var v2= var(l2);

 	 			assert(vars[v2].clauseID>=0);
 	 			int clauseID = vars[v2].clauseID;
 	 			if(clauseID==1){
 	 				int a=1;
 	 			}
 	 			PbClause & pbclause = clauses[clauseID];
 	 			if(pbclause.isSatisfied)
 	 				continue;

 	 			Lit rhs = pbclause.rhs.lit;
 	 			lbool rhs_val = value(rhs);
 	 			if(pbclause.oneSided && rhs_val==l_False){
 	 				//this clause is free.
 	 				pbclause.isSatisfied=true;
 	 				trail.push({false,clauseID});
 	 				continue;
 	 			}

 	 			if(v2 == var(pbclause.rhs.lit)){


 	 			}else{

					unsigned int total = pbclause.rhs.weight;
					//compute over and under approximations...
					unsigned int underApprox=0;
					unsigned int overApprox=0;
					int n_Free=0;
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
								largestUnassigned = l;
							}
							if (e.weight<=smallestUnassignedWeight){
								smallestUnassignedWeight=e.weight;
							}
						}else{
							int  a=1;
						}
					}

					overApprox=underApprox+unassignedWeight;
					if(rhs_val==l_True && underApprox>=total){
						//this is a satisfied constraint
						pbclause.isSatisfied=true;
						trail.push({false,clauseID});
					}else if(rhs_val==l_False && overApprox<total){
						//this is a satisfied constraint
						pbclause.isSatisfied=true;
						trail.push({false,clauseID});
					}else if(rhs_val==l_True && overApprox<total){
						//conflict
						conflict.push(~rhs);
						buildSumLTReason(clauseID,conflict);
						 dbg_prove(pbclause,conflict);
						toSolver(conflict);
						return false;
					}else if (rhs_val==l_False && underApprox>=total){
						//conflict
						assert(!pbclause.oneSided);
						conflict.push(rhs);

						buildSumGEReason(clauseID,conflict);
						 dbg_prove(pbclause,conflict);
						toSolver(conflict);
						return false;
					}else if (underApprox>=total && rhs_val==l_Undef){
						pbclause.isSatisfied=true;
						pbclause.isOverEq=true;
						dbg_prop(pbclause,rhs);
						//trail.push({false,clauseID});
						enqueue(rhs,pbclause.reason);
					}else if (overApprox<total && rhs_val==l_Undef){
						pbclause.isSatisfied=true;
						pbclause.isOverEq=false;
						dbg_prop(pbclause,~rhs);
						//trail.push({false,clauseID});
						enqueue(~rhs,pbclause.reason);
					}else if(rhs_val==l_False ){
						pbclause.isOverEq=false;
						//then _all_ unassigned lits whose weight would push the under approx over over the limit must be assigned to false
						for(PbElement e:pbclause.clause){
							Lit l = e.lit;
							if(value(l)==l_Undef && underApprox+e.weight>=total){
								dbg_prop(pbclause,~l);
								enqueue(~l,pbclause.reason);
							}
						}

					}else if (rhs_val==l_True){
						assert(n_Free>0);
						//if(n_Free==1){
						if(overApprox - largestUnassignedWeight <total){
							//then the largest unassigned weight is forced
							//assert(smallestUnassignedWeight==largestUnassignedWeight);
							assert(largestUnassigned!=lit_Undef);
							//assert(underApprox+largestUnassignedWeight>=total);//else the over approx would have triggered a conflict
							dbg_prop(pbclause,largestUnassigned);
							enqueue(largestUnassigned,pbclause.reason);

							//it may also be the case that the second largest weight is forced.
							if(n_Free>1){
								for(PbElement e:pbclause.clause){
									Lit l = e.lit;
									lbool val = value(l);
									if(val==l_True){

									}else if(val==l_Undef){
										if(overApprox-e.weight<total){
											dbg_prop(pbclause,e.lit);
											enqueue(e.lit,pbclause.reason);
										}
									}
								}
							}
						}
					}

 	 			}

 	 		}


 	 		double elapsed = rtime(2)-startproptime;
 	 		propagationtime+=elapsed;
 	 		dbg_fully_propped();
 	 		return true;
 	 	};

private:
 	 		void buildElementForcedFalseReason(int clauseID, Lit element,vec<Lit> & conflict){
 	 			//the reason why a lit is forced to false is the set of true literals in the clause, and the rhs
 	 			PbClause & pbclause = clauses[clauseID];
				assert(pbclause.isSatisfied);
				assert(value(pbclause.rhs.lit)==l_False);

				conflict.push(pbclause.rhs.lit);

				int underApprox = 0;
				int unassignedWeight = 0;
				int forcedWeight = 0;
				for(PbElement e:pbclause.clause){
					Lit l = e.lit;
					if(e.lit ==element){
						forcedWeight = e.weight;
						assert(value(e.lit)==l_False);
					}else{
						assert(var(element)!=var(e.lit));
					}

					if(value(l)==l_True){
						if(level(var(l))>0)
							conflict.push(~l);
						underApprox +=e.weight;
					}
				}
				assert(forcedWeight+underApprox>=pbclause.rhs.weight);
 	 		}
 	 		void buildElementForcedTrueReason(int clauseID, Lit element,vec<Lit> & conflict){
 	 			//the reason why a lit is forced to false is the set of true literals in the clause, and the rhs
 	 			PbClause & pbclause = clauses[clauseID];
				assert(pbclause.isSatisfied);
				assert(value(pbclause.rhs.lit)==l_True);
				conflict.push(~pbclause.rhs.lit);

				int overApprox = 0;
				int forcedWeight = 0;
				for(PbElement e:pbclause.clause){
					Lit l = e.lit;
					if(e.lit ==element){
						forcedWeight = e.weight;
						assert(value(e.lit)==l_True);
					}else{
						assert(var(element)!=var(e.lit));
						assert(value(l)!=l_Undef);
					}

					if(value(l)==l_False){
						if(level(var(l))>0)
							conflict.push(l);
					}else{
						overApprox +=e.weight;
					}

				}
				assert(overApprox>=pbclause.rhs.weight);
				assert(overApprox-forcedWeight <pbclause.rhs.weight);

 	 		}
 	 		void buildSumGEReason(int clauseID,vec<Lit> & conflict){
 	 			PbClause & pbclause = clauses[clauseID];
 	 			assert(!pbclause.isSatisfied);

 	 			int underApprox = 0;
 	 			int unassignedWeight = 0;
 	 			for(PbElement e:pbclause.clause){
						Lit l = e.lit;
						if(value(l)==l_True){
							if(level(var(l))>0)
								conflict.push(~l);
							underApprox +=e.weight;
						}else if(value(l)==l_Undef){
							unassignedWeight+=e.weight;
						}
					}
 	 			assert(underApprox>=pbclause.rhs.weight);
 	 		}

 	 		void buildSumLTReason(int clauseID,vec<Lit> & conflict){
 	 			PbClause & pbclause = clauses[clauseID];
 	 			assert(!pbclause.isSatisfied);

 	 			int underApprox = 0;
 	 			int unassignedWeight = 0;
 	 			for(PbElement e:pbclause.clause){
						Lit l = e.lit;
						if(value(l)==l_True){
							underApprox +=e.weight;
						}else if(value(l)==l_Undef){
							unassignedWeight+=e.weight;
						}else if(value(l)==l_False){
							if(level(var(l))>0)
								conflict.push(l);
						}
					}
 	 			assert(underApprox+unassignedWeight<pbclause.rhs.weight);
 	 		}

 	 		//Follow the normalization guidelines from the minisatp paper.
 	 		bool normalize(vec<Lit> & clause,vec<int> & weights,  Lit & rhs_lit,int &rhs, bool oneSided=false){
 	 			int i,j=0;
 	 			assert(clause.size()==weights.size());
 	 			assert(S->decisionLevel()==0);
 	 			for(i = 0;i<clause.size();i++){
 	 				Lit l = clause[i];
 	 				int w = weights[i];

 	 				if(S->value(l)==l_True){
 	 					//drop this literal, and subtract is total from the rhs
 	 					rhs-=w;
 	 				}else if (S->value(l)==l_False){
 	 					//then drop this literal
 	 				}else if(w==0){
 	 					//then drop this literal
 	 				}else if (w<0){
 	 					//then invert this literal and update the right hand side
 	 					weights[j]=-w;
 	 					clause[j++]=~l;
 	 					rhs+=-w;
 	 				}else{
 	 					//don't change
 	 					weights[j]=w;
 	 					clause[j++]=l;
 	 				}
 	 			}

 	 			clause.shrink(i-j);
 	 			weights.shrink(i-j);

 	 			if(rhs<=0){
 	 				//then this is a trivially true constraint (even if the total weights are 0, because this is a >= constraint)
 	 				clause.clear();
 	 				weights.clear();
 	 				S->addClause(rhs_lit);
 	 				return false;
 	 			}

 	 			for(int i = 0;i<weights.size();i++){
 	 				assert(weights[i]>0);
 	 				if(weights[i]>rhs){
 	 					weights[i]=rhs;
 	 				}
 	 			}

 	 			int weightTotal = 0;
 	 			bool allWeights1=true;
 	 			//as per the minisatp paper, divide the left and right hand sides by the gcd of the left hand side.
 	 			int d = gcd(weights);
 	 			for(int i = 0;i<weights.size();i++){
 	 				weights[i]/=d;
 	 				weightTotal +=weights[i];
 	 				allWeights1 &= (weights[i]==1);
 	 			}

 	 			int rhs2 = (int)ceil((double)rhs / (double)d);
 	 			assert(rhs/d == rhs2 || rhs/d==rhs2-1);
 	 			rhs=rhs2;

 	 			allWeights1 &= (rhs==1);

 	 			if(weightTotal<rhs){
 	 				//then this is a trivially unsat constraint
 	 				clause.clear();
 	 				weights.clear();
 	 				S->addClause(~rhs_lit);
 	 				return false;
 	 			}
 	 			assert(clause.size()>0);//else weightTotal is 0;

 	 			//It is always the case that at most one element of the clause must be true, OR the value must be false.
 	 			//so we can learn this clause right away
 				clause.push(~rhs_lit);
 				S->addClause(clause);
				clause.pop();

				if(!oneSided){
					//then enforce the other side of the constraint
					for(int i = 0;i<clause.size();i++){
						if(weights[i]>=rhs){
							S->addClause(rhs_lit, ~clause[i]);//If value is false, then all the lits in the clause that have >= weight as rhs must be false.
						}
					}
				}


 	 			if(allWeights1){
 	 				//this as a pure SAT constraint
 	 				clause.clear();
					weights.clear();
 	 				return false;
 	 			}else{
 	 				assert(clause.size()>1);//because we divided by gcd...
 	 			}

 	 			// the minisatp paper suggests splitting the clause into a pb part and a pure-sat part by introducing a new variable
 	 			// might want to implement that in the future...

 	 			return true;
 	 		}


 	  Lit greaterThanEq_implement(const vec<Lit> & _clause,const vec<int> &_weights,  Lit rhs_lit,int total, bool oneSided=false){
 		  _clause.copyTo(clause);
 		  _weights.copyTo(weights);
		  if(normalize(clause,weights,rhs_lit,total)){

			 for(int w:weights){
				 assert(w>0);

			 }
			  int clauseID = clauses.size();
			  assert(total>0);

			 assert(weights.size()==clause.size());

			  clauses.push();
			  CRef reason = S->newReasonMarker(this);
			  PbClause & pbclause = clauses.last();
			  rhs_lit = mkLit(newVar(var(rhs_lit),clauseID,total,true) ,sign(rhs_lit));
			  pbclause.rhs.lit = rhs_lit;
			  pbclause.rhs.weight=(unsigned int) total;
			  pbclause.oneSided = oneSided;
			  pbclause.reason = reason;
			  reasonMap.insert(reason,clauseID);



/*

			  if(vars[var(rhs_lit)].clauseID>=0){
				  //then create a new var
				  Lit l = mkLit(newVar());
				  makeEqual(l,rhs_lit);
				  rhs_lit=l;
			  }
*/

			  vars[var(rhs_lit)].clauseID = clauseID;
			  vars[var(rhs_lit)].rhs = true;
			  vars[var(rhs_lit)].weight=total;
			  for(int i = 0;i<clause.size();i++){
				  Lit l = clause[i];
				  l = mkLit(newVar(var(l),clauseID,weights[i]) ,sign(l));

				  pbclause.clause.push({l,(unsigned int)weights[i]});
				 /* if(vars[var(rhs_lit)].clauseID>=0){
					  //then create a new var
					  Lit l2 = mkLit(newVar());
					  makeEqual(l,l2);
					  l=l2;
				  }*/
	/*			  Var v= var(l);
				  vars[v].clauseID = clauseID;
				  vars[v].rhs = false;
				  vars[v].weight=total;*/
			  }
		  }
		  return rhs_lit;

	 }
 	 void dbg_fully_propped(){
#ifndef NDEBUG
 		for(PbClause & c:clauses){
 			dbg_fully_propped(c);

 		}
#endif
 	 }

 	  void dbg_fully_propped(const PbClause & c){
#ifndef NDEBUG

 		  vec<Lit> prove;
 		  for(PbElement p:c.clause){
 			  if(value(p.lit)==l_False){
 				  prove.push(~p.lit);
 			  }else if(value(p.lit)==l_True){
 				  prove.push(p.lit);
 			  }
 		  }

 		  if(value(c.rhs.lit)==l_False){
 			  prove.push(~c.rhs.lit);
		  }else if(value(c.rhs.lit)==l_True){
			  prove.push(c.rhs.lit);
		  }


 		 dbg_sat(c,prove);


		  for(PbElement p:c.clause){
			 // if(var(p.lit)==23 || toInt(p.lit)==23){


			  if(value(p.lit)==l_Undef){
				  prove.push(p.lit);
				  dbg_sat(c,prove);
				  prove.pop();
				  prove.push(~p.lit);
				  dbg_sat(c,prove);
				  prove.pop();

			  }
			// }
		  }

		  if(value(c.rhs.lit)==l_Undef){
			  prove.push(c.rhs.lit);
				  dbg_sat(c,prove);
				  prove.pop();
				  prove.push(~c.rhs.lit);
				  dbg_sat(c,prove);
				  prove.pop();
		  }
#endif
 	  }

 	  void dbg_prop(const PbClause & c,Lit e){
#ifndef NDEBUG
 		  vec<Lit> prove;
 		  for(PbElement p:c.clause){
 			  if(value(p.lit)==l_False){
 				  prove.push(p.lit);
 			  }else if(value(p.lit)==l_True){
 				  prove.push(~p.lit);
 			  }
 		  }

 		  if(value(c.rhs.lit)==l_False){
 			  prove.push(c.rhs.lit);
		  }else if(value(c.rhs.lit)==l_True){
			  prove.push(~c.rhs.lit);
		  }

 		 prove.push(e);
 		 dbg_prove(c,prove);

#endif
 	  }

 	  void dbg_prove(const PbClause & c, const vec<Lit> & clause){
#ifndef NDEBUG
 		 bool rhs_val = true;
 		 for(Lit l:clause){
 			 if(var(l)==var(c.rhs.lit)){
 				 if(l==c.rhs.lit){
 					 rhs_val=false;
 				 }else{
 					 assert(l==~c.rhs.lit);
 					 rhs_val=true;
 				 }
 			 }
 		 }
 		 FILE *f =  fopen("test.opb","w");
 		 for(PbElement & e:c.clause){
 			 fprintf(f,"%s%d x%d ",e.weight>=0 ? "+":"",e.weight,var(e.lit)+1);
 		 }

 		 if(rhs_val){
 			 fprintf(f," >= %d ;\n",c.rhs.weight);
 		 }else{
 			 fprintf(f," < %d ;\n",c.rhs.weight);
 		 }

 		 fflush(f);
 		 int r =system("minisat+ test.opb >/dev/null")>>8;
 		 assert(r==10);


 		 for(Lit l:clause){
 			 if(var(l)==var(c.rhs.lit)){

 			 }else{
 				 bool found=false;
 				 bool sign = false;
 				 for(PbElement & e:c.clause){
 					 if(e.lit == l){
 						 found=true;
 					 }else if(e.lit==~l){
 						 sign=true;
 						 found=true;
 					 }
 				 }
 				 assert(found);
 				 if(!sign){
 					 fprintf(f,"1 x%d < 1 ;\n",var(l)+1);
 				 }else{
 					 fprintf(f,"1 x%d >= 1 ;\n",var(l)+1);
 				 }

 			 }
 		 }
 		 fflush(f);
 		 fclose(f);
 		 r =system("minisat+ test.opb >/dev/null")>>8;
 		 assert(r==20);
#endif
 	  }

 	  void dbg_sat(const PbClause & c, const vec<Lit> & clause){
#ifndef NDEBUG
 		 bool rhs_val = true;
 		 for(Lit l:clause){
 			 if(var(l)==var(c.rhs.lit)){
 				 if(l==c.rhs.lit){
 					 rhs_val=false;
 				 }else{
 					 assert(l==~c.rhs.lit);
 					 rhs_val=true;
 				 }
 			 }
 		 }
 		 FILE *f =  fopen("test.opb","w");
 		 for(PbElement & e:c.clause){
 			 fprintf(f,"%s%d x%d ",e.weight>=0 ? "+":"",e.weight,var(e.lit)+1);
 		 }

 		 if(rhs_val){
 			 fprintf(f," >= %d ;\n",c.rhs.weight);
 		 }else{
 			 fprintf(f," < %d ;\n",c.rhs.weight);
 		 }

 		 fflush(f);
 		// int r =system("minisat+ test.opb >/dev/null")>>8;
 		// assert(r==10);


 		 for(Lit l:clause){
 			 if(var(l)==var(c.rhs.lit)){

 			 }else{
 				 bool found=false;
 				 bool sign = false;
 				 for(PbElement & e:c.clause){
 					 if(e.lit == l){
 						 found=true;
 					 }else if(e.lit==~l){
 						 sign=true;
 						 found=true;
 					 }
 				 }
 				 assert(found);
 				 if(!sign){
 					 fprintf(f,"1 x%d >= 1 ;\n",var(l)+1);
 				 }else{
 					 fprintf(f,"1 x%d < 1 ;\n",var(l)+1);
 				 }

 			 }
 		 }
 		 fflush(f);
 		 fclose(f);
 		 int r =system("minisat+ test.opb >/dev/null")>>8;
 		 assert(r==10);
#endif
 	  }
public:
 	 void buildReason(Lit p, vec<Lit> & reason){
 		 CRef marker = S->reason(var(toSolver(p)));
 		 assert(marker!=CRef_Undef);
 		 assert(reasonMap.has(marker));
 		 int clauseID = reasonMap[marker];
 		 PbClause & pbclause = clauses[clauseID];
 		reason.push(p);
 		 if(var(p)==var(pbclause.rhs.lit)){
 			 if(value(p)==l_False)
 				buildSumLTReason(clauseID,reason);
			 else{
				 assert(value(p)==l_True);
				 buildSumGEReason(clauseID,reason);
			 }

 		 }else{
 			 //this must be an element of the clause that was forced
 			 if(value(p)==l_False)
 				 buildElementForcedFalseReason(clauseID,p,reason);
 			 else{
 				 assert(value(p)==l_True);
 				 buildElementForcedTrueReason(clauseID,p,reason);
 			 }
 		 }
 		 dbg_prove(pbclause,reason);
 		toSolver(reason);
 	 }

	 void printStats(int detailLevel){

	}
	 void preprocess(){

	}
	 void implementConstraints(){
		 for(ConstraintToImplement &c:constraintsToImplement){
			 //convert to >= constraints
			 if(c.rhs_lit ==lit_Undef){
				 //this is an unconditional constraint
				 c.rhs_lit = mkLit(S->newVar());
				 S->addClause(c.rhs_lit);//constraint must hold.
			 }
			 if(c.op==PbType::GE){
				 greaterThanEq_implement(c.clause,c.weights,c.rhs_lit,c.total);
			 }else  if(c.op==PbType::GT){
				 greaterThanEq_implement(c.clause,c.weights,c.rhs_lit,c.total+1);
			 }else if (c.op==PbType::LT){

				 greaterThanEq_implement(c.clause,c.weights,~c.rhs_lit,c.total);
			 }else if(c.op==PbType::LE){

				 greaterThanEq_implement(c.clause,c.weights,~c.rhs_lit,c.total-1);
			 }else if(c.op==PbType::EQ){
				greaterThanEq_implement(c.clause,c.weights,~c.rhs_lit,c.total+1);
				 greaterThanEq_implement(c.clause,c.weights,c.rhs_lit,c.total+1);
			 }

		 }
		// constraintsToImplement.clear();
	 }

	 //If the rhs_lit is lit_Undef, then this is an unconditional constraint
	 void addConstraint(vec<Lit> & clause, vec<int> & weights, int rhs, Lit rhs_lit=lit_Undef , PbType op=PbType::GE, bool oneSided=false){
		 constraintsToImplement.push();
		 assert(clause.size());
		 assert(clause.size()==weights.size());
		 clause.copyTo(constraintsToImplement.last().clause);
		 weights.copyTo(constraintsToImplement.last().weights);
		 constraintsToImplement.last().rhs_lit = rhs_lit;
		 constraintsToImplement.last().total = rhs;
		 constraintsToImplement.last().op = op;
		 constraintsToImplement.last().oneSided=oneSided;
	 }

	  void backtrackUntil(int level){
		  bool changed=false;
			//need to remove and add edges in the two graphs accordingly.
			if(trail_lim.size()>level){
				int stop = trail_lim[level];
				for(int i = trail.size()-1;i>=trail_lim[level];i--){
					if(!trail[i].isLit){
						int clauseID = trail[i].lit_clauseID;
						assert(clauseID>=0);assert(clauseID<clauses.size());
						assert(clauses[clauseID].isSatisfied);
						clauses[clauseID].isSatisfied=false;
						clauses[clauseID].isOverEq=false;
					}
					changed=true;
				}
				trail.shrink(trail.size()-stop);
				trail_lim.shrink(trail_lim.size()-level);
				assert(trail_lim.size()==level);
				if(qhead>trail.size())
					qhead=trail.size();
			}
	  }
	  bool check_solved(){
		  for(int i = 0;i<clauses.size();i++){
			  PbClause & c = clauses[i];
			  Lit rhs_lit = c.rhs.lit;
			  lbool rhs_val = value(rhs_lit);
			  if(c.oneSided && rhs_val==l_Undef){
				 continue;//this constraint is unenforced
			  }

			  long lhs = 0;
			  for(PbElement & e:c.clause){
				  if(value(e.lit)==l_True){
					  lhs+=e.weight;
				  }
			  }
			  if(rhs_val==l_True && lhs<c.rhs.weight){
				  assert(false);
				  return false;
			  }else if(rhs_val==l_False && lhs>=c.rhs.weight){
				  assert(false);
				  return false;
			  }

		  }

		  for(ConstraintToImplement & c:constraintsToImplement){

			  Lit rhs_lit = c.rhs_lit;
			  lbool rhs_val = S->value(rhs_lit);
			  if(c.oneSided && rhs_val==l_False){
				 continue;//this constraint is unenforced
			  }

			  long lhs = 0;
			  for(int i = 0;i<c.clause.size();i++){
				  Lit l = c.clause[i];
				  if(S->value(l)==l_True){
					  lhs+=c.weights[i];
				  }
			  }
			  if(c.op==PbType::GE){
				  if(rhs_val==l_True && lhs<c.total){
					  assert(false);
					  return false;
				  }else if(rhs_val==l_False && lhs>=c.total){
					  assert(false);
					  return false;
				  }
			  }else if(c.op==PbType::GT){
				  if(rhs_val==l_True && lhs<=c.total){
					  assert(false);
					  return false;
				  }else if(rhs_val==l_False && lhs>c.total){
					  assert(false);
					  return false;
				  }
			  }else if(c.op==PbType::LE){
				  if(rhs_val==l_True && lhs>c.total){
					  assert(false);
					  return false;
				  }else if(rhs_val==l_False && lhs<=c.total){
					  assert(false);
					  return false;
				  }
			  }else if(c.op==PbType::LT){
				  if(rhs_val==l_True && lhs>=c.total){
					  assert(false);
					  return false;
				  }else if(rhs_val==l_False && lhs<c.total){
					  assert(false);
					  return false;
				  }
			  }else if(c.op==PbType::EQ){
				  if(rhs_val==l_True && lhs!=c.total){
					  assert(false);
					  return false;
				  }else if(rhs_val==l_False && lhs==c.total){
					  assert(false);
					  return false;
				  }
			  }
		  }

		  return true;

	  }

	 bool solveTheory(vec<Lit> & conflict){
		 return propagateTheory(conflict);
	 }

	 static unsigned  int gcd (unsigned int a, unsigned int b) {
	   //Adapted from lingeling
	   assert (a>0);
	   assert (b>0);
	   if (a < b)
		  std::swap (a, b);
	   while (b>0) {
		  unsigned int tmp = b;
		  b = a % b;
		  a = tmp;
	   }
	   return a;
	 }

	 static unsigned int gcd(vec<int> & of ){
		if(of.size()==0){
			return 0;
		}else if(of.size()==1){
			return(of[0]);
		}

		int curgcd = gcd(of[0],of[1]);
		for(int i = 2;i<of.size();i++){
			curgcd = gcd(curgcd,of[2]);
		}
		return curgcd;
	 }

};


};


#endif /* PBTHEORY_H_ */
