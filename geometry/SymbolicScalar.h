/*
 * SymbolicScalar.h
 *
 *  Created on: Sep 5, 2014
 *      Author: sam
 */

#ifndef SYMBOLICSCALAR_H_
#define SYMBOLICSCALAR_H_
#include "Polygon.h"
#include "GeometryDetector.h"
#include "GeometryTypes.h"
#include "core/SolverTypes.h"
#include "PointSet.h"
#include "GeometryDetector.h"
#include "core/Config.h"
#include "core/TheorySolver.h"
#include <vector>

//very simple symbolic rational scalar class.
//should be swapped out for (or interfaced to) a proper LRA solver.

template<class T>
class SymbolicScalar:public GeometryDetector{
protected:

	double rnd_seed;
	TheorySolver * S;
	CRef value_less_than=CRef_Undef;
	CRef value_less_than_equal_to=CRef_Undef;
	CRef value_greater_than_equal_to=CRef_Undef;
	CRef value_greater_than=CRef_Undef;



public:
	SymbolicScalar(int detectorID,TheorySolver * solver,  double seed=1):GeometryDetector(detectorID),S(solver),rnd_seed(seed){

	}

	virtual ~SymbolicScalar(){

	}
	virtual bool isConstant(){
		return false;
	}
	virtual T & getOverApprox()=0;
	virtual T & getUnderApprox()=0;
};

template<class T>
class ConstantScalar:public SymbolicScalar<T>{
	T value;

public:
	ConstantScalar(int detectorID,TheorySolver * solver, T value , double seed=1): SymbolicScalar<T>(detectorID,solver,seed),value(value){

	}
	bool isConstant(){
		return true;
	}

	void buildValueLessThanReason(const T & v, vec<Lit> & conflict, bool inclusive){
		if(inclusive){
			assert(value<=v);
		}else{
			assert(value<v);
		}
		//no justification set required
	}

	void buildValueGreaterThanReason(const T & v, vec<Lit> & conflict, bool inclusive){
		if(inclusive){
			assert(value>=v);
		}else{
			assert(value>v);
		}
		//no justification set required
	}

	T & getOverApprox(){
		return value;
	}
	T & getUnderApprox(){
		return value;
	}

};

template<class T>
class ConditionalScalar:public SymbolicScalar<T>{
	SymbolicScalar<T> & thn;
	SymbolicScalar<T> & els;
	Lit condition;

	T over_approx;
	T under_approx;
public:

	ConditionalScalar(Lit condition, SymbolicScalar<T> & thn, SymbolicScalar<T> & els,int detectorID,TheorySolver * solver,  double seed=1): SymbolicScalar<T>(detectorID,solver,seed), thn(thn),els(els){
		condition= mkLit( this->solver->newVar(var(condition),this->getID()),sign(condition));
	}
	bool isConstant(){
		if (this->solver->isConstant(condition)){
			if(this->solver->value(condition)==l_True){
				return thn.isConstant();
			}else{
				return els.isConstant();
			}
		}
		return false;
	}
	T & getOverApprox(){
		lbool val = this->solver->value(condition);
		if(val==l_True){
			return thn.getOverApprox();
		}else if (val==l_False){
			return els.getOverApprox();
		}

		return max(*thn.getOverApprox(),els.getOverApprox());
	}
	T & getUnderApprox(){
		lbool val = this->solver->value(condition);
		if(val==l_True){
			return thn.getUnderApprox();
		}else if (val==l_False){
			return els.getUnderApprox();
		}
		return min(*thn.getOverApprox(),els.getOverApprox());
	}


	void buildValueLessThanReason(const T & v, vec<Lit> & conflict, bool inclusive){
		if(inclusive){
			assert(getOverApprox() <= v);
		}else{
			assert(getOverApprox() < v);
		}
		lbool val = this->solver->value(condition);
		if(val==l_True){
			thn.buildValueLessThanReason(v,conflict,inclusive);
			conflict.push(~condition);
		}else if (val==l_False){
			els.buildValueLessThanReason(v,conflict,inclusive);
			conflict.push(condition);
		}else{
			assert(val==l_Undef);
			if(inclusive){
				if(thn.getOverApprox() <= v){
					thn.buildValueLessThanReason(v, conflict,inclusive);
				}else{
					assert(els.getOverApprox() <= v);
					els.buildValueLessThanReason(v, conflict,inclusive);
				}
			}else{
				if(thn.getOverApprox() < v){
					thn.buildValueLessThanReason(v, conflict,inclusive);
				}else{
					assert(els.getOverApprox() <v);
					els.buildValueLessThanReason(v, conflict,inclusive);
				}
			}

		}
	}

	void buildValueGreaterThanReason(const T & v, vec<Lit> & conflict, bool inclusive){
		if(inclusive){
			assert(getUnderApprox() >= v);
		}else{
			assert(getUnderApprox() > v);
		}
		lbool val = this->solver->value(condition);
		if(val==l_True){
			thn.buildValueGreaterThanReason(v,conflict,inclusive);
			conflict.push(~condition);
		}else if (val==l_False){
			els.buildValueGreaterThanReason(v,conflict,inclusive);
			conflict.push(condition);
		}else{
			assert(val==l_Undef);
			if(thn.getUnderApprox().contains(v)){
				thn.buildValueGreaterThanReason(v, conflict,inclusive);
			}else{
				assert(els.getUnderApprox().contains(v));
				els.buildValueGreaterThanReason(v, conflict,inclusive);
			}
		}
	}


};


#endif /* SYMBOLICSCALAR_H_ */
