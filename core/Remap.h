/*
 * Remap.h
 *
 *  Created on: Oct 14, 2015
 *      Author: sam
 */

#ifndef REMAP_H_
#define REMAP_H_
#include "core/SolverTypes.h"
#include "mtl/Vec.h"
#include "core/Solver.h"
#include "bv/BVTheorySolver.h"
//Classes for handling mapping to and from input formula numberings
namespace Monosat{


class DimacsMap{
	bool remap_vars;
	vec<Var> var_map;//Map of external variables to internal varialbes
	vec<Var> var_reverse_map; //Map of internal variables to internal variables
public:
	DimacsMap(bool remap_vars=true):remap_vars(remap_vars){
	}
	virtual ~DimacsMap() {
	}
	inline Var mapVar(Solver & S, Var var){
		if(!remap_vars){
			while (var >= S.nVars())
				S.newVar();
			return var;
		}
		if (inVarMap(var)){
			return getVarFromExternalVar(var);
		}else{
			Var v =S.newVar();
			addVarToMap(var,v);
			return v;
		}
	}

	inline bool inVarMap(Var externalVar){
		if(externalVar==var_Undef){
			return false;
		}
		return externalVar< var_map.size() && var_map[externalVar] !=var_Undef;
	}
	inline	Var unmap(Var internalVar){
			if(!remap_vars){
				return internalVar;
			}
			if (internalVar< var_map.size() && var_reverse_map[internalVar]!=var_Undef){
				return var_reverse_map[internalVar];
			}
			return var_Undef;
		}
	inline	Lit unmap(Lit internalLit){
		if(!remap_vars){
			return internalLit;
		}
		Var internalVar = var(internalLit);
		if (internalVar< var_map.size() && var_reverse_map[internalVar]!=var_Undef){
			return mkLit(var_reverse_map[internalVar],sign(internalLit));
		}
		return lit_Undef;
	}
	inline	void addVarToMap(Var v, Var map_to){
			if(!remap_vars){
				return;
			}
			if(inVarMap(v)){
				throw std::runtime_error("Variables can only be mapped at most once!");
			}
			var_map.growTo(v+1,var_Undef);
			var_map[v]=map_to;
			var_reverse_map.growTo(map_to+1,var_Undef);
			var_reverse_map[map_to]=v;
		}
	inline	bool hasMappedVar(Var internalVar){
			if(internalVar==var_Undef){
				return false;
			}
			return internalVar< var_reverse_map.size() && var_reverse_map[internalVar] !=var_Undef;
		}
	inline	Var getVarFromExternalVar(Var externalVar){
		if(!remap_vars){
					return externalVar;
				}
			if (inVarMap(externalVar)){
				return var_map[externalVar];
			}
			return var_Undef;
		}

	inline	int dimacs(Var v){
				return dimacs(mkLit(v));
			}
	inline	int dimacs(Lit l){
				return sign(l) ? -(unmap(var(l)) + 1) : (unmap(var(l)) + 1);
			}

};

class BVMap{
	bool remap_vars;
	vec<int> bv_map;//Map of external variables to internal varialbes
	vec<int> bv_reverse_map; //Map of internal variables to internal variables
public:
	BVMap(bool remap_vars=true):remap_vars(remap_vars){
	}
	virtual ~BVMap() {
	}
	inline int mapBV(Solver & S, int bv){
		BVTheorySolver<long>* theory = (BVTheorySolver<long>* )S.getBVTheory();
		if(!remap_vars)
			return bv;
		if(!inBVMap(bv)){
			parse_errorf("Undefined bitvector bv%d",bv);
		}
		return bv_map[bv];
	}

	inline bool inBVMap(int externalBV){
		return externalBV>=0 && externalBV< bv_map.size() && bv_map[externalBV] !=-1;
	}

	inline	void addBVToMap(int bvID, int map_to){
		if(!remap_vars)
			return;
		if(inBVMap(bvID)){
			parse_errorf("Bitvector %d mapped multiple times!",bvID);
		}
		bv_map.growTo(bvID+1,-1);
		bv_map[bvID]=map_to;
		bv_reverse_map.growTo(map_to+1,-1);
		bv_reverse_map[map_to]=bvID;
	}
	inline	bool hasMappedBV(int internalBV){
		if(!remap_vars)
			return false;
		if(internalBV==-1){
			return false;
		}
		return internalBV>=0 && internalBV< bv_reverse_map.size() && bv_reverse_map[internalBV] !=-1;
	}
	inline	int getBVFromExternalBV(int externalBV){
		if(!remap_vars)
			return externalBV;
		if (inBVMap(externalBV)){
			return bv_map[externalBV];
		}
		return -1;
	}
	inline int unmapBV(int bv){
		if(!remap_vars)
			return bv;
		if(bv>=0 && bv< bv_reverse_map.size() && bv_reverse_map[bv] !=-1){
			return bv_reverse_map[bv];
		}else{
			return bv;
		}
	}

};

};


#endif /* REMAP_H_ */
