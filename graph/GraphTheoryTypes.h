/*
 * GraphTheoryTypes.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef GRAPHTHEORYTYPES_H_
#define GRAPHTHEORYTYPES_H_

#include "core/SolverTypes.h"
namespace Minisat{

struct ReachabilityConstraint{
	int from;
	int to;
	int distance;
	Var reach_var;
};

struct DefaultEdgeStatus{

	vec<char> status;

	bool operator [] (int index) const {
		return status[index];
	}
	void setStatus(int index, bool value){
		status[index]=value;
	}

	int size(){
		return status.size();
	}
	void growTo(int size) {
		status.growTo(size);
	}

};

typedef DefaultEdgeStatus PositiveEdgeStatus;
typedef DefaultEdgeStatus NegativeEdgeStatus;
typedef DefaultEdgeStatus CutEdgeStatus;

	struct Assignment{
		bool isEdge:1;
		bool assign:1;
		int from:30;
		int to;
		Var var;
		Assignment(bool _isEdge,bool _assign,int _from, int _to, Var v):isEdge(_isEdge),assign(_assign),from(_from),to(_to),var(v){

		}
		Assignment():isEdge(false),assign(false),from(0),to(0),var(var_Undef){

		}
	};
	struct Edge{
		Var v;
		int from;
		int to;
		int edgeID;
		int weight;
	};

	// Returns a random float 0 <= x < 1. Seed must never be 0.
	   static inline double drand(double& seed) {
	   	assert(seed!=0);
	       seed *= 1389796;
	       int q = (int)(seed / 2147483647);
	       seed -= (double)q * 2147483647;
	       return seed / 2147483647; }

	   // Returns a random integer 0 <= x < size. Seed must never be 0.
	   static inline int irand(double& seed, int size) {
	       return (int)(drand(seed) * size); }
};

#endif /* GRAPHTHEORYTYPES_H_ */
