/*
 * GraphTheoryTypes.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef GRAPHTHEORYTYPES_H_
#define GRAPHTHEORYTYPES_H_

#include "core/SolverTypes.h"
#include "mtl/Rnd.h"
namespace Minisat{

struct ReachabilityConstraint{
	int from;
	int to;
	int distance;
	Var reach_var;
};

struct ConnectivityConstraint{
	int from;
	int to;
	int distance;
	Var connect_var;
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
//extern DefaultEdgeStatus defaultStatus;
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
		Var outerVar;
		int from;
		int to;
		int edgeID;
		int weight;
		Edge(Var v, Var outerVar,int from,int to, int edgeID, int weight=0):v(v),outerVar(outerVar),from(from),to(to),edgeID(edgeID),weight(weight){

		}
		Edge():v(var_Undef),outerVar(var_Undef),from(-1),to(-1),edgeID(-1),weight(0){

		}
	};

};

#endif /* GRAPHTHEORYTYPES_H_ */
