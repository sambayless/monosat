//adapted from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/euler.js
//untested at this point.

#ifndef EULER_TREE_H
#define EULER_TREE_H


#include "TreapCustom.h"
#include "mtl/Vec.h"


using namespace Minisat;

class EulerTree{
	struct EulerHalfEdge;
/*	struct NodeData{
		bool vertex;
		int count()const{
			return vertex;
		}
	};*/
	TreapCustom t;
	struct EulerHalfEdge{
		//Value value;
		//NodeData d;
		int value;
		int s;
		int t;
		TreapCustom::Node * node;//node in the treap
		int opposite;
	};

	struct EulerVertex{
		//Value value;
		  TreapCustom::Node * node;
	};

	vec<EulerHalfEdge*> forward_edges;
	vec<EulerHalfEdge*> backward_edges;
public:
	//vec<EulerVertex*> nodes;
	void cut(int edgeID){

		EulerHalfEdge * edge = forward_edges[edgeID];
		EulerHalfEdge *opposite = edge->opposite;
		assert(opposite==backward_edges[edgeID]);
		TreapCustom::Node * a = edge->node;
		TreapCustom::Node * b= t.split(a);
		TreapCustom::Node * c = opposite->node;
		TreapCustom::Node * d = t.split(c);

		if(d && t.findRoot(a) != t.findRoot(d)){
			t.concat(a,d);
		}else if (b && t.findRoot(c)!= t.findRoot(b)){
			t.concat(c,b);
		}
		t.remove(edge->node);
		t.remove(opposite->node);
	}


	void setFlag (EulerVertex* node, bool f) {
	  t.setFlag(node->node,f);
	}

	bool path(EulerVertex*  from, EulerVertex*  to){
		return t.findRoot(from->node)==t.findRoot(from->node);
	}

	void makeRoot(EulerVertex*  node){
		TreapCustom::Node * a = node->node;
		TreapCustom::Node * b= t.split(a);
		if(b){
			t.concat(b,a);
		}
	}
	int link(EulerVertex*  node, EulerVertex*  otherNode, int edgeID,int value) {
		  //Move both vertices to root
		  makeRoot(node);
		  makeRoot(otherNode);

		  //Create half edges and link them to each other
		  forward_edges.growTo(edgeID+1);
		  backward_edges.growTo(edgeID+1);

		  if(forward_edges[edgeID]==NULL){
			  forward_edges[edgeID] = new EulerHalfEdge();
			  backward_edges[edgeID] = new EulerHalfEdge();
		  }

		  forward_edges[edgeID]->value=value;
		  forward_edges[edgeID]->s=node;
		  forward_edges[edgeID]->t = otherNode;
		  forward_edges[edgeID]->opposite = edgeID;

		  backward_edges[edgeID]->value=value;
		  backward_edges[edgeID]->s=otherNode;
		  backward_edges[edgeID]->t = node;
		  backward_edges[edgeID]->opposite = edgeID;

		  //var st = new EulerHalfEdge(value, this, other, null, null)
		  //var ts = new EulerHalfEdge(value, other, this, null, st)
		  //st.opposite = ts
		  forward_edges[edgeID]->node = t.insert( node->node,forward_edges[edgeID],1);
		  backward_edges[edgeID]->node = t.insert(otherNode->node,backward_edges[edgeID],1);


		  //Link tours together
		  t.concat(node,otherNode);

		  //Return half edge
		  return  edgeID;
		}


	EulerVertex * createVertex() {
		EulerVertex *  v = new EulerVertex( );
	  v->node = t.createNode(v);
	  return v;
	}
};

#endif



