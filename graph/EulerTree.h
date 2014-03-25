//adapted from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/euler.js
//untested at this point.

#ifndef EULER_TREE_H
#define EULER_TREE_H


#include "TreapCustom.h"
#include "mtl/Vec.h"


using namespace Minisat;

class EulerTree{
	struct EulerHalfEdge;
	struct EulerVertex;
	typedef TreapCustom<int> Treap;
/*	struct NodeData{
		bool vertex;
		int count()const{
			return vertex;
		}
	};*/
	Treap t;

	struct EulerHalfEdge{
		//Value value;
		//NodeData d;
		int value;
		int index;
		int from;
		int to;
		Treap::Node * node;//node in the treap

	};



	struct EulerVertex{
		//Value value;
		Treap::Node * first;
		Treap::Node * last;
		int index;
#ifndef NDEBUG

		EulerVertex * left;
		EulerVertex*right;
		EulerVertex*parent;
#endif
		EulerVertex(Treap::Node * first,Treap::Node * last):first(first),last(last){
			index=0;
#ifndef NDEBUG

		 left=NULL;
		 right=NULL;
		 parent=NULL;
#endif
		}
	};
	EulerVertex * root;
	vec<EulerVertex*> vertices;
	vec<EulerHalfEdge*> forward_edges;
	vec<EulerHalfEdge*> backward_edges;

	void dbg_tour(){
#ifndef NDEBUG
		//build the tour
		vec<int> tour;
		{
			vec<bool> seenLeft;
			vec<bool> seenRight;
			seenLeft.growTo(vertices.size());
			seenRight.growTo(vertices.size());

			EulerVertex * r = root;
			while(r){
				tour.push(r->index);
				//visit left
				if(!r->left){
					seenLeft[r->index]=true;
				}
				if(!r->right){
					seenRight[r->index]=true;
				}

				if(! seenLeft[r->index]){
					seenLeft[r->index]=true;
					r=r->left;
				}else if (!seenRight[r->index]){
					seenRight[r->index]=true;
					r=r->right;
				}else{
					r= r->parent;
				}
			}
		}

		vec<bool> visited;
		visited.growTo(vertices.size());

		//traverse the treap
		Treap::Node * r = root->first;
		assert(r==t.first(r));//the root of the euler tour tree must be the start of the euler tour, which must be the first element in the bst
		int i = 0;
		assert(i==r->value);
		//ok, visit in order
		 vec<Treap::Node*> stack;
		  while (stack.size()){
		    if (r){
		    	stack.push(r);
		    	r = r->left;
		    }else{
		      r = stack.last();stack.pop();
		      visited[r->value]=true;
		      int p = tour[i++];
		      assert(r->value==p);
		      r = r->right;
		    }
		  }

		  for(bool v : visited){
			  assert(v);
		  }
#endif
	}

public:

	//vec<EulerVertex*> nodes;
	void cut(int edgeID){

		EulerHalfEdge * edge = forward_edges[edgeID];
		cut (edge);
	}

	void cut(EulerHalfEdge * edge){
		EulerHalfEdge *opposite = backward_edges[edge->opposite];

		Treap::Node * a = edge->node;
		Treap::Node * b= t.split(a);
		Treap::Node * c = opposite->node;
		Treap::Node * d = t.split(c);

		if(d && t.findRoot(a) != t.findRoot(d)){
			t.concat(a,d);
		}else if (b && t.findRoot(c)!= t.findRoot(b)){
			t.concat(c,b);
		}
		t.remove(edge->node);
		t.remove(opposite->node);
	}

	bool connected(EulerVertex*  from, EulerVertex*  to){
		return t.findRoot(from->node)==t.findRoot(from->node);
	}

	void makeRoot(EulerVertex*  node){
		Treap::Node * a = node->node;
		Treap::Node * b= t.split(a);
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
			  forward_edges[edgeID]->index=edgeID;
			  backward_edges[edgeID]->index=edgeID;

			  forward_edges[edgeID]->value=value;
			  forward_edges[edgeID]->from=node->index;
			  forward_edges[edgeID]->to = otherNode->index;


			  backward_edges[edgeID]->value=value;
			  backward_edges[edgeID]->from=otherNode->index;
			  backward_edges[edgeID]->to = node->index;
		  }


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
	  Treap::Node * n1= t.createNode(vertices.size());
	  Treap::Node * n2= t.createNode(vertices.size());
	  vertices.push(new EulerVertex{n1,n2});
	  vertices.last()->index = vertices.size()-1;
	  return vertices.last();
	}
};

#endif



