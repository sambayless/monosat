//adapted from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/euler.js
//untested at this point.

#ifndef EULER_TREE_H
#define EULER_TREE_H


#include "TreapCustom.h"
#include "mtl/Vec.h"


using namespace Minisat;

class EulerTree{
public:
	struct EulerHalfEdge;
	struct EulerVertex;
private:
	typedef TreapCustom<EulerHalfEdge*> Treap;

	//make this non-static later
	static Treap t;
	static EulerVertex * root;
	static vec<EulerVertex*> vertices;
	static vec<EulerHalfEdge*> forward_edges;
	static vec<EulerHalfEdge*> backward_edges;
public:
	struct EulerHalfEdge{
		//Value value;
		//NodeData d;
		int value;
		int index;
		EulerVertex * from;
		EulerVertex * to;
		Treap::Node * node;//node in the treap

	};



	struct EulerVertex{
		//Value value;

		//Treap::Node * first;
		//Treap::Node * last;
		EulerHalfEdge * left_out;//the half edge leading out of the vertex. This half edge has the first occurrence of the euler vertex in the tour.
		EulerHalfEdge * right_in;//the half edge returning to the vertex.  This half edge has the last occurrence of the euler vertex in the tour.
		int index;
#ifndef NDEBUG

		EulerVertex * dbg_left;
		EulerVertex*dbg_right;
		EulerVertex*dbg_parent;
#endif
		EulerVertex():left_out(nullptr),right_in(nullptr){
			index=0;
#ifndef NDEBUG

			dbg_left=NULL;
			dbg_right=NULL;
			dbg_parent=NULL;
#endif
		}

		Treap::Node * first(){
			if(left_out){
				return left_out->node;
			}else{
				return nullptr;
			}
		}

		Treap::Node * last(){
			if(right_in){
				return right_in->node;
			}else{
				return nullptr;
			}
		}
		 void dbg_remove(){
				if(!dbg_left && !dbg_right){
					if(dbg_parent){
						if(dbg_parent->dbg_left==this){
							dbg_parent->dbg_left=nullptr;
						}else {
							assert(dbg_parent->dbg_right==this);
							dbg_parent->dbg_right=nullptr;
						}

					}else{
						//do nothing
					}
				}else if (dbg_left && ! dbg_right){
					if(dbg_parent){
						if(dbg_parent->dbg_left==this){
							dbg_parent->dbg_left=dbg_left;
						}else {
							assert(dbg_parent->dbg_right==this);
							dbg_parent->dbg_right=dbg_left;
						}
					}else{
						dbg_left->dbg_parent=nullptr;

					}
				}else if (dbg_right && ! dbg_left){
					if(dbg_parent){
						if(dbg_parent->dbg_left==this){
							dbg_parent->dbg_left=dbg_right;
						}else {
							assert(dbg_parent->dbg_right==this);
							dbg_parent->dbg_right=dbg_right;
						}

					}else{
						dbg_right->dbg_parent=nullptr;

					}

				}else{
					//not handled yet...
					assert(false);
					/*EulerVertex * next = dbg_right->findMin();
					this->key=next->key;
					assert(!(next->dbg_right || next->dbg_left));
					next->remove();
					checkTree();*/
				}
				dbg_parent=nullptr;
			}

		void dbg_insert(EulerVertex* node) {
#ifndef NDEBUG
		     if (!dbg_left) {
		    	 dbg_left = node;
		    	 dbg_left->dbg_parent = this;
		      }else if (!dbg_right) {
		    	  dbg_right = node;
				  dbg_right->dbg_parent = this;
			  } else {
				  dbg_right->dbg_insert(node);
			  }
#endif
		  }
		void dbg_tour(){
		#ifndef NDEBUG
				//build the tour
				vec<int> tour;
				{
					vec<bool> seenLeft;
					vec<bool> seenRight;
					seenLeft.growTo(vertices.size());
					seenRight.growTo(vertices.size());

					EulerVertex * r = this;
					while(r){
						tour.push(r->index);
						//visit dbg_left
						if(!r->dbg_left){
							seenLeft[r->index]=true;
						}
						if(!r->dbg_right){
							seenRight[r->index]=true;
						}

						if(! seenLeft[r->index]){
							seenLeft[r->index]=true;
							r=r->dbg_left;
						}else if (!seenRight[r->index]){
							seenRight[r->index]=true;
							r=r->dbg_right;
						}else{
							r= r->dbg_parent;
						}
					}
				}

				vec<bool> visited;
				visited.growTo(vertices.size());

				//traverse the treap
				Treap::Node * r = this->left_out->node;
				assert(r==t.first(r));//the root of the euler tour tree must be the start of the euler tour, which must be the first element in the bst
				int i = 0;
				assert(i==r->value->from->index);
				//ok, visit in order
				 vec<Treap::Node*> stack;
				  while (stack.size()){
				    if (r){
				    	stack.push(r);
				    	r = r->left;
				    }else{
				      r = stack.last();stack.pop();
				      int nodeindex =r->value->from->index;
				      visited[nodeindex]=true;
				      int p = tour[i++];
				      assert(nodeindex==p);
				      r = r->right;
				    }
				  }

				  for(bool v : visited){
					  assert(v);
				  }
		#endif
			}
	};




public:

	//vec<EulerVertex*> nodes;
/*	void cut(int edgeID){

		EulerHalfEdge * edge = forward_edges[edgeID];
		cut (edge);
	}*/
	void cut(EulerVertex * v){
		if(root==v){
			assert(false);
			return;//do nothing...
		}
		root->dbg_tour();
		v->dbg_tour();
		v->dbg_remove();

		Treap::Node * prev = v->first()->prev;
		Treap::Node * next = v->last()->next;
		assert(v->first()->value->from==v);
		assert(v->last()->value->to==v);
		assert(prev->value->to==v);
		assert(next->value->from==v);
#ifndef NDEBUG
		{
			EulerVertex* parent = v->dbg_parent;
			assert(prev->value->from==parent);
			assert(next->value->to==parent);
		}
#endif

		t.split(v->first());
		t.split(v->last());

		v->dbg_tour();

		//remove the two half edges connected the vertex to the tree

		Treap::Node * prev_kept=prev->prev;
		Treap::Node * next_kept=next->next;

		t.remove(next);
		t.remove(prev);

		next->value->node=nullptr;
		prev->value->node=nullptr;

		t.concat(prev_kept,next_kept);
		root->dbg_tour();

	}

	bool connected(EulerVertex*  from, EulerVertex*  to){
		return t.findRoot(from->first())==t.findRoot(from->first());
	}

	void makeRoot(EulerVertex*  node){
		Treap::Node * a = node->first();
		Treap::Node * b= t.split(a);
		if(b){
			t.concat(b,a);
		}
	}

	//Make othernode a child of node.
	int link(EulerVertex*  node, EulerVertex*  otherNode, int edgeID) {
		  //Move both vertices to root
		  //Is this really required?
		  makeRoot(node);
		  makeRoot(otherNode);

		  node->dbg_tour();
		  otherNode->dbg_tour();
		  root->dbg_tour();

		  node->dbg_insert(otherNode);
		  //Create half edges and link them to each other
		  forward_edges.growTo(edgeID+1);
		  backward_edges.growTo(edgeID+1);

		  if(forward_edges[edgeID]==NULL){
			  forward_edges[edgeID] = new EulerHalfEdge();
			  backward_edges[edgeID] = new EulerHalfEdge();
			  forward_edges[edgeID]->index=edgeID;
			  backward_edges[edgeID]->index=edgeID;


			  forward_edges[edgeID]->from=node;
			  forward_edges[edgeID]->to = otherNode;



			  backward_edges[edgeID]->from=otherNode;
			  backward_edges[edgeID]->to = node;
		  }else{
			  assert(forward_edges[edgeID]->index==edgeID);
			  assert(backward_edges[edgeID]->index==edgeID);


			  assert(forward_edges[edgeID]->from==node);
			  assert(forward_edges[edgeID]->to == otherNode);


			  assert(backward_edges[edgeID]->from==otherNode);
			  assert(backward_edges[edgeID]->to == node);
		  }

		  //forward_edges[edgeID]->value==value;
		  //backward_edges[edgeID]->value==value;

		  forward_edges[edgeID]->node = t.insert( node->last(),forward_edges[edgeID],1);
		  backward_edges[edgeID]->node = t.insert(otherNode->last(),backward_edges[edgeID],1);

		  //Link tours together
		  t.concat(forward_edges[edgeID]->node ,otherNode->first());
		  t.concat(otherNode->last() ,backward_edges[edgeID]->node);

		  node->dbg_tour();
		  otherNode->dbg_tour();
		  root->dbg_tour();
		  //Return half edge
		  return  edgeID;
		}

	void link(int u, int v, int edgeID) {
		link(vertices[u],vertices[v],edgeID);
	}
	void cut(int u) {
		cut(vertices[u]);
	}
	bool connected(int u, int v){
		return connected(vertices[u],vertices[v]);
	}
	EulerVertex* getVertex(int n){
		return vertices[n];
	}

	EulerVertex * createVertex() {
	  vertices.push(new EulerVertex());
	  vertices.last()->index = vertices.size()-1;
	  return vertices.last();
	}
};

#endif



