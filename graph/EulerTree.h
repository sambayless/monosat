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
	int nComponents;
	//make this non-static later
	static Treap t;
	//static EulerVertex * root;
	static vec<EulerVertex*> vertices;
	static vec<EulerHalfEdge*> forward_edges;
	static vec<EulerHalfEdge*> backward_edges;
public:
	struct EulerHalfEdge{
		//Value value;
		//NodeData d;
		int value;
		int index;
		bool isforward;//forward edges go from a parent to a child

		EulerVertex * from;
		EulerVertex * to;
		Treap::Node * node;//node in the treap

#ifndef NDEBUG
		int rank;//this is the position of this edge in the euler tour. This information is only maintained implicitly as the edges position in the underlying binary search tree.
#endif
	};



	struct EulerVertex{
		//Value value;

		//Treap::Node * first;
		//Treap::Node * last;

		EulerHalfEdge * left_out;//the half edge leading out of the vertex. This half edge has the first occurrence of the euler vertex in the tour.
		EulerHalfEdge * right_in;//the half edge returning to the vertex.  This half edge has the last occurrence of the euler vertex in the tour.
		int index;
		int subtree_size;
		bool subtree_has_incident_edges;
		bool has_incident_edges;
#ifndef NDEBUG
		//Note: in an euler-tour tree representation, these values are not explicitly maintained
		EulerVertex * dbg_left;
		EulerVertex*dbg_right;
		EulerVertex*dbg_parent;
#endif
		EulerVertex():left_out(nullptr),right_in(nullptr){
			index=0;
			subtree_size=1;
			has_incident_edges=false;
			subtree_has_incident_edges=false;
#ifndef NDEBUG

			dbg_left=NULL;
			dbg_right=NULL;
			dbg_parent=NULL;
#endif
		}

		int getSize(){
			if(!left_out){
				assert(!right_in);
				return 1;//size is just a single node
			}
			assert(right_in);
			assert(left_out->rank<right_in->rank);
			assert(((left_out->rank - right_in->rank) %2) ==0);
			int subtree_size = (left_out->rank - right_in->rank)/2;
			assert(subtree_size==dbg_getSize());
			return subtree_size;
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
		//Get the next node in the tour.
		EulerVertex * getNext(){
			if (left_out)
				return left_out->to;
			else if(last())
				return last()->value->to;
			else
				return nullptr;
		}
		int dbg_getSize(){
#ifndef NDEBUG
			int s = 1;
			if(dbg_left){
				s+=dbg_left->dbg_getSize();
			}
			if(dbg_right){
				s+=dbg_right->dbg_getSize();
			}
			return s;
#endif
			return 0;
		}

		 void dbg_remove(){
#ifndef NDEBUG
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
#endif
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
				      assert(r->value->rank==i);
				      visited[nodeindex]=true;
				      int p = tour[i++];
				      assert(nodeindex==p);
				      r = r->right;
				    }
				  }

				  for(bool v : visited){
					  assert(v);
				  }

				  EulerVertex * v = this;
				  for(int i = 0;i<tour.size();i++){
					  int p = tour[i];
					  assert(v->index==p);
					  v =v->getNext();
				  }
				  assert(v->getNext() ==this);

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
		if(getParent(v)==nullptr){
			return;//do nothing
		}else{
			nComponents++;
		}

		v->dbg_tour();
#ifndef NDEBUG
		EulerVertex * dbg_parent = v->dbg_parent;
#endif
		v->dbg_remove();

		 int removedNodes = v->subtree_size;
		  //update subtree sizes going up the tree
		  EulerVertex * p = getParent(v);
		  while(p){
			  p->subtree_size-= removedNodes;
			  assert(p->subtree_size==p->dbg_getSize());
			  p = getParent(p);
		  }


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

		v->dbg_tour();
#ifndef NDEBUG
		dbg_parent->dbg_tour();
#endif
	}

	//Get the parent of the vertex in the euler tour representation. This is O(1).
	EulerVertex * getParent(EulerVertex * v){
		EulerVertex * parent = nullptr;
		if(!v->first()){
			//this is a root
			assert(!v->dbg_parent);
			return nullptr;
		}
		assert(v->first()->value->isforward);
		parent= v->first()->value->from;
		assert(parent==v->dbg_parent);
		return parent;
	}
	EulerVertex * getLeft(EulerVertex * v){
			if(!v->left_out){
				assert(!v->dbg_left);
				return nullptr;
			}

			EulerVertex * left =  v->left_out->to;
			assert(left==v->dbg_left);
			return left;
		}

	EulerVertex * getRight(EulerVertex * v){
			if(!v->right_in){
				assert(!v->dbg_right);
				return nullptr;
			}

			EulerVertex * right =  v->right_in->from;
			assert(right==v->dbg_right);
			return right;
		}

	//Get the next node in the tour.
	EulerVertex * getNext(EulerVertex * from){
		return from->getNext();
	}

	EulerVertex * findRoot(EulerVertex * v){
		EulerVertex * parent = getParent(v);
		while(parent){
			v = parent;
			parent = getParent(v);
		}
		return v;
	}

	//Return the size of the subtree rooted at v, including v
	int getSubtreeSize(EulerVertex * v){
		assert(v->subtree_size==v->dbg_getSize());
		return v->subtree_size;
	}

	bool hasIncidentEdges(EulerVertex* v){
		return v->has_incident_edges;
	}
	bool subtreeHasIncidentEdges(EulerVertex * v){
		return v->subtree_has_incident_edges;
	}

	void setHasIncidentEdges(EulerVertex * v, bool hasIncident){
		if(!v->has_incident_edges && hasIncident){
			v->has_incident_edges=hasIncident;

			//ok, now follow the tour up the tree and mark that the subtrees have incident edges
			EulerVertex * parent = getParent(v);
			while(parent && !parent->subtree_has_incident_edges){
				parent->subtree_has_incident_edges=true;
				parent = getParent(parent);
			}

#ifndef NDEBUG
			EulerVertex * p = v->dbg_parent;
			while(p){
				assert(p->subtree_has_incident_edges);
				p=p->dbg_parent;
			}
#endif
		}else if(v->has_incident_edges &&! hasIncident){
			v->has_incident_edges=false;//don't update parent information yet.
		}
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
		  assert(!connected(node,otherNode));
		  nComponents--;
		  //Move both vertices to root
		  //Is this really required?
		  makeRoot(node);
		  makeRoot(otherNode);

		  node->dbg_tour();
		  otherNode->dbg_tour();


		  node->dbg_insert(otherNode);

		  assert(otherNode->dbg_parent == node);

		  //Create half edges and link them to each other
		  forward_edges.growTo(edgeID+1);
		  backward_edges.growTo(edgeID+1);

		  if(forward_edges[edgeID]==NULL){
			  forward_edges[edgeID] = new EulerHalfEdge();
			  backward_edges[edgeID] = new EulerHalfEdge();
			  forward_edges[edgeID]->index=edgeID;
			  backward_edges[edgeID]->index=edgeID;
			  forward_edges[edgeID]->isforward=true;

			  forward_edges[edgeID]->from=node;
			  forward_edges[edgeID]->to = otherNode;


			  backward_edges[edgeID]->isforward=false;
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

		  assert(forward_edges[edgeID]->to->dbg_parent==forward_edges[edgeID]->from);
		  assert(backward_edges[edgeID]->from->dbg_parent==backward_edges[edgeID]->to);

		  //forward_edges[edgeID]->value==value;
		  //backward_edges[edgeID]->value==value;

		  forward_edges[edgeID]->node = t.insert( node->last(),forward_edges[edgeID],1);
		  backward_edges[edgeID]->node = t.insert(otherNode->last(),backward_edges[edgeID],1);

		  //Link tours together
		  t.concat(forward_edges[edgeID]->node ,otherNode->first());
		  t.concat(otherNode->last() ,backward_edges[edgeID]->node);
		  assert(otherNode->subtree_size==otherNode->dbg_getSize());
		  int addedNodes = otherNode->subtree_size;
		  //update subtree sizes going up the tree
		  EulerVertex * p = getParent(node);
		  while(p){
			  p->subtree_size+= addedNodes;
			  assert(p->subtree_size==p->dbg_getSize());
			  p = getParent(p);
		  }
		  node->dbg_tour();
		  otherNode->dbg_tour();

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

	//Get the parent of the vertex in the euler tour representation. This is O(1).
	int getParent(int v){
		EulerVertex * parent = getParent(vertices[v]);
		if(!parent){
			return -1;
		}else{
			return parent->index;
		}
	}



	int findRoot(int v){
		EulerVertex * parent = getParent(vertices[v]);
		while(parent){
			v = parent->index;
			parent = getParent(vertices[v]);
		}
		return v;
	}

	int numComponents(){
		return nComponents;
	}

	//Return the size of the subtree rooted at v, including v
	int getSubtreeSize(int v){
		return getSubtreeSize(vertices[v]);
	}

	bool hasIncidentEdges(int v){
		return hasIncidentEdges(vertices[v]);
	}
	bool subtreeHasIncidentEdges(int v){
		return subtreeHasIncidentEdges(vertices[v]);
	}
	void setHasIncidentEdges(int v, bool hasIncident){
		setHasIncidentEdges(vertices[v],hasIncident);
	}

	void cutEdge(int edgeID){
		assert(forward_edges[edgeID]->from);
		assert(forward_edges[edgeID]->to->dbg_parent==forward_edges[edgeID]->from);
		cut(forward_edges[edgeID]->to);
		assert(!connected(forward_edges[edgeID]->to,forward_edges[edgeID]->from));
	}

	EulerVertex* getVertex(int n){
		return vertices[n];
	}

	EulerVertex * createVertex() {
	  nComponents++;
	  vertices.push(new EulerVertex());
	  vertices.last()->index = vertices.size()-1;
	  return vertices.last();
	}
};

#endif



