//adapted from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/euler.js
//untested at this point.

#ifndef EULER_TREE_H
#define EULER_TREE_H


#include "TreapCustom.h"
#include "mtl/Vec.h"
#include <cstdio>
#include "SplayTree.h"
#include "SearchTree.h"
using namespace Minisat;

class EulerTree{
public:
	struct EulerHalfEdge;
	struct EulerVertex;
private:


	typedef SplayTree<EulerHalfEdge*> Tree;
	int nComponents;

	Tree t;
	//static EulerVertex * root;
	vec<EulerVertex*> vertices;
	vec<EulerHalfEdge*> forward_edges;
	vec<EulerHalfEdge*> backward_edges;
public:
	struct EulerHalfEdge{
		//Value value;
		//NodeData d;
		int value;
		int index;
		bool isforward;//forward edges go from a parent to a child

		EulerVertex * from;
		EulerVertex * to;
		Tree::Node * node;//node in the treap

#ifndef NDEBUG
		int rank;//this is the position of this edge in the euler tour. This information is only maintained implicitly as the edges position in the underlying binary search tree.
#endif
	};



	struct EulerVertex{
		//Value value;

		//Treap::Node * first;
		//Treap::Node * last;
		EulerHalfEdge * left_out;//the half edge leading into the vertex, from its parent. The 'to' field in this halfedge the first occurrence of the euler vertex in the tour.
								//IFF the vertex is root, then this is instead the half edge leading to its first child. In that case, the 'from' field is the first occurence of the vertex in the tour
		EulerHalfEdge * right_in;//the half edge returning to the vertex.  This half edge has the last occurrence of the euler vertex in the tour.
		int index;
		int subtree_size;
		bool subtree_has_incident_edges;
		bool has_incident_edges;
#ifndef NDEBUG
		EulerTree * owner;
		//Note: in an euler-tour tree representation, these values are not explicitly maintained
		EulerVertex * dbg_parent;
		vec<EulerVertex*>dbg_children;

#endif
		EulerVertex():left_out(nullptr),right_in(nullptr){

			index=0;
			subtree_size=1;
			has_incident_edges=false;
			subtree_has_incident_edges=false;
#ifndef NDEBUG
			owner=nullptr;
			dbg_parent=nullptr;
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

		bool isSingleton(){
			assert(!left_out == !right_in);
			return !left_out;
		}

		bool isRoot(){
			EulerVertex * parent = nullptr;
			if(isSingleton()){

				return true;
			}else if (isLeaf()){

				assert(first()->value->from==dbg_parent);
				return false;
			}else{
				if(!first()->prev()){
					assert(owner->t.findRoot(first()) ==owner->t.findRoot(last()));
					//assert(owner->t.first(owner->t.findRoot(first())) ==first());
					//assert(owner->t.last(owner->t.findRoot(last())) ==last());
				}
				return !( first()->prev());
			}


		}

/*		bool isRoot(){
			//Only a node that is both a leaf and root doesn't have half edges stored.
			if(isSingleton()){
				return true;
			}

			assert((left_out->from == this)==(right_in->to == this));
			return left_out->from == this;
		}*/

		bool isLeaf(){
			if(isSingleton()){
				return true;
			}
			//leaf vertices are special. Instead of storing the edges leading out to their first and last child, they store the (only two) edges leading into them in the euler tour.
			//For non-singleton leaf vertices only, left_out is instead left_in, and right_in is instead right_out;
			assert(!left_out || (left_out->from==this || left_out->to==this));
			assert(!right_in || (right_in->from==this || right_in->to==this));
			assert((left_out->to == this)==(right_in->from == this));
			return left_out->to == this;
		}

		Tree::Node * first(){
			if(left_out){
				return left_out->node;
			}else{
				return nullptr;
			}
		}
		Tree::Node * last(){

			if(right_in){
				return right_in->node;
			}else{
				return nullptr;
			}
		}
		void setFirst(EulerHalfEdge * n){
			assert(!n || (n->from==this || n->to==this));
			left_out=n;
		}
		void setLast(EulerHalfEdge * n){
			assert(!n || (n->from==this || n->to==this));
			right_in=n;
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

		void tour(vec<int> & tour_list){
			tour_list.clear();

			if(isSingleton() || isLeaf()){
				tour_list.push(index);
				return;
			}
			Tree::Node* n= first();
			Tree::Node* l= last();

			assert(owner->t.findRoot(first())==owner->t.findRoot(l));
			while(n!=last()){
				tour_list.push(n->value->from->index);
				assert( n->next());
				assert(n->value->to == n->next()->value->from);
				//printf("(%d,%d)\n",n->value->from->index,n->value->to->index);
				assert(n->next());
				n=n->next();
				assert(n);
			}
			tour_list.push(n->value->from->index);
			tour_list.push(n->value->to->index);
		}

		int dbg_getSize(){
#ifndef NDEBUG
			int s = 1;
			for(EulerVertex* c:dbg_children){
				s+=c->dbg_getSize();
			}
			assert(s==subtree_size);
			return s;
#endif
			return 0;
		}



		void dbg_make_root(){
#ifndef NDEBUG
			if(dbg_parent){
				EulerVertex * p = dbg_parent;
				EulerVertex *r= this;
				while(r->dbg_parent){
					r= r->dbg_parent;
				}
				p->dbg_children.remove(this);

				dbg_children.push(r);
				r->dbg_parent=this;
				dbg_parent=nullptr;
			}
#endif
		}
		 void dbg_remove(){
#ifndef NDEBUG



				if(dbg_parent){
					dbg_parent->dbg_children.remove(this);
					dbg_parent=nullptr;
				}else{
					//do nothing
				}
				/*else if (dbg_left && ! dbg_right){
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
*/
			/*	}else{
					//not handled yet...
					assert(false);
					EulerVertex * next = dbg_right->findMin();
					this->key=next->key;
					assert(!(next->dbg_right || next->dbg_left));
					next->remove();
					checkTree();
				}*/

#endif
			}

		void dbg_insert(EulerVertex* node) {
#ifndef NDEBUG
			assert(!node->dbg_parent);
			dbg_children.push(node);
			node->dbg_parent=this;
		/*     if (!dbg_left) {
		    	 dbg_left = node;
		    	 dbg_left->dbg_parent = this;
		      }else if (!dbg_right) {
		    	  dbg_right = node;
				  dbg_right->dbg_parent = this;
			  } else {
				  dbg_right->dbg_insert(node);
			  }*/
#endif
		  }

		void dbg_build_tour(EulerVertex * r, vec<int> & tour){
			tour.push(r->index);
			for(EulerVertex * c:r->dbg_children){
				dbg_build_tour(c,tour);
				tour.push(r->index);
			}
		}

		void dbg_tour(){
		#ifndef NDEBUG
			if(!left_out){
				assert(!right_in);
				return;
			}
				//build the tour
				vec<int> dbg_tour;
				dbg_build_tour(this,dbg_tour);
				vec<int> real_tour;
				tour(real_tour);
				assert(dbg_tour.size()==real_tour.size());
				for(int i = 0;i<real_tour.size();i++){
					assert(dbg_tour[i]==real_tour[i]);
				}

				for(EulerVertex * c: dbg_children)
					c->dbg_tour();

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


		Tree::Node * prev = v->first()->prev();
		Tree::Node * next = v->last()->next();
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

		Tree::Node * T2 = t.splitBefore(v->first());
		Tree::Node * T3 = t.splitAfter(v->last());
		//possibly have to delete a redundant edge here in some cases, watch out!
		t.concat(T2,T3);


		v->dbg_tour();

		//remove the two half edges connected the vertex to the tree

		Tree::Node * prev_kept=prev->prev();
		Tree::Node * next_kept=next->next();

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

/*
	void cut(int edgeID){

			EulerHalfEdge * forwardEdge = forward_edges[edgeID];
			EulerHalfEdge * backwardEdge = backward_edges[edgeID];

			EulerVertex * v = forwardEdge->to;

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


			Tree::Node * prev = v->first()->prev();
			Tree::Node * next = v->last()->next();
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

			Tree::Node * T2 = t.split(forwardEdge->node);
			Tree::Node * T3 = t.split(backwardEdge->node);



			v->dbg_tour();

			//remove the two half edges connected the vertex to the tree

			Tree::Node * prev_kept=prev->prev();
			Tree::Node * next_kept=next->next();

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
*/

	//Get the parent of the vertex in the euler tour representation. This is O(1).
	EulerVertex * getParent(EulerVertex * v){
		EulerVertex * parent = nullptr;
		if(v->isSingleton()){
			assert(!v->dbg_parent);
			return nullptr;
		}else if (v->isLeaf()){
			assert(v->first()->value->from==v->dbg_parent);
			assert( v->first()->value->from==v->dbg_parent);
			return v->first()->value->from;
		}else{
			if(v->first()->prev() ==nullptr){
				assert(!v->dbg_parent);
				return nullptr;
			}
			assert(v->first()->prev()->value->to==v);
			assert(v->first()->prev()->value->from==v->dbg_parent);

			assert( v->first()->value->from==v);
			return v->first()->prev()->value->from;
		}


/*		if(!v->first() ){
			//this is wrong... is it? even the root should still have first and last pointers... shouldn't it?
			//this is a root
			assert(!v->dbg_parent);
			return nullptr;
		}
		assert(v->first()->value->isforward);
		parent= v->first()->value->from;
		assert(parent==v->dbg_parent);*/
		return parent;
	}
/*	EulerVertex * getLeft(EulerVertex * v){
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
		}*/

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
		if(from==to){
			return true;
		}
		if(from->first()==nullptr || to->first()==nullptr){
			return false;
		}
		return t.findRoot(from->first())==t.findRoot(to->first());
	}

	void makeRoot(EulerVertex*  node){
		node->dbg_tour();
	/*	if(!node->isRoot()){
			assert(node->right_in);
			Tree::Node * a = node->first();
			Tree::Node * b= t.splitBefore(a);
			assert(a==t.findRoot(a));
			assert(b==t.findRoot(b));

			assert(t.findMin(a)==a);//a is the first element in this tree now
			assert(t.findMax(b)==b);
			//now reverse the order of these
			if(b){
				t.concat(a,b);
			}
			assert(t.findMin(t.findRoot(a))==node->first());

			assert(t.findMax(t.findRoot(a))==node->last());
		}*/
		if(!node->isRoot()){
			//adjust subtree sizes
			EulerVertex * p = getParent(node);
			EulerVertex * orginalParent = p;
			int rootSize = 0;
			while(p){
				p->subtree_size -=node->subtree_size;
				rootSize = p->subtree_size;
				p=getParent(p);
			}

			if(node->isLeaf()){
				Tree::Node * left= node->first();

				assert(left->next()==node->last());

				Tree::Node * right = t.splitAfter(left);
				assert(right==node->last());

				//ok, now join these in reverse order
				t.concat(right,left);

				node->setLast(left->value);
				node->setFirst(right->value);
				assert( node->subtree_size==1);

			}else{
				Tree::Node * a = node->first();
				Tree::Node * r1 = t.findMin(t.findRoot(a));
				if(r1==a){
					bool b = node->isRoot();
				}
				assert(r1!=a);
				Tree::Node * prev= t.splitBefore(a);
				assert(prev);
				assert(t.findRoot(a)==a);
				assert(t.findRoot(prev)==prev);
				Tree::Node * rr = t.findMax(prev);
				if(prev!=r1){
					assert(t.findRoot(r1)==a);
					Tree::Node * Tl = t.splitBefore(r1);
					t.concat(prev,Tl);
				}
				assert(t.findMin(t.findRoot(prev))==prev);
				Tree::Node * e = t.findMax(t.findRoot(prev));
				t.concat(node->last(),prev);//is this correct?
				node->setLast(e->value);
			}
			node->dbg_make_root();
			//also, the previous parent of this node _may_ now be a leaf, in which case it needs to be adjusted...
			if(orginalParent->first()->next()==orginalParent->last()){
				assert(orginalParent->dbg_children.size()==0);
				orginalParent->setFirst(orginalParent->first()->prev()->value);
				orginalParent->setLast(orginalParent->last()->next()->value);
			}else{
				assert(orginalParent->dbg_children.size());
			}


			node->subtree_size+=rootSize;
			assert(node->subtree_size==node->dbg_getSize());
			assert(t.findMin(t.findRoot(node->first()))==node->first());
			assert(t.findMax(t.findRoot(node->first()))==node->last());
			assert(t.findRoot( node->first()) ==t.findRoot(node->last()));
			assert(node->isRoot());

		}

		node->dbg_tour();
	}

	//Make othernode a child of node.
	int link(EulerVertex*  node, EulerVertex*  otherNode, int edgeID) {
		  assert (node !=otherNode);
		  assert(!connected(node,otherNode));
		  nComponents--;
		  //Move both vertices to root
		  //Is this really required?



		  //makeRoot(node);
		  //makeRoot(otherNode);

		  node->dbg_tour();
		  otherNode->dbg_tour();

		  dbg_printTour(node);
		  dbg_printTour(otherNode);
		  assert(findRoot(node)->isRoot());
		  assert(findRoot(otherNode)->isRoot());
		  assert(otherNode->isRoot());

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

			  forward_edges[edgeID]->node = t.createNode(forward_edges[edgeID]);
			  backward_edges[edgeID]->node = t.createNode(backward_edges[edgeID]);
		  }else{
			  assert(forward_edges[edgeID]->index==edgeID);
			  assert(backward_edges[edgeID]->index==edgeID);


			  assert(forward_edges[edgeID]->from==node);
			  assert(forward_edges[edgeID]->to == otherNode);


			  assert(backward_edges[edgeID]->from==otherNode);
			  assert(backward_edges[edgeID]->to == node);
		  }



		  makeRoot(node);//this should be avoidable
		  makeRoot(otherNode);
		  dbg_printTour(node);
		  dbg_printTour(otherNode);
		  node->dbg_tour();
		  otherNode->dbg_tour();

		  Tree::Node * tleft = nullptr;
		  if(node->first())
			  tleft=t.splitBefore(node->first());
		  if(tleft)
			  assert(t.findRoot(tleft)==tleft);
		  if(node->first())
			  assert(t.findRoot(node->first())==node->first());
		  if(node->first())
			  t.concat(node->first(),forward_edges[edgeID]->node);
		  else{
			  node->setFirst(forward_edges[edgeID]);
		  }
		  if(otherNode->first()){
			  assert(t.findRoot( otherNode->first()) ==t.findRoot(otherNode->last()));

			  t.concat(forward_edges[edgeID]->node,otherNode->first());
			  assert(t.findRoot( otherNode->first()) ==t.findRoot(otherNode->last()));
			  assert(t.findRoot( forward_edges[edgeID]->node) ==t.findRoot(otherNode->last()));
		  }else
			  otherNode->setFirst(forward_edges[edgeID]);

		  if(otherNode->last()){
			  assert(t.findRoot( otherNode->first()) ==t.findRoot(otherNode->last()));
			  t.concat(otherNode->last(),backward_edges[edgeID]->node);
			  assert(t.findRoot( otherNode->first()) ==t.findRoot(backward_edges[edgeID]->node));
			  assert(t.findRoot( node->first()) ==t.findRoot(backward_edges[edgeID]->node));
		  }else{
			  otherNode->setLast(backward_edges[edgeID]);
			  t.concat( node->first(), backward_edges[edgeID]->node);
		  }
		  if(tleft){
			  t.concat(backward_edges[edgeID]->node,tleft);
			  node->setLast(tleft->value);
		  }else{
			  node->setLast(backward_edges[edgeID]);
		  }
		  node->dbg_insert(otherNode);
		  dbg_printTour(node);
		  assert(t.findMin(t.findRoot(node->first())) == node->first() );
		  Tree::Node * max = t.findMax(t.findRoot(node->first()));
		  Tree::Node * l = node->last();
		  assert(t.findMax(t.findRoot(node->first())) == node->last() );

		  assert(t.findRoot( forward_edges[edgeID]->node) == t.findRoot(backward_edges[edgeID]->node));

		  node->dbg_tour();
		  otherNode->dbg_tour();
		  dbg_printTour(findRoot(node));
		  findRoot(node)->dbg_tour();

		  assert(otherNode->subtree_size==otherNode->dbg_getSize());
		  int addedNodes = otherNode->subtree_size;
		  //update subtree sizes going up the tree
		  EulerVertex * p = node;
		  while(p){
			  p->subtree_size+= addedNodes;
			  int dbg = p->dbg_getSize();
			  assert(p->subtree_size==p->dbg_getSize());
			  p = getParent(p);
		  }
		  assert(node->subtree_size==node->dbg_getSize());

		  assert(connected(node,otherNode));
/*		  EulerVertex * r = findRoot(node);
		  findRoot(node)->dbg_tour();
		  Treap::Node * last = r->last();
		  Treap::Node * reallast = t.last(t.findRoot(r->last()));*/
		  assert(findRoot(node)->isRoot());
		  //Return half edge
		  return  edgeID;
		}


	void tour(EulerVertex* v, vec<int> & tour_out){
		v->tour(tour_out);
	}

	void dbg_printTour(EulerVertex * v){
		dbg_printDbgTour(v);

		vec<int> tour_list;
		v->tour(tour_list);
		printf("tour:");
		for(int i:tour_list){
			printf("%d,",i);
		}
		printf("\n");
	}
	void dbg_printDbgTour(EulerVertex * v){
		vec<int> tour_list;
		v->dbg_build_tour(v,tour_list);
		printf("dbgtour:");
		for(int i:tour_list){
			printf("%d,",i);
		}
		printf("\n");
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


	bool isRoot(int v){
		return vertices[v]->isRoot();
	}

	EulerVertex* getVertex(int n){
		return vertices[n];
	}

	EulerVertex * createVertex() {
	  nComponents++;
	  vertices.push(new EulerVertex());
#ifndef NDEBUG
	  vertices.last()->owner = this;
#endif
	  vertices.last()->index = vertices.size()-1;
	  return vertices.last();
	}
};

#endif



