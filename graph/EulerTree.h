//adapted from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/euler.js
//untested at this point.

#ifndef EULER_TREE_H
#define EULER_TREE_H


#include "TreapCustom.h"
#include "mtl/Vec.h"
#include <cstdio>
#include <algorithm>
#include "AugmentedSplayTree.h"
#include "SearchTree.h"
using namespace Minisat;

class EulerTree{
public:
	struct EulerHalfEdge;
	struct EulerVertex;
private:


	typedef AugmentedSplayTree<EulerHalfEdge*> Tree;
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
		bool contains(EulerVertex * v){
			return from==v||to==v;
		}
	};



	struct EulerVertex{
		//Value value;

		//Treap::Node * first;
		//Treap::Node * last;
		EulerHalfEdge * left_out;//the half edge leading into the vertex, from its parent. The 'to' field in this halfedge the first occurrence of the euler vertex in the tour.
								//IFF the vertex is root, then this is instead the half edge leading to its first child. In that case, the 'from' field is the first occurence of the vertex in the tour
		EulerHalfEdge * right_in;//the half edge returning to the vertex.  This half edge has the last occurrence of the euler vertex in the tour.
		int index;
/*		int subtree_size;
		bool subtree_has_incident_edges;
		bool has_incident_edges;*/
#ifndef NDEBUG
		EulerTree * owner;
		//Note: in an euler-tour tree representation, these values are not explicitly maintained
		EulerVertex * dbg_parent;
		vec<EulerVertex*>dbg_children;
		bool dbg_visited;
		vec<EulerVertex*> dbg_children_t;
#endif
		EulerVertex():left_out(nullptr),right_in(nullptr){

			index=0;
/*			subtree_size=1;
			has_incident_edges=false;
			subtree_has_incident_edges=false;*/
#ifndef NDEBUG
			dbg_visited=false;
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

		/*bool isRoot(){
			EulerVertex * parent = nullptr;
			if(isSingleton()){

				return true;
			}else if (isLeaf()){

				assert(incidentEdgeA()->value->from==dbg_parent);
				return false;
			}else{
				if(!incidentEdgeA()->prev()){
					assert(owner->t.findRoot(incidentEdgeA()) ==owner->t.findRoot(incidentEdgeB()));
					//assert(owner->t.first(owner->t.findRoot(first())) ==first());
					//assert(owner->t.last(owner->t.findRoot(last())) ==last());
				}
				return !( incidentEdgeA()->prev());
			}


		}*/

/*		bool isRoot(){
			//Only a node that is both a leaf and root doesn't have half edges stored.
			if(isSingleton()){
				return true;
			}

			assert((left_out->from == this)==(right_in->to == this));
			return left_out->from == this;
		}*/

/*		bool isLeaf(){
			if(isSingleton()){
				return true;
			}
			//leaf vertices are special. Instead of storing the edges leading out to their first and last child, they store the (only two) edges leading into them in the euler tour.
			//For non-singleton leaf vertices only, left_out is instead left_in, and right_in is instead right_out;
			assert(!left_out || (left_out->from==this || left_out->to==this));
			assert(!right_in || (right_in->from==this || right_in->to==this));
			assert((left_out->to == this)==(right_in->from == this));
			return left_out->to == this;
		}*/
		//This is an _ARBITRARY_ incident edge.
		Tree::Node * incidentEdgeA(){
			if(left_out){
				return left_out->node;
			}else{
				return nullptr;
			}
		}
		//This is an _ARBITRARY_ incident edge.
		Tree::Node * incidentEdgeB(){

			if(right_in){
				return right_in->node;
			}else{
				return nullptr;
			}
		}

	/*	Tree::Node * getIncomingEdgeA(){
			if(!left_out){
				return nullptr;
			}else{
				if (left_out->to==this){
					return left_out;
				}else{
					assert(left_out->node->next());
					assert(left_out->from == this);
					return left_out->node->prev();
				}
			}
		}*/

		void setIncidentEdgeA(EulerHalfEdge * n){
			assert(!n || (n->from==this || n->to==this));
			left_out=n;
		}
		void setIncidentEdgeB(EulerHalfEdge * n){
			assert(!n || (n->from==this || n->to==this));
			right_in=n;
		}


		//Get the next node in the tour.
		EulerVertex * getNext(){
			if (left_out)
				return left_out->to;
			else if(incidentEdgeB())
				return incidentEdgeB()->value->to;
			else
				return nullptr;
		}

		//
		void tour(vec<int> & tour_list){
			tour_list.clear();

			if(isSingleton() ){
				tour_list.push(index);
				return;
			}

			Tree::Node* n= owner->t.findMin(owner->t.findRoot(incidentEdgeA()));
			/*Tree::Node* l= owner->t.findMax(owner->t.findRoot(incidentEdgeA()));
			Tree::Node * f = n;
			Tree::Node * b = l;
			//ok, now find the first and last occurences of this vertex in the tour...
			while(n->value->to != this && n->value->from != this){
				n=n->next();
			}
			while(l->value->to != this && l->value->from != this){
				l=l->prev();
			}*/
	/*		if(n->next()->value->from == this ||  n->next()->value->to == this)
				n=n->next();

			if(l->prev()->value->from == this ||  l->prev()->value->to == this)
				l=l->prev();
*/


			while(n->next()){
				printf("(%d -> %d)",n->value->from->index,n->value->to->index);
				tour_list.push(n->value->from->index);
				assert( n->next());

				//printf("(%d,%d)\n",n->value->from->index,n->value->to->index);
				assert(n->next());
				n=n->next();
				assert(n);
			}
			printf("(%d -> %d)\n",n->value->from->index,n->value->to->index);
			tour_list.push(n->value->from->index);
			tour_list.push(n->value->to->index);

		}

		int dbg_getSize(){
#ifndef NDEBUG
			int s = 1;
			for(EulerVertex* c:dbg_children){
				s+=c->dbg_getSize();
			}
			//assert(s==subtree_size);
			return s;
#endif
			return 0;
		}

		void dbg_make_parent(EulerVertex * new_parent){
		#ifndef NDEBUG
					if(new_parent){
						assert(dbg_children.contains(new_parent));
						dbg_children.remove(new_parent);
						EulerVertex * p = dbg_parent;
						dbg_parent=new_parent;
						if(p){
							dbg_children.push(p);

							p->dbg_make_parent(this);

						}
					}
		#endif
				}

		void dbg_make_root(){
#ifndef NDEBUG
			if(dbg_parent){
				dbg_children.push(dbg_parent);
				dbg_parent->dbg_make_parent(this);

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
		void dbg_build_tour_helper(EulerVertex * r, vec<int> & tour){
			tour.push(r->index);
					for(EulerVertex * c:r->dbg_children){
						dbg_build_tour_helper(c,tour);
						tour.push(r->index);
					}
		}
		void dbg_build_tour(EulerVertex * r, vec<int> & tour){
			while(r->dbg_parent)
				r=r->dbg_parent;
			dbg_build_tour_helper(r,tour);
		}

		void dbg_clear(){
			assert(dbg_visited);
			dbg_visited=false;
			assert(dbg_children.size()==dbg_children_t.size());
			dbg_children_t.copyTo(dbg_children);
			dbg_children_t.clear();
			for(EulerVertex * v:dbg_children)
				v->dbg_clear();
		}

		void dbg_tour(){
		#ifndef NDEBUG
			if(!left_out){
				assert(!right_in);
				return;
			}

			EulerHalfEdge * from = owner->t.findMin(owner->t.findRoot(incidentEdgeA()))->value;

			EulerVertex * dbgFrom = this;
			while(dbgFrom->dbg_parent)
				dbgFrom=dbgFrom->dbg_parent;

			assert(from->contains(dbgFrom));
			vec<EulerVertex*> dbg_stack;
			dbgFrom->dbg_visited=true;
			dbg_stack.push(dbgFrom);
			while(from && dbg_stack.size()){

				EulerVertex* p = dbg_stack.last();
				bool found =false;
				assert(p->incidentEdgeA());
				assert(p->incidentEdgeB());
				for(EulerVertex* t:p->dbg_children){
					if(!t->dbg_visited && from->contains(t)){
						t->dbg_visited=true;
						found=true;
						dbg_stack.push(t);
						p->dbg_children_t.push(t);
						from=from->node->next()->value;
						break;
					}
				}
				if(!found){
					dbg_stack.pop();
					if(dbg_stack.size())
						assert(from->contains(p->dbg_parent));
					else{
						assert(!p->dbg_parent);
					}
					if(from->node->next())
						from=from->node->next()->value;
					else
						from=nullptr;
				}
			}
			assert(!from);
			assert(dbg_stack.size()==1);
			dbgFrom->dbg_clear();



				//build the tour
				vec<int> dbg_tour;
				dbg_build_tour(this,dbg_tour);
				vec<int> real_tour;
				tour(real_tour);
				assert(dbg_tour.size()==real_tour.size());
				/*for(int i = 0;i<real_tour.size();i++){
					//assert(dbg_tour[i]==real_tour[i]);
				}*/



		#endif
			}
	};




public:

	//vec<EulerVertex*> nodes;
/*	void cut(int edgeID){

		EulerHalfEdge * edge = forward_edges[edgeID];
		cut (edge);
	}*/
	//Cut an edge in the tree, splitting it into two
	void cut(int edgeID){

		nComponents++;

		EulerHalfEdge * f = forward_edges[edgeID];
		EulerHalfEdge * b = backward_edges[edgeID];
		EulerVertex * from = f->from;
		EulerVertex * to = f->to;
		dbg_printTour(from);
		dbg_printTour(to);

		f->from->dbg_tour();
		f->to->dbg_tour();

		assert(connected(f->from, f->to));
		assert(t.findRoot(f->node)==t.findRoot(b->node));

		if(t.compare(f->node, b->node)>0){
			std::swap(f,b);
		}

		//ok, f is before b in the tour now
		Tree::Node * p = f->node->prev();
		Tree::Node * n = b->node->next();
		Tree::Node * pn = f->node->next();
		Tree::Node * nn = b->node->prev();

		Tree::Node * t1 = t.splitBefore(f->node);
		assert(!t1|| t.findMax(t1)==p);
		assert(t.findRoot(t1)!=t.findRoot(f->node));
		Tree::Node * t2 = t.splitAfter(b->node);
		assert(!t2|| t.findMin(t2)==n);

		assert(t.findRoot(t1)!=t.findRoot(f->node));
		assert(t.findRoot(t2)!=t.findRoot(b->node));
		assert(t.findRoot(t1)!=t.findRoot(b->node));
		assert(t.findRoot(t2)!=t.findRoot(f->node));

		dbg_printTour(from);
		dbg_printTour(to);

		if(t1 && t2)
			t.concat(t1,t2);//rejoin the two ends of the outer tour

		//dbg_printTour(from);
		//dbg_printTour(to);

		assert(t.findRoot(t1)!=t.findRoot(f->node));
		assert(t.findRoot(t2)!=t.findRoot(b->node));
		assert(t.findRoot(t1)!=t.findRoot(b->node));
		assert(t.findRoot(t2)!=t.findRoot(f->node));
		dbg_printTour(from);
		dbg_printTour(to);
		//ok, now we need to pick new incident edges for both vertices
		//these are the previous and next edges

		assert(pn);assert(nn);

		if(pn->value->contains(from) && pn->value->contains(to)){
			//at least one of from, to is a leaf
			if((n||p) && ! (n&&p)){
				if(!n){
					n = t.findMin(t.findRoot(p));
				}else{
					p = t.findMax(t.findRoot(n));
				}
				assert(n!=p);
			}
			if(n){
				assert(p);
				dbg_printTour(from);
				dbg_printTour(to);
				assert(n->value->contains(from) ||n->value->contains(to) );
				assert(p->value->contains(from) ||p->value->contains(to) );

				if(n->value->contains(from)){
					assert(p->value->contains(from));
					from->setIncidentEdgeA(n->value);
					from->setIncidentEdgeB(p->value);
					to->setIncidentEdgeA(nullptr);
					to->setIncidentEdgeB(nullptr);
				}else{
					assert(n->value->contains(to));
					assert(p->value->contains(to));
					to->setIncidentEdgeA(n->value);
					to->setIncidentEdgeB(p->value);
					from->setIncidentEdgeA(nullptr);
					from->setIncidentEdgeB(nullptr);
				}

			}else{
				from->setIncidentEdgeA(nullptr);
				from->setIncidentEdgeB(nullptr);
				to->setIncidentEdgeA(nullptr);
				to->setIncidentEdgeB(nullptr);

				//both of these vertices are now singletons.
			}

		}else if(pn->value->contains(from)){
			assert(nn->value->contains(from));
			from->setIncidentEdgeA(pn->value);
			from->setIncidentEdgeB(nn->value);
			if((n||p) && ! (n&&p)){
				if(!n){
					n = t.findMin(t.findRoot(p));
				}else{
					p = t.findMax(t.findRoot(n));
				}
				assert(n!=p);
			}

			if(n){
				assert(p);
				assert(n->value->contains(to));
				assert(p->value->contains(to));
				to->setIncidentEdgeA(n->value);
				to->setIncidentEdgeB(p->value);
			}else{
				//to is now a singleton.
				to->setIncidentEdgeA(nullptr);
				to->setIncidentEdgeB(nullptr);
			}
		}else if(pn->value->contains(to)){
			assert(nn->value->contains(to));
			to->setIncidentEdgeA(pn->value);
			to->setIncidentEdgeB(nn->value);
			if((n||p) && ! (n&&p)){
				if(!n){
					n = t.findMin(t.findRoot(p));
				}else{
					p = t.findMax(t.findRoot(n));
				}
				assert(n!=p);
			}

			if(n){
				assert(p);
				assert(n->value->contains(from));
				assert(p->value->contains(from));
				from->setIncidentEdgeA(n->value);
				from->setIncidentEdgeB(p->value);
			}else{
				//to is now a singleton.
				from->setIncidentEdgeA(nullptr);
				from->setIncidentEdgeB(nullptr);
			}
		}

		/*else{
			assert(nn->value->contains(to));
			to->setIncidentEdgeA(pn->value);
			to->setIncidentEdgeB(nn->value);
			if(n && p){
				assert(n->value->contains(from));
				assert(p->value->contains(from));
				from->setIncidentEdgeA(n->value);
				from->setIncidentEdgeB(p->value);
			}else if (n){
				assert(n->value->contains(from));
				from->setIncidentEdgeA(n->value);
				p=t.findMax( t.findRoot(n));
				assert(p!=n);
				assert(p->value->contains(from));
				from->setIncidentEdgeB(p->value);
			}else if (p){
				assert(p->value->contains(from));
				from->setIncidentEdgeA(p->value);
				n=t.findMax( t.findRoot(p));
				assert(p!=n);
				assert(n->value->contains(from));
				from->setIncidentEdgeB(n->value);
			}else{
				//from is now a singleton.
				from->setIncidentEdgeA(nullptr);
				from->setIncidentEdgeB(nullptr);
			}
		}*/


		/*if(!t1 && ! t2){
			//one or both of the nodes is now a singleton
			f->from->setIncidentEdgeA(nullptr);
			f->from->setIncidentEdgeB(nullptr);
			f->to->setIncidentEdgeA(nullptr);
			f->to->setIncidentEdgeB(nullptr);

		}else if (!t1 || ! t2){
			Tree::Node * t3 = p?p:n;
			//at most one vertex is now a singleton... and we need to figure out which one, if either.

			assert(t3->value->contains(from) || t3->value->contains(to));
			EulerHalfEdge * m = t.findMin(t3)->value;
			if(m->contains(from)){
				from->setIncidentEdgeA(m);
				assert(t.findMax(t3)->value->contains(from));
				from->setIncidentEdgeB(t.findMax(t3)->value);
				//to->setIncidentEdgeA(nullptr);
				//to->setIncidentEdgeB(nullptr);
			}else{
				assert(m->contains(to));
				to->setIncidentEdgeA(m);
				assert(t.findMax(t3)->value->contains(to));
				to->setIncidentEdgeB(t.findMax(t3)->value);
				//from->setIncidentEdgeA(nullptr);
				//from->setIncidentEdgeB(nullptr);
			}
		}else{


			assert(p->value->contains(to) || p->value->contains(from));
			assert(n->value->contains(to) || n->value->contains(from));
			if(p->value->contains(from)){

				assert(!p->value->contains(to));
				assert(!t.findMin(n)->value->contains(to));
				assert(t.findMin(n)->value->contains(from));
				f->from->setIncidentEdgeA(p->value);
				f->from->setIncidentEdgeB(n->value);

				Tree::Node * fn = f->node->next();
				Tree::Node * bp = b->node->prev();

				assert(fn->value->to == f->to ||fn->value->from == f->to);
				assert(bp->value->to == f->to ||bp->value->from == f->to);

				f->to->setIncidentEdgeA(fn->value);
				f->to->setIncidentEdgeB(bp->value);

			}else{
				assert(p->value->contains(to));
				assert(n->value->contains(to));
				f->to->setIncidentEdgeA(p->value);
				f->to->setIncidentEdgeB(n->value);

				Tree::Node * fn = f->node->next();
				Tree::Node * bp = b->node->prev();

				assert(fn->value->to == f->from ||fn->value->from == f->from);
				assert(bp->value->to == f->from ||bp->value->from == f->from);

				f->from->setIncidentEdgeA(fn->value);
				f->from->setIncidentEdgeB(bp->value);
			}
		}*/
		dbg_printTour(from);
		dbg_printTour(to);
		//finally, remove these half edges from the inner tour
		t.splitAfter(f->node);
		t.splitBefore(b->node);
		//check if either of these has become a singleton
		if(from->incidentEdgeA() && t.size( t.findRoot(from->incidentEdgeA()) ) ==1){
			from->setIncidentEdgeA(nullptr);
			from->setIncidentEdgeB(nullptr);
		}
		if(to->incidentEdgeA() && t.size( t.findRoot(to->incidentEdgeA()) ) ==1){
			to->setIncidentEdgeA(nullptr);
			to->setIncidentEdgeB(nullptr);
		}
		assert(t.size(f->node)==1);
		assert(t.size(b->node)==1);

		t.setIncident(f->node,0);
		t.setIncident(b->node,0);
		assert(!connected(from,to));
#ifndef NDEBUG
		if (f->to->dbg_parent==f->from)
			f->to->dbg_remove();
		else
			f->from->dbg_remove();
#endif

		dbg_printTour(from);
		dbg_printTour(to);
		from->dbg_tour();
		to->dbg_tour();


		f->from->dbg_tour();
		f->to->dbg_tour();


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
	/*EulerVertex * getParent(EulerVertex * v){
		EulerVertex * parent = nullptr;
		if(v->isSingleton()){
			assert(!v->dbg_parent);
			return nullptr;
		}else if (v->isLeaf()){
			assert(v->incidentEdgeA()->value->from==v->dbg_parent);
			assert( v->incidentEdgeA()->value->from==v->dbg_parent);
			return v->incidentEdgeA()->value->from;
		}else{
			if(v->incidentEdgeA()->prev() ==nullptr){
				assert(!v->dbg_parent);
				return nullptr;
			}
			assert(v->incidentEdgeA()->prev()->value->to==v);
			assert(v->incidentEdgeA()->prev()->value->from==v->dbg_parent);

			assert( v->incidentEdgeA()->value->from==v);
			return v->incidentEdgeA()->prev()->value->from;
		}


		if(!v->first() ){
			//this is wrong... is it? even the root should still have first and last pointers... shouldn't it?
			//this is a root
			assert(!v->dbg_parent);
			return nullptr;
		}
		assert(v->first()->value->isforward);
		parent= v->first()->value->from;
		assert(parent==v->dbg_parent);
		return parent;
	}*/
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

	//Get the next node in the tour. This is log time using splay trees!
	EulerVertex * getNext(EulerVertex * from){
		return from->getNext();
	}

/*
	EulerVertex * findRoot(EulerVertex * v){
		EulerVertex * parent = getParent(v);
		while(parent){
			v = parent;
			parent = getParent(v);
		}
		return v;
	}
*/

	//Return the size of the subtree rooted at v, including v
/*	int getSubtreeSize(EulerVertex * v){
		assert(v->subtree_size==v->dbg_getSize());
		return v->subtree_size;
	}*/

	//return the size of the complete tree that v is an element of
	int getFullTreeSize(EulerVertex * v){
		if(v->isSingleton())
			return 1;
		else
			return t.size( t.findRoot(v->incidentEdgeA()));
	}

	/*bool hasIncidentEdges(EulerVertex* v){
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
	}*/

	bool connected(EulerVertex*  from, EulerVertex*  to){
		if(from==to){
			return true;
		}
		if(from->incidentEdgeA()==nullptr || to->incidentEdgeA()==nullptr){
			return false;
		}
		return t.findRoot(from->incidentEdgeA())==t.findRoot(to->incidentEdgeA());
	}

	void makeRoot(EulerVertex*  node){
		node->dbg_tour();

		if(!node->isSingleton()){
			if(t.size( t.findRoot(node->incidentEdgeA()))==2){
				//then both nodes are equivalent to being root.
				node->dbg_make_root();
				assert(t.findRoot( node->incidentEdgeA()) ==t.findRoot(node->incidentEdgeB()));

				node->dbg_tour();
				return;
			}

			node->dbg_tour();
			dbg_printTour(node);
			//split the tree at the two incident edges
			EulerHalfEdge * f = node->incidentEdgeA()->value;
			EulerHalfEdge * b =  node->incidentEdgeB()->value;
			assert(f!=b);
			if(t.compare(f->node, b->node)>0){
				std::swap(f,b);
			}
			assert(t.size(t.findRoot(f->node))%2==0);//full tree is always even

			//we have to do something slightly tricky here, which is to figure out whether f, or one if its neighbours, is the edge we should be splitting on.
			//this is tricky because we aren't keeping track of this information explicitly with extra 'vertex' nodes in the binary search tree, t, which is the usual solution

			//t.splay(f->node);
			//assert(f->node->right);//because f is less than the other incident edge, so it can't be the rightmost edge.

			//if there is a node to the right, then there is also a successor node (which may or may not be f->node->right).
			EulerVertex * other = (f->to==node) ? f->from:f->to; assert(other!=node);
			EulerHalfEdge * next = f->node->next()->value;
			//the idea is to look at up to two adjacent nodes to 'node'. This will give us enough information to figure out where to split the tour.

			if(!next->contains(node)){
				//1) next does NOT contain node. in that case, we need to split on the previous node (if it exists)
				Tree::Node * fn =f->node->prev();
				if(!fn){
					//node is ALREADY the root, don't do anything.
					node->dbg_make_root();
					assert(t.findRoot( node->incidentEdgeA()) ==t.findRoot(node->incidentEdgeB()));

					node->dbg_tour();
					return;
				}else
					f = fn->value;
			}else if(next->contains(other)){
				//2) next connects the exact same nodes as f, in which case one of them is a leaf.
				//we need to figure out which one is the leaf and which one is the node, to figure out where to split.
				//The node is the only one that will be an element of a third edge
				Tree::Node * next_next = next->node->next();
				if(!next_next)
					next_next=f->node->prev();

				if(next_next->value->contains(node)){
					assert(other->dbg_parent==node);
					f = next;//We want to cut AFTER returning to node.
				}else{
					//node is the leaf, so f is the correct place to cut.
				}

			}else{
				//3) next contains node, but not the other vertex of f.
				//then f is the right place to cut
			}

			//ok, f is before b in the tour. Now we are going to split the tour after f, and then append the section that ends with f  to the right hand side of the tour.
			Tree::Node * right = t.splitAfter(f->node);

			dbg_printEdge(f->node);


			if(right){
				//assert(t.size(right)%2==0);
				assert(t.findMax(t.findRoot(f->node))==f->node);
				//assert(t.findMin(t.findRoot(right))->value->to == node || t.findMin(t.findRoot(right))->value->from == node );
				assert(t.findRoot(right)!=t.findRoot(f->node));
				t.concat(right,f->node);//'right' was previously the right hand side of the tour. After this operation, it is the left hand side
			}else{
				//this already is root
			}
			assert(t.findMax(t.findRoot(f->node))==f->node);

			//assert(t.findMin(t.findRoot(f->node))->value->to ==node || t.findMin(t.findRoot(f->node))->value->from ==node );

			node->dbg_make_root();
			assert(t.findRoot( node->incidentEdgeA()) ==t.findRoot(node->incidentEdgeB()));
			dbg_printTour(node);
			node->dbg_tour();
		}


	}

	//Make othernode a child of node.
	int link(EulerVertex*  node, EulerVertex*  otherNode, int edgeID) {
		  assert (node !=otherNode);
		  assert(!connected(node,otherNode));
		  nComponents--;




		  node->dbg_tour();
		  otherNode->dbg_tour();

		  dbg_printTour(node);
		  dbg_printTour(otherNode);


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



		  makeRoot(node);
		  makeRoot(otherNode);
		  dbg_printTour(node);
		  dbg_printTour(otherNode);
		  node->dbg_tour();
		  otherNode->dbg_tour();

		  Tree::Node* f =nullptr;
		  if(node->incidentEdgeA())
			  f= t.findMin(t.findRoot(node->incidentEdgeA()));

		 // Tree::Node* b =t.findMax(t.findRoot(node->incidentEdgeA()));
/*
		  if(f ){
			  dbg_printEdge(f);
			  dbg_printEdge(b);
			  Tree::Node * nf = f->next();

			  if(nf && (nf->value->from==node || nf->value->to==node)){
				  f=nf;
			  }
			  Tree::Node * pb = b->prev();
			  if(pb && (pb->value->from==node || pb->value->to==node)){
				  b=pb;
			  }
		  }*/
		  dbg_printEdge(f);


		  //dbg_printEdge(b);

		  //ok, f is before b in the tour


		  Tree::Node * tleft = nullptr;
	/*	  if(f)
			  tleft=t.splitBefore(f);
		  if(tleft)
			  assert(t.findRoot(tleft)==tleft);
		  if(f)
			  assert(t.findRoot(f)==f);
		  if(f)
			  t.concat(f,forward_edges[edgeID]->node);
		  else{
			  node->setIncidentEdgeA(forward_edges[edgeID]);
		  }*/
		  if(f)
			  t.concat(f,forward_edges[edgeID]->node);
		  else
			  node->setIncidentEdgeA(forward_edges[edgeID]);

		  if(otherNode->incidentEdgeA()){
			  assert(t.findRoot( otherNode->incidentEdgeA()) ==t.findRoot(otherNode->incidentEdgeB()));

			  t.concat(forward_edges[edgeID]->node,otherNode->incidentEdgeA());
			  assert(t.findRoot( otherNode->incidentEdgeA()) ==t.findRoot(otherNode->incidentEdgeB()));
			  assert(t.findRoot( forward_edges[edgeID]->node) ==t.findRoot(otherNode->incidentEdgeB()));
		  }else
			  otherNode->setIncidentEdgeA(forward_edges[edgeID]);

		  if(otherNode->incidentEdgeB()){
			  assert(t.findRoot( otherNode->incidentEdgeA()) ==t.findRoot(otherNode->incidentEdgeB()));
			  t.concat(otherNode->incidentEdgeB(),backward_edges[edgeID]->node);
			  assert(t.findRoot( otherNode->incidentEdgeA()) ==t.findRoot(backward_edges[edgeID]->node));
			  assert(t.findRoot( node->incidentEdgeA()) ==t.findRoot(backward_edges[edgeID]->node));
		  }else{
			  otherNode->setIncidentEdgeB(backward_edges[edgeID]);
			  t.concat(node->incidentEdgeA() , backward_edges[edgeID]->node);
		  }
		  if(tleft){
			  t.concat(backward_edges[edgeID]->node,tleft);
			  node->setIncidentEdgeB(tleft->value);
		  }else{
			  node->setIncidentEdgeB(backward_edges[edgeID]);
		  }
		  node->dbg_insert(otherNode);
		  dbg_printTour(node);


		  assert(t.findRoot( forward_edges[edgeID]->node) == t.findRoot(backward_edges[edgeID]->node));

		  node->dbg_tour();
		  otherNode->dbg_tour();

		  t.findMin(t.findRoot(node->incidentEdgeA()))->value->to->dbg_tour();
		  t.findMin(t.findRoot(node->incidentEdgeA()))->value->from->dbg_tour();
		  //update subtree sizes going up the tree


		  assert(connected(node,otherNode));

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

	void dbg_printEdge(Tree::Node * e){
		if(!e)
			return;
		printf("(%d -> %d)\n",e->value->from->index,e->value->to->index);
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
/*	void cut(int u) {
		cut(vertices[u]);
	}*/
	bool connected(int u, int v){
		return connected(vertices[u],vertices[v]);
	}

	//Get the parent of the vertex in the euler tour representation. This is O(1).
/*
	int getParent(int v){
		EulerVertex * parent = getParent(vertices[v]);
		if(!parent){
			return -1;
		}else{
			return parent->index;
		}
	}
*/



/*
	int findRoot(int v){
		EulerVertex * parent = getParent(vertices[v]);
		while(parent){
			v = parent->index;
			parent = getParent(vertices[v]);
		}
		return v;
	}
*/

	int numComponents(){
		return nComponents;
	}

	//Return the size of the subtree rooted at v, including v
/*	int getSubtreeSize(int v){
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
	}*/

/*	void cutEdge(int edgeID){
		assert(forward_edges[edgeID]->from);
		assert(forward_edges[edgeID]->to->dbg_parent==forward_edges[edgeID]->from);
		cut(forward_edges[edgeID]->to);
		assert(!connected(forward_edges[edgeID]->to,forward_edges[edgeID]->from));
	}*/


	/*bool isRoot(int v){
		return vertices[v]->isRoot();
	}*/

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



