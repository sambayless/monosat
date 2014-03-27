//translated from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/treap.js
//untested at this point.
/**
The MIT License (MIT)

Copyright (c) 2013 Mikola Lysenko

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
**/
#ifndef TREAPCUSTOM_H_
#define TREAPCUSTOM_H_
#include <stdlib.h>
#include "mtl/Rnd.h"
#include "BST.h"
using namespace Minisat;
template<class Value=int>
class TreapCustom{
public:
	double seed;
	struct Node{

		  Value value;

		  bool flag;
		  int flagAggregate;
		  int count;
		   int priority;
		   Node * parent;
		   Node *  left;
		   Node *   right;
		   Node *  next;
		   Node *  prev;
		   Node(Value _value, bool _flag,bool _flagAggregate, int valueCount, int _priority, Node * _parent, Node * _left, Node * _right, Node * _next, Node * _prev):
		   value(_value),flag(_flag),flagAggregate(_flagAggregate),count(valueCount),priority(_priority),parent(_parent),left(_left),right(_right),next(_next),prev(_prev)
		   {


		   }
	};

	Node * createNode(Value value){
		int rnd = irand(seed,10000);
		return new Node(value,false,false,0,rnd,NULL,NULL,NULL,NULL,NULL);
	}

void bubbleUp(Node * node) {
  while(true) {
    Node* p = node->parent;
    if(!p || p->priority < node->priority) {
      break;
    }
    if(node == p->left) {
      Node* b = node->right;
      p->left = b;
      if(b) {
        b->parent = p;
      }
      node->right = p;
    } else {
      Node* b = node->left;
      p->right = b;
      if(b) {
        b->parent = p;
      }
      node->left = p;
    }
    update(p);
	update(node);
    Node* gp = p->parent;
    p->parent = node;
    node->parent = gp;
    if(gp) {
      if(gp->left == p) {
        gp->left = node;
      } else {
        gp->right = node;
      }
    }
  }
  Node * p = node->parent;
  while(p) {
    update(p);
    p = p->parent;
  }
}

Node* findRoot(Node* node) {
  Node * n = node;
  while(n->parent) {
    n = n->parent;
  }
  return n;
}
Node * first(Node * node){
	Node * l = findRoot(node);
	  while(l->left) {
	    l = l->left;
	  }
	  return l;
}
Node * last (Node * node){
	Node * r = findRoot(node);
	  while(r->right) {
	    r = r->right;
	  }
	  return r;
}




Node * insert(Node* node, Value  value, int valueCount){
	 if(!node->right) {
		 int rnd = irand(seed,10000);
		 Node * nn = node->right = new Node(value, false, false, valueCount, rnd, node, NULL, NULL, node->next, node);
	    if(node->next) {
	      node->next->prev = nn;
	    }
	    node->next = nn;
	    bubbleUp(nn);
	    return nn;;
	  }
	 int rnd = irand(seed,10000);
	 Node * v = node->next;
	 v->left = new Node(value, false, false,valueCount, rnd, v, NULL, NULL, v, node);
	 Node * nn = v->left;
	  v->prev = nn;
	  node->next = nn;
	  bubbleUp(nn);
	  return nn;
	}

void swapNodes(Node* a, Node*b) {
  int p = a->priority;
  a->priority = b->priority;
  b->priority = p;
  Node* t = a->parent;
  a->parent = b->parent;
  if(b->parent) {
    if(b->parent->left == b) {
      b->parent->left = a;
    } else {
      b->parent->right = a;
    }
  }
  b->parent = t;
  if(t) {
    if(t->left == a) {
      t->left = b;
    } else {
      t->right = b;
    }
  }
  t = a->left;
  a->left = b->left;
  if(b->left) {
    b->left->parent = a;
  }
  b->left = t;
  if(t) {
    t->parent = b;
  }
  t = a->right;
  a->right = b->right;
  if(b->right) {
    b->right->parent = a;
  }
  b->right = t;
  if(t) {
    t->parent = b;
  }
  t = a->next;
  a->next = b->next;
  if(b->next) {
    b->next->prev = a;
  }
  b->next = t;
  if(t) {
    t->prev = b;
  }
  t = a->prev;
  a->prev = b->prev;
  if(b->prev) {
    b->prev->next = a;
  }
  b->prev = t;
  if(t) {
    t->next = b;
  }
  int c = a->count;
  a->count = b->count;
  b->count = c;
  bool f = a->flag;
  a->flag = b->flag;
  b->flag = f;
  f = a->flagAggregate;
  a->flagAggregate = b->flagAggregate;
  b->flagAggregate = f;
}

void update(Node * node){

	  int c = 0;// node->value;//countOfValue(node->value);
	  bool f = node->flag;
	  if(node->left) {
	    c += node->left->count;
	    f = f || node->left->flagAggregate;
	  }
	  if(node->right) {
	    c += node->right->count;
	    f = f || node->right->flagAggregate;
	  }
	  node->count = c;
	  node->flagAggregate = f;

}

void setFlag(Node * node, bool f) {
	node->flag = f;
 for(Node* v=node; v; v=v->parent) {
   bool pstate = v->flagAggregate;
   update(v);
   if(pstate == v->flagAggregate) {
     break;
   }
 }
}

void remove(Node * node) {

  if(node->left && node->right) {
    Node* other = node->next;
    swapNodes(other, node);
  }
  if(node->next) {
    node->next->prev = node->prev;
  }
  if(node->prev) {
    node->prev->next = node->next;
  }
  Node* r = NULL;
  if(node->left) {
    r = node->left;
  } else {
    r = node->right;
  }
  if(r) {
    r->parent = node->parent;
  }
  if(node->parent) {
    if(node->parent->left == node) {
      node->parent->left = r;
    } else {
      node->parent->right = r;
    }
    //Update all ancestor counts
    Node* p = node->parent;
    while(p) {
      update(p);
      p = p->parent;
    }
  }
  //Remove all pointers from detached node
  node->parent = node->left = node->right = node->prev = node->next = NULL;
  node->count = 1;
}
Node * split(Node * node) {

 Node * s = insert(node,Value(),0);
 s->priority = -100000;
 bubbleUp(s);
 Node * l = s->left;
 Node * r = s->right;
 if(l) {
   l->parent = NULL;
 }
 if(r) {
   r->parent = NULL;
 }
 if(s->prev) {
   s->prev->next = NULL;
 }
 if(s->next) {
   s->next->prev = NULL;
 }
 return r;
}

Node * concatRecurse(Node *a,Node * b) {
 if(a == NULL) {
   return b;
 } else if(b == NULL) {
   return a;
 } else if(a->priority < b->priority) {
   a->right = concatRecurse(a->right, b);
   a->right->parent = a;
   update(a);
   return a;
 } else {
   b->left = concatRecurse(a, b->left);
   b->left->parent = b;
   update(b);
   return b;
 }
}

 Node * concat(Node * node,Node * other){

	   if(!other) {
	     return NULL;
	   }
	   Node* ra = findRoot(node);
		Node* ta = ra;
	   while(ta->right) {
	     ta = ta->right;
	   }
	   Node* rb = findRoot(other);
	   assert(rb!=ra);
		Node* sb = rb;
	   while(sb->left) {
	     sb = sb->left;
	   }
	   ta->next = sb;
	   sb->prev = ta;
	   Node* r = concatRecurse(ra, rb);
	   r->parent = NULL;
	   return r;

 }



};
#endif /* TREAPCUSTOM_H_ */
