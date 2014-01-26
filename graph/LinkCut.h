#include <cstddef>
#include <cassert>
#include "mtl/Vec.h"
/**
15-451 Algorithms
Fall 2012
D. Sleator

Link/Cut trees    October 9, 2012

Full paper at: http://www.cs.cmu.edu/~sleator/papers/self-adjusting.pdf

-------------------------------------------------------------------------

The goal of link/cut trees is to represent a forest of trees under these
kinds of operations:

  link(p,q)   p is a tree root, and q is any node in another tree.
              make p a child of q.

  cut(p)      p is not a tree root.  Delete the edge from p to its
              parent, thus separating the tree in two

  path(p)     This just means to do something to the path from p to the
              root of p's tree.  This could be counting its length,
              finding the minimum cost edge on this path, adding a
              constant to all the costs on this path, etc.  (Or, in this
              case, it's going to be to toggle whether or not each edge
              on this path is currently turned on.)  The operations
              supportable are those you can do on a subsequence of nodes
              in a binary search tree.

  findroot(p) return the root of the tree containing p

All of these operations are supported by link/cut trees in O(log n) time
where n is the size of the tree.

I'll use the terminology from our paper.  The "real tree" is the tree
that the data structure is trying to represent.  It's the tree that the
API presents to the user.  The "virtual tree" is the actual tree
represented in the memory of the computer.  (We should probably should
have reversed these two terms, to better fit the analogy with the real
and virtual images in optics.)  The virtual tree has exactly the same
set of nodes as the real tree, they're just linked together in a
different way.

The virtual tree has dashed and solid edges.  The connected components
of solid edges are binary trees (I'll call them splay trees).  The splay
tree is represented with the usual parent, left, and right and child
pointers.  The root of each splay tree is linked via a dashed edge to
some other node in the virtual tree.

The relationship between the real tree and the virtual tree that
represents it can be made clear by explaining how to transform the
virtual tree into the real tree.  To do this, convert each splay tree
into a path by traversing it in left to right order.  This is now a path
of nodes up the real tree.  (In other words, if nodes a and b are in the
same splay tree, and b is the successor of a, then a is a child of b in
the real tree.)  The rightmost node of this path is linked to the node
that is attached to the root of that splay tree by a dashed edge.

Note that the real tree has been partitioned into solid paths (which
always go up the tree from child to parent) and dashed edges.  This
partitioning is determined by the algorithm and will change with time.
It is not under the control of the client.

Every node in the virtual tree has the usual parent, left, and right
child pointers.  How do we tell if a node is a root of a splay tree?
Easy.  It's a root if (1) its parent pointer is NULL or (2) if its
parent's left and right children ARN'T it.  This is computed by the
isroot() function below.

Now that we know how to identify the root, it's easy to implement the
splay(p) operation that takes a node p and splays it to the root of its
splay tree.

And using splay() we can implement expose(p), which transforms the
virtual tree in a way that puts node p at the root of the virtual
tree.  (The real tree, of course, does not change.)  Expose works by
doing a sequence of splays and relinking operations, called splices.
A splice converts the left solid edge down from a node to dashed and
converts one of the dashed edges down from that node to solid.
(Picture is really required to explain it.  Check out the code, which
is quite short.)  After expose(p) the path from p to the root of the
real tree consists of all the nodes to the right of p in its splay
tree.  (Actually, the way the code below works, after expose(p), p is
the root and also the leftmost node in its splay tree.)

[The expose algorithm is the weakest part of these notes, that must be
explained on the blackboard.]

This design does not allow us to walk down the tree along the dashed
edges.  But walking down is not necessary for the basic operations of
expose, link, cut, path, findroot, etc.

Theorem:
  A sequence of m link/cut operations on starting from a set of n
  separate nodes is O(m log n).

Proof:
  We need to assign a weight to each node so that we can carry out the
  analysis using the access lemma for splay trees.  Recall from the
  splay tree analysis that the size of a node is the total weight of
  all nodes in the subtree rooted there.  We will define the weights
  of the nodes so that the size of a node is the number of nodes in
  the virtual tree rooted there.

  To make this be the case, we simply assign the weight of a node to
  be 1 (for itself) plus the sizes of all the subtrees connected to it
  via dashed edges.  It's easy to see that with this definition of
  weight, a splice operation does not change the sizes of any nodes,
  and thus does not change the potential of the tree.

  Recall that the rank of a node is the binary log of the size of the
  node.  Let the potential function be (Sum rank(x)) for all nodes x.
  Here is the access lemma:

  Access Lemma: The amortized number of rotations done when splaying a
  node x (of rank r(x)) in a tree rooted at t (of rank r(t)) is
  bounded by 3(r(t) - r(x)) + 1.

  We will focus on the cost of the expose() operation.  (The other
  operations are easy to analyze using the bound (derived below) for
  expose() combined with an analysis of the potential function changes
  involved.)

  Let the cost of expose() be measured as the number of edges on the
  path from the exposed node to the root of the virtual tree.  (This
  is exactly the same as the number of rotations done in the expose.)

  We will prove that with this potential function, the amortized cost
  of expose() is at most 12m(log n) + 2m.

  Let's look at one expose().  Say that there are k splay operations
  until p gets to the root (see the code below).  These splays are
  called phase1 splays.  After this we splay q, this splay is a phase2
  splay.

  The cost of the phase1 splays telescope, because the starting size
  of one splay is greater than the ending size of the previous splay.
  So the cost of the phase1 splays is at most 3(log n) + k.

  The phase2 splay costs amortized 3(log n)+1.  So the total amortized
  cost of the expose is 6(log n)+k+1.  If we add this together for all
  m expose() operations we can write:

      total cost <= 6 m (log n) + m + (#splays in phase1)

  Where the latter term comes from summing all the "k"s from each
  expose() operation.  It remains to bound this term.

  Consider a modified analysis of splaying where we do not count the
  rotation done in the zig case.  In this case the access lemma bound
  becomes 3(r(t) - r(x)).  The total modified cost of splaying in a
  sequence of exposes is then at most 6 m (log n).  But the number of
  phase1 splays is at most the true cost of the phase2 splays.  (This
  is because if there are k phase1 splays in an expose() then the
  corresponding phase2 splay does k rotations.)  So the true cost of
  the phase2 splays is at most m more than the modified cost we
  defined in this paragraph (there are at most m zigs that have to be
  accounted for).  Thus:

              (#splays in phase1) <= 6 m (log n) + m.

  Putting this together gives us:

                total cost <= 12 m (log n) + 2m

  We should note in conclusion that the initial potential is 0 (all
  the nodes are separate, so the sizes are all 1).  And the final
  potential is positive.

  QED.


It should be fairly easy to modify the code below to support various
different operations.  Most of the work should be in setting up the node
class, and adjusting normalize() and update().
**/

/*
class LinkCut {


	struct Node {
	    int s;

	    int  id;

	    Node *l,* r,* p;

	    Node ( int i) {
			id = i;

			l = r = p = NULL;

	    }

	    bool isroot() {
	    	return p==NULL || (p->l != this && p->r != this);
	    }

	     The tree structure has changed in the vicinity of this node
	       (for example, if this node is linked to a different left
	       child in a rotation).  This function fixes up the data fields
	       in the node to maintain invariants.
	    void update() {

			if (l != NULL) {
				s += l->s;

			}
			if (r != NULL) {
				s += r->s;

			}
	    }
	};
	int nodes;
	Node* addNode(){
		return new Node(nodes++);
	}

     void rotR (Node* p) {
	Node* q = p->p;
	Node* r = q->p;

	if ((q->l=p->r) != NULL) q->l->p = q;
	p->r = q;
	q->p = p;
	if ((p->p=r) != NULL) {
	    if (r->l == q) r->l = p;
	    else if (r->r == q) r->r = p;
	}
	q->update();
    }

     void rotL (Node* p) {
	Node * q = p->p;
	Node * r = q->p;

	if ((q->r=p->l) != NULL) q->r->p = q;
	p->l = q;
	q->p = p;
	if ((p->p=r) != NULL) {
	    if (r->l == q) r->l = p;
	    else if (r->r == q) r->r = p;
	}
	q->update();
    }

     void splay(Node *p) {
	while (!p->isroot()) {
	    Node* q = p->p;
	    if (q->isroot()) {
		if (q->l == p) rotR(p); else rotL(p);
	    } else {
		Node* r = q->p;
		if (r->l == q) {
		    if (q->l == p) {rotR(q); rotR(p);}
		    else {rotL(p); rotR(p);}
		} else {
		    if (q->r == p) {rotL(q); rotL(p);}
		    else {rotR(p); rotL(p);}
		}
	    }
	}

	p->update();    // only useful if p was not already a root
    }
public:
     This makes node q the root of the virtual tree, and also q is the
       leftmost node in its splay tree
     void expose(Node* q) {
	Node* r = NULL;
	for (Node* p=q; p != NULL; p=p->p) {
	    splay(p);
	    p->l = r;
	    p->update();
	    r = p;
	};
	splay(q);
    }

     assuming p and q are nodes in different trees and
       that p is a root of its tree, this links p to q
     void link(Node *p, Node *q)  {
	expose(p);
	assert(p->r==NULL);
	//if (p.r != NULL) throw new Exception("non-root link");
	p->p = q;
    }

     this returns the id of the node that is the root of the tree containing p
     int rootid(Node* p) {
	expose(p);
	while(p->r != NULL) p = p->r;
	splay(p);
	return p->id;
    }
      p is not a tree root.  Delete the edge from p to its
                 parent, thus separating the tree in two
                 void cut(Node *p)  {
             	expose(p);
             	assert(p->r==NULL);
             	//if (p.r != NULL) throw new Exception("non-root link");
             	p->p = q;
                 }


};*/
//Watch out - LinkCut cannot by itself implement a dynamic disjoint set datastructure - you need more for that.

 class LinkCut {
	 int setCount;
  struct Node {
	int id;
    Node* left;
    Node* right;
    Node *parent;
    Node(int _id):id(_id), left(NULL),right(NULL),parent(NULL){};
  };

  vec<Node*> nodes;
  // Whether x is a root of a splay tree
  bool isRoot(Node *x) {
    return x->parent == NULL || (x->parent->left != x && x->parent->right != x);
  }

  void connect(Node* ch, Node* p, bool leftChild) {
    if (leftChild)
      p->left = ch;
    else
      p->right = ch;
    if (ch != NULL)
      ch->parent = p;
  }

  void rotR (Node* p) {
	Node* q = p->parent;
	Node* r = q->parent;

	if ((q->left=p->right) != NULL) q->left->parent = q;
	p->right = q;
	q->parent = p;
	if ((p->parent=r) != NULL) {
	    if (r->left == q) r->left = p;
	    else if (r->right == q) r->right = p;
	}

 }

  void rotL (Node* p) {
	Node * q = p->parent;
	Node * r = q->parent;

	if ((q->right=p->left) != NULL) q->right->parent = q;
	p->left = q;
	q->parent = p;
	if ((p->parent=r) != NULL) {
	    if (r->left == q) r->left = p;
	    else if (r->right == q) r->right = p;
	}

 }

  void splay(Node *p) {
 	while (!isRoot(p)) {
 	    Node* q = p->parent;
 	    if (isRoot(q)) {
 		if (q->left == p) rotR(p); else rotL(p);
 	    } else {
 		Node* r = q->parent;
 		if (r->left == q) {
 		    if (q->left == p) {rotR(q); rotR(p);}
 		    else {rotL(p); rotR(p);}
 		} else {
 		    if (q->right == p) {rotL(q); rotL(p);}
 		    else {rotR(p); rotL(p);}
 		}
 	    }
 	}
     }
  // Makes node x the root of the virtual tree, and also x is the leftmost
  // node in its splay tree
   Node *expose(Node* x) {
    Node *last = NULL;
    for (Node* y = x; y != NULL; y = y->parent) {
      splay(y);
      y->left = last;
      last = y;
    }
    splay(x);
    return last;
  }


   Node * _findRoot(Node * x) {
    expose(x);
    while (x->right != NULL) {
      x = x->right;
      splay(x);
    }
    return x;
  }

   bool dbgSetCount(){
	   int count = 0;
	   for(int i = 0;i<nodes.size();i++){
		   Node * n = nodes[i];
		   Node* r = _findRoot(n);
		   if(r==n){
			   count++;
		   }
	   }
	   return count==setCount;
   }

  // prerequisite: x and y are in distinct trees
    void _link(Node* x, Node* y) {
    //assert (_findRoot(x) != _findRoot(y));
#ifndef NDEBUG
	Node* sY = _findRoot(y);
    Node* sX = _findRoot(x);
    assert(sY!=sX);//else this is a bug
#endif

    setCount--;
    expose(x);

    x->parent = y;
    assert(dbgSetCount());
  }

    bool _connected(Node* x, Node *y) {
    if (x == y)
      return true;
    expose(x);
    expose(y);
    return x->parent != NULL;
  }

    void _cut(Node *x, Node *y) {
    expose(x);
    expose(y);
    if( x->parent != NULL){
    	setCount++;
    }
    assert(! (y->right != x || x->left != NULL || x->right != NULL));

    y->right->parent = NULL;
    y->right = NULL;
    assert(dbgSetCount());
  }

public:
    LinkCut():setCount(0){

    }

   int addNode(){
	   //return new Node();
	   setCount++;
	   nodes.push(new Node(nodes.size()));
	   return nodes.size()-1;
   }
   int nNodes(){
	   return nodes.size();
   }

   int findRoot(int x) {
    return _findRoot(nodes[x])->id;
  }

  // prerequisite: x and y are in distinct trees
    void link(int x, int y) {
    	if(x==y)
    		return;
    	Node * xnode = nodes[x];
    	Node * ynode = nodes[y];
    	_link(xnode,ynode);
  }

    bool connected(int x, int y) {
    if (x == y)
      return true;
     Node * xnode = nodes[x];
 	Node * ynode = nodes[y];
     expose(xnode);
     expose(ynode);
#ifndef NDEBUG
     int s1 = findRoot(x);
     int s2 = findRoot(y);
     bool dbg_connected = s1==s2;
     if(dbg_connected){
    	 assert(xnode->parent);
     }else{
    	 assert(xnode->parent==NULL);
     }

#endif
     return xnode->parent != NULL;




  }

    void cut(int x, int y) {
        Node * xnode = nodes[x];
    	Node * ynode = nodes[y];
    	_cut(xnode,ynode);
  }

    int numRoots(){
    	assert(dbgSetCount());
    	return setCount;
    }

    void reset(){
    	for(int i = 0;i<nodes.size();i++){
    		nodes[i]->parent=NULL;
    		nodes[i]->left=NULL;
    		nodes[i]->right =NULL;
    	}
    	setCount= nodes.size();
    }

};
