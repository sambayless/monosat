#ifndef LINK_CUT
#define LINK_CUT

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


};




import java.io.*;
import java.util.*;

/* Solution to http://www.codeforces.com/contest/117/problem/E
   using link/cut trees.  See detailed documentation below.

   D. Sleator, November 2011
*/
/*

public class ESubmit {
    int n, m, q;
    int[] ss, ff, first, next, qs, qe;

    boolean[] on_cycle;
    int[] parent;
    int cycle_vertex = -1;

    void dfs(int p, int v) {
	if(parent[v]>=0) {
	    if (cycle_vertex == -1) {
		cycle_vertex = v;
		parent[v] = p;
	    }
	    return;
	}
	parent[v] = p;
	for (int e = first[v]; e>=0; e = next[e]) {
	    int w = ff[e];
	    if (w==p) continue;
	    dfs(v,w);
	}
    }

    void treedfs(int p, int v) {
	parent[v] = p;
	for (int e = first[v]; e>=0; e = next[e]) {
	    int w = ff[e];
	    if (w==p || on_cycle[w]) continue;
	    treedfs(v,w);
	}
    }

    public void read_input_and_solve() throws IOException, MyException {
	m = n = nextInt();
	q = nextInt();

	 We use the integer array representation of an undirected graph.
	   The vertices are numbered 0..n-1.  There are two versions of
	   each edge, one in each direction numbered 0 to 2m-1.  Given an
	   edge e, the two end vertices of it are ss[e] and ff[e].  The
	   edges out of a vertex are linked together in a linked list.  The
	   first edge in the list for vertex v is first[v].  The next one
	   is next[first[v]] etc.

	   The advantage this has over keeping an adjacency array (of
	   neighboring vertices) for each vertex is that given an edge we
	   can find its reverse immediately.  If an edge has index j, then
	   its reverse has index j+m mod 2*m.  There's an alternative
	   design where the reverse of j is j^1.  (This kind of thing is
	   standard in competition programming.)


	ss = new int[2*m];   one end of an edge
	ff = new int[2*m];   the other end of it
	qs = new int[q];
	qe = new int[q];

	first = new int[n];    the first edge out of a vertex
	next = new int[2*m];   the next edge out of the same vertex
	Arrays.fill(first, -1);

	for (int i=0; i<m; i++) {
	    ss[i] = ff[i+m] = nextInt()-1;
	    ff[i] = ss[i+m] = nextInt()-1;
	}
	for (int i=0; i<2*m; i++) {
	    next[i] = first[ss[i]];
	    first[ss[i]] = i;
	}
	// Done reading in the graph

	for (int i=0; i<q; i++) {
	    qs[i] = nextInt()-1;
	    qe[i] = nextInt()-1;
	}
	// Done reading in the queries

	parent = new int[n];
	Arrays.fill(parent, -1);
	on_cycle = new boolean[n];
	int nc=0;  // the length of the cycle

	dfs(n,0);
	if (cycle_vertex < 0) throw new MyException ("no cycle");
	for (int v = parent[cycle_vertex]; true; v = parent[v]) {
	    on_cycle[v] = true;
	    nc++;
	    if (v==cycle_vertex) break;
	}

	int[] nodemap = new int[n];
	    // nodemap[i] is the index into the cyclenode array where i is stored

	Node[] cyclenodes = new Node[nc];
	for (int i=0, v = parent[cycle_vertex]; true; v = parent[v], i++) {
	    cyclenodes[i] = new Node(1, v);
	    if (i>0) LinkCut.link(cyclenodes[i-1], cyclenodes[i]);
	    nodemap[v] = i;
	    if (v==cycle_vertex) break;
	}

	Arrays.fill(parent, -1);
	for (int v=0; v<n; v++) if (on_cycle[v]) treedfs(-1,v);

	Node[] treenodes = new Node[n];
	for (int i=0; i<n; i++)
	    treenodes[i] = new Node((on_cycle[i])?0:1, i);
	for (int i=0; i<n; i++)
	    if (!on_cycle[i]) LinkCut.link(treenodes[i], treenodes[parent[i]]);

	int cycle_edges = 0;
	int tree_edges = 0;
	for (int i=0; i<q; i++) {
	    tree_edges += LinkCut.toggle(treenodes[qs[i]]);
	    int place1 = nodemap[LinkCut.rootid(treenodes[qs[i]])];
	    tree_edges += LinkCut.toggle(treenodes[qe[i]]);
	    int place2 = nodemap[LinkCut.rootid(treenodes[qe[i]])];

	    // Now we have to figure out which way around the cycle it goes
	    // place1 and place2 are positions around the cycle.
	    boolean path_right;
	    int right_length = (place2-place1+nc) % nc;
	    if (2*right_length < nc) path_right = true;
	    else if (2*right_length > nc) path_right = false;
	    else {
		int nf = (place1 + 1) % nc;
		int nb = (place1 - 1 + nc) % nc;
		int idf = cyclenodes[nf].id;
		int idb = cyclenodes[nb].id;
		path_right = (idf < idb);
	    }
	    cycle_edges += LinkCut.toggle(cyclenodes[place1]);
	    cycle_edges += LinkCut.toggle(cyclenodes[place2]);
	    if (path_right) {
		if (place1 > place2) cycle_edges += LinkCut.toggle(cyclenodes[0]);
	    } else {
		if (place1 < place2) cycle_edges += LinkCut.toggle(cyclenodes[0]);
	    }
	    int components = n - tree_edges - cycle_edges;
	    if (cycle_edges == nc) components++;
	    System.out.println(components);
	}
    }

     Boilerplate stolen from Petr, because scanner is too slow
    public static void main(String[] args) {
	(new ESubmit()).run();
    }

    BufferedReader reader;
    StringTokenizer tokenizer;
    PrintWriter writer;

    public void run() {
	try {
	    reader = new BufferedReader(new InputStreamReader(System.in));
	    tokenizer = null;
	    writer = new PrintWriter(System.out);
	    read_input_and_solve();
	    reader.close();
	    writer.close();
	} catch (Exception e) {
	    e.printStackTrace();
	    System.exit(1);
	}
    }

    int nextInt() throws IOException {
	return Integer.parseInt(nextToken());
    }

    long nextLong() throws IOException {
	return Long.parseLong(nextToken());
    }

    double nextDouble() throws IOException {
	return Double.parseDouble(nextToken());
    }

    String nextToken() throws IOException {
	while (tokenizer == null || !tokenizer.hasMoreTokens()) {
	    tokenizer = new StringTokenizer(reader.readLine());
	}
	return tokenizer.nextToken();
    }
}

class MyException extends Exception {
    public String s;
    public MyException(String prob) {
	super(prob);
	s = prob;
    }
}
*/

//http://codeforces.com/contest/117/submission/860934
/* The following two classes implement the Link/Cut trees of Sleator and Tarjan,
   as described in http://www.cs.cmu.edu/~sleator/papers/self-adjusting.pdf

   I first saw this elegant approach to implementing them in a submission
   by Fan Haoqiang, who used it to solve this problem:
   http://www.ioi2011.or.th/hsc/tasks/EN/elephants.pdf His code is here:
   http://builtinclz.abcz8.com/showcode.php?id=2011/IOI2011/elephants.cpp

   The goal of link/cut trees is to represent a forest of trees under these
   "kinds" of operations:

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

   Every node in the virtual tree has the usual parent, left, and right
   child pointers.  How do we tell if a node is a root of a splay tree?
   Easy.  It's a root if (1) its parent pointer is null or (2) if its
   parent's left and right children ARN'T it.  This is computed by the
   isroot() function below.

   Now that we know how to identify the root, it's easy to implement the
   splay(p) operation that takes a node p and splays it to the root of its
   splay tree.

   And using splay() we can implement expose(p), which transforms the
   virtual tree in a way that puts node p at the root of the virtual tree.
   (The real tree, of course, does not change.)  Expose works by doing a
   sequence of splays and relinking operations. (Picture is really required
   to explain it.  Check out the code, which is quite short.)  After
   expose(p) the path from p to the root of the real tree consists of all
   the nodes to the right of p in its splay tree.  (Actually, the way the
   code below works, after expose(p), p is the root and also the leftmost
   node in its splay tree.)

   This design does not allow us to walk down the tree along the dashed
   edges.  But walking down is not necessary for the basic operations of
   expose, link, cut, path, findroot, etc.

   The code below is customized to support the operations needed for this
   specific problem.  Each node in the virtual tree stores the information
   about the edge from that node in the real tree to its parent in the real
   tree.  A node that is a root of its real tree does not correspond to any
   edge.  We capture this by giving each node an individual weight my_s.
   my_s is 1 for all nodes except a node that is the root of the real tree.
   The s field maintains the sum of all the my_s fields in the subtree of
   the splay tree rooted here.

   Each node also has a boolean value flip.  Toggling this bit implicitly
   changes the on/off state of all the edges in the subtree (of the splay
   tree) rooted here.  Each node has my_flip, which tells if the current
   edge is flipped or not.  And each node has an "on" value which keeps
   the total weight of all the on edges in the subtree rooted here.

   Using this representation it's very simple to flip the state of all
   the edges from a node p to the root of the real tree.  To see how this is
   done look at the code for toggle(p).

   It should be fairly easy to modify the code below to support various
   different operations.  Most of the work should be in setting up the node
   class, and adjusting normalize() and update().

                              Danny Sleator
*/
/*

class Node {
    int s, my_s, on, id;
    boolean flip, my_flip;
    Node l, r, p;

    Node (int c, int i) {
	id = i;
	s = my_s = c;
	on = 0;
	l = r = p = null;
	flip = my_flip = false;
    }

    boolean isroot() {
	return p==null || (p.l != this && p.r != this);
    }

     If this node is flipped, we unflip it, and push the change
       down the tree, so that it represents the same thing.
    void normalize() {
	if (flip) {
	    flip = false;
	    on = s-on;
	    my_flip = !my_flip;
	    if (l != null) l.flip = !l.flip;
	    if (r != null) r.flip = !r.flip;
	}
    }

     The tree structure has changed in the vicinity of this node
       (for example, if this node is linked to a different left
       child in a rotation).  This function fixes up the data fields
       in the node to maintain invariants.
    void update() {
	s = my_s;
	on = (my_flip)?my_s:0;
	if (l != null) {
	    s += l.s;
	    if (l.flip) on += l.s-l.on; else on += l.on;
	}
	if (r != null) {
	    s += r.s;
	    if (r.flip) on += r.s-r.on; else on += r.on;
	}
    }
}

class LinkCut {
    static void rotR (Node p) {
	Node q = p.p;
	Node r = q.p;
	q.normalize();
	p.normalize();
	if ((q.l=p.r) != null) q.l.p = q;
	p.r = q;
	q.p = p;
	if ((p.p=r) != null) {
	    if (r.l == q) r.l = p;
	    else if (r.r == q) r.r = p;
	}
	q.update();
    }

    static void rotL (Node p) {
	Node q = p.p;
	Node r = q.p;
	q.normalize();
	p.normalize();
	if ((q.r=p.l) != null) q.r.p = q;
	p.l = q;
	q.p = p;
	if ((p.p=r) != null) {
	    if (r.l == q) r.l = p;
	    else if (r.r == q) r.r = p;
	}
	q.update();
    }

    static void splay(Node p) {
	while (!p.isroot()) {
	    Node q = p.p;
	    if (q.isroot()) {
		if (q.l == p) rotR(p); else rotL(p);
	    } else {
		Node r = q.p;
		if (r.l == q) {
		    if (q.l == p) {rotR(q); rotR(p);}
		    else {rotL(p); rotR(p);}
		} else {
		    if (q.r == p) {rotL(q); rotL(p);}
		    else {rotR(p); rotL(p);}
		}
	    }
	}
	p.normalize(); // only useful if p was already a root.
	p.update();    // only useful if p was not already a root
    }

     This makes node q the root of the virtual tree, and also q is the
       leftmost node in its splay tree
    static void expose(Node q) {
	Node r = null;
	for (Node p=q; p != null; p=p.p) {
	    splay(p);
	    p.l = r;
	    p.update();
	    r = p;
	};
	splay(q);
    }

     assuming p and q are nodes in different trees and
       that p is a root of its tree, this links p to q
    static void link(Node p, Node q) throws MyException {
	expose(p);
	if (p.r != null) throw new MyException("non-root link");
	p.p = q;
    }

     Toggle all the edges on the path from p to the root
       return the count after - count before
    static int toggle(Node p) {
	expose(p);
	int before = p.on;
	p.flip = !p.flip;
	p.normalize();
	int after = p.on;
	return after - before;
    }

     this returns the id of the node that is the root of the tree containing p
    static int rootid(Node p) {
	expose(p);
	while(p.r != null) p = p.r;
	splay(p);
	return p.id;
    }
}
*/

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
    }
    //splay(x);
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
    assert(!x->parent);
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
   // and that p is a root of its tree, this links p to q
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
#endif
