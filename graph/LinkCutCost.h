#ifndef LINK_CUT_COST
#define LINK_CUT_COST

#include <cstddef>
#include <cassert>
#include <climits>
#include <algorithm>
#include "mtl/Vec.h"
using namespace Minisat;
#ifdef NDEBUG
#define NDEBUG_LINKCUT
#endif
//This implementation of link/cut trees is based of a combination of D. Sleator's implementation (//http://codeforces.com/contest/117/submission/860934),
//and the (extremely well presented!) version described by Klein and Mozes, in the Appendix of Planarity.

//This link cut tree has been augmented to support operations for computing Maximum flows.
//Specifically, each edge in the abstract tree (ie, NOT the underlying splay tree) has an associated cost.
//For each node, we can compute in log time the cost of the edge connecting it to its parent, and also in log time compute
//the minimum cost of any edge on the path from the node to its parent.

//Note that is probably not the minimum cost that you might naturally expect the tree to compute: its is NOT the minimum cost of any node in the abstract tree beneath the node.
//This is especially a point of potential confusion, given that the way this minimum path cost is computed is by finding the minimum cost of any node in the concrete splay tree beneath the node in question (after exposing it).

 class LinkCutCost {
	static const int INF = INT_MAX/2;
	 int setCount;
  struct Node {
	int id;
	int netcost=INF;//netcost=grosscost(v)-grossmin(v)   cost of the edge connecting this node to its parent.
	int netmin=0;//netmin(v)=grossmin(v) if v is root, or grossmin(v)-grossmin(parent(v)) else   minimum cost of any edge below this one.


#ifndef NDEBUG_LINKCUT
	int cost=INF;
	int min=INF;
	int dbg_min=INF;
	Node *dbg_parent=nullptr;
#endif
	Node* left;
    Node* right;
    Node *parent;
    Node(int _id):id(_id),left(NULL),right(NULL),parent(NULL){};


  };

  int grossmin(Node * v){
	  dbg_min(v);
	  expose(v);
	  return v->netmin+v->netcost;
  }
  int grosscost(Node * v){
	  //int cost = v->netcost + grossmin(v);
	  //assert(cost==v->cost);
#ifndef NDEBUG_LINKCUT
	  int cp = v->cost;
#endif
	  expose(v);
	  dbg_print_forest();
	  //int it = iter;
	  int cost = v->netcost;
#ifndef NDEBUG_LINKCUT
	  assert(cost==cp);
#endif
	  return cost;
	  //return v->cost;
  }
  void update(Node * v) {
#ifndef NDEBUG_LINKCUT
     	int mincost = v->cost;

 		if (v->left) {
 			mincost = std::min(mincost,v->left->min);
 		}
 		if (v->right) {
 			mincost = std::min(mincost, v->right->min);
 		}
 		if(mincost==INF)
 			mincost=0;
 		v->min=mincost;

 		{
			v->dbg_min =v->cost;

			if (v->left) {
				v->dbg_min = std::min(v->dbg_min, v->left->dbg_min);
			}
			if (v->right) {
				v->dbg_min = std::min(v->dbg_min, v->right->dbg_min);
			}
			if(v->dbg_min==INF)
				v->dbg_min=0;
			assert(v->dbg_min==v->min);
 		}
 #endif
    }

  //Update the costs on every edge in the path of the tree from v to root v by delta.
  //has NO EFFECT if v is a root in the abstract tree (as roots have no edges on the path to themselves).
  void updatePathCost(Node * v, int delta){
	  expose(v);
	  v->netcost+=delta;
#ifndef NDEBUG_LINKCUT

	  while(v){
		  v->cost+=delta;
		  if(v->cost<v->min){
			  v->min=v->cost;
			  v->dbg_min=v->cost;
		  }

		  v=v->dbg_parent;
	  }
#endif
	  dbg_min(v);

  }


  void dbg_min(Node * v){
#ifndef NDEBUG_LINKCUT
	  if(!v)
		  return;
	  //This debug code only checks correctly if v is exposed!
	  if(v->parent){
		  return;
	  }

	  int m = v->netmin;
	  int dbg_min = v->cost;
	  Node * p =v;

	  int c = v->netcost;

	  while(p->parent){
		  if(p==p->parent->left || p==p->parent->right){
			  p=p->parent;

			  //m+=p->netmin;
			  c+=p->netcost;

		  }else{
			  break;
		  }
	  }

	  m+=c;//minw = delta_minw+w;

	  p =v;
	  dbg_min = std::min(dbg_min,p->cost);

	  while(p->dbg_parent){

			  p=p->dbg_parent;

			  dbg_min = std::min(dbg_min,p->cost);

	  }

	  int minGrossCost =v->cost;

	  vec<Node*> Q;
	  Q.push(v);
	  while(Q.size()){
		  Node * w = Q.last();
		  Q.pop();
		  minGrossCost = std::min(minGrossCost,w->cost);
		  if(w->right)
			  Q.push(w->right);

		  if(w->left)
			  Q.push(w->left);
	  }
	 // if(minGrossCost==INF)
	//	  minGrossCost=0;
	  assert(m==dbg_min);
	  assert(m==minGrossCost);






/*	  while (v->right != NULL) {
	       v = v->right;
	     }
	  dbg_isGrossMin(v->min,v);*/
#endif
  }
  bool dbg_is_ancestor(Node * v, Node * p){
	  if (v==p){
		  return true;
	  }
	  while(v->parent){
		  if(v->parent==p)
			  return true;
		  v=v->parent;
	  }
	  return false;
  }

public:
  void dbg_print_forest(bool force = false){
#ifndef NDEBUG_LINKCUT
	  int iter = 0;
	  if(!force)
		  return;
	 /* if(++iter<= 1415550){
	 		 return;
	 	  }
	  if(iter== 1415555){
		  int a =1;
	  }*/

		printf("digraph{\n");
		for(int i = 0;i<nodes.size();i++){
			printf("n%d [label=\"%d: c%d, m%d\"]\n", i,i,nodes[i]->netcost,nodes[i]->min);
		}

		for(int i = 0;i<nodes.size();i++){

			Node * n = nodes[i];
			if(n->parent){
				const char * s = "black";
				if(n==n->parent->left){
					s="blue";
					assert(n!=n->parent->right);
				}

				if(n==n->parent->right){
					s="red";
				}

				/*if(value(e.v)==l_True)
					s="blue";
				else if (value(e.v)==l_False)
					s="red";*/
				printf("n%d -> n%d [label=\"%d: %d\",color=\"%s\"]\n", i,n->parent->id, i, n->cost, s);
			}

		}

		printf("}\n");

#endif
  }

  void dbg_cost(Node * v){
#ifndef NDEBUG_LINKCUT
	  if(!v)
		  return;
	  int c = v->netcost;
	  Node * p =v;
	  while(p->parent){
		  if(p==p->parent->left || p==p->parent->right){
			  p=p->parent;
			  c+=p->netcost;
		  }else{
			  break;
		  }
	  }
	  assert(c==v->cost);
#endif
  }

  void dbg_all(){
#ifndef NDEBUG_LINKCUT
	  for(Node * v:nodes){
		  expose(v);
		  dbg_cost(v);
		  dbg_min(v);
	  }
#endif
  }

  void dbg_isGrossMin(int min,Node * v){
#ifndef NDEBUG_LINKCUT
	  static int iter = 0;
	  if(++iter==22349){
		  int a=1;
	  }
	 // dbg_print_forest();
	  int minGrossCost =v->cost;

	  vec<Node*> Q;
	 /* for(int i = 0;i<nodes.size();i++){
		  Node * y = nodes[i];
		  if(dbg_is_ancestor(y,v)){
			  minGrossCost = std::min(minGrossCost, grosscost(y));
			  assert(minGrossCost>=v->min);
		  }
	  }*/
	  Q.push(v);
	  while(Q.size()){
		  Node * w = Q.last();
		  Q.pop();
		  minGrossCost = std::min(minGrossCost, w->cost);
		  if(w->right)
			  Q.push(w->right);

		  if(w->left)
			  Q.push(w->left);
	  }
	  if(minGrossCost==INF)
		  minGrossCost=0;
	  assert(minGrossCost == min);
#endif
  }


  vec<Node*> nodes;
  // Whether x is a root of a splay tree
  bool isRoot(Node *x) {
    return x->parent == NULL || (x->parent->left != x && x->parent->right != x);
  }

 /* void connect(Node* ch, Node* p, bool leftChild) {
    if (leftChild)
      p->left = ch;
    else
      p->right = ch;
    if (ch != NULL)
      ch->parent = p;
  }
*/

  void rotR (Node* p) {
	  static int riter=0;
	if(++riter==23){
		int a=1;
	}
	Node* q = p->parent;
	Node* r = q->parent;
#ifndef NDEBUG_LINKCUT
	int cp = p->cost;
	int cq = q->cost;
	//int cr = r? grosscost(r):0;
	int cb = p->right? p->right->cost:0;

	if((p->id==60 && q->id==50) || (r && r->id==60 && p->id==50 ) ){
		int a=1;
	}
	dbg_cost(q);dbg_cost(p); dbg_cost(r); dbg_cost(q->right);

#endif

	//delta weight/minweight update rules from Klein and Mozes's Planarity, chapter 17
	int deltap = p->netcost;
	int deltaq = q->netcost;
	int deltab = p->right?p->right->netcost:0;

	if ((q->left=p->right) != nullptr){
		q->left->parent = q;
		//q->left->netcost+=p->netcost;
		//q->left->netcost-=q->netcost;
	}
	p->right = q;
	q->parent = p;



	p->netcost+=deltaq;
	q->netcost= -deltap;
	if(q->left){
		q->left->netcost+=deltap;
	}

	if ((p->parent=r) != nullptr) {
	    if (r->left == q){
	    	r->left = p;
	    }else if (r->right == q){
	    	r->right = p;
	    }
	   // p->netcost-=r->netcost;
	}


	update(q);
	childChange(q);
	childChange(p);
	if(r){
		//childChange(r);//this may be unneeded, because the sum of r's children can't have changed...
	}
	if(q->netcost==27 || p->netcost==27){
			int a=1;
		}
	dbg_cost(q);dbg_cost(p); dbg_cost(r); dbg_cost(q->right);
	//dbg_min(q);dbg_min(p);dbg_min(r); dbg_min(q->right);
 }

  void rotL (Node* p) {
	  static int liter=0;
	if(++liter==16498){//16512
		int a=1;
	}
	Node * q = p->parent;
	Node * r = q->parent;
#ifndef NDEBUG_LINKCUT
	int cp = p->cost;
	int cq = q->cost;
	//int cr = r? grosscost(r):0;
	int cb = p->left? p->left->cost:0;
	if((p->id==60 && q->id==50) || (r && r->id==60 && p->id==50 ) ){
		int a=1;
	}
#endif

	int deltap = p->netcost;
	int deltaq = q->netcost;
	int deltab = p->left?p->left->netcost:0;

	if ((q->right=p->left) != nullptr){
		q->right->parent = q;
		//q->right->netcost+=p->netcost;
		//q->right->netcost-=q->netcost;
	}
	p->left = q;
	q->parent = p;

	if ((p->parent=r) != nullptr) {
	    if (r->left == q){
	    	r->left = p;
	    }else if (r->right == q){
	    	r->right = p;
	    }
	    //p->netcost-=r->netcost;
	}

	p->netcost+=deltaq;
	q->netcost= -deltap;
	if(q->right){
		q->right->netcost+=deltap;
	}


	update(q);
	childChange(q);
	childChange(p);
	if(r){
		//childChange(r);
	}
	if(q->netcost==27 || p->netcost==27){
		int a=1;
	}
	dbg_cost(q);dbg_cost(p); dbg_cost(r); dbg_cost(q->right);
	//dbg_min(q);dbg_min(p);dbg_min(r); dbg_min(q->right);
 }

  void splay(Node *p) {
#ifndef NDEBUG_LINKCUT
	  dbg_print_forest();
	  int cp = p->cost;
#endif
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

 	update(p);

     }
  void childChange(Node * x){
	  x->netmin = std::min(0,x->left? (x->left->netmin + x->left->netcost) :0);
	  x->netmin = std::min(x->netmin,x->right? (x->right->netmin + x->right->netcost) :0);
	  if(x->netcost==27 && x->id==50){
		  int a=1;
	  }
  }
  	/* void dashedToSolid(Node * x,Node * p){
  		 assert(x->parent==p);
  		 assert(x!=p->right);
  		 assert(x!=p->left);
  		 assert(!p->parent);
  		 x->netcost-= p->netcost;
  	 }

  	 void solidToDashed(Node * x){
  		 assert(x->parent);
  		 assert(!x->parent->parent);
  		 x->netcost+=x->parent->netcost;
  	 }
  	 */
  // Makes node x the root of the virtual tree, and also x is the leftmost
  // node in its splay tree
   Node *expose(Node* x) {
    Node *last = nullptr;
    for (Node* y = x; y != NULL; y = y->parent) {
    	dbg_cost(y);
      splay(y);
      if(y->left){
    	  y->left->netcost+=y->netcost;
      }
      y->left = last;
      if(last){
    	  last->netcost-=y->netcost;
      }
      childChange(y);
      update(y);
      last = y;
      if(y->netcost==27){
    	  int a=1;
      }
    }
    dbg_cost(x);
    splay(x);
    if(x->netcost==27){
     	  int a=1;
       }
    dbg_cost(x);
    dbg_min(x);
   // dbg_all();
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

  // prerequisite: x and y are in distinct trees, AND x is the root of its tree
    void _link(Node* x, Node* y, int cost=0) {
    //assert (_findRoot(x) != _findRoot(y));
#ifndef NDEBUG_LINKCUT
	Node* sY = _findRoot(y);
    Node* sX = _findRoot(x);
    assert(sY!=sX);//else these are already linked.
    assert(x==sX);

#endif

    setCount--;
    expose(x);
    //expose(y);//need to expose y as well, otherwise netcost will get computed incorrectly below... nope, thats not true, because x doesn't become an attached child of y below
    int x_gross_min = grossmin(x);
    int parent_min = grossmin(y);

    assert(!x->parent);
#ifndef NDEBUG_LINKCUT
    assert(x->cost==INF);
#endif
    x->parent = y;
    x->netcost= cost;// - y->netcost;
#ifndef NDEBUG_LINKCUT

        if(x->left==nullptr && x->right==nullptr){
        	assert(x->netmin==0);
        	x->min=cost;
        }else{
        	x->min=std::min(x->min,cost);
        }
#endif
    x->netmin = std::min(x->netmin,cost);

#ifndef NDEBUG_LINKCUT
    x->cost=cost;
	x->dbg_parent = y;
#endif
    assert(grosscost(x)==cost);
    assert(dbgSetCount());

    dbg_min(y);
    dbg_all();
  }

    bool _connected(Node* x, Node *y) {
    if (x == y)
      return true;
    expose(x);
    expose(y);
    return x->parent != NULL;
  }
    void _cut(Node * x){
    	expose(x);
    	assert(x->right);//else, x is a root node
    	x->right->netcost+=x->netcost;
    	x->right->parent=nullptr;
    	update(x->right);
    	//childChange(x->right);
    	x->right=nullptr;
    	x->netcost=INF;
#ifndef NDEBUG_LINKCUT
    	x->cost=INF;
    	x->dbg_parent = nullptr;
#endif
    	childChange(x);
    	update(x);
    	setCount++;
    	dbg_all();
    }

    //u is the node to search from; a is a weight greater to or equal to the delta_minw(s)
	//returns the rightmost (or leftmost, if leftdir is true) solid descendent v of u such that w(v) <= a+w(u)
	Node* solidFind(Node * u, int alpha, bool leftdir){
		assert(alpha>=u->netmin);
		Node * u1 = u->left;
		Node * u2 = u->right;
		if(!leftdir){
			std::swap(u1,u2);
		}
		if(u1 && (alpha-u1->netcost) >= u1->netmin){
			return solidFind(u1,alpha-u1->netcost,leftdir);
		}
		if (alpha>=0){
			splay(u);
			return u;
		}
		return solidFind(u2,alpha-u2->netcost,leftdir);
	}

	//Find the ancestor of s with the lowest weight.
	Node * ancecstorFindWeight(Node * u,int alpha, bool leftdir){
		expose(u);
		bool u_is_candidate = u->netcost <=alpha;
		bool candidate_on_right = u->right && (u->right->netmin+u->right->netcost+u->netcost <= alpha);
		if (u_is_candidate && (leftdir || !candidate_on_right)){
			return u;
		}
		if(!candidate_on_right)
			return nullptr;
		return solidFind(u->right,alpha-u->right->netcost - u->netcost,leftdir);
	}

	Node * ancecstorFindMin(Node * u,bool leftdir){
		expose(u);
		int alpha = u->netcost+ u->netmin;
		bool u_is_candidate = u->netcost <=alpha;
		bool candidate_on_right = u->right && (u->right->netmin+u->right->netcost+u->netcost <= alpha);
		if (u_is_candidate && (leftdir || !candidate_on_right)){
			return u;
		}
		if(!candidate_on_right)
			return nullptr;
		return solidFind(u->right,alpha-u->right->netcost - u->netcost,leftdir);
	}

public:
    LinkCutCost():setCount(0){

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
    int r= _findRoot(nodes[x])->id;
    dbg_min(nodes[x]);
    return r;
  }

  // prerequisite: x and y are in distinct trees. y becomes the parent of x.
    void link(int x, int y, int cost = 0) {
    	if(x==y)
    		return;
    	if(x==39){
    	    		int a=1;
    	    	}
    	Node * xnode = nodes[x];
    	Node * ynode = nodes[y];
    	_link(xnode,ynode,cost);

  }

    bool connected(int x, int y) {
    if (x == y)
      return true;
     Node * xnode = nodes[x];
 	 Node * ynode = nodes[y];
     expose(xnode);
     expose(ynode);
#ifndef NDEBUG_LINKCUT
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

    //Returns the cost of the edge connecting x to the tree; x must not be a root.
    int getCost(int x){

    	return grosscost(nodes[x]);
    }

    //Returns the lowest cost in the tree rooted at x.
    //x must be a root!
    int minCost(int x){

    	return grossmin(nodes[x]);
    }

    //Find the minimum weight ancestor of s that is either closest to root(s), or closest to s.
    int ancecstorFindMin(int s, bool closestToRoot=true){
    	Node * n = ancecstorFindMin(nodes[s],!closestToRoot);
    	return n?  n->id:s;
    }

    //updates the cost of each element in x's path to root.
    void updateCostOfPathToRoot(int x, int delta){
    	updatePathCost(nodes[x],delta);
    }

    void cut(int x){
    	if(x==39){
    		int a=1;
    	}
    	_cut(nodes[x]);
    	dbg_min(nodes[x]);
    }


    int numRoots(){
    	assert(dbgSetCount());
    	return setCount;
    }

    void reset(){
    	for(int i = 0;i<nodes.size();i++){
    		/*nodes[i]->parent=NULL;
    		nodes[i]->left=NULL;
    		nodes[i]->right =NULL;*/
    		assert(nodes[i]->id==i);
    		(*nodes[i]).~Node();
    		new (nodes[i]) Node(i);
    	}
    	setCount= nodes.size();
    }

};
#endif
