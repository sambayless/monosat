//Implementation modified from wikipedia

/**
 *From http://www.eli.sdsu.edu/courses/fall95/cs660/notes/splay/Splay.html
 *  Splay Operations

access(i, t): if i is in tree t return pointer to i, otherwise return null pointer

Find i, then splay tree t at i.
If i is not in tree t, then splay last node accessed looking for i



join (a, b): Return tree formed by combining tree "a", and tree "b". Assumes that every item in "a" has key less then every item in "b"

Splay largest item in "a", then add "b" as a right child of root of "a"



split (i, t): Split tree t, containing item i, into two trees: "a", containing all items with key less or equal to "i"; and "b", containing all items with key greater than "i"

Perform access(i, t) then split tree at root


insert(i, t): insert i in tree t

Perform split (i, t) then make i the root of the two trees returned by split


delete(i, t): delete i from tree t

Perform access(i, t) then perform join on t's subtrees
 */

#ifndef SPLAY_TREE
#define SPLAY_TREE
#include <functional>
#include "SearchTree.h"

template< typename T>
struct _node {
_node<T> *left, *right;
_node<T> *parent;
T value;
_node( const T& init = T( ) ) : left( 0 ), right( 0 ), parent( 0 ), value( init ) { }

_node * prev(){
	if(left)
		return left->findMax();
	else {
		  _node*n = this;
		  _node *p = n->parent;
		  while(p  && n == p->left)
		  {
		     n = p;
		     p = p->parent;
		  }
		  return p;
	}
}
_node * next(){
	if(right)
		return right->findMin();
	else {
		  _node*n = this;
		  _node *p = n->parent;
		  while(p  && n == p->right)
		  {
		     n = p;
		     p = p->parent;
		  }
		  return p;
	}

}

private:

_node* findMin(  ) {
  _node*u=this;
  while( u->left ) u = u->left;
  return u;
}

_node* findMax(  ) {
  _node*u=this;
  while( u->right ) u = u->right;
  return u;
}

} ;

template< typename T >
class SplayTree:public SearchTree<_node<T>> {
public:
  //Comp comp;
  unsigned long p_size;
  typedef _node<T> Node;
private:
  //node * root;
  void left_rotate( Node *x ) {
    Node *y = x->right;
    x->right = y->left;
    if( y->left ) y->left->parent = x;
    y->parent = x->parent;
    if( !x->parent ){
    	//root = y;
    }else if( x == x->parent->left ) x->parent->left = y;
    else x->parent->right = y;
    y->left = x;
    x->parent = y;
  }

  void right_rotate( Node *x ) {
    Node *y = x->left;
    x->left = y->right;
    if( y->right ) y->right->parent = x;
    y->parent = x->parent;
    if( !x->parent ){
    	//root = y;
    }else if( x == x->parent->left ) x->parent->left = y;
    else x->parent->right = y;
    y->right = x;
    x->parent = y;
  }

  void splay( Node *x ) {
    while( x->parent ) {
      if( !x->parent->parent ) {
        if( x->parent->left == x ) right_rotate( x->parent );
        else left_rotate( x->parent );
      } else if( x->parent->left == x && x->parent->parent->left == x->parent ) {
        right_rotate( x->parent->parent );
        right_rotate( x->parent );
      } else if( x->parent->right == x && x->parent->parent->right == x->parent ) {
        left_rotate( x->parent->parent );
        left_rotate( x->parent );
      } else if( x->parent->left == x && x->parent->parent->right == x->parent ) {
        right_rotate( x->parent );
        left_rotate( x->parent );
      } else {
        left_rotate( x->parent );
        right_rotate( x->parent );
      }
    }
  }

  void replace( Node *u, Node *v ) {
    if( !u->parent ){
    	//root = v;
    }else if( u == u->parent->left ) u->parent->left = v;
    else u->parent->right = v;
    if( v ) v->parent = u->parent;
  }

  Node* subtree_minimum( Node *u ) {
    while( u->left ) u = u->left;
    return u;
  }

  Node* subtree_maximum( Node *u ) {
    while( u->right ) u = u->right;
    return u;
  }
public:
  SplayTree( ) : p_size( 0 ) { }

  Node * createNode(const T &value){
	  return new Node(value);
  }

  void insertAfter(Node * insertAt, Node * toInsert ) {
    Node *z = insertAt;
    Node *p = nullptr;

    while( z ) {
      p = z;
      z = z->right;
      //else z = z->left;
    }

    z = toInsert;
    z->parent = p;

    if( !p ){
    	//root = z;
    }else
    	p->right = z;

    splay( z );
    p_size++;
  }

  void insertBefore(Node * insertAt, Node * toInsert ) {
    Node *z = insertAt;
    Node *p = nullptr;

    while( z ) {
      p = z;
      z = z->left;
      //else z = z->left;
    }

    z = toInsert;
    z->parent = p;

    if( !p ){
    	//root = z;
    }else
    	p->left = z;

    splay( z );
    p_size++;
  }

/*  node* find( const T &key ) {
    node *z = root;
    while( z ) {
      if( comp( z->key, key ) ) z = z->right;
      else if( comp( key, z->key ) ) z = z->left;
      else return z;
    }
    return 0;
  }*/

  //Splits and returns the node to the right of splitAt; splitAt remains connected to the nodes left of it.
  Node*split(Node * splitAt){
	  splay(splitAt);
	  assert(!splitAt->parent);

	  Node * r = splitAt->right;
	  splitAt->right = nullptr;
	  if(r){
		  r->parent=nullptr;
	  }
	  return r;
  }

  Node *concat(Node * left, Node*right){
	   Node* maxLeft = findMax(left);
	   splay(maxLeft);
	   assert(!maxLeft->parent);
	   assert(!maxLeft->right);
	   maxLeft->right=right;
	   right->parent=maxLeft;
	   return left;
  }

  void remove( Node * toRemove ) {
    Node *z =toRemove;
    if( !z ) return;

    splay( z );

    if( !z->left ) replace( z, z->right );
    else if( !z->right ) replace( z, z->left );
    else {
      Node *y = subtree_minimum( z->right );
      if( y->parent != z ) {
        replace( y, y->right );
        y->right = z->right;
        y->right->parent = y;
      }
      replace( z, y );
      y->left = z->left;
      y->left->parent = y;
    }

    delete z;
    p_size--;
  }
  Node* findRoot(Node* of ) {
	  while(of->parent){
		  of=of->parent;
	  }
	  return of;
  }
  Node* findMin(Node* of ) { return subtree_minimum( of ); }
  Node* findMax(Node* of ) { return subtree_maximum( of ); }

  //bool empty( ) const { return root == 0; }
  unsigned long size( ) const { return p_size; }
};

#endif // SPLAY_TREE
