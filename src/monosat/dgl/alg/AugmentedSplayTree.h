/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

//Implementation modified from wikipedia's
//This splay tree is augmented to keep track of the sizes of subtrees
#ifndef SPLAY_TREE
#define SPLAY_TREE

#include <functional>
#include "SearchTree.h"

#define linked

template<typename T>
class AugmentedSplayTree;

template<typename T>
struct _node {
    friend class AugmentedSplayTree<_node<T>>;

    _node<T>* left, * right;
    _node<T>* parent;
//link to the successor and predessor nodes, for quick traversal

#ifdef linked
    _node<T>* _next;
    _node<T>* _prev;

    void setPrev(_node<T>* v){
        _prev = v;
    }

    void setNext(_node<T>* v){
        _next = v;
    }

#else
    void setPrev(_node<T> * v) {

    }

    void setNext(_node<T> * v) {

    }
#endif

    int subtree_size;
public:

    T value;

    _node(const T& init = T()) :
            left(0), right(0), parent(0), subtree_size(1), value(init){
#ifdef linked
        _next = 0;
        _prev = 0;
#endif
    }

public:
    _node* prev(){
#ifdef linked
#ifdef DEBUG_DGL
        _node*exp = prev_slow();
#endif
        assert(_prev == prev_slow());
        return _prev;
#else
        return prev_slow();
#endif
    }

    _node* next(){
#ifdef linked
#ifdef DEBUG_DGL
        _node*exp = next_slow();
#endif
        assert(_next == next_slow());
        return _next;
#else
        return next_slow();
#endif
    }

private:
//it should be possible to add prev and next links to this node and maintain them in O(1) time to ensure these operations are also O(1) instead of O(log(n)), as they currently are.
    _node* prev_slow(){
        if(left)
            return left->findMax();
        else{
            _node* n = this;
            _node* p = n->parent;
            while(p && n == p->left){
                n = p;
                p = p->parent;
            }
            return p;
        }
    }

    _node* next_slow(){
        if(right)
            return right->findMin();
        else{
            _node* n = this;
            _node* p = n->parent;
            while(p && n == p->right){
                n = p;
                p = p->parent;
            }
            return p;
        }

    }

public:
    _node* findRoot(){
        _node* u = this;
        while(u->parent)
            u = u->parent;
        return u;
    }

    _node* findMin(){
        _node* u = this;
        while(u->left)
            u = u->left;
        return u;
    }

    _node* findMax(){
        _node* u = this;
        while(u->right)
            u = u->right;
        return u;
    }

};

template<typename T>
class AugmentedSplayTree : public AugmentedSearchTree<_node<T>> {
public:
    //Comp comp;
    unsigned long long p_size;
    typedef _node<T> Node;
private:
    //node * root;
    /**
     *
     * http://www.cs.cmu.edu/~avrim/451f12/lectures/lect0927.txt:
     *
     * Let's denote by n.l and n.r the left and right children
     of a node n.  We can use the following fact to help us
     keep the sizes:

     n.size = n.l.size + n.r.size + 1

     (The null nodes have size 0.)

     When we do a rotation, the only nodes whose sizes change
     are the two involved in the rotation.

     y             x
     right rotation:          /     ====>     \
                           x                 y

     To update the sizes after the rotation, we first fix y
     (by applying the above formula) and then we fix x.

     This allows us to splay while maintaining the size fields.
     */
    void left_rotate(Node* x){
        Node* y = x->right;

        // int x_incident = getIncident(x);
        // int y_incident = getIncident(y);

        x->right = y->left;
        if(y->left)
            y->left->parent = x;
        y->parent = x->parent;
        if(!x->parent){
            //root = y;
        }else if(x == x->parent->left)
            x->parent->left = y;
        else
            x->parent->right = y;
        y->left = x;
        x->parent = y;

        //update x's subtree size (have to do this first, because x is a child of y)
        x->subtree_size = 1;
        // x->n_incident=x_incident;
        if(x->left){
            x->subtree_size += x->left->subtree_size;
            //x->n_incident+=x->left->n_incident;
        }
        if(x->right){
            x->subtree_size += x->right->subtree_size;
            //x->n_incident+=x->right->n_incident;
        }

        //update y's subtree size
        y->subtree_size = 1;
        //y->n_incident=y_incident;
        if(y->left){
            y->subtree_size += y->left->subtree_size;
            //y->n_incident+=y->left->n_incident;
        }
        if(y->right){
            y->subtree_size += y->right->subtree_size;
            //y->n_incident+=y->right->n_incident;
        }
        dbg_checkSubtreeSize(y);
    }

    void right_rotate(Node* x){
        Node* y = x->left;

        //int x_incident = getIncident(x);
        // int y_incident = getIncident(y);

        x->left = y->right;
        if(y->right)
            y->right->parent = x;
        y->parent = x->parent;
        if(!x->parent){
            //root = y;
        }else if(x == x->parent->left)
            x->parent->left = y;
        else
            x->parent->right = y;
        y->right = x;
        x->parent = y;

        //update x's subtree size (have to do this first, because x is a child of y)
        x->subtree_size = 1;
        //x->n_incident=x_incident;
        if(x->left){
            x->subtree_size += x->left->subtree_size;
            //x->n_incident+=x->left->n_incident;
        }
        if(x->right){
            x->subtree_size += x->right->subtree_size;
            //x->n_incident+=x->right->n_incident;
        }

        //update y's subtree size
        y->subtree_size = 1;
        // y->n_incident=y_incident;
        if(y->left){
            y->subtree_size += y->left->subtree_size;
            //y->n_incident+=y->left->n_incident;
        }
        if(y->right){
            y->subtree_size += y->right->subtree_size;
            //y->n_incident+=y->right->n_incident;
        }
        dbg_checkSubtreeSize(y);
    }

public:
    int dbg_checkSubtreeSize(Node* x){
#ifdef DEBUG_DGL
        if (!x)
            return 0;

        int size = 1 + dbg_checkSubtreeSize(x->left) + dbg_checkSubtreeSize(x->right);
        assert(x->subtree_size == size);
        return size;
#endif
        return 0;
    }

    void dbg_checkNextSucc(Node* x){
#ifdef DEBUG_DGL
        if (!x)
            return;
        x->next();
        x->prev();
        if (x->next())
            x->next()->prev();
        if (x->prev())
            x->prev()->next();

#endif

    }

    void splay(Node* x){
        while(x->parent){
            if(!x->parent->parent){
                if(x->parent->left == x)
                    right_rotate(x->parent);
                else
                    left_rotate(x->parent);
            }else if(x->parent->left == x && x->parent->parent->left == x->parent){
                right_rotate(x->parent->parent);
                right_rotate(x->parent);
            }else if(x->parent->right == x && x->parent->parent->right == x->parent){
                left_rotate(x->parent->parent);
                left_rotate(x->parent);
            }else if(x->parent->left == x && x->parent->parent->right == x->parent){
                right_rotate(x->parent);
                left_rotate(x->parent);
            }else{
                left_rotate(x->parent);
                right_rotate(x->parent);
            }
        }
        assert(findRoot(x) == x);
        dbg_checkSubtreeSize(x);
    }

    void replace(Node* u, Node* v){
        if(!u->parent){
            //root = v;
        }else if(u == u->parent->left)
            u->parent->left = v;
        else
            u->parent->right = v;
        if(v){
#ifdef linked
            assert(!v->prev());
            assert(!v->parent);

            if(u->prev()){
                v->setPrev(u->prev());
                assert(u->prev()->next() == u);
                u->prev()->setNext(v);
            }
            if(u->next()){
                v->setNext(u->next());
                assert(u->next()->prev() == u);
                u->next()->setPrev(v);
            }
#endif
            v->parent = u->parent;
        }
        dbg_checkSubtreeSize(u);
        dbg_checkSubtreeSize(v);
        dbg_checkNextSucc(v);
        dbg_checkNextSucc(u);
    }

    Node* subtree_minimum(Node* u){
        while(u->left)
            u = u->left;
        return u;
    }

    Node* subtree_maximum(Node* u){
        while(u->right)
            u = u->right;
        return u;
    }

public:
    AugmentedSplayTree() :
            p_size(0){
    }

    Node* createNode(const T& value){
        return new Node(value);
    }

    void insertAfter(Node* insertAt, Node* toInsert){
        Node* z = insertAt;
        Node* p = nullptr;

        while(z){
            p = z;
            z = z->right;
            //else z = z->left;
        }

        z = toInsert;
        z->parent = p;

        if(!p){
            //root = z;
        }else{
            assert(!p->right);
#ifdef linked

            assert(!p->next());
            Node* m = z->findMin();
            p->setNext(m);
            m->setPrev(p);
#endif
            p->right = z;
            p->subtree_size += z->subtree_size;
        }
        splay(z);   //why?
        p_size++;
        dbg_checkSubtreeSize(z);
        dbg_checkNextSucc(z);

    }

    void insertBefore(Node* insertAt, Node* toInsert){
        Node* z = insertAt;
        Node* p = nullptr;

        while(z){
            p = z;
            z = z->left;
            //else z = z->left;
        }

        z = toInsert;
        z->parent = p;

        if(!p){
            //root = z;
        }else{

#ifdef linked

            assert(!p->prev());
            Node* m = z->findMax();
            p->setPrev(m);
            m->setNext(p);
#endif
            p->left = z;
            p->subtree_size += z->subtree_size;
        }
        splay(z);
        p_size++;
        dbg_checkSubtreeSize(z);
        dbg_checkNextSucc(z);
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
    //This is LOG time
    Node* splitAfter(Node* splitAt){
        splay(splitAt);
        assert(!splitAt->parent);
        assert(findRoot(splitAt) == splitAt);
        Node* r = splitAt->right;
        if(r){
#ifdef linked
            assert(splitAt->next());
            splitAt->next()->setPrev(nullptr);
            splitAt->setNext(nullptr);
#endif
            r->parent = nullptr;
        }
        splitAt->right = nullptr;

        assert(findMax(splitAt) == splitAt);
        assert(findMax(findRoot(splitAt)) == splitAt);

        if(r){
            assert(findRoot(r) == r);

            assert(findRoot(splitAt) != findRoot(r));

        }
        if(r){
            Node* p = splitAt;
            while(p){
                p->subtree_size -= r->subtree_size;
                dbg_checkSubtreeSize(p);
                p = p->parent;
            }
        }
        dbg_checkNextSucc(splitAt);
        dbg_checkNextSucc(r);
        return r;
    }

    //This is LOG time
    Node* splitBefore(Node* splitAt){
        splay(splitAt);
        assert(!splitAt->parent);

        Node* l = splitAt->left;

        if(l){
#ifdef linked
            assert(splitAt->prev());
            splitAt->prev()->setNext(nullptr);
            splitAt->setPrev(nullptr);
#endif
            l->parent = nullptr;
        }

        splitAt->left = nullptr;

        if(l){
            Node* p = splitAt;
            while(p){
                p->subtree_size -= l->subtree_size;
                dbg_checkSubtreeSize(p);
                p = p->parent;
            }
        }
        dbg_checkNextSucc(splitAt);
        dbg_checkNextSucc(l);
        return l;
    }

    //This is LOG time
    Node* concat(Node* left, Node* right){
        left = findRoot(left);
        right = findRoot(right);
        assert(findRoot(left) == left);
        assert(findRoot(right) == right);

        Node* maxLeft = findMax(left);
        //only needed if we are maintaining next/prev pointers
#ifdef linked
        Node* minRight = findMin(right);
        assert(!maxLeft->next());
        assert(!minRight->prev());
        maxLeft->setNext(minRight);
        minRight->setPrev(maxLeft);
#endif
        splay(maxLeft);
        assert(!maxLeft->parent);
        assert(!maxLeft->right);
        maxLeft->right = right;
        right->parent = maxLeft;
        maxLeft->subtree_size += right->subtree_size;
        dbg_checkSubtreeSize(maxLeft);
        dbg_checkNextSucc(maxLeft);
        dbg_checkNextSucc(right);
        dbg_checkNextSucc(left);
        return left;
    }
    /*  //This is LOG time
     void remove( Node * toRemove ) {
     Node *z =toRemove;
     if( !z ) return;

     splay( z );
     //need to maintain subtree sizes here...
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
     dbg_checkSubtreeSize(y);
     }

     delete z;
     p_size--;
     }*/
    //This is LOG time
    Node* findRoot(Node* of){

        while(of && of->parent){
            of = of->parent;
        }
        return of;
    }

    //This is LOG time
    Node* findMin(Node* of){
        return subtree_minimum(of);
    }

    //This is LOG time
    Node* findMax(Node* of){
        return subtree_maximum(of);
    }

    //This is LOG time
    int depth(Node* a){
        int depth = 0;
        while(a->parent){
            depth++;
            a = a->parent;
        }
        assert(depth >= 0);
        return depth;
    }

    //Return -1 if a is before b, 1 if a is after b, 0 if either they are the same node, or not in the same tree.
    //This calculation is LOG time
    int compare(Node* a, Node* b){
        if(a == b){
            return 0;
        }
        assert(findRoot(a) == findRoot(b));
        //can these be avoided? well, yes, if you are willing to leave markers behind on the visited nodes...
        int depthA = depth(a);
        int depthB = depth(b);

        while(depthA > depthB){
            if(a->parent == b){
                if(a == b->left){
                    return -1;
                }else{
                    assert(a == b->right);
                    return 1;
                }
            }
            a = a->parent;
            depthA--;

        }
        while(depthB > depthA){
            if(b->parent == a){
                if(b == a->left){
                    return 1;
                }else{
                    assert(b == a->right);
                    return -1;
                }
            }
            b = b->parent;
            depthB--;
        }
        assert(depthA == depthB);
        int depth = depthA;
        while(depth){
            assert(a != b);
            if(a->parent == b->parent){
                if(a == a->parent->left){
                    assert(b == a->parent->right);
                    return -1;
                }else{
                    assert(b == a->parent->left);
                    assert(a == a->parent->right);
                    return 1;
                }
            }
            a = a->parent;
            b = b->parent;
            depth--;
        }
        assert(a && b);
        assert(!a->parent);
        assert(!b->parent);
        return 0;
    }

    void clear(Node* n){
        n->parent = nullptr;
        n->left = nullptr;
        n->right = nullptr;
        n->subtree_size = 1;
        n->setNext(nullptr);
        n->setPrev(nullptr);
    }

    //This is CONSTANT time
    int size(Node* of){
        return of->subtree_size;
    }
    /*	void incrementIncident(Node * x){
     x->n_incident++;
     }
     void decrementIncident(Node * x){
     x->n_incident--;
     assert(x->n_incident>=0);
     assert(getLocalIncident(x)>=0);
     }*/

    //This is LOG time
    /*	void addToIncident(Node * x, int add){
     while(x){
     x->n_incident+=add;
     assert(x->n_incident>=0);
     assert(getIncident(x)>=0);
     x=x->parent;
     }
     }

     //This is LOG time
     void setIncident(Node * x, int n_incident){
     int n = getIncident(x);
     int diff = n_incident-n;
     addToIncident(x,diff);
     assert(getIncident(x)==n_incident);
     }

     //This is CONSTANT time
     int getIncident(Node * x){
     int local_incident= x->n_incident;
     if(x->left)
     local_incident-=x->left->n_incident;
     if(x->right)
     local_incident-=x->right->n_incident;
     assert(local_incident>=0);
     return local_incident;
     }

     //This is CONSTANT time
     int getSubtreeIncident(Node * of){
     return of->n_incident;
     }
     */
    //bool empty( ) const { return root == 0; }
    unsigned long long size() const{
        return p_size;
    }
};

#endif // SPLAY_TREE
