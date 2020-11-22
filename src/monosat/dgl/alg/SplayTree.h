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
//Implementation modified from wikipedia
#ifndef SPLAY_TREE
#define SPLAY_TREE

#include <functional>
#include "SearchTree.h"

#define linked

template<typename T>
class SplayTree : public SearchTree<_node < T>>


template<typename T>
struct _node {
    friend SplayTree<_node<T>>;
    _node<T>* left, * right;
    _node<T>* parent;
//link to the successor and predessor nodes, for quick traversal
private:
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
public:
    T value;

    _node(const T& init = T()) :
            left(0), right(0), parent(0), _next(0), _prev(0), value(init){
    }

public:
    _node* prev(){
#ifdef linked
        assert(_prev == prev_slow());
        return _prev;
#else
        return prev_slow();
#endif
    }

    _node* next(){
#ifdef linked
        assert(_next == next_slow());
        return _next;
#else
        return next_slow();
#endif
    }

private:
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
class SplayTree : public SearchTree<_node<T>> {
public:
    //Comp comp;
    uint64_t p_size;
    typedef _node<T> Node;
private:
    //node * root;
    void left_rotate(Node* x){
        Node* y = x->right;
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
    }

    void right_rotate(Node* x){
        Node* y = x->left;
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
    SplayTree() :
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
        }else
            p->right = z;

        splay(z);
        p_size++;
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
        }else
            p->left = z;

        splay(z);
        p_size++;
    }

    //Splits and returns the node to the right of splitAt; splitAt remains connected to the nodes left of it.
    Node* splitAfter(Node* splitAt){
        splay(splitAt);
        assert(!splitAt->parent);

        Node* r = splitAt->right;

        splitAt->right = nullptr;
        if(r){
#ifdef linked
            assert(r->prev());
            r->prev()->setNext(nullptr);
            r->setPrev(nullptr);
#endif
            r->parent = nullptr;
        }
        return r;
    }

    Node* splitBefore(Node* splitAt){
        splay(splitAt);
        assert(!splitAt->parent);

        Node* l = splitAt->left;
        splitAt->left = nullptr;
        if(l){
#ifdef linked
            assert(l->next());
            l->next()->setPrev(nullptr);
            l->setNext(nullptr);
#endif
            l->parent = nullptr;
        }
        return l;
    }

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

        return left;
    }

    void remove(Node* toRemove){
        Node* z = toRemove;
        if(!z)
            return;

        splay(z);

        if(!z->left)
            replace(z, z->right);
        else if(!z->right)
            replace(z, z->left);
        else{
            Node* y = subtree_minimum(z->right);
            if(y->parent != z){
                replace(y, y->right);
                y->right = z->right;
                y->right->parent = y;
            }
            replace(z, y);
            y->left = z->left;
            y->left->parent = y;
        }

        delete z;
        p_size--;
    }

    Node* findRoot(Node* of){
        while(of->parent){
            of = of->parent;
        }
        return of;
    }

    Node* findMin(Node* of){
        return subtree_minimum(of);
    }

    Node* findMax(Node* of){
        return subtree_maximum(of);
    }

    //bool empty( ) const { return root == 0; }
    uint64_t size() const{
        return p_size;
    }
};

#endif // SPLAY_TREE
