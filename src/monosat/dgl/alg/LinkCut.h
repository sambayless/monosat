/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless
 Copyright (c) 2011, Daniel Sleator

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

#ifndef LINK_CUT
#define LINK_CUT

#include <cstddef>
#include <cassert>
#include <vector>

//Link-Cut Tree,
//This is a C++ adaptation of Danny Sleator's (public domain) java-implementation (see http://codeforces.com/contest/117/submission/860934), used with permission from the author.

//Watch out - LinkCut cannot by itself implement a dynamic disjoint-set data structure - you need more for that.
class LinkCut {
    int setCount;

    struct Node {
        int id;
        Node* left;
        Node* right;
        Node* parent;

        Node(int _id) :
                id(_id), left(nullptr), right(nullptr), parent(nullptr){
        };
    };

    std::vector<Node*> nodes;

    // Whether x is a root of a splay tree
    bool isRoot(Node* x){
        return x->parent == nullptr || (x->parent->left != x && x->parent->right != x);
    }

    void connect(Node* ch, Node* p, bool leftChild){
        if(leftChild)
            p->left = ch;
        else
            p->right = ch;
        if(ch != nullptr)
            ch->parent = p;
    }

    void rotR(Node* p){
        Node* q = p->parent;
        Node* r = q->parent;

        if((q->left = p->right) != nullptr)
            q->left->parent = q;
        p->right = q;
        q->parent = p;
        if((p->parent = r) != nullptr){
            if(r->left == q)
                r->left = p;
            else if(r->right == q)
                r->right = p;
        }

    }

    void rotL(Node* p){
        Node* q = p->parent;
        Node* r = q->parent;

        if((q->right = p->left) != nullptr)
            q->right->parent = q;
        p->left = q;
        q->parent = p;
        if((p->parent = r) != nullptr){
            if(r->left == q)
                r->left = p;
            else if(r->right == q)
                r->right = p;
        }

    }

    void splay(Node* p){
        while(!isRoot(p)){
            Node* q = p->parent;
            if(isRoot(q)){
                if(q->left == p)
                    rotR(p);
                else
                    rotL(p);
            }else{
                Node* r = q->parent;
                if(r->left == q){
                    if(q->left == p){
                        rotR(q);
                        rotR(p);
                    }else{
                        rotL(p);
                        rotR(p);
                    }
                }else{
                    if(q->right == p){
                        rotL(q);
                        rotL(p);
                    }else{
                        rotR(p);
                        rotL(p);
                    }
                }
            }
        }
    }

    // Makes node x the root of the virtual tree, and also x is the leftmost
    // node in its splay tree
    Node* expose(Node* x){
        Node* last = nullptr;
        for(Node* y = x; y != nullptr; y = y->parent){
            splay(y);
            y->left = last;
            last = y;
        }
        splay(x);
        return last;
    }

    Node* _findRoot(Node* x){
        expose(x);
        while(x->right != nullptr){
            x = x->right;
        }
        //splay(x);
        return x;
    }

    bool dbgSetCount(){
        int count = 0;
        for(int i = 0; i < nodes.size(); i++){
            Node* n = nodes[i];
            Node* r = _findRoot(n);
            if(r == n){
                count++;
            }
        }
        return count == setCount;
    }

    // prerequisite: x and y are in distinct trees
    void _link(Node* x, Node* y){
        //assert (_findRoot(x) != _findRoot(y));
#ifdef DEBUG_DGL
        Node* sY = _findRoot(y);
        Node* sX = _findRoot(x);
        assert(sY != sX);  //else this is a bug
#endif

        setCount--;
        expose(x);
        assert(!x->parent);
        x->parent = y;
        assert(dbgSetCount());
    }

    bool _connected(Node* x, Node* y){
        if(x == y)
            return true;
        expose(x);
        expose(y);
        return x->parent != nullptr;
    }

    void _cut(Node* x, Node* y){
        expose(x);
        expose(y);
        if(x->parent != nullptr){
            setCount++;
        }
        assert(!(y->right != x || x->left != nullptr || x->right != nullptr));

        y->right->parent = nullptr;
        y->right = nullptr;
        assert(dbgSetCount());
    }

public:
    LinkCut() :
            setCount(0){

    }

    int addNode(){
        //return new Node();
        setCount++;
        nodes.push_back(new Node(nodes.size()));
        return nodes.size() - 1;
    }

    int nNodes(){
        return nodes.size();
    }

    int findRoot(int x){
        return _findRoot(nodes[x])->id;
    }

    // prerequisite: x and y are in distinct trees
    // and that p is a root of its tree, this links p to q
    void link(int x, int y){
        if(x == y)
            return;
        Node* xnode = nodes[x];
        Node* ynode = nodes[y];
        _link(xnode, ynode);
    }

    bool connected(int x, int y){
        if(x == y)
            return true;
        Node* xnode = nodes[x];
        Node* ynode = nodes[y];
        expose(xnode);
        expose(ynode);
#ifdef DEBUG_DGL
        int s1 = findRoot(x);
        int s2 = findRoot(y);
        bool dbg_connected = s1 == s2;
        if (dbg_connected) {
            assert(xnode->parent);
        } else {
            assert(xnode->parent==NULL);
        }

#endif
        return xnode->parent != nullptr;

    }

    void cut(int x, int y){
        Node* xnode = nodes[x];
        Node* ynode = nodes[y];
        _cut(xnode, ynode);
    }

    int numRoots(){
        assert(dbgSetCount());
        return setCount;
    }

    void reset(){
        for(int i = 0; i < nodes.size(); i++){
            nodes[i]->parent = nullptr;
            nodes[i]->left = nullptr;
            nodes[i]->right = nullptr;
        }
        setCount = nodes.size();
    }

};

#endif
