/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless
 Copyright (c) 2013 Mikola Lysenko

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

//adapted from https://github.com/mikolalysenko/dynamic-forest/blob/master/lib/euler.js
//untested at this point.
#ifndef EULER_TREE_H
#define EULER_TREE_H

#include "TreapCustom.h"
#include <vector>
#include <cstdio>
#include <algorithm>
#include "AugmentedSplayTree.h"
#include "SearchTree.h"

class EulerTree {

private:
    struct EulerHalfEdge;
    struct EulerVertex;

    typedef AugmentedSplayTree<EulerHalfEdge*> Tree;
    int nComponents;

    Tree t;
    //static EulerVertex * root;
    std::vector<EulerVertex*> vertices;
    std::vector<EulerHalfEdge*> forward_edges;
    std::vector<EulerHalfEdge*> backward_edges;

    struct EulerHalfEdge {
        //Value value;
        //NodeData d;
        int value;
        int index;
        bool isforward;        //forward edges go from a parent to a child

        EulerVertex* from;
        EulerVertex* to;
        Tree::Node* node;        //node in the treap

#ifdef DEBUG_DGL
        int rank;//this is the position of this edge in the euler tour. This information is only maintained implicitly as the edges position in the underlying binary search tree.
#endif

        bool contains(EulerVertex* v){
            return from == v || to == v;
        }

        bool contains(int v){
            return from->index == v || to->index == v;
        }
    };

    struct EulerVertex {
        //Value value;

        //Treap::Node * first;
        //Treap::Node * last;
        EulerHalfEdge* left_out;//the half edge leading into the vertex, from its parent. The 'to' field in this halfedge the first occurrence of the euler vertex in the tour.
        //IFF the vertex is root, then this is instead the half edge leading to its first child. In that case, the 'from' field is the first occurence of the vertex in the tour
        EulerHalfEdge* right_in;//the half edge returning to the vertex.  This half edge has the last occurrence of the euler vertex in the tour.
        int index;
        bool visited;

#ifdef DEBUG_DGL
        EulerTree * owner;
        //Note: in an euler-tour tree representation, these values are not explicitly maintained
        EulerVertex * dbg_parent;
        std::vector<EulerVertex*> dbg_children;
        bool dbg_visited;
        std::vector<EulerVertex*> dbg_children_t;
#endif

        EulerVertex() :
                left_out(nullptr), right_in(nullptr), visited(false){

            index = 0;
#ifdef DEBUG_DGL
            dbg_visited = false;
            owner = nullptr;
            dbg_parent = nullptr;
#endif
        }

        void clear(){
            left_out = nullptr;
            right_in = nullptr;
            visited = false;
#ifdef DEBUG_DGL
            dbg_visited = false;
            owner = nullptr;
            dbg_parent = nullptr;
            dbg_children.clear();
            dbg_visited = false;
#endif
        }

        bool isSingleton(){
            assert(!left_out == !right_in);
            return !left_out;
        }

        //This is an _ARBITRARY_ incident edge.
        Tree::Node* incidentEdgeA(){
            if(left_out){
                return left_out->node;
            }else{
                return nullptr;
            }
        }

        //This is an _ARBITRARY_ incident edge.
        Tree::Node* incidentEdgeB(){

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

        void setIncidentEdgeA(EulerHalfEdge* n){
            assert(!n || (n->from == this || n->to == this));
            left_out = n;
        }

        void setIncidentEdgeB(EulerHalfEdge* n){
            assert(!n || (n->from == this || n->to == this));
            right_in = n;
        }

        /*		//Get the next node in the tour.
         EulerVertex * getNext(){
         if (left_out)
         return left_out->to;
         else if(incidentEdgeB())
         return incidentEdgeB()->value->to;
         else
         return nullptr;
         }*/

        //
        void dbg_real_tour(std::vector<int>& tour_list){
#ifdef DEBUG_DGL
            tour_list.clear();

            if (isSingleton()) {
                tour_list.push_back(index);
                return;
            }

            Tree::Node* n = owner->t.findMin(owner->t.findRoot(incidentEdgeA()));

            while (n->next()) {
#ifdef dbg_print
                printf("(%d -> %d)",n->value->from->index,n->value->to->index);
#endif
                tour_list.push_back(n->value->from->index);
                assert(n->next());

                //printf("(%d,%d)\n",n->value->from->index,n->value->to->index);
                assert(n->next());
                n = n->next();
                assert(n);
            }
#ifdef dbg_print
            printf("(%d -> %d)\n",n->value->from->index,n->value->to->index);
#endif
            tour_list.push_back(n->value->from->index);
            tour_list.push_back(n->value->to->index);
#endif
        }

        int dbg_getSize(){
#ifdef DEBUG_DGL
            int s = 1;
            for (EulerVertex* c : dbg_children) {
                s += c->dbg_getSize();
            }
            //assert(s==subtree_size);
            return s;
#endif
            return 0;
        }

        void dbg_make_parent(EulerVertex* new_parent){
#ifdef DEBUG_DGL
            if (new_parent) {
                assert(std::count(dbg_children.begin(), dbg_children.end(), new_parent));

                remove(dbg_children.begin(), dbg_children.end(), new_parent);
                EulerVertex * p = dbg_parent;
                dbg_parent = new_parent;
                if (p) {
                    dbg_children.push_back(p);

                    p->dbg_make_parent(this);

                }
            }
#endif
        }

        void dbg_make_root(){
#ifdef DEBUG_DGL
            if (dbg_parent) {
                dbg_children.push_back(dbg_parent);
                dbg_parent->dbg_make_parent(this);

                dbg_parent = nullptr;
            }
#endif
        }

        void dbg_remove(){
#ifdef DEBUG_DGL

            if (dbg_parent) {
                std::remove(dbg_parent->dbg_children.begin(), dbg_parent->dbg_children.end(), this);
                //dbg_parent->dbg_children.remove(this);
                dbg_parent = nullptr;
            } else {
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

        void dbg_insert(EulerVertex* node){
#ifdef DEBUG_DGL
            assert(!node->dbg_parent);
            dbg_children.push_back(node);
            node->dbg_parent = this;

#endif
        }

        void dbg_build_tour_helper(EulerVertex* r, std::vector<int>& tour){
#ifdef DEBUG_DGL
            tour.push_back(r->index);
            for (EulerVertex * c : r->dbg_children) {
                dbg_build_tour_helper(c, tour);
                tour.push_back(r->index);

            }
#endif
        }

        void dbg_build_tour(EulerVertex* r, std::vector<int>& tour){
#ifdef DEBUG_DGL
            while (r->dbg_parent)
                r = r->dbg_parent;
            dbg_build_tour_helper(r, tour);
#endif
        }

        void dbg_clear(){
#ifdef DEBUG_DGL
            assert(dbg_visited);
            dbg_visited = false;
            assert(dbg_children.size() == dbg_children_t.size());
            dbg_children.clear();
            dbg_children = dbg_children_t;
            dbg_children_t.clear();
            for (EulerVertex * v : dbg_children)
                v->dbg_clear();
#endif
        }

        void dbg_tour(){
#ifdef DEBUG_DGL
            if (!left_out) {
                assert(!right_in);
                return;
            }

            EulerHalfEdge * from = owner->t.findMin(owner->t.findRoot(incidentEdgeA()))->value;

            EulerVertex * dbgFrom = this;
            while (dbgFrom->dbg_parent)
                dbgFrom = dbgFrom->dbg_parent;

            assert(from->contains(dbgFrom));
            std::vector<EulerVertex*> dbg_stack;
            dbgFrom->dbg_visited = true;
            dbg_stack.push_back(dbgFrom);
            while (from && dbg_stack.size()) {

                EulerVertex* p = dbg_stack.back();
                bool found = false;
                assert(p->incidentEdgeA());
                assert(p->incidentEdgeB());
                for (EulerVertex* t : p->dbg_children) {
                    if (!t->dbg_visited && from->contains(t)) {
                        t->dbg_visited = true;
                        found = true;
                        dbg_stack.push_back(t);
                        p->dbg_children_t.push_back(t);
                        from = from->node->next()->value;
                        break;
                    }
                }
                if (!found) {
                    dbg_stack.pop_back();
                    if (dbg_stack.size())
                        assert(from->contains(p->dbg_parent));
                    else {
                        assert(!p->dbg_parent);
                    }
                    if (from->node->next())
                        from = from->node->next()->value;
                    else
                        from = nullptr;
                }
            }
            assert(!from);
            assert(dbg_stack.size() == 1);
            dbgFrom->dbg_clear();

            //build the tour
            std::vector<int> dbg_tour;
            dbg_build_tour(this, dbg_tour);
            std::vector<int> real_tour;
            dbg_real_tour(real_tour);
            assert(dbg_tour.size() == real_tour.size());
            /*for(int i = 0;i<real_tour.size();i++){
             //assert(dbg_tour[i]==real_tour[i]);
             }*/

#endif
        }
    };

public:

    int getFullTreeSize(EulerVertex* v){
        if(v->isSingleton())
            return 1;
        else{
            t.dbg_checkSubtreeSize(t.findRoot(v->incidentEdgeA()));
            return t.size(t.findRoot(v->incidentEdgeA())) / 2 +
                   1;//divide by two, to go from half edges to full edges, then add one, to get the number of nodes in the tree.
        }
    }

    bool connected(EulerVertex* from, EulerVertex* to){
        if(from == to){
            return true;
        }
        if(from->incidentEdgeA() == nullptr || to->incidentEdgeA() == nullptr){
            return false;
        }
        return t.findRoot(from->incidentEdgeA()) == t.findRoot(to->incidentEdgeA());
    }

    void makeRoot(EulerVertex* node){
        node->dbg_tour();

        if(!node->isSingleton()){
            if(t.size(t.findRoot(node->incidentEdgeA())) == 2){
                //then both nodes are equivalent to being root.
                node->dbg_make_root();
                assert(t.findRoot(node->incidentEdgeA()) == t.findRoot(node->incidentEdgeB()));

                node->dbg_tour();
                return;
            }

            node->dbg_tour();
            dbg_printTour(node);
            //split the tree at the two incident edges
            EulerHalfEdge* f = node->incidentEdgeA()->value;
            EulerHalfEdge* b = node->incidentEdgeB()->value;
            assert(f != b);
            if(t.compare(f->node, b->node) > 0){
                std::swap(f, b);
            }
            assert(t.size(t.findRoot(f->node)) % 2 == 0);                                //full tree is always even

            //we have to do something slightly tricky here, which is to figure out whether f, or one if its neighbours, is the edge we should be splitting on.
            //this is tricky because we aren't keeping track of this information explicitly with extra 'vertex' nodes in the binary search tree, t, which is the usual solution

            //t.splay(f->node);
            //assert(f->node->right);//because f is less than the other incident edge, so it can't be the rightmost edge.

            //if there is a node to the right, then there is also a successor node (which may or may not be f->node->right).
            EulerVertex* other = (f->to == node) ? f->from : f->to;
            assert(other != node);
            EulerHalfEdge* next = f->node->next()->value;
            //the idea is to look at up to two adjacent nodes to 'node'. This will give us enough information to figure out where to split the tour.

            if(!next->contains(node)){
                //1) next does NOT contain node. in that case, we need to split on the previous node (if it exists)
                Tree::Node* fn = f->node->prev();
                if(!fn){
                    //node is ALREADY the root, don't do anything.
                    node->dbg_make_root();
                    assert(t.findRoot(node->incidentEdgeA()) == t.findRoot(node->incidentEdgeB()));

                    node->dbg_tour();
                    return;
                }else
                    f = fn->value;
            }else if(next->contains(other)){
                //2) next connects the exact same nodes as f, in which case one of them is a leaf.
                //we need to figure out which one is the leaf and which one is the node, to figure out where to split.
                //The node is the only one that will be an element of a third edge
                Tree::Node* next_next = next->node->next();
                if(!next_next)
                    next_next = f->node->prev();

                if(next_next->value->contains(node)){
#ifdef DEBUG_DGL
                    assert(other->dbg_parent == node);
#endif
                    f = next;                //We want to cut AFTER returning to node.
                }else{
                    //node is the leaf, so f is the correct place to cut.
                }

            }else{
                //3) next contains node, but not the other vertex of f.
                //then f is the right place to cut
            }

            //ok, f is before b in the tour. Now we are going to split the tour after f, and then append the section that ends with f  to the right hand side of the tour.
            Tree::Node* right = t.splitAfter(f->node);

            dbg_printEdge(f->node);

            if(right){
                //assert(t.size(right)%2==0);
                assert(t.findMax(t.findRoot(f->node)) == f->node);
                //assert(t.findMin(t.findRoot(right))->value->to == node || t.findMin(t.findRoot(right))->value->from == node );
                assert(t.findRoot(right) != t.findRoot(f->node));
                t.concat(right,
                         f->node);//'right' was previously the right hand side of the tour. After this operation, it is the left hand side
            }else{
                //this already is root
            }
            assert(t.findMax(t.findRoot(f->node)) == f->node);

            //assert(t.findMin(t.findRoot(f->node))->value->to ==node || t.findMin(t.findRoot(f->node))->value->from ==node );

            node->dbg_make_root();
            assert(t.findRoot(node->incidentEdgeA()) == t.findRoot(node->incidentEdgeB()));
            dbg_printTour(node);
            node->dbg_tour();
        }

    }

    //Make othernode a child of node.
    int link(EulerVertex* node, EulerVertex* otherNode, int edgeID){
        assert(node != otherNode);
        assert(!connected(node, otherNode));
        nComponents--;

        node->dbg_tour();
        otherNode->dbg_tour();

        dbg_printTour(node);
        dbg_printTour(otherNode);

        //Create half edges and link them to each other
        if(forward_edges.size() <= edgeID){
            assert(backward_edges.size() <= edgeID);
            forward_edges.resize(edgeID + 1);
            backward_edges.resize(edgeID + 1);
        }
        if(forward_edges[edgeID] == NULL){
            forward_edges[edgeID] = new EulerHalfEdge();
            backward_edges[edgeID] = new EulerHalfEdge();
            forward_edges[edgeID]->index = edgeID;
            backward_edges[edgeID]->index = edgeID;
            forward_edges[edgeID]->isforward = true;

            forward_edges[edgeID]->from = node;
            forward_edges[edgeID]->to = otherNode;

            backward_edges[edgeID]->isforward = false;
            backward_edges[edgeID]->from = otherNode;
            backward_edges[edgeID]->to = node;

            forward_edges[edgeID]->node = t.createNode(forward_edges[edgeID]);
            backward_edges[edgeID]->node = t.createNode(backward_edges[edgeID]);
        }else{
            assert(forward_edges[edgeID]->index == edgeID);
            assert(backward_edges[edgeID]->index == edgeID);

            assert(forward_edges[edgeID]->from == node || forward_edges[edgeID]->to == node);
            assert(forward_edges[edgeID]->to == node || forward_edges[edgeID]->from == node);

        }

        makeRoot(node);
        makeRoot(otherNode);
        dbg_printTour(node);
        dbg_printTour(otherNode);
        node->dbg_tour();
        otherNode->dbg_tour();

        Tree::Node* f = nullptr;
        if(node->incidentEdgeA())
            f = t.findMin(t.findRoot(node->incidentEdgeA()));

        dbg_printEdge(f);

        //ok, f is before b in the tour

        //Tree::Node * tleft = nullptr;

        if(f)
            t.concat(f, forward_edges[edgeID]->node);
        else
            node->setIncidentEdgeA(forward_edges[edgeID]);

        if(otherNode->incidentEdgeA()){
            assert(t.findRoot(otherNode->incidentEdgeA()) == t.findRoot(otherNode->incidentEdgeB()));

            t.concat(forward_edges[edgeID]->node, otherNode->incidentEdgeA());
            assert(t.findRoot(otherNode->incidentEdgeA()) == t.findRoot(otherNode->incidentEdgeB()));
            assert(t.findRoot(forward_edges[edgeID]->node) == t.findRoot(otherNode->incidentEdgeB()));
        }else
            otherNode->setIncidentEdgeA(forward_edges[edgeID]);

        if(otherNode->incidentEdgeB()){
            assert(t.findRoot(otherNode->incidentEdgeA()) == t.findRoot(otherNode->incidentEdgeB()));
            t.concat(otherNode->incidentEdgeB(), backward_edges[edgeID]->node);
            assert(t.findRoot(otherNode->incidentEdgeA()) == t.findRoot(backward_edges[edgeID]->node));
            assert(t.findRoot(node->incidentEdgeA()) == t.findRoot(backward_edges[edgeID]->node));
        }else{
            otherNode->setIncidentEdgeB(backward_edges[edgeID]);
            t.concat(node->incidentEdgeA(), backward_edges[edgeID]->node);
        }
        //Is this a bug?
        /*if (tleft) {
            t.concat(backward_edges[edgeID]->node, tleft);
            node->setIncidentEdgeB(tleft->value);
        } else {*/
        node->setIncidentEdgeB(backward_edges[edgeID]);
        //}
        node->dbg_insert(otherNode);
        dbg_printTour(node);

        assert(t.findRoot(forward_edges[edgeID]->node) == t.findRoot(backward_edges[edgeID]->node));

        node->dbg_tour();
        otherNode->dbg_tour();

        t.findMin(t.findRoot(node->incidentEdgeA()))->value->to->dbg_tour();
        t.findMin(t.findRoot(node->incidentEdgeA()))->value->from->dbg_tour();
        //update subtree sizes going up the tree

        assert(connected(node, otherNode));

        //Return half edge
        return edgeID;
    }

public:
    //Cut an edge in the tree, splitting it into two
    //Returns true if the two components have in fact split.
    void cut(int edgeID){

        nComponents++;

        EulerHalfEdge* f = forward_edges[edgeID];
        EulerHalfEdge* b = backward_edges[edgeID];
        EulerVertex* from = f->from;
        EulerVertex* to = f->to;
        dbg_printTour(from);
        dbg_printTour(to);

        f->from->dbg_tour();
        f->to->dbg_tour();

        assert(connected(f->from, f->to));
        assert(t.findRoot(f->node) == t.findRoot(b->node));

        if(t.compare(f->node, b->node) > 0){
            std::swap(f, b);
        }

        //ok, f is before b in the tour now
        Tree::Node* p = f->node->prev();
        Tree::Node* n = b->node->next();
        Tree::Node* pn = f->node->next();
        Tree::Node* nn = b->node->prev();

        Tree::Node* t1 = t.splitBefore(f->node);
        assert(!t1 || t.findMax(t1) == p);
        assert(t.findRoot(t1) != t.findRoot(f->node));
        Tree::Node* t2 = t.splitAfter(b->node);
        assert(!t2 || t.findMin(t2) == n);

        assert(t.findRoot(t1) != t.findRoot(f->node));
        assert(t.findRoot(t2) != t.findRoot(b->node));
        assert(t.findRoot(t1) != t.findRoot(b->node));
        assert(t.findRoot(t2) != t.findRoot(f->node));

        dbg_printTour(from);
        dbg_printTour(to);

        if(t1 && t2)
            t.concat(t1, t2);    //rejoin the two ends of the outer tour

        //dbg_printTour(from);
        //dbg_printTour(to);

        assert(t.findRoot(t1) != t.findRoot(f->node));
        assert(t.findRoot(t2) != t.findRoot(b->node));
        assert(t.findRoot(t1) != t.findRoot(b->node));
        assert(t.findRoot(t2) != t.findRoot(f->node));
        dbg_printTour(from);
        dbg_printTour(to);
        //ok, now we need to pick new incident edges for both vertices
        //these are the previous and next edges

        assert(pn);
        assert(nn);

        if(pn->value->contains(from) && pn->value->contains(to)){
            //at least one of from, to is a leaf
            if((n || p) && !(n && p)){
                if(!n){
                    n = t.findMin(t.findRoot(p));
                }else{
                    p = t.findMax(t.findRoot(n));
                }
                assert(n != p);
            }
            if(n){
                assert(p);
                dbg_printTour(from);
                dbg_printTour(to);
                assert(n->value->contains(from) || n->value->contains(to));
                assert(p->value->contains(from) || p->value->contains(to));

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
            if((n || p) && !(n && p)){
                if(!n){
                    n = t.findMin(t.findRoot(p));
                }else{
                    p = t.findMax(t.findRoot(n));
                }
                assert(n != p);
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
            if((n || p) && !(n && p)){
                if(!n){
                    n = t.findMin(t.findRoot(p));
                }else{
                    p = t.findMax(t.findRoot(n));
                }
                assert(n != p);
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

        dbg_printTour(from);
        dbg_printTour(to);
        //finally, remove these half edges from the inner tour
        t.splitAfter(f->node);
        t.splitBefore(b->node);
        //check if either of these has become a singleton
        if(from->incidentEdgeA() && t.size(t.findRoot(from->incidentEdgeA())) == 1){
            from->setIncidentEdgeA(nullptr);
            from->setIncidentEdgeB(nullptr);
        }
        if(to->incidentEdgeA() && t.size(t.findRoot(to->incidentEdgeA())) == 1){
            to->setIncidentEdgeA(nullptr);
            to->setIncidentEdgeB(nullptr);
        }
        assert(t.size(f->node) == 1);
        assert(t.size(b->node) == 1);

        /*	t.setIncident(f->node,0);
         t.setIncident(b->node,0);*/
        assert(!connected(from, to));
#ifdef DEBUG_DGL
        if (f->to->dbg_parent == f->from)
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

    void dbg_real_tour(EulerVertex* v, std::vector<int>& tour_out){
        v->dbg_real_tour(tour_out);
    }

    void dbg_printTour(EulerVertex* v){
        dbg_printDbgTour(v);

        std::vector<int> tour_list;
        v->dbg_real_tour(tour_list);
#ifdef dbg_print
        printf("tour:");
        for(int i:tour_list) {
            printf("%d,",i);
        }
        printf("\n");
#endif
    }

    void dbg_printEdge(Tree::Node* e){
        if(!e)
            return;
#ifdef dbg_print
        printf("(%d -> %d)\n",e->value->from->index,e->value->to->index);
#endif
    }

    void dbg_printDbgTour(EulerVertex* v){
        std::vector<int> tour_list;
        v->dbg_build_tour(v, tour_list);
#ifdef dbg_print
        printf("dbgtour:");
        for(int i:tour_list) {
            printf("%d,",i);
        }
        printf("\n");
#endif
    }

    void link(int u, int v, int edgeID){
        link(vertices[u], vertices[v], edgeID);
    }

    bool edgeInTree(int edgeID){
        return forward_edges[edgeID]->node->parent || forward_edges[edgeID]->node->left
               || forward_edges[edgeID]->node->right;
    }

    /*	void cut(int u) {
     cut(vertices[u]);
     }*/
    bool connected(int u, int v){
        return connected(vertices[u], vertices[v]);
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

    struct iterator {
        EulerHalfEdge* n;
        EulerVertex* start;

        bool backward;

        //static std::vector<Tree::Node *>  iterator_stack;
        iterator(EulerVertex* singleton) :
                n(nullptr), start(singleton), backward(false){

        }

        iterator(EulerVertex* first, EulerHalfEdge* n) :
                n(n), start(first), backward(false){

        }

        iterator() :
                n(nullptr), start(nullptr), backward(false){

        }

        bool operator!=(const iterator& other) const{
            return !(other.n == n && other.start == start);
        }

        bool operator==(const iterator& other) const{
            return other.n == n && other.start == start;
        }

        iterator& operator++(){

            if(!n){
                start = nullptr;
                return *this;
            }
            if(n->from->visited != true){
                n->from->visited = true;
                if(n->to->visited != true){
                    return *this;
                }
            }

            n->to->visited = true;
            //first traverse forwards to the end of the bst from the start node
            if(!backward){
                while(n->node->next()){
                    n = n->node->next()->value;
                    if(!n->from->visited || !n->to->visited){
                        return *this;
                    }
                }
                backward = true;        //now switch to traversing backward from the start node
                n = start->incidentEdgeA()->value;
            }
            assert(backward);
            while(n->node->prev()){
                n = n->node->prev()->value;
                if(!n->from->visited || !n->to->visited){

                    return *this;
                }
            }
            n = nullptr;
            start = nullptr;

            backward = false;
            return *this;
        }

        int operator*() const{

            if(!n){
                return start->index;
            }else{

                assert(n->from->visited != true || n->to->visited != true);
                if(n->from->visited){
                    return n->to->index;
                }else{
                    return n->from->index;
                }
            }

        }

    };

    struct tour_iterator {
        EulerHalfEdge* n;
        EulerHalfEdge* start;
        bool backward;
        bool strict;
        bool restarted;

        tour_iterator(EulerHalfEdge* start, bool strict = false) :
                n(start), start(start), backward(false), strict(strict), restarted(false){
#ifdef DEBUG_DGL
            if (start) {
                //check that the nodes are well structured...
                Tree::Node * t = start->node;
                int total = 1;
                int size = t->findRoot()->subtree_size;
                while (t) {
                    Tree::Node *t2 = t->next();
                    assert(t != t2);
                    if (t2) {
                        assert(t2->prev() == t);
                    }
                    if (t2) {
                        total++;
                        assert(total <= size);

                    }
                    t = t2;
                }
            }
#endif
        }

        bool operator!=(const tour_iterator& other) const{
            return !(other.n == n);
        }

        bool operator==(const tour_iterator& other) const{
            return other.n == n;
        }

        tour_iterator& operator++(){
            if(!n)
                return *this;

            //first traverse forwards to the end of the bst from the start node
            if(!backward){
                if(n->node->next()){
                    if(n->node->next()->value == start){
                        //we are done. this can only happen in strict mode.
                        n = nullptr;
                        start = nullptr;
                        backward = false;
                        return *this;
                    }
                    assert(n->node != n->node->next());
                    n = n->node->next()->value;
                    return *this;
                }

                if(strict && !restarted){
                    restarted = true;
                    //ok, now continue the tour from the left hand side
                    //this costs log time, so this isn't the default behaviour
                    n = start->node->findRoot()->findMin()->value;
                    if(n == start){
                        n = nullptr;
                        start = nullptr;
                        backward = false;                        //then we are done.
                    }
                    return *this;
                }else if(strict && restarted){
                    n = nullptr;
                    start = nullptr;
                    backward = false;
                    return *this;
                }else{
                    backward = true;                        //now switch to traversing backward from the start node
                    n = start;
                }
            }
            assert(backward);
            if(n->node->prev()){
                n = n->node->prev()->value;
                return *this;
            }
            n = nullptr;
            start = nullptr;

            backward = false;
            return *this;
        }

        int operator*() const{

            assert(n);
            return n->index;

        }

    };

    //Traverse a full tour starting from this node, in _arbitrary_ order.
    //If strict_tour is set, then the tour will be in 'proper' tour order, starting and ending at 'fromVertex'.
    tour_iterator begin_half_edge_tour(int fromVertex, bool strict_tour = false){
        if(vertices[fromVertex]->isSingleton()){
            return tour_iterator(nullptr);
        }else{
            EulerHalfEdge* h = vertices[fromVertex]->incidentEdgeA()->value;
            //ensure that we start on an outgoing, not an incoming, edge in the tour
            if(strict_tour){
                makeRoot(vertices[fromVertex]);
                h = t.findMin(t.findRoot(h->node))->value;
            }
            return tour_iterator(h, strict_tour);
        }

    }

    tour_iterator end_half_edge_tour(){
        return tour_iterator(nullptr);
    }

    //Traverse the full tree starting from this node, in arbitrary order.
    iterator begin(int fromVertex){
        if(vertices[fromVertex]->isSingleton()){
            return iterator(vertices[fromVertex]);
        }else{
            //Tree::Node * rootNode = t.findRoot(vertices[fromVertex]->incidentEdgeA());
            //EulerVertex * first = rootNode->value->from;
            //first, clear any visited vertices...
            //this is really inefficient! Get rid of this!
            Tree::Node* n = vertices[fromVertex]->left_out->node;
            while(n){
                n->value->from->visited = false;
                n->value->to->visited = false;
                n = n->next();
            }
            n = vertices[fromVertex]->left_out->node;
            while(n){
                n->value->from->visited = false;
                n->value->to->visited = false;
                n = n->prev();
            }
            return iterator(vertices[fromVertex], vertices[fromVertex]->incidentEdgeA()->value);
        }
        return iterator();
    }

    iterator end(){
        return iterator();
    }

    //find an _arbitrary_, but consistent, root node for the tree this node is connected to.
    int findRoot(int node){
        if(vertices[node]->isSingleton()){
            return node;
        }else{
            return t.findRoot(vertices[node]->incidentEdgeA())->value->from->index;        //arbitrarily pick this node.
        }
    }

    int nVertices(){
        return vertices.size();
    }

    void clear(){
        nComponents = vertices.size();
        for(int i = 0; i < vertices.size(); i++){
            vertices[i]->clear();
        }
        for(int i = 0; i < forward_edges.size(); i++){
            if(forward_edges[i] && forward_edges[i]->node)
                t.clear(forward_edges[i]->node);
        }
        for(int i = 0; i < backward_edges.size(); i++){
            if(backward_edges[i] && backward_edges[i]->node)
                t.clear(backward_edges[i]->node);
        }
    }

    //return the size of the complete tree that v is an element of
    int getFullTreeSize(int v){
        return getFullTreeSize(vertices[v]);
    }

    void createVertex(){
        nComponents++;
        vertices.push_back(new EulerVertex());
#ifdef DEBUG_DGL
        vertices.back()->owner = this;
#endif
        vertices.back()->index = vertices.size() - 1;
        //return vertices.back();
    }

    EulerTree(){
        nComponents = 0;

    }
};

#endif

