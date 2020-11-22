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

#ifndef SEARCH_TREE_H
#define SEARCH_TREE_H

template<typename Node>
class SearchTree {

public:


    SearchTree(){

    }

    virtual ~SearchTree(){

    }

    virtual Node* findRoot(Node* of) = 0;

    virtual Node* findMin(Node* root) = 0;

    virtual Node* findMax(Node* root) = 0;

    //virtual void remove(Node * toRemove)=0;
    virtual Node* splitAfter(Node* splitAt) = 0;

    virtual Node* splitBefore(Node* splitAt) = 0;

    //All nodes in right must be (strictly) greater than all nodes in left
    virtual Node* concat(Node* left, Node* right) = 0;
};

template<typename Node>
class AugmentedSearchTree : public SearchTree<Node> {

public:

    AugmentedSearchTree(){

    }

    virtual ~AugmentedSearchTree(){

    }

    //Returns the number of nodes in the subtree rooted at n (including n)
    virtual int size(Node* n) = 0;

    virtual Node* findRoot(Node* of) = 0;

    virtual Node* findMin(Node* root) = 0;

    virtual Node* findMax(Node* root) = 0;

    virtual int depth(Node* a) = 0;

    //Return -1 if a< b, 1 if a>b, 0 if either they are the same node, or not in the same tree.
    virtual int compare(Node* a, Node* b) = 0;

    //virtual void remove(Node * toRemove)=0;
    virtual Node* splitAfter(Node* splitAt) = 0;

    virtual Node* splitBefore(Node* splitAt) = 0;

    //All nodes in right must be (strictly) greater than all nodes in left
    virtual Node* concat(Node* left, Node* right) = 0;
};

#endif
