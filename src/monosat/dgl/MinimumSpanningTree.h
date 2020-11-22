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

#ifndef MINIMUMSPANNINGTREE_H_
#define MINIMUMSPANNINGTREE_H_

#include <vector>

namespace dgl {
template<typename Weight = int>
class MinimumSpanningTree {
public:

    struct NullStatus {
        void setMinimumSpanningTree(Weight& min_weight, bool connected){

        }

        void inMinimumSpanningTree(int edge, bool in_tree){

        }
    };

    static NullStatus nullStatus;

    virtual ~MinimumSpanningTree(){
    };

    virtual void printStats(){

    }

    virtual int numUpdates() const = 0;

    virtual void update() = 0;

    virtual bool dbg_uptodate() = 0;

    virtual bool dbg_mst() = 0;

    //Total weight of the mst (or infinite, if the graph is disconnected)
    virtual Weight& weight() = 0;

    //Sum of the weight of the mst of each tree in the forest
    virtual Weight& forestWeight() = 0;

    virtual std::vector<int>& getSpanningTree() = 0;

    virtual int getParent(int node) = 0;

    virtual int getParentEdge(int node) = 0;

    virtual bool edgeInTree(int edgeid) = 0;

    virtual int numComponents() = 0;

    virtual int getComponent(int node) = 0;

    virtual int getRoot(int component = 0) = 0;
};

template<typename Weight>
typename MinimumSpanningTree<Weight>::NullStatus MinimumSpanningTree<Weight>::nullStatus;
};

#endif /* MINIMUMSPANNINGTREE_H_ */
