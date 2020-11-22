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

#ifndef STEINERTREE_H_
#define STEINERTREE_H_

#include <vector>

namespace dgl {

template<typename Weight = int>
class SteinerTree {
public:

    struct NullStatus {

        void setMinimumSteinerTree(Weight min_weight){

        }

    };

    static NullStatus nullStatus;

    virtual ~SteinerTree(){
    };

    virtual void printStats(){

    }

    virtual void update() = 0;

    //Total weight of the steiner tree (or infinite, if the graph is disconnected)
    virtual Weight& weight() = 0;

    virtual bool disconnected() = 0;

    virtual void getSteinerTree(std::vector<int>& edges) = 0;
};

template<typename Weight>
typename SteinerTree<Weight>::NullStatus SteinerTree<Weight>::nullStatus;
};
#endif
