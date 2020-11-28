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

#ifndef ALLPAIRS_H_
#define ALLPAIRS_H_

#include <vector>

namespace dgl {
class AllPairs {
public:

    struct NullStatus {
        void setReachable(int u, int v, bool reachable){

        }

        void setMininumDistance(int u, int v, bool reachable, int distance){

        }

    };

    static NullStatus nullStatus;

    virtual ~AllPairs(){
    };

    virtual int numUpdates() const = 0;

    virtual void addSource(int s) = 0;

    virtual void update() = 0;

    virtual bool connected_unsafe(int from, int t) = 0;

    virtual bool connected_unchecked(int from, int t) = 0;

    virtual bool connected(int from, int t) = 0;

    virtual int distance(int from, int t) = 0;

    virtual int distance_unsafe(int from, int t) = 0;

    //Return a path of edges from source to 'to'
    virtual void getPath(int source, int to, std::vector<int>& path_store) = 0;
};
};

#endif
