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

#ifndef DYNAMIC_CONNECTIVITY_IMPL_H_
#define DYNAMIC_CONNECTIVITY_IMPL_H_
namespace dgl {
class DynamicConnectivityImpl {
public:

    DynamicConnectivityImpl(){

    }

    virtual ~DynamicConnectivityImpl(){

    }

    virtual int numComponents() = 0;

    virtual bool connected(int u, int v) = 0;

    virtual void addNode() = 0;

    virtual void addEdge(int from, int to, int edgeID) = 0;

    virtual bool edgeEnabled(int edgeid) const = 0;

    virtual void dbg_print() = 0;

    virtual bool setEdgeEnabled(int from, int to, int edgeid, bool enabled) = 0;
};
};
#endif
