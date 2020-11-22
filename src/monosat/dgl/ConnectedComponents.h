#ifndef CONNECTED_COMPONENTS_H_
#define CONNECTED_COMPONENTS_H_

#include <stdexcept>
#include <vector>

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

namespace dgl {
class ConnectedComponents {
public:

    struct NullConnectedComponentsStatus {
        void setConnected(int u, int v, bool connected){

        }

        void setComponents(int components){

        }
    };

    static NullConnectedComponentsStatus nullConnectedComponentsStatus;

    virtual ~ConnectedComponents(){
    };

    virtual void printStats(){

    }

    virtual void update() = 0;

    virtual void addConnectedCheck(int u, int v){
        throw std::runtime_error("Not implemented");
    }

    virtual int numComponents() = 0;

    //Get an arbitrary element from the given set
    virtual int getElement(int set){
        throw std::runtime_error("Not implemented");
    }

    //Get the component this element belongs to
    virtual int getComponent(int node) = 0;

    virtual bool connected(int from, int to) = 0;
};
};
#endif
