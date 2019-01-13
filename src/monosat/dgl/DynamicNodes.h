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
#ifndef DYNAMIC_NODES
#define DYNAMIC_NODES

#include <vector>
#include <cassert>

class DynamicNodes {
    std::vector<bool> nodeStatus;
    int n_enabled = 0;
    int n_actualNodes = 0;
public:
    void addNode(int nodeID){
        assert(nodeID >= 0);
        if(nodeStatus.size() <= nodeID)
            nodeStatus.resize(nodeID + 1);
        nodeStatus[nodeID] = true;
        n_enabled++;
        n_actualNodes++;
    }

    int numEnabled() const{
        return n_enabled;
    }

    int nodes() const{
        return nodeStatus.size();
    }

    bool nodeEnabled(int n) const{
        return nodeStatus[n];
    }

    void setNodeEnabled(int n, bool enabled){
        if(nodeStatus[n] != enabled){
            nodeStatus[n] = enabled;
            if(enabled){
                n_enabled++;
            }else{
                n_enabled--;
            }
        }
        assert(n_enabled >= 0);
        assert(n_enabled <= n_actualNodes);
    }
};

#endif
