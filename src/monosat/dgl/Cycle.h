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

#ifndef CYCLE_H_
#define CYCLE_H_

#include <vector>

namespace dgl {
class Cycle {
public:

    bool marked;

    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    int64_t stats_history_clears = 0;
    double mod_percentage;

    double stats_full_update_time;
    double stats_fast_update_time;

    virtual ~Cycle(){
    }

    //hint to the algorithm that any discovered cycles will be removed.
    virtual void forceDAG(){

    }

    virtual bool hasDirectedCycle() = 0;

    virtual bool hasUndirectedCycle() = 0;

    virtual void update() = 0;

    virtual std::vector<int>& getUndirectedCycle() = 0;

    virtual std::vector<int>& getDirectedCycle() = 0;
};
};
#endif /* REACH_H_ */
