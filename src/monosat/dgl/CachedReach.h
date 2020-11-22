/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2017, Sam Bayless

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


#ifndef MONOSAT_CACHED_REACH_H
#define MONOSAT_CACHED_REACH_H


#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "Graph.h"
#include "DynamicGraph.h"

#include "monosat/core/Config.h"
#include "Reach.h"
#include "alg/IntMap.h"
#include "alg/Rnd.h"

namespace dgl {

template<typename Weight, typename Graph = DynamicGraph<Weight>, class Status = typename Distance<Weight>::NullStatus, bool undirected = false>
class CachedReach : public Reach {
public:

    Graph& g;
    Status& status;
    int last_modification = -1;
    int last_addition = -1;
    int last_deletion = -1;
    int history_qhead = 0;
    int reportPolarity = 0;
    Reach* reach;
    int last_history_clear = 0;


    int64_t stats_full_updates = 0;
    int64_t stats_fast_updates = 0;
    int64_t stats_fast_failed_updates = 0;
    int64_t stats_skip_deletes = 0;
    int64_t stats_skipped_updates = 0;
    int64_t stats_num_skipable_deletions = 0;
    int64_t stats_random_shortest_paths = 0;
    int64_t stats_random_shortest_edges = 0;
    int64_t stats_n_recomputes = 0;
    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;
    double random_seed = 0;
    double randomShortestPathFrequency = 0;//With this probability, select a random (but still shortest) path
    double randomShortestEdgeFrequency = 0;
public:


    CachedReach(Reach* reach, Graph& graph, Status& status, int reportPolarity = 0,
                double randomShortestPathFrequency = 0, double randomShortestEdgeFrequency = 0, double random_seed = 0)
            :
            g(graph), status(status), reach(reach), reportPolarity(reportPolarity),
            randomShortestPathFrequency(randomShortestPathFrequency),
            randomShortestEdgeFrequency(randomShortestEdgeFrequency), random_seed(random_seed){

    }

    CachedReach(Reach* reach, Graph& graph, int reportPolarity = 0, double randomShortestPathFrequency = 0,
                double randomShortestEdgeFrequency = 0, double random_seed = 0) :
            g(graph), status(Distance<Weight>::nullStatus), reach(reach), reportPolarity(reportPolarity),
            randomShortestPathFrequency(randomShortestPathFrequency),
            randomShortestEdgeFrequency(randomShortestEdgeFrequency), random_seed(random_seed){

    }

    void setSource(int s) override{
        if(s != reach->getSource()){
            needs_recompute = true;
            last_modification = -1;
            last_addition = -1;
            last_deletion = -1;
            reach->setSource(s);
        }
    }

    int getSource() override{
        return reach->getSource();
    }


    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }


    alg::IntSet<int> destinations;
    alg::IntSet<int> edge_in_path;
    std::vector<bool> has_path_to;
    std::vector<int> previous_edge;
    bool has_non_reach_destinations = false;
    bool needs_recompute = true;


    void clearCache() override{
        needs_recompute = true;
    }

    void addDestination(int node) override{
        if(!destinations.has(node)){
            destinations.insert(node);
            has_non_reach_destinations = true;
        }
    }

    void removeDestination(int node) override{

    }

    void printStats() override{
        reach->printStats();
        printf("Cached Reach Recomputes: %" PRId64 "\n", stats_n_recomputes);
        if(randomShortestPathFrequency > 0){
            printf("Random shortest paths (%f freq): %" PRId64 "\n", randomShortestPathFrequency,
                   stats_random_shortest_paths);
        }
        if(randomShortestEdgeFrequency > 0){
            printf("Random shortest edges (%f freq): %" PRId64 "\n", randomShortestEdgeFrequency,
                   stats_random_shortest_edges);
        }

    }

    void recompute(){
        needs_recompute = false;
        stats_n_recomputes++;
        reach->update();
        //do not need to reset previous_edge vector here; instead we allow it to contain incorrect values on the assumption they will be corrected before being accessed
        edge_in_path.clear();//clear and rebuild the path tree
        int source = getSource();
        assert(previous_edge[source] == -1);
        has_non_reach_destinations = false;

        bool randomShortestPath = alg::drand(random_seed) < randomShortestPathFrequency;
        if(randomShortestPath){
            stats_random_shortest_paths++;
        }

        for(int d:destinations){

            if(d == getSource()){
                has_path_to[d] = true;
                continue;
            }

            if(reach->connected(d)){
                has_path_to[d] = true;
                if(reportPolarity >= 0){
                    status.setReachable(d, true);
                }
                //extract path
                bool randomShortestEdge = randomShortestPath || (alg::drand(random_seed) < randomShortestEdgeFrequency);
                if(randomShortestEdge && !randomShortestPath){
                    stats_random_shortest_edges++;
                }
                int edgeID = randomShortestEdge ? reach->randomIncomingEdge(d, random_seed) : reach->incomingEdge(d);
                assert(edgeID >= 0);
                int p = d;
                while(!edge_in_path.has(edgeID)){
                    assert(edgeID >= 0);
                    previous_edge[p] = edgeID;//do not need to reset previous_edge vector here; instead we allow it to contain incorrect values on the assumption they will be corrected before being accessed

                    edge_in_path.insert(edgeID);
                    p = g.getEdge(edgeID).from;
                    if(p == source){
                        edgeID = -1;
                        break;
                    }
                    randomShortestEdge = randomShortestPath || (alg::drand(random_seed) < randomShortestEdgeFrequency);
                    if(randomShortestEdge && !randomShortestPath){
                        stats_random_shortest_edges++;
                    }
                    edgeID = randomShortestEdge ? reach->randomIncomingEdge(p, random_seed) : reach->incomingEdge(p);
                    assert(edgeID >= 0);
                }
                assert(randomShortestPath || previous_edge[p] == edgeID);
            }else{
                has_path_to[d] = false;
                previous_edge[d] = -1;
                has_non_reach_destinations = true;
                if(reportPolarity <= 0){
                    status.setReachable(d, false);
                }
            }
        }

        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();

        last_history_clear = g.nHistoryClears();
    }

    void update() override{

        if(!needs_recompute && last_modification > 0 && g.getCurrentHistory() == last_modification){
            return;
        }


        if(last_modification <= 0 || g.changed()){//Note for the future: there is probably room to improve this further.
            edge_in_path.clear();

            has_path_to.clear();
            has_path_to.resize(g.nodes(), false);
            previous_edge.clear();
            previous_edge.resize(g.nodes(), -1);
            needs_recompute = true;

        }

        if(!needs_recompute && last_history_clear != g.nHistoryClears()){
            if(!g.changed() && last_history_clear >= 0 && last_history_clear == g.nHistoryClears() - 1 &&
               history_qhead == g.getPreviousHistorySize()){
                //no information was lost in the history clear
                history_qhead = 0;
                last_history_clear = g.nHistoryClears();
            }else{
                history_qhead = g.historySize();
                last_history_clear = g.nHistoryClears();
                for(int edgeid = 0; edgeid < g.edges(); edgeid++){
                    if(g.edgeEnabled(edgeid)){
                        if(has_non_reach_destinations){
                            //needs recompute
                            needs_recompute = true;
                            break;
                        }
                    }else if(edge_in_path.has(edgeid)){
                        //needs recompute
                        needs_recompute = true;
                        break;
                    }
                }
            }
        }
        if(!needs_recompute){
            for(int i = history_qhead; i < g.historySize(); i++){
                int edgeid = g.getChange(i).id;

                if(g.getChange(i).addition && g.edgeEnabled(edgeid)){
                    int v = g.getEdge(edgeid).to;
                    if(has_non_reach_destinations){
                        //needs recompute
                        needs_recompute = true;
                        break;
                    }
                }else if(!g.getChange(i).addition && !g.edgeEnabled(edgeid)){
                    int v = g.getEdge(edgeid).to;
                    if(edge_in_path.has(edgeid)){
                        needs_recompute = true;
                        break;
                    }
                }
            }
        }

        if(needs_recompute){
            recompute();

        }


        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();

        last_history_clear = g.nHistoryClears();

    }

    bool dbg_manual_uptodate() override{
        return reach->dbg_manual_uptodate();
    }

    bool connected_unsafe(int t) override{
        return has_path_to[t];
    }

    bool connected_unchecked(int t) override{
        assert(last_modification == g.getCurrentHistory() && !needs_recompute);
        return connected_unsafe(t);
    }

    bool connected(int t) override{
        update();
        return has_path_to[t];
    }

    int distance(int t){
        if(connected(t))
            return 1;
        else
            return -1;
    }

    int distance_unsafe(int t){
        if(connected_unsafe(t))
            return 1;
        else
            return -1;
    }

    int incomingEdge(int t) override{
        assert(last_modification == g.getCurrentHistory() && !needs_recompute);
        assert(previous_edge[t] >= 0);
        return previous_edge[t];
    }

    int previous(int t) override{
        assert(last_modification == g.getCurrentHistory() && !needs_recompute);
        if(previous_edge[t] < 0){
            return -1;
        }
        int edgeID = incomingEdge(t);
        assert(edgeID >= 0);
        assert(g.edgeEnabled(edgeID));
        assert(g.getEdge(edgeID).to == t);
        return g.getEdge(edgeID).from;
    }

};
};

#endif //MONOSAT_CACHED_REACH_H
