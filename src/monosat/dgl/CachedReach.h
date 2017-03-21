//
// Created by sam on 20/03/17.
//

#ifndef MONOSAT_CACHED_REACH_H
#define MONOSAT_CACHED_REACH_H


#include <vector>
#include "monosat/dgl/alg/Heap.h"
#include "DynamicGraph.h"
#include "monosat/core/Config.h"
#include "Reach.h"
#include "alg/IntMap.h"

namespace dgl {

template<typename Weight, bool undirected = false>
class CachedReach: public Reach {
public:

    DynamicGraph<Weight> & g;

    int last_modification=-1;
    int last_addition=-1;
    int last_deletion=-1;
    int history_qhead=0;
    Reach * reach;
    int last_history_clear=0;



    long stats_full_updates=0;
    long stats_fast_updates=0;
    long stats_fast_failed_updates=0;
    long stats_skip_deletes=0;
    long stats_skipped_updates=0;
    long stats_num_skipable_deletions=0;


    double stats_full_update_time=0;
    double stats_fast_update_time=0;


public:


    CachedReach(Reach * reach,DynamicGraph<Weight> & graph) :
            g(graph), reach(reach) {


    }

    //Connectivity(const Connectivity& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),mod_percentage(0.2),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};

    void setSource(int s) {

        last_modification = -1;
        last_addition = -1;
        last_deletion = -1;
        reach->setSource(s);
    }
    int getSource() {
        return reach->getSource();
    }


    long num_updates = 0;
    int numUpdates() const {
        return num_updates;
    }


    alg::IntSet<int> destinations;
    alg::IntSet<int> edge_in_path;
    //std::vector<bool> edge_in_path;
    std::vector<bool> has_path_to;

    bool has_non_reach_destinations=false;

     void addDestination(int node) override{
         if(!destinations.has(node)) {
             destinations.insert(node);
             has_non_reach_destinations=true;
         }
    }
     void removeDestination(int node) override{

    }
    void recompute(){
        reach->update();
        edge_in_path.clear();//clear and rebuild the path tree
        int source = getSource();
        has_non_reach_destinations=false;
        for(int d:destinations){

            if(d==getSource()){
                has_path_to[d]=true;
                continue;
            }

            if(reach->connected(d)){
                has_path_to[d]=true;
                //extract path
                int edgeID = reach->incomingEdge(d);
                assert(edgeID>=0);
                int p = d;
                while(!edge_in_path.has(edgeID)){
                    assert(edgeID>=0);

                    edge_in_path.insert(edgeID);
                    p = reach->previous(p);
                    if(p==source){
                        break;
                    }
                    edgeID = reach->incomingEdge(p);
                    assert(edgeID>=0);
                }
            }else{
                has_path_to[d]=false;
                has_non_reach_destinations=true;
            }
        }

        num_updates++;
        last_modification = g.modifications;
        last_deletion = g.deletions;
        last_addition = g.additions;

        history_qhead = g.historySize();

        last_history_clear = g.historyclears;
    }
    void update() override{
        static int iteration = 0;
        int local_it = ++iteration;
        if (last_modification > 0 && g.modifications == last_modification){
            return;
        }
        bool needs_recompute=false;

        if (last_modification <= 0 || g.changed()) {//Note for the future: there is probably room to improve this further.
            edge_in_path.clear();

            //edge_in_path.resize(g.edges(),false);
            has_path_to.clear();
            has_path_to.resize(g.nodes(),false);
           // prev_edge.clear();
            //prev_edge.resize(g.nodes(),-1);
            needs_recompute=true;

        }

        if (last_history_clear != g.historyclears) {
            if (!g.changed() && last_history_clear>=0 && last_history_clear == g.historyclears-1 && history_qhead==g.getPreviousHistorySize()){
                //no information was lost in the history clear
                history_qhead = 0;
                last_history_clear = g.historyclears;
            }else {
                history_qhead = g.historySize();
                last_history_clear = g.historyclears;
                for (int edgeid = 0; edgeid < g.edges(); edgeid++) {
                    if (g.edgeEnabled(edgeid)) {
                        if(has_non_reach_destinations){
                            //needs recompute
                            needs_recompute=true;
                            break;
                        }
                    } else if (edge_in_path.has(edgeid)){
                        //needs recompute
                        needs_recompute=true;
                        break;
                    }
                }
            }
        }

        for (int i = history_qhead; i < g.historySize(); i++) {
            int edgeid = g.getChange(i).id;

            if (g.getChange(i).addition && g.edgeEnabled(edgeid)) {
                int v = g.getEdge(edgeid).to;
                if(has_non_reach_destinations){
                    //needs recompute
                    needs_recompute=true;
                    break;
                }
            } else if (!g.getChange(i).addition && !g.edgeEnabled(edgeid)) {
                int v = g.getEdge(edgeid).to;
                if(edge_in_path.has(edgeid)){
                    needs_recompute=true;
                    break;
                }
            }
        }


        if(needs_recompute){
            recompute();

        }


        num_updates++;
        last_modification = g.modifications;
        last_deletion = g.deletions;
        last_addition = g.additions;

        history_qhead = g.historySize();

        last_history_clear = g.historyclears;

    }

    bool connected_unsafe(int t) {
        return has_path_to[t];
    }
    bool connected_unchecked(int t) {
        assert(last_modification == g.modifications);
        return connected_unsafe(t);
    }
    bool connected(int t) {
        if (last_modification != g.modifications)
            update();
        return has_path_to[t];
    }
    int distance(int t) {
        if (connected(t))
            return 1;
        else
            return -1;
    }
    int distance_unsafe(int t) {
        if (connected_unsafe(t))
            return 1;
        else
            return -1;
    }
    int incomingEdge(int t) {
        return reach->incomingEdge(t);
    }
    int previous(int t) {
        return reach->previous(t);
    }

};
}
;

#endif //MONOSAT_CACHED_REACH_H
