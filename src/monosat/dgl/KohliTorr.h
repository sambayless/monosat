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

#ifndef MAXFLOW_KOHLI_TORR_H
#define MAXFLOW_KOHLI_TORR_H

//Wrapper around Kohli and Torr's  Dynamic Graph Cuts algorithm (version 2)
//Note that the Kohli Torr implementation itself is under the GPL (Version 2); only this wrapper code is MIT licensed.

#include "Graph.h"
#include "DynamicGraph.h"
#include "MaxFlow.h"
#include <vector>

#include "EdmondsKarpDynamic.h"
#include <algorithm>
#include <limits>

#ifdef DEBUG_DGL
//#define DEBUG_MAXFLOW2
#endif

#ifdef LINK_GPL
//IF GPL sources are included, link them here.
//Else, the 'KohliTorr' class will instead just be the much sloder EdmondsKarp algorithm
#include "monosat/dgl/alg/dyncut/graph.h"

namespace dgl {
template<typename Weight>
class KohliTorr : public MaxFlow<Weight>, public DynamicGraphAlgorithm {
    Weight f = 0;
    Weight static_maxflow = -1;//the maximum possible flow in the graph.
    Weight source_max = 0;
    Weight dest_max = 0;
    Graph<Weight>& g;
    DynamicGraph<Weight>* flow_graph = nullptr;//optional graph that can be used to record the edges in the actual flow
    /**
     * Note: The Kohli Torr implementation does _not_ support multiple edges between the same nodes.
     */
    kohli_torr::Graph<Weight, Weight, Weight>* kt = nullptr;

    struct LocalEdge {
        int from;
        int id;
        bool backward = false;

        LocalEdge(int from = -1, int id = -1, bool backward = false) :
                from(from), id(id), backward(backward){

        }
    };

    int source = -1;
    int sink = -1;
    bool dynamic = true;
    Weight curflow;
    int last_modification;
    int last_deletion;
    int last_addition;
    std::vector<int> tmp_edges;
    std::vector<std::vector<int>> tmp_edge_map;
    int history_qhead;
    int last_history_clear;
    //bool backward_maxflow=false;
    bool kt_preserve_order = false;
    //std::vector<int> all_multi_edges;
    //std::vector<std::vector<int>> multi_edges;

    typedef typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id arc;


    std::vector<int> arc_map; //map from edgeids to arcs (note: this may be many to one)
    //A brief explanation here:
    //The implementation of Kohli and Torr's maxflow algorithm assumes at most one (directed) arc between any two nodes (and a separate reverse arc for every forward arc).
    //In contrast, the DGL allows for multiple (directed) edges between the same two nodes, and for forward edges without matching backward edges.
    //This makes translating between KT and the DGL a pain; the vectors below form a one-to-many map of KT arc ids to all directed dgl edge ids with the same start and end point.
    std::vector<std::pair<int, int>> edge_map_ind;
    std::vector<int> edge_map;

    Weight INF;

    static const auto KT_SOURCE = kohli_torr::Graph<Weight, Weight, Weight>::SOURCE;
    static const auto KT_SINK = kohli_torr::Graph<Weight, Weight, Weight>::SINK;


    Weight max_capacity = 1; //start this at 1, not 0
    Weight sum_of_edge_capacities = 0;
    std::vector<Weight> local_weights;
    std::vector<bool> edge_enabled;
    std::vector<int> changed_edges;
    std::vector<int> changed_partition;
    bool flow_needs_recalc = true;
    bool needs_recompute = true;
    int alg_id = -1;

    inline typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id getArc(int edgeID){
        assert(edgeID < arc_map.size());
        assert(arc_map[edgeID] != -1);
        return kt->get_arc(arc_map[edgeID]);
    }

    std::vector<int>::iterator getArcEdges(int arcID){
        assert(arcID < edge_map_ind.size());
        return edge_map.begin() + edge_map_ind[arcID].first;
    }

    std::vector<int>::iterator getArcEdgesEnd(int arcID){
        assert(arcID < edge_map_ind.size());
        return edge_map.begin() + edge_map_ind[arcID].second;
    }

public:
    double stats_calc_time = 0;
    double stats_flow_time = 0;
    double stats_init_time = 0;
    int64_t stats_flow_calcs = 0;
    int64_t stats_inits = 0;
    int64_t stats_reinits = 0;
    bool same_source_sink = false;//if the sink and source are the same node, flow is infinite, and lots of logic is bypassed to avoid putting KT in a bad state.
    KohliTorr(Graph<Weight>& g, int source, int sink, bool kt_preserve_order = false) :
            g(g), source(source), sink(sink), kt_preserve_order(kt_preserve_order), INF(0xF0F0F0){
        curflow = 0;
        last_modification = -1;
        last_deletion = -1;
        last_addition = -1;
        needs_recompute = true;
        history_qhead = 0;
        last_history_clear = -1;
        if(source == sink){
            if(source == 0)
                this->sink = 1; //but what if the graph has only 1 node...?
            else
                this->sink = 0;

            same_source_sink = true;
        }
        alg_id = g.addDynamicAlgorithm(this);

    }

    std::string getName() override{
        return "KohliTorr(" + std::to_string(getSource()) + ", " + std::to_string(getSink()) + ")";
    }

    int getSource() const override{
        return source;
    }

    int getSink() const override{
        return sink;
    }

    Graph<Weight>* getFlowGraph() override{
        if(!flow_graph){
            flow_graph = new DynamicGraph<Weight>();
            while(flow_graph->nodes() < g.nodes()){
                flow_graph->addNode();
            }
            while(flow_graph->edges() < g.edges()){
                int edgeID = flow_graph->edges();
                int f = g.getEdge(edgeID).from;
                int t = g.getEdge(edgeID).to;
                flow_graph->addEdge(f, t, edgeID, 0);
                flow_graph->disableEdge(edgeID);
            }
        }

        return flow_graph;
    }

    void updateMaxCapacity(Weight new_max_capacity){

        if(same_source_sink)
            return;

        max_capacity_increased = false;
        if(new_max_capacity > max_capacity){
            {
                if(new_max_capacity < 1){
                    new_max_capacity = 1;
                }
                assert(new_max_capacity >= max_capacity);
                if(kt){
                    if(dynamic){
                        Weight dif = new_max_capacity - max_capacity;
                        Weight curweight = kt->getTweight(source);
                        Weight newweight = curweight + dif;
                        assert(newweight >= 0);
                        kt->edit_tweights(source, newweight, 0);
                    }else{
                        Weight dif = new_max_capacity - max_capacity;
                        Weight curweight = kt->getTweight(source);
                        Weight newweight = curweight + dif;
                        assert(newweight >= 0);
                        kt->edit_tweights_wt(source, newweight, 0);

                    }
                }
            }
            {

                if(kt){
                    if(dynamic){
                        Weight dif = new_max_capacity - max_capacity;
                        Weight curweight = -kt->getTweight(sink);
                        assert(curweight >= 0);
                        Weight newweight = curweight + dif;
                        assert(newweight >= curweight);
                        assert(newweight >= 0);
                        kt->edit_tweights(sink, 0, newweight);

                    }else{
                        Weight dif = new_max_capacity - max_capacity;
                        Weight curweight = -kt->getTweight(sink);
                        assert(curweight >= 0);
                        Weight newweight = curweight + dif;
                        assert(newweight >= curweight);
                        assert(newweight >= 0);
                        kt->edit_tweights_wt(sink, 0, newweight);
                    }
                }
            }
            max_capacity = new_max_capacity;
            last_modification = g.getCurrentHistory() - 1;
        }
        flow_needs_recalc = true;
        needs_recompute = true;
    }

    void setSource(int s) override{
        same_source_sink = false;
        assert(source != sink);
        if(source == s){
            return;
        }
        if(s == sink){
            same_source_sink = true;
            return;
        }
        if(kt){

            if(dynamic){
                Weight curweight = kt->getTweight(source);
                Weight newWeight = curweight - max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights(source, newWeight, 0);
                }else{
                    kt->edit_tweights(source, 0, -newWeight);
                }
                //if(!backward_maxflow){
                curweight = kt->getTweight(s);
                newWeight = curweight + max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights(s, newWeight, 0);
                }else{
                    kt->edit_tweights(s, 0, -newWeight);
                }

                /*}else{
                 kt->edit_tweights(s,0,max_capacity);
                 }*/
            }else{
                Weight curweight = kt->getTweight(source);
                Weight newWeight = curweight - max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights_wt(source, newWeight, 0);
                }else{
                    kt->edit_tweights_wt(source, 0, -newWeight);
                }
                //if(!backward_maxflow){
                curweight = kt->getTweight(s);
                newWeight = curweight + max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights_wt(s, newWeight, 0);
                }else{
                    kt->edit_tweights_wt(s, 0, -newWeight);
                }
            }
        }
        needs_recompute = true;
        source = s;
        last_modification = g.getCurrentHistory() - 1;
        flow_needs_recalc = true;
        updateSourceSink();
    }

    void updateSourceSink(){
        static_maxflow = -1;
        if(source != sink){
            source_max = 0;
            for(int i = 0; i < g.nIncident(getSource()); i++){
                int edgeID = g.incident(getSource(), i).id;
                int t = g.incident(getSource(), i).node;
                if(t != getSource()){
                    Weight w = g.getWeight(edgeID);
                    assert(w >= 0);
                    source_max += w;
                }
            }
            dest_max = 0;
            for(int i = 0; i < g.nIncoming(getSink()); i++){
                int edgeID = g.incoming(getSink(), i).id;
                int t = g.incoming(getSink(), i).node;
                if(t != getSink()){
                    Weight w = g.getWeight(edgeID);
                    assert(w >= 0);
                    dest_max += w;
                }
            }
            if(source_max < dest_max){
                static_maxflow = source_max;
            }else{
                static_maxflow = dest_max;
            }
        }
    }

    void setSink(int t) override{
        same_source_sink = false;
        assert(source != sink);
        if(sink == t){
            return;
        }
        if(t == source){
            same_source_sink = true;
            return;
        }
        if(kt){
            if(dynamic){
                Weight curweight = -kt->getTweight(sink);
                Weight newWeight = curweight - max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights(sink, 0, newWeight);
                }else{
                    kt->edit_tweights(sink, -newWeight, 0);
                }
                curweight = -kt->getTweight(t);
                newWeight = curweight + max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights(t, 0, newWeight);
                }else{
                    kt->edit_tweights(t, -newWeight, 0);
                }
            }else{
                Weight curweight = -kt->getTweight(sink);
                Weight newWeight = curweight - max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights_wt(sink, 0, newWeight);
                }else{
                    kt->edit_tweights_wt(sink, -newWeight, 0);
                }
                //if(!backward_maxflow){
                curweight = -kt->getTweight(t);
                newWeight = curweight + max_capacity;
                if(newWeight >= 0){
                    kt->edit_tweights_wt(t, 0, newWeight);
                }else{
                    kt->edit_tweights_wt(t, -newWeight, 0);
                }

            }
        }
        sink = t;
        needs_recompute = true;
        last_modification = std::min(last_modification, g.getCurrentHistory() - 1);
        flow_needs_recalc = true;
        updateSourceSink();
    }

    void setCapacity(int u, int w, Weight c){

    }

    void setAllEdgeCapacities(Weight c){

    }

    int64_t num_updates = 0;

    int numUpdates() const override{
        return num_updates;
    }

    void updateHistory() override{
        update();
    }


    void initKT(){
        if(same_source_sink)
            return;
        stats_inits++;
        needs_recompute = true;
        clearChangedEdges();
        edge_enabled.clear();
        flow_needs_recalc = true;
        if(kt){
            delete (kt);
            kt = nullptr;
        }
        if(!kt){
            kt = new kohli_torr::Graph<Weight, Weight, Weight>(g.nodes(), g.edges());
            kt->preserve_backward_order = kt_preserve_order;
            kt->maxflow(false); //just to initialize things.
        }

        while(kt->get_node_num() < g.nodes()){
            int node_id = kt->add_node();
            assert(node_id == kt->get_node_num() - 1);
        }
        max_capacity = 0;
        sum_of_edge_capacities = 0;
        local_weights.clear();
        local_weights.resize(g.edges(), 0);
        edge_enabled.resize(g.edges(), false);
        //multi_edges.clear();//fix this
        //multi_edges.resize(g.edges());

        arc_map.clear();
        arc_map.resize(g.edges(), -1);
        edge_map_ind.clear();
        edge_map_ind.resize(kt->get_arc_num(), {-1, -1});

        edge_map.clear();//.resize(kt->get_arc_num());
        for(int i = 0; i < tmp_edge_map.size(); i++)
            tmp_edge_map[i].clear();


        tmp_edges.clear();
        tmp_edges.resize(g.nodes(), -1);
        for(int n = 0; n < g.nodes(); n++){
            int first_arc_id = kt->get_arc_num();
            for(int i = 0; i < g.nIncident(n, false); i++){

                int edgeID = g.incident(n, i, false).id;

                if(g.selfLoop(edgeID))
                    continue;

                if(arc_map[edgeID] != -1){
                    continue;
                }

                int from = g.getEdge(edgeID).from;
                int to = g.getEdge(edgeID).to;
                assert(from == n);
                assert(from != to);

                max_capacity += g.getWeight(edgeID);
                edge_enabled[edgeID] = false;

                if(tmp_edges[to] > -1){
                    int arc_id = tmp_edges[to];
                    assert(g.getEdge(tmp_edge_map[arc_id - first_arc_id].back()).to == to);
                    assert(g.getEdge(tmp_edge_map[arc_id - first_arc_id].back()).from == from);
                    tmp_edge_map[arc_id - first_arc_id].push_back(edgeID);
                    arc_map[edgeID] = arc_id;
                }else{

                    assert(!kt->has_edge(from, to));
                    assert(arc_map[edgeID] == -1);
                    int arc_id = kt->add_edge(from, to, 0, 0);

                    tmp_edges[to] = arc_id;
                    if(tmp_edge_map.size() < arc_id - first_arc_id + 2)
                        tmp_edge_map.resize((arc_id - first_arc_id) + 2);
                    tmp_edge_map[arc_id - first_arc_id].push_back(edgeID);
                    arc_map[edgeID] = arc_id;
                }
            }
            //set the corresponding arc for each other  from-to edge
            for(int i = 0; i < g.nIncoming(n, false); i++){
                int edgeID = g.incoming(n, i, false).id;
                if(g.selfLoop(edgeID))
                    continue;
                int from = g.getEdge(edgeID).from;
                int to = g.getEdge(edgeID).to;
                assert(to == n);

                if(tmp_edges[from] != -1 && arc_map[edgeID] == -1){
                    //this is a backward edge corresponding to a forward edge
                    int arc_id = tmp_edges[from] + 1;// the reverse arc is always stored right after the forward arc
                    assert(arc_id > -1);
                    tmp_edge_map[arc_id - first_arc_id].push_back(edgeID);
                    arc_map[edgeID] = arc_id;

                    max_capacity += g.getWeight(edgeID);
                    edge_enabled[edgeID] = false;
                }
            }
            edge_map_ind.resize(kt->get_arc_num(), {-1, -1});
            for(int i = 0; i < tmp_edge_map.size(); i++){
                int arc_id = i + first_arc_id;
                if(arc_id >= kt->get_arc_num())
                    break;
                int first_ind = edge_map.size();
                for(int edgeID:tmp_edge_map[i]){
                    edge_map.push_back(edgeID);
                }
                int last_ind = edge_map.size();
                assert(arc_id < kt->get_arc_num());
                edge_map_ind[arc_id] = {first_ind, last_ind};


#ifdef DEBUG_DGL
                for(int edgeID:tmp_edge_map[i]){
                    int arcID=arc_map[edgeID];
                    assert(arcID==arc_id);
                    bool found=false;
                    for (auto  it = getArcEdges(arcID);it!=getArcEdgesEnd(arcID);++it){
                        int eID = *it;
                        if(eID==edgeID){
                            found=true;
                            break;
                        }
                    }
                    assert(found);
                }
#endif
                tmp_edge_map[i].clear();
            }

            for(int i = 0; i < g.nIncident(n, false); i++){
                int edgeID = g.incident(n, i, false).id;
                int to = g.getEdge(edgeID).to;
                tmp_edges[to] = -1;
            }

#ifdef DEBUG_DGL
            for(int i = 0;i<tmp_edges.size();i++)
                assert(tmp_edges[i]==-1);
#endif

        }

#ifdef DEBUG_DGL
        for (int edgeID = 0; edgeID < g.edges(); edgeID++) {
                    if (!g.hasEdge(edgeID) || g.selfLoop(edgeID))
                            continue;
                    int arcID = arc_map[edgeID];
                    assert(arcID>-1);
                    bool found=false;
                    for (auto  it = getArcEdges(arcID);it!=getArcEdgesEnd(arcID);++it){
                        int eID = *it;
                        if(eID==edgeID){
                            found=true;
                            break;
                        }
                    }
                    assert(found);
                }
#endif

        for(int edgeID = 0; edgeID < g.edges(); edgeID++){
            if(!g.hasEdge(edgeID) || g.selfLoop(edgeID))
                continue;
            assert(arc_map[edgeID] > -1);
            if(g.edgeEnabled(edgeID)){
                int from = g.getEdge(edgeID).from;
                int to = g.getEdge(edgeID).to;
                edge_enabled[edgeID] = true;
                set_local_weight(edgeID, g.getWeight(edgeID));
                kt->edit_edge_inc(from, to, g.getWeight(edgeID), 0, getArc(edgeID));
            }else{
                set_local_weight(edgeID, 0);
            }
        }
        int s = source;
        int t = sink;
        if(source != sink){
            assert(max_capacity >= 0);
            if(max_capacity < 1){
                max_capacity = 1;//if max capacity is 0, then kt allows the cut to be between the actual kt source and s, which can lead to incorrect results here.
            }
            if(dynamic){
                //if(!backward_maxflow){
                kt->edit_tweights(s, max_capacity, 0);
                kt->edit_tweights(t, 0, max_capacity);
                /*	}else{
                     kt->edit_tweights(t,max_capacity,0);
                     kt->edit_tweights(s,0,max_capacity);
                     }*/
            }else{
                //if(!backward_maxflow){
                kt->edit_tweights_wt(s, max_capacity, 0);
                kt->edit_tweights_wt(t, 0, max_capacity);
                /*			}else{
                         kt->edit_tweights_wt(t,max_capacity,0);
                         kt->edit_tweights_wt(s,0,max_capacity);
                         }*/
            }
        }
        updateSourceSink();
    }

    bool max_capacity_increased = true;

    const Weight update() override{
        int s = source;
        int t = sink;
        if(same_source_sink)
            return INF;
        //see http://cstheory.stackexchange.com/a/10186

        if(g.outfile()){
            fprintf(g.outfile(), "f %d %d\n", s, t);
            fflush(g.outfile());
        }

        //C.resize(g.nodes());
        bool any_changed = false;
        if(!needs_recompute && last_modification > 0 && g.getCurrentHistory() == last_modification){
#ifdef DEBUG_MAXFLOW2
            EdmondsKarpDynamic<Weight> ek(g,source,sink);
            Weight expected_flow =ek.maxFlow(source,sink);
            assert(curflow==expected_flow);
            bassert(curflow == expected_flow);
#endif
            return curflow;
        }else if(!kt || last_modification <= 0 || kt->get_node_num() != g.nodes()
                 || edge_enabled.size() != g.edges()){
            initKT();
        }else if(!g.changed() && last_history_clear >= 0 && last_history_clear == g.nHistoryClears() - 1 &&
                 history_qhead == g.getPreviousHistorySize()){
            //no information was lost in the history clear
            history_qhead = 0;
            last_history_clear = g.nHistoryClears();
        }else if(g.nHistoryClears() != last_history_clear || g.changed()){
            stats_reinits++;

            for(int edgeid = 0; edgeid < g.edges(); edgeid++){
                if(!g.hasEdge(edgeid) || g.selfLoop(edgeid))
                    continue;
                int from = g.getEdge(edgeid).from;
                int to = g.getEdge(edgeid).to;
                if(g.edgeEnabled(edgeid) && !edge_enabled[edgeid]){


                    assert(local_weight(edgeid) == 0);
                    edge_enabled[edgeid] = true;
                    set_local_weight(edgeid, g.getWeight(edgeid));

                    kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, g.getWeight(edgeid), 0,
                                      getArc(edgeid));
                    if(curflow < static_maxflow){

                        needs_recompute = true;
                    }
                }else if(g.edgeEnabled(edgeid) && edge_enabled[edgeid] && g.getWeight(edgeid) != local_weight(edgeid)){
                    Weight dif = g.getWeight(edgeid) - local_weight(edgeid);
                    typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id arcID = getArc(edgeid);

                    if(dif < 0 && kt->get_flow(arcID) > 0){
                        flow_needs_recalc = true;
                        needs_recompute = true;
                    }


                    kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, dif, 0, getArc(edgeid));
                    set_local_weight(edgeid, g.getWeight(edgeid));
                    if(dif > 0 && curflow < static_maxflow){

                        needs_recompute = true;
                    }
                }else if(!g.edgeEnabled(edgeid) && edge_enabled[edgeid]){
                    typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id arcID = getArc(edgeid);
                    assert(edge_enabled[edgeid]);
                    edge_enabled[edgeid] = false;
                    if(kt->get_flow(arcID) > 0){

                        needs_recompute = true;
                    }
                    kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, -local_weight(edgeid), 0,
                                      getArc(edgeid));
                    set_local_weight(edgeid, 0);

                }
            }
            history_qhead = g.historySize();
        }

        assert(kt);

        for(int i = history_qhead; i < g.historySize(); i++){
            int edgeid = g.getChange(i).id;
            int from = g.getEdge(edgeid).from;
            int to = g.getEdge(edgeid).to;
            if(g.selfLoop(edgeid))
                continue; //skip self loops
            if(g.getChange(i).addition && g.edgeEnabled(edgeid) && !edge_enabled[edgeid]){


                assert(local_weight(edgeid) == 0);
                edge_enabled[edgeid] = true;
                set_local_weight(edgeid, g.getWeight(edgeid));

                kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, g.getWeight(edgeid), 0, getArc(edgeid));
                if(curflow < static_maxflow){

                    needs_recompute = true;
                }
            }else if((g.getChange(i).weight_increase || g.getChange(i).weight_decrease) && g.edgeEnabled(edgeid) &&
                     edge_enabled[edgeid] && g.getWeight(edgeid) != local_weight(edgeid)){
                typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id arcID = getArc(edgeid);
                Weight dif = g.getWeight(edgeid) - local_weight(edgeid);

                if(dif < 0 && kt->get_flow(arcID) > 0){

                    needs_recompute = true;
                }
                kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, dif, 0, getArc(edgeid));
                set_local_weight(edgeid, g.getWeight(edgeid));
                if(dif > 0 && curflow < static_maxflow){

                    needs_recompute = true;
                }
            }else if(g.getChange(i).deletion && !g.edgeEnabled(edgeid) && edge_enabled[edgeid]){
                typename kohli_torr::Graph<Weight, Weight, Weight>::arc_id arcID = getArc(edgeid);
                assert(edge_enabled[edgeid]);
                edge_enabled[edgeid] = false;

                if(kt->get_flow(arcID) > 0){

                    needs_recompute = true;
                }
                kt->edit_edge_inc(g.getEdge(edgeid).from, g.getEdge(edgeid).to, -local_weight(edgeid), 0, arcID);
                set_local_weight(edgeid, 0);
            }
        }
        if(sum_of_edge_capacities > max_capacity || max_capacity_increased){
            updateMaxCapacity(sum_of_edge_capacities);
        }
        if(needs_recompute){
            f = kt->maxflow(dynamic);
            needs_recompute = false;
        }
#ifdef DEBUG_MAXFLOW2
        EdmondsKarpDynamic<Weight> ek(g,source,sink);
        Weight expected_flow =ek.maxFlow(source,sink);
        bassert(f == expected_flow);
#endif
        //dbg_print_graph(s,t,true);
        //g.drawFull(true);
#ifdef DEBUG_DGL

        for (int i = 0; i < g.edges(); i++) {
            if (g.getEdge(i).from == g.getEdge(i).to)
                continue; //skip self edges.
            assert(edge_enabled[i] == g.edgeEnabled(i));
        }

#endif

        dbg_check_flow(s, t);
        curflow = f;
        num_updates++;
        last_modification = g.getCurrentHistory();
        last_deletion = g.nDeletions();
        last_addition = g.nAdditions();

        history_qhead = g.historySize();
        g.updateAlgorithmHistory(this, alg_id, history_qhead);
        last_history_clear = g.nHistoryClears();

        if(flow_graph){
            getChangedEdges();//keep the flow graph up to date, if it is defined
        }

        return f;
    }

    std::vector<int>& getChangedEdges() override{
        calc_flow();
        while(kt->changed_edges.size()){
            int arc = kt->changed_edges.back();
            kt->changed_edges.pop_back();
            kt->unmarkFlowEdge(kt->get_arc(arc));
            for(auto it = getArcEdges(arc); it != getArcEdgesEnd(arc); ++it){
                int edgeID = *it;
                assert(arc_map[edgeID] == arc);
                //for (int edgeID =  : edge_map[arc]) {
                changed_edges.push_back(edgeID);
                if(flow_graph){
                    Weight f = getEdgeFlow(edgeID);
                    if(f > 0){
                        flow_graph->enableEdge(edgeID);
                        flow_graph->setEdgeWeight(edgeID, f);
                    }else{
                        flow_graph->disableEdge(edgeID);
                        flow_graph->setEdgeWeight(edgeID, 0);
                    }
                }
            }
        }
        return changed_edges;
    }

    void clearChangedEdges() override{
        changed_edges.clear();
    }

    std::vector<int>& getChangedPartition() override{
        while(kt->nodes_partition_changed.size()){
            int node = kt->nodes_partition_changed.back();
            kt->nodes_partition_changed.pop_back();
            changed_partition.push_back(node);
        }
        return changed_partition;
    }

    void clearChangedPartition() override{
        while(changed_partition.size()){
            int node_id = changed_partition.back();
            changed_partition.pop_back();
            kt->unmarkPartitionAssignment(node_id);
        }

    }

private:
    inline Weight& local_weight(int edgeid){
        return local_weights[edgeid];
    }

    inline void set_local_weight(int edgeid, Weight w){
        flow_needs_recalc = true;
        Weight dif = w - local_weight(edgeid);
        Weight old_max_cap = sum_of_edge_capacities;
        if(dif != 0){
            sum_of_edge_capacities -= local_weights[edgeid];
            local_weights[edgeid] = w;

            sum_of_edge_capacities += w;

            if(sum_of_edge_capacities <
               0){//this is a total hack... how should overflows be dealt with, properly, here?
                sum_of_edge_capacities = std::numeric_limits<Weight>::max() / 2;
            }
            int from = g.getEdge(edgeid).from;
            int to = g.getEdge(edgeid).to;
            if(from != to){
                if(from == getSource()){
                    source_max += dif;


                    if(source_max < dest_max){
                        if(static_maxflow != source_max){
                            static_maxflow = source_max;
                            //needs_recompute = true;
                        }
                    }else{
                        if(static_maxflow != dest_max){
                            static_maxflow = dest_max;
                            //needs_recompute = true;
                        }
                    }

                }else if(to == getSink()){
                    dest_max += dif;
                    if(source_max < dest_max){
                        if(static_maxflow != source_max){
                            static_maxflow = source_max;
                            //needs_recompute = true;
                        }
                    }else{
                        if(static_maxflow != dest_max){
                            static_maxflow = dest_max;
                            //needs_recompute = true;
                        }
                    }
                }
            }
        }
        if(sum_of_edge_capacities > old_max_cap){
            max_capacity_increased = true;
        }
    }

/*	inline void collect_multi_edges(int for_edge) {
		//this has to go - it is too slow!
		if (multi_edges[for_edge].size() == 0) {
			int from = g.getEdge(for_edge).from;
			int to = g.getEdge(for_edge).to;
			//for(int i = g.nIncident(from,true)-1;i>=0;i--){
			for (int i = 0; i < g.nIncident(from, true); i++) {
				int edgeid = g.incident(from, i, true).id;
				if ((g.getEdge(edgeid).from == from && g.getEdge(edgeid).to == to)) {
					multi_edges[for_edge].push_back(edgeid);
				}
			}
		}
	}*/

    void dbg_print_graph(int from, int to, bool only_flow = false){
#ifdef DEBUG_DGL

        if (edge_enabled.size() < g.edges())
            return;
        printf("Graph %d\n", it);
        printf("digraph{\n");
        for (int i = 0; i < g.nodes(); i++) {
            if (i == from) {
                printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
            } else if (i == to) {
                printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
            }
            if (kt->what_segment(i, kohli_torr::Graph<Weight, Weight, Weight>::SOURCE)
                    == kohli_torr::Graph<Weight, Weight, Weight>::SOURCE) {
                printf("n%d [label=\"%d\", style=filled, fillcolor=blue]\n", i, i);
            } else if (kt->what_segment(i, kohli_torr::Graph<Weight, Weight, Weight>::SOURCE)
                    == kohli_torr::Graph<Weight, Weight, Weight>::SINK) {
                printf("n%d [label=\"%d\", style=filled, fillcolor=red]\n", i, i);
            } else
                printf("n%d\n", i);
        }
        printf("outer_source\n");
        printf("outer_sink\n");
        for (int i = 0; i < g.edges(); i++) {
            if (edge_enabled[i]) {
                auto & e = g.getEdge(i);
                const char * s = "black";
                if (only_flow && getEdgeFlow(e.id) == 0)
                    continue;
                std::cout << "n" << e.from << " -> n" << e.to << " [label=\"" << i << ": " << getEdgeFlow(e.id) << "/"
                        << g.getWeight(i) << "\" color=\"" << s << "\"]\n";
                //printf("n%d -> n%d [label=\"%d: %d/%d\",color=\"%s\"]\n", e.from,e.to, i, F[i],capacity[i] , s);
            }
        }
        for (int i = 0; i < g.nodes(); i++) {
            Weight w = kt->get_trcap(i);
            //Weight c = kt->t_cap(i);
            /**
             * 	// if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
             // otherwise         -tr_cap is residual capacity of the arc node->SINK
             */
            /*if(w<0){
             std::cout<<"n" << from <<" -> n" << i << " [label=\"" << w <<  "\" color=\"blue" <<"\"]\n";
             }else if(w>0){
             std::cout<<"n" << i <<" -> n" << to << " [label=\""<< w  <<  "\" color=\"red" << "\"]\n";
             }*/
            if (w < 0) {
                std::cout << "outer_source" << " -> n" << i << " [label=\"" << w << "\" color=\"blue" << "\"]\n";
            } else if (w > 0) {
                std::cout << "n" << i << " -> outer_sink " << " [label=\"" << w << "\" color=\"red" << "\"]\n";
            }
        }
        printf("}\n");
#endif
    }

    void bassert(bool c){
#ifdef DEBUG_DGL
        assert(c);
        if (!c) {
            throw std::logic_error("Assertion Fail!");
        }
#endif
    }


    void dbg_check_flow(int s, int t){
/*#ifdef DEBUG_DGL
		//check that the flow is legal
		for (int i = 0; i < g.edges(); i++) {
			Weight flow = getEdgeFlow(i);
			bassert(flow >= 0);
			bassert(flow <= g.getWeight(i));
			if (flow != 0) {
				bassert(g.edgeEnabled(i));
			}
		}

		for (int n = 0; n < g.nodes(); n++) {
			Weight flow_in = 0;
			Weight flow_out = 0;
			for (int i = 0; i < g.nIncident(n); i++) {
				int edge = g.incident(n, i).id;
				if (g.edgeEnabled(edge)) {
					Weight flow = getEdgeFlow(edge);
					bassert(flow >= 0);
					flow_out += flow;
				}
			}
			for (int i = 0; i < g.nIncoming(n); i++) {
				int edge = g.incoming(n, i).id;
				if (g.edgeEnabled(edge)) {
					Weight flow = getEdgeFlow(edge);
					bassert(flow >= 0);
					flow_in += flow;
				}
			}
			if (n == s) {
				//if(!backward_maxflow){
				Weight excess = flow_out-flow_in;
				bassert(excess==f);
				//bassert(flow_in == 0);//this doesn't have to be the case - their can be spurious flow loops...
				//bassert(flow_out == f);
				*//*}else{
				 bassert(flow_in==0);
				 bassert(flow_out==f);
				 }*//*
			} else if (n == t) {
				//if(!backward_maxflow){
				Weight excess = flow_in-flow_out;
				bassert(excess==f);
				//bassert(flow_out == 0);
				//bassert(flow_in == f);
				*//*}else{
				 bassert(flow_out==0);
				 bassert(flow_in==f);
				 }*//*
			} else {
				bassert(flow_in == flow_out);
			}

		}

#endif*/
    }

    std::vector<bool> seen;
    std::vector<bool> visited;

    inline void calc_flow(){
        if(same_source_sink)
            return;

        //if(needs_recompute){
        update();
        //}
        if(!flow_needs_recalc)
            return;

        //double startflowtime = Monosat::rtime(0);

        stats_flow_calcs++;


        Weight maxflow = kt->maxflow(true, nullptr);
        needs_recompute = false;
        flow_needs_recalc = false;
        kt->clear_t_edges(source, sink);

        dbg_check_flow(source, sink);
        //}

        assert(kt->maxflow(true, nullptr) == maxflow);
        //stats_calc_time+=  Monosat::rtime(0)-startcalctime;
        //stats_flow_time+=  startcalctime-startflowtime;
        //printf("flow calc time %f %f\n", stats_flow_time,stats_calc_time);
    }

    std::vector<int> Q;

public:

    const Weight minCut(std::vector<MaxFlowEdge>& cut) override{
        if(same_source_sink)
            return INF;
        Weight f = this->maxFlow();
        calc_flow();//it is critical to ensure that the flow assignment is up to date, or else the identified cut may be inaccurate
        int s = source;
        int t = sink;

        if(g.outfile()){
            fprintf(g.outfile(), "m %d %d\n", s, t);
            fflush(g.outfile());
        }

        cut.clear();
        //dbg_print_graph(s, t);
        /*   	if(f==0)
         return 0;*/

        auto SOURCE = kohli_torr::Graph<Weight, Weight, Weight>::SOURCE; // backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SINK : kohli_torr::Graph<Weight,Weight,Weight>::SOURCE;
        auto SINK = kohli_torr::Graph<Weight, Weight, Weight>::SINK; //backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SOURCE :  kohli_torr::Graph<Weight,Weight,Weight>::SINK;

        //KT allows for all the nodes to be source or sink, if the cut is placed at the final source->'outer source' or 'outer sink' -> sink edge.
        //Ideally, this shouldn't happen here, because those should have infinite weight edges.
        if(kt->what_segment(s, SOURCE) == SINK){
            throw std::logic_error("Error in mincut analysis");
        }else if(kt->what_segment(t, SOURCE) == SOURCE){
            throw std::logic_error("Error in mincut analysis");
        }

        for(int n = 0; n < g.nodes(); n++){
            auto t = kt->what_segment(n, SINK);
            if(t == SINK){
                //check to see if any neighbouring edges are in the source segment.
                for(int i = 0; i < g.nIncoming(n); i++){
                    //any edge that crosses from the source side into the sink side is on the cut.
                    auto& e = g.incoming(n, i);
                    if(g.edgeEnabled(e.id)){
                        auto segment = kt->what_segment(e.node, SINK);
                        if(segment == SOURCE){
                            Weight fe = g.getWeight(e.id);
                            //then this edge is on the cut
                            cut.push_back(MaxFlowEdge{n, e.node, e.id});
                        }
                    }
                }
            }
        }

        //kohli-torr is unusual in that it maintains the _mincut_, explicitly.

        /*calc_flow();
         cut.clear();
         Q.clear();
         Q.push_back(s);
         seen.clear();
         seen.resize(g.nodes());
         seen[s]=true;
         //
         //explore the residual graph
         for(int j = 0;j<Q.size();j++){
         int u = Q[j];

         for(int i = 0;i<g.nIncident(u);i++){
         if(!g.edgeEnabled(g.incident(u,i).id))
         continue;
         int v = g.incident(u,i).node;
         int id = g.incident(u,i).id;
         if(getEdgeResidualCapacity(id) == 0){
         cut.push_back(MaxFlowEdge{u,v,id});//potential element of the cut
         }else if(!seen[v]){
         Q.push_back(v);
         seen[v]=true;
         }
         }
         for(int i = 0;i<g.nIncoming(u);i++){
         if(!g.edgeEnabled(g.incoming(u,i).id))
         continue;
         int v = g.incoming(u,i).node;
         int id = g.incoming(u,i).id;
         if(getEdgeFlow(id) == 0){

         }else if(!seen[v]){
         Q.push_back(v);
         seen[v]=true;
         }
         }
         }
         //Now keep only the edges from a seen vertex to an unseen vertex
         int i, j = 0;
         for( i = 0;i<cut.size();i++){
         if(!seen[cut[i].v] && seen[cut[i].u]){
         cut[j++]=cut[i];
         }
         }
         cut.resize(j);
         */
#ifdef DEBUG_DGL
        Weight dbg_sum = 0;
        for (int i = 0; i < cut.size(); i++) {
            int id = cut[i].id;
            bassert(getEdgeFlow(id) == g.getWeight(id));
            dbg_sum += getEdgeFlow(id);
        }
        bassert(dbg_sum == f);
#endif

        return f;
    }

    const Weight getEdgeCapacity(int id) override{
        assert(g.edgeEnabled(id));
        return local_weight(id);
    }

public:
    const bool inSourcePartition(int node) override{
        if(same_source_sink)
            return node == source;
        update();
        auto SOURCE = kohli_torr::Graph<Weight, Weight, Weight>::SOURCE; // backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SINK : kohli_torr::Graph<Weight,Weight,Weight>::SOURCE;
        auto SINK = kohli_torr::Graph<Weight, Weight, Weight>::SINK; //backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SOURCE :  kohli_torr::Graph<Weight,Weight,Weight>::SINK;
        return kt->what_segment(node, SOURCE) == SOURCE;
    }

    const bool isOnCut(int edgeID) override{
        if(same_source_sink)
            return false;
        if(g.getEdge(edgeID).from == g.getEdge(edgeID).to)
            return false;        //self edges are never on the cut

        update();

        auto SOURCE = kohli_torr::Graph<Weight, Weight, Weight>::SOURCE; // backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SINK : kohli_torr::Graph<Weight,Weight,Weight>::SOURCE;
        auto SINK = kohli_torr::Graph<Weight, Weight, Weight>::SINK; //backward_maxflow? kohli_torr::Graph<Weight,Weight,Weight>::SOURCE :  kohli_torr::Graph<Weight,Weight,Weight>::SINK;

        int u = g.getEdge(edgeID).from;
        int v = g.getEdge(edgeID).to;
        return kt->what_segment(u, SOURCE) == SOURCE && kt->what_segment(v, SOURCE) == SINK;
    }

    const Weight getEdgeFlow(int flow_edge) override{
        if(same_source_sink)
            return 0;
        if(g.getEdge(flow_edge).from == g.getEdge(flow_edge).to)
            return 0;        //self edges have no flow.

        calc_flow();

        //we need to pick, possibly arbitrarily (but deterministically), which of the edges have flow
        int arc_id = arc_map[flow_edge];
        assert(arc_id >= 0);
        arc a;
        a = kt->get_arc(arc_id);

        Weight start_cap = kt->get_ecap(a);
        Weight end_cap = kt->get_rcap(a);
        Weight remaining_flow = kt->get_ecap(a) - kt->get_rcap(a);

        if(remaining_flow <= 0)
            return 0;

        //for (int edgeid : edge_map [arc_id]) {
        for(auto it = getArcEdges(arc_id); it != getArcEdgesEnd(arc_id); ++it){
            int edgeid = *it;
            assert(g.getEdge(edgeid).from == g.getEdge(flow_edge).from);
            assert(g.getEdge(edgeid).to == g.getEdge(flow_edge).to);

            if(g.edgeEnabled(edgeid)){
                assert(arc_map[edgeid] == arc_id);
                //assert(local_weight(edgeid)==g.getWeight(edgeid)); //these might not be in sync, if this is called during backtracking.
                Weight edge_cap = local_weight(edgeid);
                if(edgeid == flow_edge){
                    if(remaining_flow >= edge_cap)
                        return edge_cap;
                    else
                        return remaining_flow;
                }
                if(remaining_flow >= edge_cap){
                    remaining_flow -= edge_cap;
                }else{
                    assert(remaining_flow < edge_cap);
                    return 0;//all the flow goes into this edge, which is not the edge being queried, so the returned flow is 0.
                }
                if(remaining_flow <= 0){
                    return 0;
                }
            }
        }

        return 0;
    }

    const Weight getEdgeResidualCapacity(int id) override{
        return getEdgeCapacity(id) - getEdgeFlow(id);
    }
};
}
#else
//Simple wrapper around EdmondsKarpDynamic, to be used if GPL sources are not linked.
//Notably, in this case, the 'KohliTorr' class does NOT implement Kohli and Torr's algorithm at all.
namespace dgl {
template<typename Weight>
class KohliTorr : public EdmondsKarpDynamic<Weight> {
    static bool warning_issued;
public:
    double stats_calc_time = 0;
    double stats_flow_time = 0;
    double stats_init_time = 0;
    int64_t stats_flow_calcs = 0;
    int64_t stats_inits=0;
    int64_t stats_reinits=0;
    KohliTorr(Graph<Weight>& g,  int source, int sink,bool kt_preserve_order) :EdmondsKarpDynamic<Weight>(g,source,sink){

    }
    std::string getName() override {
        return "EdmondsKarpDynamicKT()";
    }
    const Weight update() override {
        if (!warning_issued){
            warning_issued=true;
            fprintf(stderr,"Warning: MonoSAT was built without GPL sources, which will greatly decrease performance of the solver on maximum flow predicates.\n");
        }
        return EdmondsKarpDynamic<Weight>::update();
    }
};
template<typename Weight>
bool KohliTorr<Weight>::warning_issued = false;
}
#endif

#endif

