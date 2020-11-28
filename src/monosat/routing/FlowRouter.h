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

#ifndef MONOSAT_FLOWROUTER_H
#define MONOSAT_FLOWROUTER_H

#include "monosat/graph/GraphTheory.h"
#include "monosat/graph/ReachDetector.h"
#include "monosat/mtl/Vec.h"
#include "monosat/mtl/IntMap.h"
#include "monosat/mtl/Heap.h"
#include "monosat/mtl/Alg.h"
#include "monosat/mtl/Rnd.h"
#include "monosat/utils/Options.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/core/Theory.h"
#include "monosat/core/TheorySolver.h"
#include "monosat/core/Config.h"
#include "monosat/dgl/RamalReps.h"

using namespace dgl;
namespace Monosat {
template<typename Weight>
class FlowRouter : public Theory, public MaxflowDetector<Weight>::FlowListener {
    IntSet<> edgevars;
    int theory_index = -1;
public:
    class OpportunisticHeuristic;

private:
    struct Net {
        //vec<int> members;
        vec<ReachDetector<Weight>*> detectors;
        vec<OpportunisticHeuristic*> heuristics;
        vec<Lit> reach_lits;
        //vec<int> dest_edges;
        vec<Lit> dest_edgelits;
        //int disconnectedEdge=-1;
        Lit disconnectedEdgeLit = lit_Undef;
        int cur_routing_dest = -1;
        int n_satisifed = 0;
        int netID = -1;

        int size() const{
            return reach_lits.size();
        }

        int nUnsat() const{
            return size() - n_satisifed;
        }
    };

    Solver* S;
    GraphTheorySolver<Weight>* g_theory = nullptr;
    Theory* maxFlowTheory = nullptr;
    MaxflowDetector<Weight>* maxflow_detector = nullptr;
    vec<Net> nets;

    int sourceNode;
    int destNode;
    Lit flowLit;
    int routerID = -1;
    vec<Lit> enabled_routing_lits;
    vec<Lit> disabled_routing_lits;
    vec<Lit> inner_conflict;
    int64_t stats_flow_decisions = 0;
    int64_t stats_flow_conflicts = 0;

    struct NetHeurisitc {
        //vec<Reach*> reaches;
        Reach* reach;
        int source;
        int outer_dest;
        IntSet<> destinations;

    };
    vec<NetHeurisitc> dest_sets;
    DynamicGraph<Weight> heuristic_graph;
    //UnweightedRamalReps<Weight> * reach=nullptr;
public:
    FlowRouter(Solver* S, GraphTheorySolver<Weight>* g, int sourceNode, int destNode, Lit maxflowLit);

    const char* getTheoryType() override{
        return "FlowRouter";
    }

    int getRouterID(){
        return routerID;
    }

    void addNet(Lit disabledEdge, vec<Lit>& dest_edges, vec<Lit>& reach_lits);


    void edgeFlowChange(int edgeID, const Weight& flow) override{
        if(flow > 0){
            heuristic_graph.enableEdge(edgeID);
        }else{
            heuristic_graph.disableEdge(edgeID);
        }
    }


    void printStats(int detailLevel = 0) override{
        printf("Flow Router:\n");
        printf("\tFlow conflicts: %" PRId64 "\n", stats_flow_conflicts);
        printf("\tFlow decisions: %" PRId64 "\n", stats_flow_decisions);
    }

    bool propagateTheory(vec<Lit>& conflict) override{
        return propagateTheory(conflict, false);
    }

    void drawGrid(DynamicGraph<Weight>& g, int net_to_draw = -1, bool only_path = false);

    bool propagateTheory(vec<Lit>& conflict, bool solve);

    bool solveTheory(vec<Lit>& conflict) override{
        return propagateTheory(conflict, true);
    }

    inline int getTheoryIndex() const override{
        return theory_index;
    }

    inline void setTheoryIndex(int id) override{
        theory_index = id;
    }

    inline void newDecisionLevel() override{

    }

    inline void backtrackUntil(int untilLevel) override{

    }

    inline int decisionLevel(){
        return S->decisionLevel();
    }

    inline void undecideTheory(Lit l) override{

    }

    void enqueueTheory(Lit l) override{

    }


    Lit decideTheory(CRef& decision_reason) override;

    class OpportunisticHeuristic : public Heuristic {
        FlowRouter* router;
        GraphTheorySolver<Weight>* outer;
        ReachDetector<Weight>* r;
        LSet to_decide;


        Lit reach_lit = lit_Undef;
        int last_over_modification = -1;
        int last_over_addition = -1;
        int last_over_deletion = -1;
        int over_history_qhead = 0;
        int last_over_history_clear = 0;

        int last_under_history_clear = 0;
        int under_history_qhead = 0;
        int last_under_modification = -1;
        int last_under_addition = -1;
        int last_under_deletion = -1;

        DynamicGraph<Weight>& g_over;
        DynamicGraph<Weight>& g_under;
        DynamicGraph<Weight>& g_h;
        IntSet<int> path_edges;
        bool path_is_cut = false;
        int dest_node = -1;

        int current_path_dest = -1;
    public:

        int netID;
        int dest;

        OpportunisticHeuristic(FlowRouter* router, int netNum, ReachDetector<Weight>* rd, int dest, Lit reach_lit) :
                router(router), outer(rd->outer), netID(netNum), dest(dest), r(rd),
                g_over(rd->g_over), g_under(rd->g_under), g_h(router->heuristic_graph), reach_lit(reach_lit){

        }


        void computePath(int to){
            path_is_cut = false;
            current_path_dest = to;
            Reach* reach = router->dest_sets[netID].reach;
            auto* over_reach = reach;
            auto* over_path = reach;

            assert(over_path);
            assert(over_reach);

            if(opt_graph_use_cache_for_decisions == 1){
                over_path->clearCache();
            }

            to_decide.clear();
            path_edges.clear();

            Lit l = reach_lit;
            assert (l != lit_Undef);

            assert(outer->getSolver()->value(l) == l_True);
            if(outer->getSolver()->value(l) == l_True){

                if(over_path->connected(to)){

                    to_decide.clear();


                    assert(over_path->connected(to));            //Else, we would already be in conflict
                    int p = to;
                    int last_edge = -1;
                    int last = to;
                    over_path->update();


                    while(p != r->source){
                        last = p;
                        assert(p != r->source);
                        last_edge = over_path->incomingEdge(p);
                        if(g_over.hasEdge(last_edge)){
                            path_edges.insert(last_edge);
                            Var edge_var = outer->getEdgeVar(last_edge);
                            if(outer->value(edge_var) == l_Undef){
                                to_decide.push(mkLit(edge_var, false));
                            }
                        }
                        int prev = over_path->previous(p);
                        p = prev;

                    }
                }

            }

        }

        bool needsRecompute(int to){
            if(to != current_path_dest)
                return true;
            if(outer->getSolver()->value(reach_lit) != l_False && path_is_cut){
                return true;
            }else if(outer->getSolver()->value(reach_lit) != l_True && !path_is_cut){
                return true;
            }


            //check if any edges on the path have been removed since the last update
            if(last_over_modification > 0 && g_over.modifications == last_over_modification &&
               last_under_modification > 0 && g_under.modifications == last_under_modification){
                return false;
            }
            if(last_over_modification <= 0 || g_over.changed() || last_under_modification <= 0 ||
               g_under.changed()){//Note for the future: there is probably room to improve this further.
                return true;
            }

            if(!path_edges.size()){
                return true;
            }

            if(last_over_history_clear != g_over.historyclears || last_under_history_clear != g_under.historyclears){
                over_history_qhead = g_over.historySize();
                last_over_history_clear = g_over.historyclears;
                under_history_qhead = g_under.historySize();
                last_under_history_clear = g_under.historyclears;
                to_decide.clear();
                for(int edgeID:path_edges){
                    if(!path_is_cut){
                        if(!g_over.edgeEnabled(edgeID))
                            return true;
                        Var edge_var = outer->getEdgeVar(edgeID);
                        Lit l = mkLit(edge_var, false);

                        if(outer->value(edge_var) == l_Undef){

                            to_decide.push(l);
                        }
                    }else{
                        if(g_under.edgeEnabled(edgeID))
                            return true;
                        Var edge_var = outer->getEdgeVar(edgeID);
                        Lit l = mkLit(edge_var, true);
                        if(outer->value(edge_var) == l_Undef){

                            to_decide.push(l);
                        }
                    }
                }
            }

            for(int i = over_history_qhead; i < g_over.historySize(); i++){
                int edgeid = g_over.getChange(i).id;

                if(g_over.getChange(i).addition && g_over.edgeEnabled(edgeid)){
                    if(!path_is_cut){

                    }else{
                        Var edge_var = outer->getEdgeVar(edgeid);
                        Lit l = mkLit(edge_var, true);
                        if(outer->value(edge_var) == l_Undef){

                            to_decide.push(mkLit(edge_var, true));
                        }
                    }
                }else if(!g_over.getChange(i).addition && !g_over.edgeEnabled(edgeid)){
                    if(!path_is_cut){
                        if(path_edges.has(edgeid)){
                            return true;
                        }
                    }else{

                    }
                }
            }

            for(int i = under_history_qhead; i < g_under.historySize(); i++){
                int edgeid = g_under.getChange(i).id;

                if(path_edges.has(edgeid)){
                    if(!g_under.getChange(i).addition && !g_under.edgeEnabled(edgeid)){
                        if(!path_is_cut){
                            Var edge_var = outer->getEdgeVar(edgeid);
                            Lit l = mkLit(edge_var, false);
                            if(outer->value(edge_var) == l_Undef){

                                to_decide.push(mkLit(edge_var, false));
                            }
                        }else{

                        }
                    }else if(g_under.getChange(i).addition && g_under.edgeEnabled(edgeid)){
                        if(!path_is_cut){

                        }else{
                            if(path_edges.has(edgeid)){
                                return true;
                            }
                        }
                    }
                }
            }

            last_over_modification = g_over.modifications;
            last_over_deletion = g_over.deletions;
            last_over_addition = g_over.additions;
            over_history_qhead = g_over.historySize();
            last_over_history_clear = g_over.historyclears;

            last_under_modification = g_under.modifications;
            last_under_deletion = g_under.deletions;
            last_under_addition = g_under.additions;
            under_history_qhead = g_under.historySize();
            last_under_history_clear = g_under.historyclears;
            return false;
        }


        Lit decideTheory(CRef& decision_reason) override{
            decision_reason = CRef_Undef;
            NetHeurisitc& n = router->dest_sets[netID];
            Reach* reach = router->dest_sets[netID].reach;

            /*
             * Possible routing heuristic policies:
             * 1) route if the current net happens to connect to the right destination node
             * 2) route if the current net happens to connect to the right destination node, and it is not the last unrouted net
             * 3) route only if _all_ nets happen to connect to all the right destination nodes
             * 4) route the shortest path as would have been chosen by the normal reach heuristic,
             * but after removing from the space of options any paths in the flow solution that connect valid targets for other nets             *
             */

            if((opt_flow_router_heuristic == 2 || opt_flow_router_heuristic == 5) && router->nets[netID].nUnsat() <= 1){
                //suppress flow based decisions if only 1 path is left to route.
                return lit_Undef;
            }
            if(opt_flow_router_heuristic == 3 || opt_flow_router_heuristic == 5){
                //check if _all_ of the paths happen to connect valid targets
                bool all_connected = true;
                for(int i = 0; i < router->nets.size(); i++){
                    Net& n2 = router->nets[i];
                    assert(n2.netID >= 0);
                    int current_dest = n2.cur_routing_dest;
                    if(current_dest >= 0){
                        Reach* reach = router->dest_sets[n2.netID].reach;
                        NetHeurisitc& n2b = router->dest_sets[n2.netID];
                        assert(reach);
                        if(!reach->connected(n2b.outer_dest)){
                            all_connected = false;
                            break;
                        }
                    }
                }
                if(!all_connected){
                    if(opt_verb >= 4){ //&& (path_edges.size() || to_decide.size())) {
                        // printf("r over graph:\n");
                        //drawGrid(g_over,dest);
                        /*printf("r under graph:\n");
                        drawGrid(g_under,dest);*/
                        printf("r flow graph:\n");
                        drawGrid(g_h, dest);
                    }
                    //suppress flow based decisions if any net is unsatisfied by the current flow
                    return lit_Undef;
                }
            }

            /*if(opt_verb>2){
                g_h.drawFull(false,true);
            }*/
            if(opt_verb >= 4){ //&& (path_edges.size() || to_decide.size())) {
                //printf("r over graph:\n");
                //drawGrid(g_over,dest);
                /*printf("r under graph:\n");
                drawGrid(g_under,dest);*/
                printf("r flow graph:\n");
                drawGrid(g_h, dest);
            }
            if(reach->connected(n.outer_dest)){
                //step back from outer dest
                int prev = reach->previous(n.outer_dest);
                assert(prev >= 0);
                if(n.destinations.has(prev)){
                    //this is a valid path, that we can extract out of the network flow solution

                    int to = prev;


                    if(outer->getSolver()->value(reach_lit) == l_Undef){
                        assert(false);
                        return lit_Undef;//if the reach lit is unassigned, do not make any decisions here
                    }

                    if(needsRecompute(to)){
                        r->stats_heuristic_recomputes++;
                        computePath(to);

                        last_over_modification = g_over.modifications;
                        last_over_deletion = g_over.deletions;
                        last_over_addition = g_over.additions;
                        over_history_qhead = g_over.historySize();
                        last_over_history_clear = g_over.historyclears;

                        last_under_modification = g_under.modifications;
                        last_under_deletion = g_under.deletions;
                        last_under_addition = g_under.additions;
                        under_history_qhead = g_under.historySize();
                        last_under_history_clear = g_under.historyclears;

                    }
                    if(opt_verb >= 4){ //&& (path_edges.size() || to_decide.size())) {
                        printf("r over graph:\n");
                        drawGrid(g_over);
                        printf("r under graph:\n");
                        drawGrid(g_under);
                    }
#ifdef DEBUG_GRAPH
                    for(Lit l:to_decide){
                assert(outer->value(l)!=l_False);
            }

            if(!path_is_cut){
                for(int edgeID:path_edges){
                    assert(g_over.edgeEnabled(edgeID));
                    Lit l = mkLit(outer->getEdgeVar(edgeID));
                    if(outer->value(l)==l_Undef){
                        assert(to_decide.contains(l));
                    }
                }

            }else{
                for(int edgeID:path_edges){
                    assert(!g_under.edgeEnabled(edgeID));
                    Lit l = mkLit(outer->getEdgeVar(edgeID));
                    if(outer->value(l)==l_Undef){
                        assert(to_decide.contains(~l));
                    }
                }
            }
#endif

                    if(to_decide.size()){//  && last_decision_status == over_path->numUpdates() the numUpdates() monitoring strategy doesn't work if the constraints always force edges to be assigned false when other edges
                        // are assigned true, as the over approx graph will always register as having been updated.
                        //instead, check the history of the dynamic graph to see if any of the edges on the path have been assigned false, and only recompute in that case.
                        while(to_decide.size()){
                            Lit l = to_decide.last();

                            if(outer->value(l) == l_Undef){
                                //stats_decide_time += rtime(2) - startdecidetime;
                                router->stats_flow_decisions++;
                                return l;
                            }else if(outer->value(l) == l_True){
                                //only pop a decision that was actually made
                                to_decide.pop();
                            }else if(outer->value(l) == l_False){
                                //is this even a reachable state? it probably shouldn't be.
                                to_decide.clear();
                                path_edges.clear();
                                //needs recompute!
                                return decideTheory(decision_reason);
                            }
                        }
                    }
                }
            }
            return lit_Undef;
        }

        void drawGrid(DynamicGraph<Weight>& g, int highlight_dest = -1, bool only_path = false){
            int width = 60;
            int height = 60;

            vec<int> startsX;
            vec<int> endsX;
            vec<int> startsY;
            vec<int> endsY;

            vec<std::pair<int, int>> starts;
            vec<std::pair<int, int>> ends;


            bool found_source = false;
            bool found_dest = false;
            printf("\n");
            for(int y = 0; y < height; y++){
                for(int x = 0; x < width; x++){
                    int f_node = (y * height + x) * 2 + 1;
                    std::pair<int, int> p = std::make_pair(x, y);

                    if(f_node == r->source || f_node - 1 == r->source){
                        //assert(starts.contains(p));
                        found_source = true;
                        if(f_node == highlight_dest || f_node - 1 == highlight_dest){
                            printf("!");
                        }else{
                            printf("*");
                        }
                    }else if(f_node == dest || f_node - 1 == dest){
                        found_dest = true;
                        // assert(ends.contains(p));
                        if(f_node == highlight_dest || f_node - 1 == highlight_dest){
                            printf("!");
                        }else{
                            printf("@");
                        }
                    }else{

                        if(starts.contains(p)){
                            int id = starts.indexOf(p) % 10;
                            printf("%d", id);
                        }else if(ends.contains(p)){
                            //printf("%%");
                            int id = ends.indexOf(p) % 10;
                            printf("%d", id);
                        }else{
                            printf("+");
                        }
                    }
                    if(x < width - 1){
                        int t_node = (y * height + x + 1) * 2;
                        //assert(g.hasEdge(f_node,t_node));
                        if(g.hasEdge(f_node, t_node) || g.hasEdge(t_node + 1, f_node - 1) ||
                           g.hasEdge(t_node, f_node) || g.hasEdge(f_node - 1, t_node + 1)){
                            if(only_path){
                                int edgea = g.getEdge(f_node, t_node);
                                int edgeb = g.getEdge(t_node, f_node);
                                if(path_edges.has(edgea) || path_edges.has(edgeb)){
                                    printf("-");
                                }else{
                                    printf(" ");
                                }
                            }else{
                                //int to_edge = g.getEdge(f_node,t_node);
                                int edgea = g.getEdge(f_node, t_node);
                                int edgeb = g.getEdge(t_node, f_node);
                                printf("-");
                                /*if((edgea>=0 && g.edgeEnabled(edgea)) || (edgeb>=0 &&  g.edgeEnabled(edgeb))) {
                                    printf("-");
                                }else{
                                    printf(" ");
                                }*/
                            }
                        }else{
                            printf(" ");
                        }
                    }
                }
                printf("\n");
                for(int x = 0; x < width; x++){
                    int f_node = (y * height + x) * 2 + 1;
                    if(y < height - 1){
                        int t_node = ((y + 1) * height + x) * 2;
                        //assert(g.hasEdge(f_node, t_node));
                        if(only_path){
                            int edgea = g.getEdge(f_node, t_node);
                            int edgeb = g.getEdge(t_node, f_node);

                            if(path_edges.has(edgea)){
                                assert(g.hasEdge(edgea));
                                printf("| ");
                            }else if(path_edges.has(edgeb)){
                                assert(g.hasEdge(edgeb));
                                printf("| ");
                            }else{
                                printf("  ");
                            }
                        }else if(g.hasEdge(f_node, t_node) || g.hasEdge(t_node + 1, f_node - 1) ||
                                 g.hasEdge(t_node, f_node) || g.hasEdge(f_node - 1, t_node + 1)){
                            int edgea = g.getEdge(f_node, t_node);
                            int edgeb = g.getEdge(t_node, f_node);
                            //if((edgea >=0 && g.edgeEnabled(edgea)) || (edgeb>=0 && g.edgeEnabled(edgeb))) {
                            printf("| ");
                            /*}else{
                                printf("  ");
                            }*/
                        }else{
                            printf("  ");
                        }
                    }
                }
                printf("\n");
            }
            printf("\n");
            assert(found_source);
            assert(found_dest);
        }


    };
};


template<typename Weight>
FlowRouter<Weight>::FlowRouter(Solver* S, GraphTheorySolver<Weight>* g, int sourceNode, int destNode, Lit maxflowLit):S(
        S), g_theory(g), sourceNode(sourceNode), destNode(destNode), flowLit(maxflowLit){

    S->addTheory(this);

    routerID = getTheoryIndex();
    maxFlowTheory = S->theories[S->getTheoryID(maxflowLit, 0)];
    //this is really a mess...
    //and should be done more safely in the future...
    //This code is recovering the reachability detector associated with each net member
    GraphTheorySolver<Weight>* g2 = (GraphTheorySolver<Weight>*) maxFlowTheory;
    Detector* d = g2->detectors[g2->getDetector(var(S->getTheoryLit(maxflowLit, maxFlowTheory)))];
    assert(d);
    assert(d->getName().compare("Max-flow Detector") == 0);
    MaxflowDetector<Weight>* rd = (MaxflowDetector<Weight>*) d;
    this->maxflow_detector = rd;

    if(opt_flow_router_heuristic){
        while(heuristic_graph.nodes() < maxflow_detector->g_over.nodes()){
            heuristic_graph.addNode();
        }
        for(int i = 0; i < maxflow_detector->g_over.edges(); i++){
            int from = maxflow_detector->g_over.getEdge(i).from;
            int to = maxflow_detector->g_over.getEdge(i).to;
            int edgeID = heuristic_graph.addEdge(from, to);
            heuristic_graph.setEdgeEnabled(edgeID, false);
        }

        rd->attachFlowListener(this);

    }
}

template<typename Weight>
void FlowRouter<Weight>::addNet(Lit disabledEdge, vec<Lit>& dest_edges, vec<Lit>& reach_lits){
    if(!S->okay())
        return;
    /*
     *   struct Net{
        //vec<int> members;
        ReachDetector<Weight> * det;
        vec<Lit> member_lits;
        //vec<int> dest_edges;
        vec<Lit> dest_edgelits;
        //int disconnectedEdge=-1;
        Lit disconnectedEdgeLit=lit_Undef;
    };
     */
    nets.push();
    dest_sets.push();
    int common_source = -1;
    nets.last().netID = nets.size() - 1;

    nets.last().disconnectedEdgeLit = disabledEdge;
    dest_edges.copyTo(nets.last().dest_edgelits);
    for(Lit l:dest_edges){
        S->setDecisionVar(var(l), false);
    }
    S->setDecisionVar(var(disabledEdge), false);

    if(opt_flow_router_heuristic > 0){
        int outer_dest = heuristic_graph.addNode();
        dest_sets.last().outer_dest = outer_dest;

    }

    for(Lit l:reach_lits){
        nets.last().reach_lits.push(l);

        assert(S->hasTheory(l));
        Theory* t = S->theories[S->getTheoryID(l, 0)];
        //this is really a mess...
        //and should be done more safely in the future...
        //This code is recovering the reachability detector associated with each net member
        GraphTheorySolver<Weight>* g = (GraphTheorySolver<Weight>*) t;
        Detector* d = g->detectors[g->getDetector(var(S->getTheoryLit(l, g)))];
        assert(d);
        assert(d->getName().compare("Reachability Detector") == 0);
        ReachDetector<Weight>* rd = (ReachDetector<Weight>*) d;
        nets.last().detectors.push(rd);

        int source = rd->source;
        if(common_source < 0){
            common_source = source;
        }else{
            assert(source == common_source);
        }
        int dest = rd->getReachNode(S->getTheoryLit(l, t));

        assert(dest >= 0);
        if(opt_flow_router_heuristic > 0){
            int outer_dest = dest_sets.last().outer_dest;
            heuristic_graph.addEdge(dest, outer_dest);
            dest_sets.last().destinations.insert(dest);
            auto* h = new OpportunisticHeuristic(this, nets.size() - 1, rd, dest, l);
            nets.last().heuristics.push(h);
            rd->attachSubHeuristic(h, dest);
        }

    }
    if(opt_flow_router_heuristic > 0){
        dest_sets.last().source = common_source;
        dest_sets.last().reach = new UnweightedRamalReps<Weight>(common_source, heuristic_graph);
    }


}

template<typename Weight>
Lit FlowRouter<Weight>::decideTheory(CRef& decision_reason){
#ifdef DEBUG_GRAPH
    for(int i = 0;i<maxflow_detector->g_over.edges();i++){
        Weight w = maxflow_detector->overapprox_conflict_detector->getEdgeFlow(i);
        if(w>0){
            assert(heuristic_graph.edgeEnabled(i));
        }else{
            assert(!heuristic_graph.edgeEnabled(i));
        }
    }
#endif
    return lit_Undef;
}

template<typename Weight>
bool FlowRouter<Weight>::propagateTheory(vec<Lit>& conflict, bool solve){
    //for each net to be routed, pick one unrouted endpoint (if any).
    //connect it to destination in g.
    if(!maxflow_detector->propagate(conflict)){
        return false;
    }
    if(opt_flow_router_policy == 0)
        return true;

    bool has_level_decision = false;
    int lev = S->decisionLevel();
    enabled_routing_lits.clear();
    disabled_routing_lits.clear();
    inner_conflict.clear();

    DynamicGraph<Weight>& g = g_theory->getOverApproximationGraph();

    /*
    * Possible net choice policies
    * 1) for each net, pick the first unrouted destination
    * 2) for each net, pick the heuristically hardest unrouted destination
    * 3) for each net, pick the unrouted destination that is next to be decided in the RUC heuristic
    * 4) for each net, pick the unrouted destination that is last to be decided in the RUC heuristic
    */

    //nets that are completely routed are connected directly to the dest node; if all are connected then do nothing
    for(int j = 0; j < nets.size(); j++){
        Net& net = nets[j];
        net.cur_routing_dest = -1;
        net.n_satisifed = 0;
        //check if any members of this detector are unrouted
        Lit unrouted = lit_Undef;
        int unrouted_n = -1;
        if(opt_flow_router_policy == 2 && net.dest_edgelits.size() > 1){
            int best_heuristic_value = -1;
            int best_heuristic_net = -1;
            for(int i = 0; i < net.dest_edgelits.size(); i++){
                ReachDetector<Weight>* r = net.detectors[i];
                Lit l = net.reach_lits[i];
                if(S->value(l) == l_True){
                    Lit edge_lit = net.dest_edgelits[i];
                    //find an endpoint that is not yet connected in the under approx graph, or arbitrarily use the last one if all of them are connected
                    if((!r->isConnected(S->getTheoryLit(l, g_theory), false))){
                        OpportunisticHeuristic* h = net.heuristics[i];

                        int heuristic_value = h->getParentHeuristic()->getHeuristicOrder();
                        assert(heuristic_value != best_heuristic_value);//all heuristic values should be unique
                        if(best_heuristic_value < 0 || heuristic_value < best_heuristic_value){
                            best_heuristic_net = i;
                            best_heuristic_value = heuristic_value;
                        }
                    }

                }
            }
            if(best_heuristic_net >= 0){
                unrouted = net.dest_edgelits[best_heuristic_net];
                unrouted_n = best_heuristic_net;
            }
        }
        for(int i = 0; i < net.dest_edgelits.size(); i++){

            ReachDetector<Weight>* r = net.detectors[i];
            Lit l = net.reach_lits[i];
            if(S->value(l) == l_True){
                Lit edge_lit = net.dest_edgelits[i];
                //find an endpoint that is not yet connected in the under approx graph, or arbitrarily use the last one if all of them are connected
                if((!r->isConnected(S->getTheoryLit(l, g_theory), false))){
                    if(unrouted == lit_Undef || unrouted == edge_lit){
                        unrouted = edge_lit;
                        unrouted_n = i;

                        //assert(S->value(net.disconnectedEdgeLit)!=l_True);
                        if(S->value(edge_lit) == l_Undef){
                            assert(S->value(edge_lit) != l_False);
                            enabled_routing_lits.push(edge_lit);
                            if(!has_level_decision){
                                has_level_decision = true;
                                S->newDecisionLevel();
                            }
                            S->enqueue((edge_lit), CRef_Undef);
                        }

                    }else{
                        if(i == net.dest_edgelits.size() - 1)
                            if(S->value(edge_lit) == l_Undef){
                                disabled_routing_lits.push(~edge_lit);
                                //g_theory->enqueue(~S->getTheoryLit(edge_lit),CRef_Undef);
                                if(!has_level_decision){
                                    has_level_decision = true;
                                    S->newDecisionLevel();
                                }
                                S->enqueue(~(edge_lit), CRef_Undef);
                            }
                    }
                }else{
                    net.n_satisifed++;
                    if(unrouted == lit_Undef && i == net.dest_edgelits.size() - 1){
                        //if all lits are routed, route an arbitrary lit.
                        unrouted = edge_lit;
                        unrouted_n = i;
                        assert(S->value(edge_lit) != l_False);
                        if(S->value(edge_lit) == l_Undef){
                            assert(S->value(edge_lit) != l_False);
                            enabled_routing_lits.push(edge_lit);
                            if(!has_level_decision){
                                has_level_decision = true;
                                S->newDecisionLevel();
                            }
                            S->enqueue((edge_lit), CRef_Undef);
                        }
                    }else{
                        if(S->value(edge_lit) == l_Undef){
                            disabled_routing_lits.push(~edge_lit);
                            //g_theory->enqueue(~S->getTheoryLit(edge_lit),CRef_Undef);
                            if(!has_level_decision){
                                has_level_decision = true;
                                S->newDecisionLevel();
                            }
                            S->enqueue(~(edge_lit), CRef_Undef);
                        }
                    }
                }
            }
        }


        assert(unrouted != lit_Undef);
        net.cur_routing_dest = unrouted_n;
        /*
        if(unrouted==lit_Undef){
            //connect source to dest directly
            //g.enableEdge(net.disconnectedEdge);
        	if(S->value(net.disconnectedEdgeLit)==l_Undef){
                if(!has_level_decision){
                    has_level_decision=true;
                    S->newDecisionLevel();
                }
        		S->enqueue((net.disconnectedEdgeLit),CRef_Undef);
        		enabled_routing_lits.push(net.disconnectedEdgeLit);
        	}else{
                assert(S->value(net.disconnectedEdgeLit)==l_True);
            }
            //routing_edges.push(net.disconnectedEdgeLit);
        }*//*else{
            //connect unrouted to dest
            //int edge_ID = net.dest_edges[unrouted_n];
            //g.enableEdge(edgeID);
            //Lit l = net.dest_edgelits[unrouted_n];
            //routing_edges.push(edge_ID);
        	if(S->value(net.disconnectedEdgeLit)==l_Undef){
                if(!has_level_decision){
                    has_level_decision=true;
                    S->newDecisionLevel();
                }
        		S->enqueue(~(net.disconnectedEdgeLit),CRef_Undef);
        		disabled_routing_lits.push(~net.disconnectedEdgeLit);
        	}else{
                assert(S->value(net.disconnectedEdgeLit)==l_True);
            }
        }*/
    }
    inner_conflict.clear();
    conflict.clear();
    //assert(lev>=1);
    if(!maxflow_detector->propagate(inner_conflict)){
        //if there is a conflict, remove any of the edge lits that were assigned here
        for(int i = 0; i < inner_conflict.size(); i++){
            Lit l = inner_conflict[i];
            assert(l != lit_Undef);
            if(edgevars.has(var(l))){
                //drop this literal
            }else if(S->level(var(l)) > lev){
                //drop this literal.
                //Note that in addition to edge literals assigned temporarily above, it may be possible for the solver to
                //assign other literals through bcp.
                //I _believe_ that any such literals are safe to remove from the learnt clause.
            }else{
                conflict.push(l);
            }
        }
        S->cancelUntil(lev);
        for(Lit l : enabled_routing_lits){
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            //S->unsafeUnassign(l);
            assert(S->value(l) == l_Undef);
        }

        for(Lit l : disabled_routing_lits){
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            //S->unsafeUnassign(l);
            assert(S->value(l) == l_Undef);
        }
        for(Lit l:conflict){
            assert(S->value(l) != l_Undef);
        }
        stats_flow_conflicts++;
        return false;
    }
    maxflow_detector->collectChangedEdges();//for edge collection here to keep the heuristics up to date
    /*if(opt_verb>=2){
        heuristic_graph.drawFull(false,true);
    }*/
    if(opt_verb > 2 && opt_flow_router_heuristic > 0){
        printf("Flow graph");
        drawGrid(heuristic_graph);
    }

    if(!solve){
        S->cancelUntil(lev);

        for(Lit l : enabled_routing_lits){
            //S->unsafeUnassign(l);
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            assert(S->value(l) == l_Undef);
        }

        for(Lit l : disabled_routing_lits){
            //S->unsafeUnassign(l);
            //g_theory->backtrackAssign(S->getTheoryLit(l));
            assert(S->value(l) == l_Undef);
        }
    }
    return true;

}

template<typename Weight>
void FlowRouter<Weight>::drawGrid(DynamicGraph<Weight>& g, int net_to_draw, bool only_path){
    int width = 60;
    int height = 60;

    vec<int> startsX;
    vec<int> endsX;
    vec<int> startsY;
    vec<int> endsY;

    vec<std::pair<int, int>> starts;
    vec<std::pair<int, int>> ends;

    vec<int> all_starts;
    vec<int> all_dests;
    for(int i = 0; i < nets.size(); i++){
        if(net_to_draw >= 0 && i != net_to_draw)
            continue;
        Net& net = nets[i];
        ReachDetector<Weight>* r = net.detectors[0];
        all_starts.push(r->source);
        for(int j = 0; j < net.reach_lits.size(); j++){
            Lit l = net.reach_lits[j];
            int node = net.detectors[j]->getNode(var(S->getTheoryLit(l, g_theory)));
            all_dests.push(node);
        }
    }

    bool found_source = false;
    bool found_dest = false;
    printf("\n");
    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++){
            int f_node = (y * height + x) * 2 + 1;
            std::pair<int, int> p = std::make_pair(x, y);
            if(all_starts.contains(f_node) || all_starts.contains(f_node - 1)){
                //assert(starts.contains(p));
                found_source = true;
                printf("*");
            }else if(all_dests.contains(f_node) || all_dests.contains(f_node - 1)){
                found_dest = true;
                // assert(ends.contains(p));
                printf("@");
            }else{

                if(starts.contains(p)){
                    int id = starts.indexOf(p) % 10;
                    printf("%d", id);
                }else if(ends.contains(p)){
                    //printf("%%");
                    int id = ends.indexOf(p) % 10;
                    printf("%d", id);
                }else{
                    printf("+");
                }
            }
            if(x < width - 1){
                int t_node = (y * height + x + 1) * 2;
                //assert(g.hasEdge(f_node,t_node));
                if(g.hasEdge(f_node, t_node) || g.hasEdge(t_node + 1, f_node - 1) || g.hasEdge(t_node, f_node) ||
                   g.hasEdge(f_node - 1, t_node + 1)){
                    if(only_path){
                        int edgea = g.getEdge(f_node, t_node);
                        int edgeb = g.getEdge(t_node, f_node);
                        if(g.edgeEnabled(edgea) || g.edgeEnabled(edgeb)){
                            printf("-");
                        }else{
                            printf(" ");
                        }
                        /*if(path_edges.has(edgea) || path_edges.has(edgeb)){
                            printf("-");
                        }else{
                            printf(" ");
                        }*/
                    }else{
                        //int to_edge = g.getEdge(f_node,t_node);
                        //if(g.edgeEnabled(to_edge)){
                        printf("-");
                    }
                }else{
                    printf(" ");
                }
            }
        }
        printf("\n");
        for(int x = 0; x < width; x++){
            int f_node = (y * height + x) * 2 + 1;
            if(y < height - 1){
                int t_node = ((y + 1) * height + x) * 2;
                //assert(g.hasEdge(f_node, t_node));
                if(only_path){
                    int edgea = g.getEdge(f_node, t_node);
                    int edgeb = g.getEdge(t_node, f_node);
                    if(g.edgeEnabled(edgea) || g.edgeEnabled(edgeb)){
                        printf("| ");
                    }else{
                        printf("  ");
                    }

                }else if(g.hasEdge(f_node, t_node) || g.hasEdge(t_node + 1, f_node - 1) || g.hasEdge(t_node, f_node) ||
                         g.hasEdge(f_node - 1, t_node + 1)){
                    //int to_edge = g.getEdge(f_node, t_node);
                    //if(g.edgeEnabled(to_edge)){
                    //int to_edge = g.getEdge(f_node,t_node);
                    //if(g.edgeEnabled(to_edge)){
                    printf("| ");
                }else{
                    printf("  ");
                }
            }
        }
        printf("\n");
    }
    printf("\n");
    //assert(found_source);
    //assert(found_dest);
}

};

#endif //MONOSAT_FLOWROUTER_H
