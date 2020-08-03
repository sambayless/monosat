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
#include "monosat/mtl/Vec.h"
#include "FSMGeneratorAcceptorDetector.h"
#include "FSMTheory.h"

using namespace Monosat;

FSMGeneratorAcceptorDetector::FSMGeneratorAcceptorDetector(int detectorID, FSMTheorySolver* outer, DynamicFSM& g_under,
                                                           DynamicFSM& g_over, DynamicFSM& acceptor_under,
                                                           DynamicFSM& acceptor_over, int gen_source, int accept_source,
                                                           double seed) :
        FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), acceptor_under(acceptor_under),
        acceptor_over(acceptor_over), gen_source(gen_source), accept_source(accept_source), rnd_seed(seed){

    if(!opt_fsm_as_graph){
        underReachStatus = new FSMGeneratorAcceptorDetector::AcceptStatus(*this, true);
        overReachStatus = new FSMGeneratorAcceptorDetector::AcceptStatus(*this, false);

        underapprox_detector = new NFALinearGeneratorAcceptor<FSMGeneratorAcceptorDetector::AcceptStatus>(g_under,
                                                                                                          acceptor_under,
                                                                                                          gen_source,
                                                                                                          accept_source,
                                                                                                          *underReachStatus);
        overapprox_detector = new NFALinearGeneratorAcceptor<FSMGeneratorAcceptorDetector::AcceptStatus>(g_over,
                                                                                                         acceptor_over,
                                                                                                         gen_source,
                                                                                                         accept_source,
                                                                                                         *overReachStatus);
        if(opt_fsm_negate_underapprox){
            inverted_overapprox_detector = overapprox_detector;
        }

        underprop_marker = outer->newReasonMarker(getID());
        overprop_marker = outer->newReasonMarker(getID());
        forcededge_marker = outer->newReasonMarker(getID());
        forced_positive_nondet_generator_transition_marker = outer->newReasonMarker(getID());
        deterministic_forcededge_marker = outer->newReasonMarker(getID());
        chokepoint_acceptor_transition_marker = outer->newReasonMarker(getID());

        cur_seen.growTo(acceptor_over.states());
        gen_cur_seen.growTo(g_over.states());

        next_seen.growTo(acceptor_over.states());
        gen_next_seen.growTo(g_over.states());
        seen_chars.growTo(g_over.outAlphabet() + 1);
    }else{
        if(graph == nullptr){
            graph = new GraphTheorySolver<int64_t>(outer->getSolver());
            //fully unroll the fsm you get by feeding the generator into the acceptor
            //(can improve on this by pruning some nodes if we know the final states already).
            constructAllPaths();
        }
    }
}

void FSMGeneratorAcceptorDetector::constructAllPaths(bool add_graph_symbols){
    assert(graph);
    //this assumes that the generator is LINEAR!
    if(!g_over.isLinear()){
        throw std::runtime_error("Only linear finite state machine generators are supported");
    }
    int gen_pos = 0;
    cur_seen.growTo(acceptor_over.states());
    gen_cur_seen.growTo(g_over.states());

    next_seen.growTo(acceptor_over.states());
    gen_next_seen.growTo(g_over.states());
    seen_chars.growTo(g_over.outAlphabet() + 1);

    nodes.growTo(g_over.states());

    for(int i = 0; i < g_over.states(); i++){
        for(int j = 0; j < acceptor_over.states(); j++){
            int n = graph->newNode();
            nodes[i].push(n);
            if(add_graph_symbols){
                std::stringstream ss;
                ss << "" << i << ":" << j;
                graph->setNodeName(n, ss.str().c_str());
            }
        }
    }


    assert(outer->decisionLevel() == 0);


    for(int s:cur){
        assert(cur_seen);
        cur_seen[s] = false;
    }
    cur.clear();
    assert(next.size() == 0);


    cur_seen[accept_source] = true;
    cur.push(accept_source);

    for(int s:gen_cur){
        assert(gen_cur_seen);
        gen_cur_seen[s] = false;
    }
    gen_cur.clear();
    assert(next.size() == 0);
    gen_cur_seen[gen_source] = true;
    gen_cur.push(gen_source);
    vec<Transition> chars;
    chars.clear();
    DynamicFSM& g = acceptor_over;
    DynamicFSM& gen = g_over;

    //initial emove pass:
    if(gen.emovesEnabled()){
        for(int i = 0; i < gen_cur.size(); i++){
            int s = gen_cur[i];
            for(int j = 0; j < gen.nIncident(s); j++){
                //now check if the label is active
                int edgeID = gen.incident(s, j).id;
                int to = gen.incident(s, j).node;
                int from = i;
                if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID, 0, 0)){
                    gen_cur_seen[to] = true;
                    gen_cur.push(to);

                }
                if(gen.transitionEnabled(edgeID, 0, 0)){
                    Var transitionVar = outer->getTransition(gen.getID(), edgeID, 0, 0).outerVar;
                    Var outerVar = outer->toSolver(outer->newAuxVar());
                    Lit edgeLit = mkLit(outerVar);
                    Lit innerLit = graph->newEdge(nodes[from][accept_source], nodes[to][accept_source], outerVar);
                    if(add_graph_symbols){
                        graph->setEdgeName(var(innerLit), "{}");
                    }
                    outer->makeEqualInSolver(mkLit(transitionVar), edgeLit);
                }
            }
        }
    }


    int gen_prev_state = gen_source;
    //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
    while(gen_pos < g_over.states()){
        int gen_state = 0;
        stepGeneratorForward(chars, seen_chars, gen_state);//get set of next strings

        if(chars.size() == 0){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
                    int to = g.incident(s, j).node;
                    int from = s;
                    if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                        cur_seen[to] = true;
                        cur.push(to);
                    }
                    if(g.transitionEnabled(edgeID, 0, 0)){
                        Var transitionVar = outer->getTransition(g.getID(), edgeID, 0, 0).outerVar;
                        Var outerVar = outer->toSolver(outer->newAuxVar());
                        Lit edgeLit = mkLit(outerVar);
                        Lit innerLit = graph->newEdge(nodes[gen_prev_state][from], nodes[gen_prev_state][to], outerVar);
                        outer->makeEqualInSolver(mkLit(transitionVar), edgeLit);
                        if(add_graph_symbols){
                            graph->setEdgeName(var(innerLit), "{}");
                        }
                    }
                }
            }
        }else{

            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
                    int to = g.incident(s, j).node;
                    int from = s;
                    if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                        cur_seen[to] = true;
                        cur.push(to);
                    }
                    if(g.transitionEnabled(edgeID, 0, 0)){
                        Var transitionVar = outer->getTransition(g.getID(), edgeID, 0, 0).outerVar;
                        Var outerVar = outer->toSolver(outer->newAuxVar());
                        Lit edgeLit = mkLit(outerVar);
                        Lit innerLit = graph->newEdge(nodes[gen_prev_state][from], nodes[gen_prev_state][to], outerVar);
                        outer->makeEqualInSolver(mkLit(transitionVar), edgeLit);
                        if(add_graph_symbols){
                            graph->setEdgeName(var(innerLit), "{}");
                        }
                    }
                    for(auto transition:chars){
                        //int genEdge = transition.edgeID;
                        int l = transition.out;
                        assert(l > 0);
                        if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0)){
                            next_seen[to] = true;
                            next.push(to);
                        }

                        if(g.transitionEnabled(edgeID, l, 0)){
                            Lit genLit = mkLit(outer->getTransition(gen.getID(), transition.edgeID, transition.in,
                                                                    transition.out).outerVar);
                            Lit acceptLit = mkLit(outer->getTransition(g.getID(), edgeID, l, 0).outerVar);
                            Var outerVar = outer->toSolver(outer->newAuxVar());
                            Lit edgeLit = mkLit(outerVar);
                            Var inner_edge = var(
                                    graph->newEdge(nodes[gen_prev_state][from], nodes[gen_state][to], outerVar));
                            //EdgeLit == (genVar && acceptVar)
                            outer->getSolver()->addClause(genLit, ~edgeLit);
                            outer->getSolver()->addClause(acceptLit, ~edgeLit);
                            outer->getSolver()->addClause(edgeLit, ~genLit, ~acceptLit);

                            if(add_graph_symbols){
                                std::stringstream ss;
                                ss << l;
                                graph->setEdgeName(inner_edge, ss.str().c_str());
                            }
                        }

                    }
                }
            }

        }

        next.swap(cur);
        next_seen.swap(cur_seen);

        for(int s:next){
            assert(next_seen[s]);
            next_seen[s] = false;
        }
        next.clear();
        if(chars.size() == 0){
            //must eventually happen because the generator is linear.
            break;
        }

        for(auto c :chars){
            assert(seen_chars[c.out]);
            seen_chars[c.out] = false;
        }
        chars.clear();
        gen_prev_state = gen_state;
        gen_pos++;
    }
}

void
FSMGeneratorAcceptorDetector::stepGeneratorForward(vec<Transition>& store, vec<bool>& store_seen, int& cur_gen_state){
    DynamicFSM& g = g_over;
    cur_gen_state = 0;

    for(int i = 0; i < gen_cur.size(); i++){
        int s = gen_cur[i];

        for(int j = 0; j < g.nIncident(s); j++){
            //now check if the label is active
            int edgeID = g.incident(s, j).id;
            int to = g.incident(s, j).node;

            if(g.transitionEnabled(edgeID, 0, 0)){
                if(!gen_cur_seen[to]){
                    gen_cur_seen[to] = true;
                    gen_cur.push(to);
                }
                cur_gen_state = to;
            }
            for(int l = 1; l < g.outAlphabet(); l++){
                if(g.transitionEnabled(edgeID, 0, l)){
                    if(!gen_next_seen[to]){
                        gen_next_seen[to] = true;
                        gen_next.push(to);
                    }
                    cur_gen_state = to;
                    if(!store_seen[l]){
                        store_seen[l] = true;
                        store.push({edgeID, 0, l});
                    }

                }
            }
        }
    }

    gen_next.swap(gen_cur);
    gen_next_seen.swap(gen_cur_seen);

    for(int s:gen_next){
        assert(gen_next_seen[s]);
        gen_next_seen[s] = false;
    }
    gen_next.clear();


}

void FSMGeneratorAcceptorDetector::preprocess(){

}

void FSMGeneratorAcceptorDetector::addAcceptLit(int generatorFinalState, int acceptorFinalState, Var outer_reach_var){
    if(gen_source < 0 || accept_source < 0
       || gen_source >= g_over.states()
       || accept_source >= acceptor_over.states()){
        throw std::runtime_error("Invalid fsm state");
    }

    if(generatorFinalState < 0 || acceptorFinalState < 0
       || generatorFinalState >= g_over.states()
       || acceptorFinalState >= acceptor_over.states()){
        throw std::runtime_error("Invalid fsm state");
    }

    if(opt_fsm_as_graph){
        assert(graph);
        if(gen_source < 0 || accept_source < 0
           || gen_source >= nodes.size()
           || accept_source >= nodes[gen_source].size()){
            throw std::runtime_error("Invalid fsm state for graph");
        }

        if(generatorFinalState < 0 || acceptorFinalState < 0
           || generatorFinalState >= nodes.size()
           || acceptorFinalState >= nodes[generatorFinalState].size()){
            throw std::runtime_error("Invalid fsm state for graph");
        }
        graph->reaches(nodes[gen_source][accept_source], nodes[generatorFinalState][acceptorFinalState],
                       outer_reach_var);

        return;
    }

    int litIndex = generatorFinalState + acceptorFinalState * g_over.states();

    lit_backward_map.growTo(g_over.states() * acceptor_over.states(), var_Undef);
    assert(litIndex < lit_backward_map.size());
    if(lit_backward_map[litIndex] != var_Undef){
        Var v = lit_backward_map[litIndex];
        outer->makeEqualInSolver(mkLit(outer->toSolver(v)), mkLit(outer_reach_var));
        return;
    }

    if(first_destination == -1)
        first_destination = generatorFinalState;

    g_under.invalidate();
    g_over.invalidate();

    Var accept_var = outer->newVar(outer_reach_var, getID());

    if(first_var == var_Undef){
        first_var = accept_var;
    }else{
        assert(accept_var >= first_var);
    }
    int index = accept_var - first_var;

    //is_changed.growTo(index+1);
    Lit acceptLit = mkLit(accept_var, false);
    all_lits.push(acceptLit);

    lit_backward_map[litIndex] = var(acceptLit);

    while(accept_lit_map.size() <= accept_var - first_var){
        accept_lit_map.push(-1);
    }
    accept_lit_map[accept_var - first_var] = all_accept_lits.size();
    all_accept_lits.push(AcceptLit(acceptLit, generatorFinalState, acceptorFinalState, nullptr));

    updatePrefixTable(generatorFinalState, acceptorFinalState);
}

void FSMGeneratorAcceptorDetector::addSuffixLit(Lit acceptorLit, int suffixStartState, int suffixEndState,
                                                Var outer_reach_var){
    AcceptLit& a = all_accept_lits[indexOf(var(acceptorLit))];
    if(!a.suffixLits)
        a.suffixLits = new vec<SuffixLit>();


    Var suffix_var = outer->newVar(outer_reach_var, getID());

    if(first_var == var_Undef){
        first_var = suffix_var;
    }else{
        assert(suffix_var >= first_var);
    }
    int index = suffix_var - first_var;

    Lit suffixLit = mkLit(suffix_var, false);

    a.suffixLits->push({suffixLit, acceptorLit, suffixStartState, suffixEndState});
}

void FSMGeneratorAcceptorDetector::AcceptStatus::accepts(int string, int state, int edgeID, int label, bool accepts){

}

bool FSMGeneratorAcceptorDetector::checkNegatedPolarity(){
    if(opt_detect_satisfied_predicates > 0 && outer->decisionLevel() == 0){
        return true;
    }else{
        return FSMDetector::checkNegatedPolarity();
    }
};

bool FSMGeneratorAcceptorDetector::propagate(vec<Lit>& conflict){


    bool skipped_positive = false;
    if(underapprox_detector && (!opt_detect_pure_theory_lits || unassigned_positives > 0)){

    }else{
        skipped_positive = true;
        //outer->stats_pure_skipped++;
        stats_skipped_under_updates++;
    }
    bool skipped_negative = false;
    if(overapprox_detector && (!opt_detect_pure_theory_lits || unassigned_negatives > 0)){

    }else{
        skipped_negative = true;
        stats_skipped_over_updates++;
    }

    if(skipped_negative && skipped_positive){
        return true;
    }

    bool global_check_negated_polarity = checkNegatedPolarity();


    for(auto& t:all_accept_lits){


        forced_edges.clear();
        chokepoint_edges.clear();
        Lit l = t.l;
        Lit lit = t.l;


        int gen_to = t.gen_to;
        int accept_to = t.accept_to;
        int failedSuffixLit = -1;


        bool check_negated_polarity = global_check_negated_polarity || (opt_detect_satisfied_predicates > 0 &&
                                                                        outer->level(var(l)) ==
                                                                        outer->decisionLevel() &&
                                                                        outer->decisionLevel() > 0);

        vec<SuffixLit>* suffixLits = t.suffixLits;

        pre_accepting_states.clear();

        if(outer->litIsRelevant(~l) && !opt_fsm_negate_underapprox &&
           (outer->value(l) != l_True || check_negated_polarity) &&
           underapprox_detector->accepts(gen_to, accept_to, false, opt_fsm_forced_edge_prop ? &forced_edges : nullptr)){
            if(outer->value(l) == l_True){
                assert(check_negated_polarity);
                outer->enqueueSat(l);
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, underprop_marker);
                outer->enqueueSat(l);
            }else{
                conflict.push(l);
                buildAcceptReason(gen_to, accept_to, conflict);
                return false;
            }
        }else if(outer->litIsRelevant(~l) && opt_fsm_negate_underapprox &&
                 (outer->value(l) != l_True || check_negated_polarity) &&
                 !inverted_overapprox_detector->accepts(gen_to, accept_to, true,
                                                        opt_fsm_forced_edge_prop ? &forced_edges : nullptr)){
            //This optimization appears to be broken. For example, if the input generator over-approximation has everything free, then this always triggers.
            //possibly the inverted overapprox acceptor has to be paired with the generator underapprox instead...
            assert(false);

            if(outer->value(l) == l_True){
                assert(check_negated_polarity);
                outer->enqueueSat(l);
            }else if(outer->value(l) == l_Undef){

                outer->enqueue(l, underprop_marker);
                outer->enqueueSat(l);
            }else{
                conflict.push(l);
                buildAcceptReason(gen_to, accept_to, conflict);
                return false;
            }
        }else if(outer->litIsRelevant(l) && (outer->value(l) != l_False || check_negated_polarity) &&
                 !overapprox_detector->accepts(gen_to, accept_to, false, opt_fsm_edge_prop ? &forced_edges : nullptr,
                                               opt_fsm_chokepoint_prop ? &chokepoint_edges : nullptr,
                                               suffixLits ? &pre_accepting_states : nullptr)){

            if(outer->value(l) == l_False){
                assert(check_negated_polarity);
                outer->enqueueSat(~l);
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(~l, overprop_marker);
                outer->enqueueSat(~l);
            }else{
                conflict.push(~l);
                buildNonAcceptReason(gen_to, accept_to, conflict);
                return false;
            }
        }else if(outer->value(l) != l_False && suffix_fsm && suffixLits && suffixLits->size()){
            //note: NOT attempting to propagate unassigned suffix lits here; it would be too expensive!
            //further, a suffix lit being assigned false doesn't mean that it cannot be accepted, but simply that it _may_ not be accepted.
            for(int i = 0; i < suffixLits->size(); i++){
                SuffixLit& t = (*suffixLits)[i];
                if(outer->value(t.l) == l_True && outer->value(t.accept_l) == l_True){

                    if(!acceptsSuffix(*suffix_fsm, acceptor_over, accept_to, t.suffix_from, t.suffix_to,
                                      pre_accepting_states)){


                        assert(!acceptsSuffix(*suffix_fsm, acceptor_over, accept_to, t.suffix_from, t.suffix_to,
                                              pre_accepting_states));
                        buildNonSuffixAcceptReason(gen_to, accept_to, pre_accepting_states, (*suffixLits)[i], conflict);
                        return false;
                    }
                }
            }

        }

        if(outer->litIsRelevant(~l) && outer->value(lit) == l_False){
            if(opt_fsm_forced_edge_prop){
                for(int s = 0; s < forced_edges.size(); s++){
                    for(auto& t:forced_edges[s]){

                        int generator_state = t.generator_state;
                        assert(generator_state == s);
                        int label = t.character;


                        for(int i = 0; i < g_over.nIncident(generator_state); i++){
                            int edgeID = g_over.incident(generator_state, i).id;
                            if(g_over.transitionEnabled(edgeID, 0, label)){
                                Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, label);
                                Lit f = mkLit(v, true);
                                //Var test = var(lit);
                                assert(outer->value(lit) == l_False);
                                //int lev = outer->level(var(lit));

                                if(outer->value(f) == l_Undef){
                                    stats_forced_edge_propagations++;
                                    setForcedVar(v, ~lit);
                                    outer->enqueue(f, forcededge_marker);
                                }else if(outer->value(f) == l_False){
                                    conflict.push(f);
                                    buildForcedEdgeReason(gen_to, accept_to, edgeID, label, conflict);
                                    return false;
                                }
                            }
                        }
                    }
                }
            }
        }else if(outer->litIsRelevant(l) && outer->value(lit) == l_True){

            if(opt_fsm_chokepoint_prop){

                //if a transition _must_ be traversed in order to accept the instance (a chokepoint),
                //then assert it true
                for(auto& t:chokepoint_edges){


                    int edgeID = t.acceptorEdgeID;
                    int from = acceptor_over.getEdge(edgeID).from;
                    int to = acceptor_over.getEdge(edgeID).to;
                    if(!acceptor_over.edgeEnabled(edgeID)){
                        throw std::runtime_error("Error in fsm chokepoint propagation");
                    }
                    int label = t.character;

                    Lit e = mkLit(outer->getTransitionVar(acceptor_over.getID(), edgeID, label, 0));
                    if(outer->value(e) == l_True){
                        //don't need to do anything
                        assert(acceptor_under.transitionEnabled(edgeID, label, 0));
                    }else if(outer->value(e) == l_False){
                        //shouldn't be possible
                        throw std::runtime_error("Error in fsm chokepoint propagation");
                    }else{
                        stats_chokepoint_edge_propagations++;
                        assert(outer->value(e) == l_Undef);
                        setForcedVar(var(e), lit);
                        outer->enqueue(e, chokepoint_acceptor_transition_marker);
                    }
                }

            }

            if(opt_fsm_edge_prop){
                //an edge in the generator can be asserted to be _false_ if the generator is known to be deterministic (in the final model, not now),
                //and the generator must follow some other set of edges in order to be accepted.
                //(alternatively, if it _must_ follow exactly one edge that edge can be asserted true, even if the generator is not deterministic).

                //to find the set of _transitions_ that the generator must pass through, run the acceptor backward first, and then forward, to get the
                //set of states that are reached on all paths to the accepting state. Then find the transitions between them that CANNOT lead to the accepting state.


                for(int s2 = 0; s2 < forced_edges.size(); s2++){
                    int n_forced = forced_edges[s2].size();
                    if(forced_edges[s2].size() == 0)
                        continue;
                    int generator_state = forced_edges[s2][0].generator_state;

                    int edgeID = -1;
                    for(int i = 0; i < g_over.nIncident(generator_state); i++){
                        edgeID = g_over.incident(generator_state, i).id;
                    }
                    if(edgeID < 0)
                        continue;

                    int remaining_transitions = 0;
                    for(int c = 0; c < g_over.outAlphabet(); c++){
                        if(g_over.transitionEnabled(edgeID, 0, c))
                            remaining_transitions++;
                    }
                    assert(n_forced < remaining_transitions);//else this would be a conflict
                    if(n_forced == remaining_transitions - 1){
                        tmp_seen_chars.clear();
                        tmp_seen_chars.growTo(g_over.outAlphabet(), false);



                        //whether or not the generator is deterministic, the remaining unforced edge must be assigned.
                        //find the unforced edge, and assign it
                        for(auto& t:forced_edges[s2]){
                            int generator_state2 = t.generator_state;
                            assert(generator_state2 == generator_state);
                            int label = t.character;
                            assert(!tmp_seen_chars[label]);
                            tmp_seen_chars[label] = true;
                        }
                        for(int c = 0; c < g_over.outAlphabet(); c++){
                            if(g_over.transitionEnabled(edgeID, 0, c) && !tmp_seen_chars[c]){
                                //this edge is forced to be true



                                Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, c);
                                Lit f = mkLit(v);
                                if(outer->value(f) == l_Undef){
                                    stats_nondet_gen_edge_propagations++;
                                    setForcedVar(v, lit);
                                    outer->enqueue(f, forced_positive_nondet_generator_transition_marker);
                                }
                                assert(outer->value(f) != l_False);//as that would have been a conflict above


                                break;
                            }
                        }
                    }else if(generator_is_deterministic){
                        //any transition that cannot lead to the accepting state must be disabled (if the generator is deterministic).
                        for(auto& t:forced_edges[s2]){


                            int generator_state = t.generator_state;

                            int label = t.character;

                            for(int i = 0; i < g_over.nIncident(generator_state); i++){
                                int edgeID = g_over.incident(generator_state, i).id;
                                if(g_over.transitionEnabled(edgeID, 0, label)){
                                    Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, label);
                                    Lit f = mkLit(v);
                                    if(outer->value(f) == l_Undef){
                                        setForcedVar(v, lit);
                                        stats_det_gen_edge_propagations++;
                                        outer->enqueue(~f, deterministic_forcededge_marker);
                                    }else if(outer->value(f) == l_True){
                                        conflict.push(~f);
                                        buildDeterministicForcedEdgeReason(gen_to, accept_to, edgeID, label, conflict);
                                        return false;
                                    }
                                }
                            }
                        }
                    }

                }
            }

        }
    }

    return true;
}


void FSMGeneratorAcceptorDetector::buildReason(Lit p, vec<Lit>& reason, CRef marker){
    if(marker == underprop_marker){
        reason.push(p);
        Var v = var(p);
        int gen_final = getGeneratorFinal(v);
        int accept_final = getAcceptorFinal(v);
        buildAcceptReason(gen_final, accept_final, reason);
    }else if(marker == overprop_marker){
        reason.push(p);
        Var v = var(p);
        int gen_final = getGeneratorFinal(v);
        int accept_final = getAcceptorFinal(v);
        buildNonAcceptReason(gen_final, accept_final, reason);
    }else if(marker == forcededge_marker){
        reason.push(p);
        Var v = var(p);

        Lit forL = getForcedVar(v);
        //lbool val = outer->value(forL);
        Var forV = var(forL);
        //lbool val2 = outer->value(forV);
        assert(outer->value(forL) == l_True);
        assert(outer->value(forV) == l_False);
        reason.push(~forL);
        int gen_final = getGeneratorFinal(forV);
        int accept_final = getAcceptorFinal(forV);
        int forcedEdge = outer->getEdgeID(v);
        int forcedLabel = outer->getOutput(v);
        buildForcedEdgeReason(gen_final, accept_final, forcedEdge, forcedLabel, reason);
    }else if(marker == chokepoint_acceptor_transition_marker){
        reason.push(p);
        Var v = var(p);
        Lit forL = getForcedVar(v);

        Var forV = var(forL);
        assert(outer->value(forL) == l_True);

        reason.push(~forL);
        int gen_final = getGeneratorFinal(forV);
        int accept_final = getAcceptorFinal(forV);
        int forcedEdge = outer->getEdgeID(v);
        int forcedLabel = outer->getInput(v);
        buildAcceptorChokepointEdgeReason(gen_final, accept_final, forcedEdge, forcedLabel, reason);
    }else if(marker == deterministic_forcededge_marker){
        assert(false);
    }else if(marker == forced_positive_nondet_generator_transition_marker){
        reason.push(p);
        Var v = var(p);
        Lit forL = getForcedVar(v);

        Var forV = var(forL);
        assert(outer->value(forL) == l_True);
        //assert(outer->value(forV)==l_False);
        reason.push(~forL);
        int gen_final = getGeneratorFinal(forV);
        int accept_final = getAcceptorFinal(forV);
        int forcedEdge = outer->getEdgeID(v);
        int forcedLabel = outer->getOutput(v);
        buildForcedNondetEdgeReason(gen_final, accept_final, forcedEdge, forcedLabel, reason);
    }else{
        assert(false);
    }
}

void FSMGeneratorAcceptorDetector::buildAcceptReason(int genFinal, int acceptFinal, vec<Lit>& conflict){

    //find a path - ideally, the one that traverses the fewest unique transitions - from source to node, learn that one of the transitions on that path must be disabled.
    if(!opt_fsm_negate_underapprox){
        static vec<NFATransition> generator_path;
        static vec<NFATransition> acceptor_path;
        generator_path.clear();
        underapprox_detector->getGeneratorAceptorPath(genFinal, acceptFinal, generator_path, acceptor_path);

        assert(underapprox_detector->accepts(genFinal, acceptFinal));
        for(auto& t:generator_path){
            int edgeID = t.edgeID;
            int input = t.input;
            assert(input == 0);
            int output = t.output;
            Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
            assert(outer->value(v) == l_True);
            conflict.push(mkLit(v, true));
        }
        for(auto& t:acceptor_path){
            int edgeID = t.edgeID;
            int input = t.input;
            int output = t.output;
            assert(output ==0);
            Var v = outer->getTransitionVar(acceptor_over.getID(), edgeID, input, 0);
            assert(outer->value(v) == l_True);
            conflict.push(mkLit(v, true));
        }
    }else{

        forced_edges.clear();
        assert(!inverted_overapprox_detector->accepts(genFinal, acceptFinal, true, &forced_edges));
        assert(!inverted_overapprox_detector->accepts(genFinal, acceptFinal, true));
        //run an NFA to find all transitions that accepting prefixes use.
        //printf("conflict %d\n",iter);
        //g_over.draw(gen_source,genFinal);
        static vec<NFATransition> generator_path;
        generator_path.clear();
        inverted_overapprox_detector->getGeneratorPath(genFinal, acceptFinal, generator_path, false, true);

        for(auto& t:generator_path){
            //printf("%d(c%d), ", t.edgeID, t.output);
            int edgeID = t.edgeID;
            int input = t.input;
            assert(input == 0);
            int output = t.output;
            Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
            assert(outer->value(v) != l_False);
            stats_forced_edges++;
            if(outer->value(v) == l_True){
                conflict.push(mkLit(v, true));
            }
        }


    }

}

bool FSMGeneratorAcceptorDetector::stepGenerator(int final, int forcedEdge, int forcedLabel, vec<int>& store,
                                                 vec<bool>& store_seen, int& cur_gen_state, vec<NFATransition>* path){
    DynamicFSM& g = g_over;


    for(int i = 0; i < gen_cur.size(); i++){
        int s = gen_cur[i];
        cur_gen_state = s;
        for(int j = 0; j < g.nIncident(s); j++){
//now check if the label is active
            int edgeID = g.incident(s, j).id;
            int to = g.incident(s, j).node;

            if(g.transitionEnabled(edgeID, 0, 0)){
                if(!gen_cur_seen[to]){
                    gen_cur_seen[to] = true;
                    gen_cur.push(to);
                }
                if(path){
                    path->push({edgeID, 0, 0});
                }
            }

            for(int l = 1; l < g.outAlphabet(); l++){
                if(g.transitionEnabled(edgeID, 0, l) || (edgeID == forcedEdge && l == forcedLabel)){
                    if(!gen_next_seen[to]){
                        gen_next_seen[to] = true;
                        gen_next.push(to);
                    }
                    if(path && !((edgeID == forcedEdge && l == forcedLabel))){
                        path->push({edgeID, 0, l});
                    }
                    if(!store_seen[l]){
                        store_seen[l] = true;
                        store.push(l);
                    }

                }
            }

        }
    }

    gen_next.swap(gen_cur);
    gen_next_seen.swap(gen_cur_seen);

    for(int s:gen_next){
        assert(gen_next_seen[s]);
        gen_next_seen[s] = false;
    }
    gen_next.clear();

    return gen_cur_seen[final];

}

bool FSMGeneratorAcceptorDetector::isAttractor(int acceptorState){
    DynamicFSM& accept = acceptor_over;
    if(acceptorState < 0){
        return true;
    }else{
        //should really fix this to work correctly for epsilon transitions between multiple acceptor states...
        for(int i = 0; i < accept.nIncident(acceptorState); i++){
            int edgeID = accept.incident(acceptorState, i).id;
            int to = accept.incident(acceptorState, i).node;
            if(to == acceptorState){
                for(int c = 1; c < accept.inAlphabet(); c++){
                    if(!accept.transitionEnabled(edgeID, c, -1)){
                        return false;
                    }
                }
                return true;
            }
        }
    }
    return false;
}

bool FSMGeneratorAcceptorDetector::find_gen_path(int gen_final, int accept_final, int forcedEdge, int forcedLabel,
                                                 vec<NFATransition>& path, bool invertAcceptance, bool all_paths,
                                                 bool under_approx, int ignoredEdge, int ignoredLabel){
    bool accepting_state_is_attractor = !invertAcceptance && isAttractor(accept_final);
    path.clear();
    for(int s:cur){
        assert(cur_seen);
        cur_seen[s] = false;
    }
    cur.clear();
    assert(next.size() == 0);
    cur_seen[accept_source] = true;
    cur.push(accept_source);


    for(int s:gen_cur){
        assert(gen_cur_seen);
        gen_cur_seen[s] = false;
    }
    gen_cur.clear();
    assert(next.size() == 0);
    gen_cur_seen[gen_source] = true;
    gen_cur.push(gen_source);
    chars.clear();
    DynamicFSM& g = under_approx ? acceptor_under : acceptor_over;
    DynamicFSM& gen = under_approx ? g_under : g_over;
    bool any_non_acceptors = accept_source != accept_final;
    bool any_non_acceptors_next = false;
    //initial emove pass:
    if(g.emovesEnabled()){
        for(int i = 0; i < cur.size(); i++){
            int s = cur[i];
            for(int j = 0; j < g.nIncident(s); j++){
                //now check if the label is active
                int edgeID = g.incident(s, j).id;
                if(edgeID == ignoredEdge && ignoredLabel == 0)
                    continue;
                int to = g.incident(s, j).node;
                if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                    cur_seen[to] = true;
                    cur.push(to);
                    if(to != accept_final)
                        any_non_acceptors = true;
                }

            }
        }
    }
    //initial emove pass:
    if(gen.emovesEnabled()){
        for(int i = 0; i < gen_cur.size(); i++){
            int s = gen_cur[i];
            for(int j = 0; j < gen.nIncident(s); j++){
                //now check if the label is active
                int edgeID = gen.incident(s, j).id;
                int to = gen.incident(s, j).node;
                //note: intentionally not skipping if (edgeID==ignoredEdge && ignoredLabel==0), because that applied to the acceptor FSM, not the generator FSM
                if(gen.transitionEnabled(edgeID, 0, 0)){
                    if(!gen_cur_seen[to]){
                        gen_cur_seen[to] = true;
                        gen_cur.push(to);
                    }
                    path.push({edgeID, 0, 0});
                }
            }
        }
    }

    bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_final];
    bool accepted = false;
    //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
    while(!accepted){
        int prev_path_size = path.size();
        int cur_gen_state = 0;
        bool accepting = stepGenerator(gen_final, forcedEdge, forcedLabel, chars, seen_chars, cur_gen_state,
                                       &path);//get set of next strings
        if(accepting_state_is_attractor){
            accepting = true;
        }

        if(chars.size() == 0){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
                    int to = g.incident(s, j).node;
                    if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0) &&
                       !(edgeID == ignoredEdge && ignoredLabel == 0)){
                        cur_seen[to] = true;
                        cur.push(to);
                        if(to != accept_final)
                            any_non_acceptors = true;

                    }
                }
            }
        }else{
            for(int l:chars){
                assert(l > 0);
                for(int i = 0; i < cur.size(); i++){
                    int s = cur[i];
                    for(int j = 0; j < g.nIncident(s); j++){
                        //now check if the label is active
                        int edgeID = g.incident(s, j).id;
                        int to = g.incident(s, j).node;


                        if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0) &&
                           !(edgeID == ignoredEdge && ignoredLabel == 0)){
                            cur_seen[to] = true;
                            cur.push(to);
                            if(to != accept_final)
                                any_non_acceptors = true;
                            //status.reaches(str,to,edgeID,0);
                        }

                        if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0) &&
                           !(edgeID == ignoredEdge && ignoredLabel == l)){
                            //status.reaches(str,to,edgeID,l);
                            next_seen[to] = true;
                            next.push(to);
                            if(to != accept_final)
                                any_non_acceptors_next = true;
                        }
                    }
                }
            }
        }
        if(!invertAcceptance && !all_paths){
            if(prev_accepting && cur_seen[accept_final]){
                accepted = true;
                path.shrink(path.size() - prev_path_size);
            }
        }else if(invertAcceptance && !all_paths){
            if(prev_accepting && any_non_acceptors){
                accepted = true;
                path.shrink(path.size() - prev_path_size);
            }
        }else if(all_paths){
            if(prev_accepting && !any_non_acceptors){
                accepted = true;
                path.shrink(path.size() - prev_path_size);
            }
        }

        next.swap(cur);
        next_seen.swap(cur_seen);

        for(int s:next){
            assert(next_seen[s]);
            next_seen[s] = false;
        }
        next.clear();
        if(chars.size() == 0){
            //must eventually happen because the generator is linear.
            break;
        }

        for(int l :chars){
            assert(seen_chars[l]);
            seen_chars[l] = false;
        }
        chars.clear();
        prev_accepting = accepting;
        any_non_acceptors = any_non_acceptors_next;
        any_non_acceptors_next = false;
    }

    return accepted;
}


void FSMGeneratorAcceptorDetector::buildForcedEdgeReason(int genFinal, int acceptFinal, int forcedEdge, int forcedLabel,
                                                         vec<Lit>& conflict){
    //the reason that an edge was forced by an _unreachability_ condition, was that if it was negated then that state would be reachable.
    //so the 'reason' is the same as the reachability reason for that state IF the forced edge was flipped.

    tmp_path.clear();
    find_gen_path(genFinal, acceptFinal, forcedEdge, forcedLabel, tmp_path, false, true);

    for(auto& t:tmp_path){
        //printf("%d(c%d), ", t.edgeID, t.output);
        int edgeID = t.edgeID;
        int input = t.input;
        assert(input == 0);
        int output = t.output;
        Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
        assert(outer->value(v) != l_False);
        stats_forced_edges++;
        if(outer->value(v) == l_True){
            conflict.push(mkLit(v, true));
        }
    }
}

void
FSMGeneratorAcceptorDetector::buildAcceptorChokepointEdgeReason(int genFinal, int acceptFinal, int forcedAcceptorEdge,
                                                                int forcedAcceptorLabel, vec<Lit>& conflict){
    if(!opt_fsm_chokepoint_prop){
        assert(false);
        throw std::logic_error("Bad fsm option");
    }else{

        tmp_path.clear();

        //add an 'ignore edge' parameter to find_gen_path, and show that if the edge is ignored, then the inverse underapprox acceptor will accept the instance.
        bool c = buildSuffixCut(genFinal, acceptFinal, conflict, false, false, -1, -1, forcedAcceptorEdge,
                                forcedAcceptorLabel);
        if(c){
            throw std::runtime_error("Internal error in fsm theory");
        }

    }
}

void FSMGeneratorAcceptorDetector::buildForcedNondetEdgeReason(int genFinal, int acceptFinal, int forcedGeneratorEdge,
                                                               int forcedGeneratorLabel, vec<Lit>& conflict){
    if(!opt_fsm_edge_prop){
        assert(false);
        throw std::logic_error("Bad fsm option");
    }else{

        tmp_path.clear();

        //add an 'ignore edge' parameter to find_gen_path, and show that if the edge is ignored, then the inverse underapprox acceptor will accept the instance.
        bool c = buildSuffixCut(genFinal, acceptFinal, conflict, false, false, forcedGeneratorEdge,
                                forcedGeneratorLabel);
        if(c){
            throw std::runtime_error("Internal error in fsm theory");
        }
        //this edge was forced to be false, because if it was true, then the _inverted_ UNDERapprox acceptor would have an accepting path

        //underapprox_detector->accepts(gen_to,accept_to,false,&forced_edges)
/*
		bool c = find_gen_path(genFinal,acceptFinal,-1,-1,tmp_path,true,false,true,forcedEdge,forcedLabel);
		if(!c){
			throw std::runtime_error("Internal error in fsm theory");
		}
		for(auto & t:tmp_path){
			//printf("%d(c%d), ", t.edgeID, t.output);
			int edgeID = t.edgeID;
			int input = t.input;
			assert(input==0);
			int output = t.output;
			Var v = outer->getTransitionVar(g_over.getID(),edgeID,0,output);
			assert(outer->value(v)==l_False);
			stats_forced_edges++;
			if(outer->value(v)==l_True){
				conflict.push(mkLit(v,true));
			}
		}*/
        //printf("\n");

    }
}

void FSMGeneratorAcceptorDetector::buildDeterministicForcedEdgeReason(int genFinal, int acceptFinal, int forcedEdge,
                                                                      int forcedLabel, vec<Lit>& conflict){
    if(!opt_fsm_edge_prop){
        assert(false);
        throw std::logic_error("Bad fsm option");
    }else{
        tmp_path.clear();

        //this edge was forced to be false, because if it was true, then the _inverted_ UNDERapprox acceptor would have an accepting path

        //underapprox_detector->accepts(gen_to,accept_to,false,&forced_edges)
        find_gen_path(genFinal, acceptFinal, forcedEdge, forcedLabel, tmp_path, true, false, true);

        for(auto& t:tmp_path){
            //printf("%d(c%d), ", t.edgeID, t.output);
            int edgeID = t.edgeID;
            int input = t.input;
            assert(input == 0);
            int output = t.output;
            Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
            assert(outer->value(v) != l_False);
            stats_forced_edges++;
            if(outer->value(v) == l_True){
                conflict.push(mkLit(v, true));
            }
        }
        //printf("\n");

    }
}

void FSMGeneratorAcceptorDetector::buildNonSuffixAcceptReason(int genFinal, int acceptor_accept_state,
                                                              vec<int>& acceptor_start_states, SuffixLit& suffixLit,
                                                              vec<Lit>& conflict){
    conflict.push(~suffixLit.l);
    conflict.push(~suffixLit.accept_l);
    printf("Acceptor start states: ");
    for(int s:acceptor_start_states){
        printf("%d ", s);
    }
    printf("\n");
    if(!acceptor_start_states.size()){
        throw std::runtime_error("Bad suffix reason");
        return;
    }
    vec<Bitset> reachable_acceptor_states;
    vec<Bitset> previous_reachable_acceptor_states;

    DynamicFSM& g_suffix = *suffix_fsm;
    DynamicFSM& acceptor = acceptor_over;
    int suffix_start_state = suffixLit.suffix_from;
    int suffix_accept_state = suffixLit.suffix_to;

    reachable_acceptor_states.growTo(g_suffix.states());
    previous_reachable_acceptor_states.growTo(g_suffix.states());

    for(int i = 0; i < reachable_acceptor_states.size(); i++){
        reachable_acceptor_states[i].growTo(acceptor.states());
        reachable_acceptor_states[i].zero();
        previous_reachable_acceptor_states[i].growTo(acceptor.states());
        previous_reachable_acceptor_states[i].ones();
    }

    vec<int> emove_acceptor_states;
    vec<bool> emove_acceptor_seen;
    vec<int> next;
    vec<int> cur;

    vec<bool> next_seen;
    vec<bool> cur_seen;

    cur_seen.growTo(g_suffix.states());
    next_seen.growTo(g_suffix.states());
    emove_acceptor_seen.growTo(acceptor.states());
    if(acceptor.emovesEnabled()){
        for(int s : acceptor_start_states){
            emove_acceptor_seen[s] = true;
        }
        for(int i = 0; i < acceptor_start_states.size(); i++){
            int s = acceptor_start_states[i];
            for(int j = 0; j < acceptor.nIncident(s); j++){
                //now check if the label is active
                int edgeID = acceptor.incident(s, j).id;
                int accept_to = acceptor.incident(s, j).node;
                if(!emove_acceptor_seen[accept_to] && acceptor.transitionEnabled(edgeID, 0, 0)){
                    emove_acceptor_seen[accept_to] = true;
                    acceptor_start_states.push(accept_to);
                }else if(!emove_acceptor_seen[accept_to] && !acceptor.transitionEnabled(edgeID, 0, 0)){
                    Var v = outer->getTransitionVar(acceptor.getID(), edgeID, 0, 0);
                    if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                        conflict.push(mkLit(v));
                    }
                }
            }
        }
        for(int s : acceptor_start_states){
            emove_acceptor_seen[s] = false;
        }
    }

    for(int i = 0; i < reachable_acceptor_states.size(); i++){
        Bitset& bt = reachable_acceptor_states[i];
        previous_reachable_acceptor_states[i].copyFrom(bt);
        bt.zero();
    }
    {
        Bitset& bt = reachable_acceptor_states[suffix_start_state];
        bt.zero();

        for(int j = 0; j < acceptor_start_states.size(); j++){
            int s = acceptor_start_states[j];
            bt.set(s);
        }
    }
    cur.push(suffix_start_state);
    cur_seen[suffix_start_state] = true;

    if(g_suffix.emovesEnabled()){
        for(int i = 0; i < cur.size(); i++){
            int s = cur[i];
            for(int j = 0; j < g_suffix.nIncident(s); j++){
                //now check if the label is active
                int edgeID = g_suffix.incident(s, j).id;
                int to = g_suffix.incident(s, j).node;
                if(!cur_seen[to] && g_suffix.transitionEnabled(edgeID, 0, 0)){
                    cur_seen[to] = true;
                    cur.push(to);
                }
            }
        }
    }

    while(cur.size()){
        for(int s = 0; s < reachable_acceptor_states.size(); s++){
            reachable_acceptor_states[s].swap(previous_reachable_acceptor_states[s]);
        }
        for(int i = 0; i < cur.size(); i++){
            int s = cur[i];
            Bitset& prev_bt = previous_reachable_acceptor_states[s];
            Bitset& cur_bt = reachable_acceptor_states[s];
            for(int j = 0; j < g_suffix.nIncident(s); j++){
                //now check if the label is active
                int edgeID = g_suffix.incident(s, j).id;
                int gen_to = g_suffix.incident(s, j).node;
                if(!cur_seen[gen_to] && g_suffix.transitionEnabled(edgeID, 0, 0)){
                    cur_seen[gen_to] = true;
                    cur.push(gen_to);
                    //status.reaches(str,to,edgeID,0);
                }

                //for each lit accepted by a next state of the NFA acceptor from one of its current states...
                for(int acceptor_state = 0; acceptor_state < prev_bt.size(); acceptor_state++){
                    if(prev_bt[acceptor_state]){
                        //this is a state of the acceptor that could be reached at this state
                        for(int k = 0; k < acceptor.nIncident(acceptor_state); k++){
                            int acceptor_edgeID = acceptor.incident(acceptor_state, k).id;
                            int acceptor_to = acceptor.incident(acceptor_state, k).node;
                            for(int l = 1; l < acceptor.inAlphabet(); l++){

                                if(g_suffix.transitionEnabled(edgeID, l, 0) &&
                                   !acceptor.transitionEnabled(acceptor_edgeID, l, 0)){
                                    //if the suffix would accept this edge, but the acceptor does NOT
                                    Var v = outer->getTransitionVar(acceptor.getID(), acceptor_edgeID, l, 0);
                                    if(v != var_Undef && outer->level(v) > 0){
                                        assert(outer->value(v) == l_False);
                                        conflict.push(mkLit(v));
                                    }
                                }else if(acceptor.transitionEnabled(acceptor_edgeID, l, 0) &&
                                         g_suffix.transitionEnabled(edgeID, l, 0)){
                                    if(!reachable_acceptor_states[gen_to][acceptor_to]){
                                        reachable_acceptor_states[gen_to].set(acceptor_to);
                                        //also need to add any 0-reachable states

                                        if(acceptor.emovesEnabled()){
                                            emove_acceptor_states.clear();
                                            emove_acceptor_states.push(acceptor_to);
                                            for(int i = 0; i < emove_acceptor_states.size(); i++){
                                                int s = emove_acceptor_states[i];
                                                assert(reachable_acceptor_states[gen_to][s]);
                                                for(int j = 0; j < acceptor.nIncident(s); j++){
                                                    //now check if the label is active
                                                    int edgeID2 = acceptor.incident(s, j).id;
                                                    int acceptor_emove_to = acceptor.incident(s, j).node;
                                                    if(!reachable_acceptor_states[gen_to][acceptor_emove_to] &&
                                                       acceptor.transitionEnabled(edgeID2, 0, 0)){
                                                        reachable_acceptor_states[gen_to].set(acceptor_emove_to);
                                                        emove_acceptor_states.push(acceptor_emove_to);
                                                    }else if(!reachable_acceptor_states[gen_to][acceptor_emove_to] &&
                                                             !acceptor.transitionEnabled(edgeID2, 0, 0)){
                                                        Var v = outer->getTransitionVar(acceptor.getID(), edgeID2, 0,
                                                                                        0);
                                                        if(v != var_Undef && outer->value(v) == l_False &&
                                                           outer->level(v) > 0){
                                                            conflict.push(mkLit(v));
                                                        }
                                                    }
                                                }
                                            }
                                        }else{
                                            reachable_acceptor_states[gen_to].set(acceptor_to);
                                        }

                                        if(!next_seen[gen_to]){
                                            next_seen[gen_to] = true;
                                            next.push(gen_to);
                                        }
                                    }/*else if ((!reachable_acceptor_states[gen_to][acceptor_to] || !next_seen[gen_to]) && !acceptor.transitionEnabled(edgeID,l,0)){
										Var v = outer->getTransitionVar(acceptor.getID(),edgeID,0,0);
										if(v!=var_Undef &&  outer->value(v)==l_False && outer->level(v)>0){
											conflict.push(mkLit(v));
										}
									}*/
                                }
                            }
                        }
                    }
                }

            }
        }

        //if(!next.size()){
        //done processing generator
        if(cur_seen[suffix_accept_state] && reachable_acceptor_states[suffix_accept_state][acceptor_accept_state]){
            throw std::runtime_error("Internal error in FSM suffix acceptor");
            return;
        }else if(!next.size()){
            /*g_over.draw(gen_source,30);
			g_suffix.draw(suffix_start_state,suffix_accept_state);
			acceptor.draw(accept_source, acceptor_accept_state);

			printf("suffix reason: ");
			for(Lit l:conflict){
				printf("%d ", dimacs(l));
			}
			printf("\n");*/
            return;
        }
        //}

        next.swap(cur);
        next_seen.swap(cur_seen);

        for(int s:next){
            assert(next_seen[s]);
            next_seen[s] = false;
        }
        next.clear();
    }

}

void FSMGeneratorAcceptorDetector::buildNonAcceptReason(int genFinal, int acceptFinal, vec<Lit>& conflict){

    static vec<NFATransition> generatorPath;
    generatorPath.clear();
    static vec<NFATransition> acceptorPath;
    acceptorPath.clear();
    assert(!overapprox_detector->accepts(genFinal, acceptFinal, false));

    static vec<bool> seen_states;
    buildSuffixCut(genFinal, acceptFinal, conflict, false, false);
}

void FSMGeneratorAcceptorDetector::updatePrefixTable(int gen_final, int accept_final){
    if(last_prefix_update < 0 && outer->decisionLevel() == 0){
        Var v = getDetectorVar(gen_final, accept_final);
        assert(v != var_Undef);
        int index = v - first_var;
        prefixTables.growTo(index + 1);
        overapprox_detector->buildPrefixTable(gen_final, accept_final, prefixTables[index], false, false);
        last_prefix_update = g_over.modifications;
    }
}

bool FSMGeneratorAcceptorDetector::buildSuffixCut(int gen_final, int accept_final, vec<Lit>& cut,
                                                  bool accepting_state_is_attractor, bool invertAcceptance,
                                                  int ignoreGeneratorEdge, int ignoreGeneratorLabel,
                                                  int ignoreAcceptorEdge, int ignoreAcceptorLabel){
    //run the nfa backwards from the end of the generator, and collect the set of reachable fsa states at each step in the (linear) generator.
    DynamicFSM& gen = g_over;
    DynamicFSM& accept = acceptor_over;
    cur_seen.growTo(accept.states());
    gen_cur_seen.growTo(gen.states());

    next_seen.growTo(accept.states());
    gen_next_seen.growTo(gen.states());
    seen_chars.growTo(gen.outAlphabet() + 1);
    vec<Bitset>& prefixTable = getPrefixTable(gen_final, accept_final);

    for(int s:cur){
        assert(cur_seen);
        cur_seen[s] = false;
    }
    cur.clear();
    assert(next.size() == 0);
    int gen_pos = gen.states() - 1;
    if(!invertAcceptance){
        cur_seen[accept_final] = true;
        cur.push(accept_final);
    }else{
        for(int i = 0; i < accept.states(); i++){
            if(i != accept_final){
                cur_seen[i] = true;
                cur.push(i);
            }
        }
    }

    for(int s:gen_cur){
        assert(gen_cur_seen);
        gen_cur_seen[s] = false;
    }
    gen_cur.clear();
    assert(next.size() == 0);
    gen_cur_seen[gen_final] = true;
    gen_cur.push(gen_final);


    chars.clear();
    DynamicFSM& g = accept;
    bool any_non_acceptors = !cur_seen[accept_final];
    bool any_non_acceptors_next = false;
    //initial emove pass:
    if(g.emovesEnabled()){
        for(int i = 0; i < cur.size(); i++){
            int s = cur[i];
            for(int j = 0; j < g.nIncoming(s); j++){
                //now check if the label is active
                int edgeID = g.incoming(s, j).id;
                int to = g.incoming(s, j).node;
                if(edgeID != ignoreAcceptorEdge || ignoreAcceptorLabel != 0){
                    if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                        cur_seen[to] = true;
                        cur.push(to);
                        if(to != accept_final)
                            any_non_acceptors = true;
                    }else if(!g.transitionEnabled(edgeID, 0, 0)){
                        Var v = outer->getTransitionVar(g.getID(), edgeID, 0, 0);
                        if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                            cut.push(mkLit(v));
                        }
                    }
                }
            }
        }
    }
    //initial emove pass:
    if(gen.emovesEnabled()){
        for(int i = 0; i < gen_cur.size(); i++){
            int s = gen_cur[i];
            for(int j = 0; j < gen.nIncoming(s); j++){
                //now check if the label is active
                int edgeID = gen.incoming(s, j).id;
                int to = gen.incoming(s, j).node;
                if(edgeID != ignoreGeneratorEdge || ignoreGeneratorLabel != 0){
                    if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID, 0, 0)){
                        gen_cur_seen[to] = true;
                        gen_cur.push(to);
                    }else if(!gen.transitionEnabled(edgeID, 0, 0)){// && prefixTable[gen_final][to]
                        Var v = outer->getTransitionVar(gen.getID(), edgeID, 0, 0);
                        if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                            cut.push(mkLit(v));
                        }
                    }
                }
            }
        }
    }

    bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_source];
    bool accepted = false;


    //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
    while(!accepted){
        bool accepting = stepGeneratorBackward(gen_final, prefixTable, cut, chars, seen_chars, nullptr,
                                               ignoreGeneratorEdge, ignoreGeneratorLabel);//get set of next strings
        if(accepting_state_is_attractor){
            accepting = true;
        }
        if(chars.size() == 0){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncoming(s); j++){
                    //now check if the label is active
                    int edgeID = g.incoming(s, j).id;
                    int to = g.incoming(s, j).node;
                    if(edgeID != ignoreAcceptorEdge || ignoreAcceptorLabel != 0){
                        if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                            cur_seen[to] = true;
                            cur.push(to);

                            if(to != accept_final)
                                any_non_acceptors = true;

                        }else if(!cur_seen[to] && !g.transitionEnabled(edgeID, 0, 0)){
                            Var v = outer->getTransitionVar(g.getID(), edgeID, 0, 0);
                            if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                                cut.push(mkLit(v));
                            }
                        }
                    }
                }
            }
        }else{
            for(int l:chars){
                assert(l > 0);
                for(int i = 0; i < cur.size(); i++){
                    int s = cur[i];
                    for(int j = 0; j < g.nIncoming(s); j++){
                        //now check if the label is active
                        int edgeID = g.incoming(s, j).id;
                        int to = g.incoming(s, j).node;
                        if(edgeID != ignoreAcceptorEdge || ignoreAcceptorLabel != 0){
                            if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                                cur_seen[to] = true;
                                cur.push(to);

                                if(to != accept_final)
                                    any_non_acceptors = true;
                                //status.reaches(str,to,edgeID,0);
                            }else if(!cur_seen[to] && !g.transitionEnabled(edgeID, 0, 0)){
                                Var v = outer->getTransitionVar(g.getID(), edgeID, 0, 0);
                                if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                                    cut.push(mkLit(v));
                                }
                            }
                        }
                        if(edgeID != ignoreAcceptorEdge || ignoreAcceptorLabel != l){
                            if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0)){

                                //status.reaches(str,to,edgeID,l);
                                next_seen[to] = true;
                                next.push(to);
                                if(to != accept_final)
                                    any_non_acceptors_next = true;
                            }else if(!next_seen[to] && !g.transitionEnabled(edgeID, l, 0)){
                                //can improve this by fully computing next_seen before adding any lits to the cut
                                Var v = outer->getTransitionVar(g.getID(), edgeID, l, 0);
                                if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                                    cut.push(mkLit(v));
                                }
                            }
                        }
                    }
                }
            }
        }
        if(prev_accepting && cur_seen[accept_source]){
            accepted = true;
        }

        next.swap(cur);
        next_seen.swap(cur_seen);

        for(int s:next){
            assert(next_seen[s]);
            next_seen[s] = false;
        }
        next.clear();
        if(chars.size() == 0){
            //must eventually happen because the generator is linear.
            break;
        }

        for(int l :chars){
            assert(seen_chars[l]);
            seen_chars[l] = false;
        }
        chars.clear();
        prev_accepting = accepting;
        any_non_acceptors = any_non_acceptors_next;
        any_non_acceptors_next = false;
        gen_pos--;
    }

    return accepted;

}

bool
FSMGeneratorAcceptorDetector::stepGeneratorBackward(int final, vec<Bitset>& prefixTable, vec<Lit>& cut, vec<int>& store,
                                                    vec<bool>& store_seen, vec<NFATransition>* path,
                                                    int ignoreGeneratorEdge, int ignoreGeneratorLabel){
    DynamicFSM& g = g_over;

    for(int i = 0; i < gen_cur.size(); i++){
        int s = gen_cur[i];
        for(int j = 0; j < g.nIncoming(s); j++){
//now check if the label is active
            int edgeID = g.incoming(s, j).id;
            int to = g.incoming(s, j).node;
            if(edgeID != ignoreGeneratorEdge || ignoreGeneratorLabel != 0){
                if(g.transitionEnabled(edgeID, 0, 0)){
                    if(!gen_cur_seen[to]){
                        gen_cur_seen[to] = true;
                        gen_cur.push(to);
                    }
                    if(path){
                        path->push({edgeID, 0, 0});
                    }
                }else{//if (prefixTable[s][to])
                    Var v = outer->getTransitionVar(g.getID(), edgeID, 0, 0);
                    if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                        cut.push(mkLit(v));
                    }
                }
            }

            int edge_assigned_true = -1;
            for(int l = 1; l < g_under.outAlphabet(); l++){
                if(edgeID != ignoreGeneratorEdge || ignoreGeneratorLabel != l){
                    if(g_under.transitionEnabled(edgeID, 0, l)){
                        edge_assigned_true = l;
                        break;
                    }
                }
            }
            if(edge_assigned_true >= 0){

//TODO: Is this only correct for deterministic generators?
//if a character has been decided true, then since exactly one character is learnt in each position, we can just learn !character, rather than learning the negation of the set of disabled characters
                Var v = outer->getTransitionVar(g.getID(), edgeID, 0, edge_assigned_true);
                if(v != var_Undef && outer->value(v) == l_True && outer->level(v) > 0){
                    cut.push(mkLit(v, true));
                }
            }

            for(int l = 1; l < g.outAlphabet(); l++){
                if(edgeID != ignoreGeneratorEdge || ignoreGeneratorLabel != l){
                    if(g.transitionEnabled(edgeID, 0, l)){
                        if(!gen_next_seen[to]){
                            gen_next_seen[to] = true;
                            gen_next.push(to);
                        }
                        if(path){
                            path->push({edgeID, 0, l});
                        }
                        if(!store_seen[l]){
                            store_seen[l] = true;
                            store.push(l);
                        }

                    }else if(edge_assigned_true == -1){//if (prefixTable[s][to])
                        Var v = outer->getTransitionVar(g.getID(), edgeID, 0, l);
                        if(v != var_Undef && outer->value(v) == l_False && outer->level(v) > 0){
                            cut.push(mkLit(v));
                        }
                    }
                }
            }

        }
    }

    gen_next.swap(gen_cur);
    gen_next_seen.swap(gen_cur_seen);

    for(int s:gen_next){
        assert(gen_next_seen[s]);
        gen_next_seen[s] = false;
    }
    gen_next.clear();

    return gen_cur_seen[gen_source];

}


Lit FSMGeneratorAcceptorDetector::decide(int level){


    double startdecidetime = rtime(2);
    for(auto& t:all_accept_lits){

        Lit l = t.l;
        if(outer->value(l) != l_False){
            int gen_to = t.gen_to;
            int accept_to = t.accept_to;
            static vec<NFATransition> generatorPath;
            generatorPath.clear();
            static vec<NFATransition> acceptorPath;
            acceptorPath.clear();
            if(overapprox_detector->getGeneratorAceptorPath(gen_to, accept_to, generatorPath, acceptorPath)){
                for(auto& t:generatorPath){
                    Var v = outer->getTransitionVar(g_over.getID(), t.edgeID, t.input, t.output);
                    if(outer->value(v) == l_Undef){
                        stats_decisions++;
                        stats_decide_time += rtime(2) - startdecidetime;
                        return mkLit(v, false);
                    }
                }
                for(auto& t:acceptorPath){
                    Var v = outer->getTransitionVar(acceptor_over.getID(), t.edgeID, t.input, t.output);
                    if(outer->value(v) == l_Undef){
                        stats_decisions++;
                        stats_decide_time += rtime(2) - startdecidetime;
                        return mkLit(v, false);
                    }
                }
            }
        }
    }
    stats_decide_time += rtime(2) - startdecidetime;
    return lit_Undef;
}

void FSMGeneratorAcceptorDetector::printSolution(std::ostream& out){
    if(graph){
        return;
    }

    for(auto& t:all_accept_lits){
        Lit l = t.l;

        int gen_to = t.gen_to;
        int accept_to = t.accept_to;

        if(outer->value(l) == l_True){
            static vec<NFATransition> generatorPath;
            generatorPath.clear();
            static vec<NFATransition> acceptorPath;
            acceptorPath.clear();
            assert(underapprox_detector->accepts(gen_to, accept_to));
            underapprox_detector->getGeneratorPath(gen_to, accept_to, generatorPath, acceptorPath);

            out << "Generated string: ";
            for(auto& t:generatorPath){
                int edgeID = t.edgeID;
                int output = t.output;
                Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
                assert(outer->value(v) == l_True);
                out << output;
            }
            out << "\n";
        }else{
            static vec<NFATransition> path;
            path.clear();
            if(g_under.generates(gen_source, gen_to, path)){
                out << "Generated string: ";
                for(auto& t:path){
                    int edgeID = t.edgeID;
                    int output = t.output;
                    Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
                    assert(outer->value(v) == l_True);
                    out << output;
                }
                out << "\n";
            }

        }
    }
}

bool FSMGeneratorAcceptorDetector::checkSatisfied(){

    NFALinearGeneratorAcceptor<> check(g_under, acceptor_under, gen_source, accept_source);

    for(auto& t:all_accept_lits){
        Lit l = t.l;

        int gen_to = t.gen_to;
        int accept_to = t.accept_to;


        if(outer->value(l) == l_Undef){
            continue;//allowing this case for now, although we probably shouldn't
        }else if(outer->value(l) == l_False && check.accepts(gen_to, accept_to)){
            return false;
        }else if(outer->value(l) == l_True && !check.accepts(gen_to, accept_to)){
            return false;
        }


    }
    vec<Lit> ignore;
    if(!propagate(ignore)){
        return false;
    }

    return true;
}


bool FSMGeneratorAcceptorDetector::acceptsSuffix(DynamicFSM& g_suffix, DynamicFSM& acceptor, int acceptor_accept_state,
                                                 int suffix_start_state, int suffix_accept_state,
                                                 vec<int>& acceptor_start_states){

    //there may or may not be a better way to do this.
    //but for now, take advantage of knowing that g_suffix_over is shallow and acyclic, with a small alphabet, to just do a dfs

    //for every time step, mark all states of the acceptor with all states that each suffix acceptor start state can lead to
    //g_suffix.draw(suffix_start_state,suffix_accept_state);
    //acceptor.draw(accept_source, acceptor_accept_state);

    //in order for the solver to accept a string, it must be possible for each suffix to generate a string that can be accepted by the same
    //sequence of acceptor states. (This is a neccessary but not a sufficient condition)
    if(!acceptor_start_states.size())
        return false;

    vec<Bitset> reachable_acceptor_states;
    vec<Bitset> previous_reachable_acceptor_states;
    reachable_acceptor_states.growTo(g_suffix.states());
    previous_reachable_acceptor_states.growTo(g_suffix.states());

    for(int i = 0; i < reachable_acceptor_states.size(); i++){
        reachable_acceptor_states[i].growTo(acceptor.states());
        reachable_acceptor_states[i].zero();
        previous_reachable_acceptor_states[i].growTo(acceptor.states());
        previous_reachable_acceptor_states[i].ones();
    }

    vec<int> emove_acceptor_states;
    vec<bool> emove_acceptor_seen;
    vec<int> next;
    vec<int> cur;

    static vec<bool> next_seen;
    static vec<bool> cur_seen;
    cur_seen.clear();
    next_seen.clear();
    cur_seen.growTo(g_suffix.states());
    next_seen.growTo(g_suffix.states());
    emove_acceptor_seen.growTo(acceptor.states());
    if(acceptor.emovesEnabled()){
        for(int s : acceptor_start_states){
            emove_acceptor_seen[s] = true;
        }
        for(int i = 0; i < acceptor_start_states.size(); i++){
            int s = acceptor_start_states[i];
            for(int j = 0; j < acceptor.nIncident(s); j++){
                //now check if the label is active
                int edgeID = acceptor.incident(s, j).id;
                int to = acceptor.incident(s, j).node;
                if(!emove_acceptor_seen[to] && acceptor.transitionEnabled(edgeID, 0, 0)){
                    emove_acceptor_seen[to] = true;
                    acceptor_start_states.push(to);
                }
            }
        }
        for(int s : acceptor_start_states){
            emove_acceptor_seen[s] = false;
        }
    }

    for(int i = 0; i < reachable_acceptor_states.size(); i++){
        Bitset& bt = reachable_acceptor_states[i];
        previous_reachable_acceptor_states[i].copyFrom(bt);
        bt.zero();
    }
    {
        Bitset& bt = reachable_acceptor_states[suffix_start_state];
        bt.zero();

        for(int j = 0; j < acceptor_start_states.size(); j++){
            int s = acceptor_start_states[j];
            bt.set(s);
        }
    }
    cur.push(suffix_start_state);
    cur_seen[suffix_start_state] = true;

    if(g_suffix.emovesEnabled()){
        for(int i = 0; i < cur.size(); i++){
            int s = cur[i];
            for(int j = 0; j < g_suffix.nIncident(s); j++){
                //now check if the label is active
                int edgeID = g_suffix.incident(s, j).id;
                int to = g_suffix.incident(s, j).node;
                if(!cur_seen[to] && g_suffix.transitionEnabled(edgeID, 0, 0)){
                    cur_seen[to] = true;
                    cur.push(to);
                }
            }
        }
    }
    bool accepts = true;
    while(cur.size()){
        for(int s = 0; s < reachable_acceptor_states.size(); s++){
            reachable_acceptor_states[s].swap(previous_reachable_acceptor_states[s]);
        }
        for(int i = 0; i < cur.size(); i++){
            int s = cur[i];
            Bitset& prev_bt = previous_reachable_acceptor_states[s];
            Bitset& cur_bt = reachable_acceptor_states[s];
            for(int j = 0; j < g_suffix.nIncident(s); j++){
                //now check if the label is active
                int edgeID = g_suffix.incident(s, j).id;
                int gen_to = g_suffix.incident(s, j).node;
                if(!cur_seen[gen_to] && g_suffix.transitionEnabled(edgeID, 0, 0)){
                    cur_seen[gen_to] = true;
                    cur.push(gen_to);
                    //status.reaches(str,to,edgeID,0);
                }

                //for each lit accepted by a next state of the NFA acceptor from one of its current states...
                for(int acceptor_state = 0; acceptor_state < prev_bt.size(); acceptor_state++){
                    if(prev_bt[acceptor_state]){
                        //this is a state of the acceptor that could be reached at this state
                        for(int k = 0; k < acceptor.nIncident(acceptor_state); k++){
                            int acceptor_edgeID = acceptor.incident(acceptor_state, k).id;
                            if(!acceptor.edgeEnabled(acceptor_edgeID))
                                continue;
                            int acceptor_to = acceptor.incident(acceptor_state, k).node;

                            for(int l = 1; l < acceptor.inAlphabet(); l++){
                                bool a_accept = acceptor.transitionEnabled(acceptor_edgeID, l, 0);
                                bool g_accept = g_suffix.transitionEnabled(edgeID, l, 0);
                                //bool g_accept2 = g_suffix.transitionEnabled(edgeID,0,l);
                                if(acceptor.transitionEnabled(acceptor_edgeID, l, 0) &&
                                   g_suffix.transitionEnabled(edgeID, l, 0)){
                                    if(!reachable_acceptor_states[gen_to][acceptor_to]){
                                        reachable_acceptor_states[gen_to].set(acceptor_to);
                                        //also need to add any 0-reachable states

                                        if(acceptor.emovesEnabled()){
                                            emove_acceptor_states.clear();
                                            emove_acceptor_states.push(acceptor_to);
                                            for(int i = 0; i < emove_acceptor_states.size(); i++){
                                                int s = emove_acceptor_states[i];
                                                assert(reachable_acceptor_states[gen_to][s]);
                                                for(int j = 0; j < acceptor.nIncident(s); j++){
                                                    //now check if the label is active
                                                    int acceptor_emove_edgeID = acceptor.incident(s, j).id;
                                                    int acceptor_emove_to = acceptor.incident(s, j).node;
                                                    if(!reachable_acceptor_states[gen_to][acceptor_emove_to] &&
                                                       acceptor.transitionEnabled(acceptor_emove_edgeID, 0, 0)){
                                                        reachable_acceptor_states[gen_to].set(acceptor_emove_to);
                                                        emove_acceptor_states.push(acceptor_emove_to);
                                                    }
                                                }
                                            }
                                        }
                                        if(!next_seen[gen_to]){
                                            next_seen[gen_to] = true;
                                            next.push(gen_to);
                                        }
                                    }
                                }
                            }

                        }
                    }
                }

            }
        }

        //done processing generator
        if(cur_seen[suffix_accept_state] && reachable_acceptor_states[suffix_accept_state][acceptor_accept_state]){
            return true;
        }else{
            accepts = false;
        }

        next.swap(cur);
        next_seen.swap(cur_seen);

        for(int s:next){
            assert(next_seen[s]);
            next_seen[s] = false;
        }
        next.clear();
    }

    return accepts;
}