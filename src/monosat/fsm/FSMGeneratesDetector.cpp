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

#include "FSMTheory.h"

using namespace Monosat;


FSMGeneratesDetector::FSMGeneratesDetector(int detectorID, FSMTheorySolver* outer, DynamicFSM& g_under,
                                           DynamicFSM& g_over, int source, vec<vec<int>>& str, double seed) :
        FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), source(source), strings(str),
        rnd_seed(seed){

    underReachStatus = new FSMGeneratesDetector::GenerateStatus(*this, true);
    overReachStatus = new FSMGeneratesDetector::GenerateStatus(*this, false);

    underapprox_detector = new NFAGenerate<FSMGeneratesDetector::GenerateStatus>(g_under, source, str,
                                                                                 *underReachStatus);
    overapprox_detector = new NFAGenerate<FSMGeneratesDetector::GenerateStatus>(g_over, source, str, *overReachStatus);

    underprop_marker = outer->newReasonMarker(getID());
    overprop_marker = outer->newReasonMarker(getID());
}

void FSMGeneratesDetector::addGeneratesLit(int strID, Var outer_reach_var){

    generate_lits.growTo(strings.size(), lit_Undef);


    if(generate_lits[strID] != lit_Undef){
        Lit r = generate_lits[strID];
        //force equality between the new lit and the old reach lit, in the SAT solver
        outer->makeEqualInSolver(outer->toSolver(r), mkLit(outer_reach_var));
        return;
    }

    g_under.invalidate();
    g_over.invalidate();

    Var accept_var = outer->newVar(g_over.getID(), outer_reach_var, getID());

    if(first_var == var_Undef){
        first_var = accept_var;
    }else{
        assert(accept_var >= first_var);
    }
    int index = accept_var - first_var;

    Lit acceptLit = mkLit(accept_var, false);
    all_lits.push(acceptLit);
    assert(generate_lits[strID] == lit_Undef);

    generate_lits[strID] = acceptLit;
    while(generate_lit_map.size() <= accept_var - first_var){
        generate_lit_map.push({-1});
    }
    generate_lit_map[accept_var - first_var] = {strID};

}

void FSMGeneratesDetector::GenerateStatus::generates(int string, bool generates){
    Lit l = detector.generate_lits[string];

    if(l != lit_Undef){
        if(!generates){
            l = ~l;
        }
        if(polarity == generates){
            lbool assign = detector.outer->value(l);

            detector.changed.push({l, string});
        }
    }
}


bool FSMGeneratesDetector::propagate(vec<Lit>& conflict){
    changed.clear();
    bool skipped_positive = false;
    if(underapprox_detector && (!opt_detect_pure_theory_lits || unassigned_positives > 0)){
        double startdreachtime = rtime(2);
        stats_under_updates++;
        underapprox_detector->update();
        double reachUpdateElapsed = rtime(2) - startdreachtime;

        stats_under_update_time += rtime(2) - startdreachtime;
    }else{
        skipped_positive = true;

        stats_skipped_under_updates++;
    }
    bool skipped_negative = false;
    if(overapprox_detector && (!opt_detect_pure_theory_lits || unassigned_negatives > 0)){
        double startunreachtime = rtime(2);
        stats_over_updates++;
        overapprox_detector->update();
        double unreachUpdateElapsed = rtime(2) - startunreachtime;
        stats_over_update_time += rtime(2) - startunreachtime;
    }else{
        skipped_negative = true;
        stats_skipped_over_updates++;
    }

    if(opt_rnd_shuffle){
        randomShuffle(rnd_seed, changed);
    }


    while(changed.size()){
        int sz = changed.size();
        Lit l = changed.last().l;


        int str = changed.last().str;

        if(sign(l)){
            assert(!overapprox_detector->generatesString(str));
        }else{
            assert(underapprox_detector->generatesString(str));
        }

        if(underapprox_detector && !skipped_positive && !sign(l)){

        }else if(overapprox_detector && !skipped_negative && sign(l)){

        }else{
            assert(sz == changed.size());

            changed.pop();
            //this can happen if the changed node's reachability status was reported before a backtrack in the solver.
            continue;
        }

        bool reach = !sign(l);
        if(outer->value(l) == l_True){
            //do nothing
        }else if(outer->value(l) == l_Undef){

            if(reach)
                outer->enqueue(l, underprop_marker);
            else
                outer->enqueue(l, overprop_marker);
        }else if(outer->value(l) == l_False){
            conflict.push(l);

            if(reach){
                buildGeneratesReason(str, conflict);
            }else{
                //The reason is a cut separating s from t
                buildNonGeneratesReason(str, conflict);
            }

            return false;
        }

        //This can be really tricky - if you are not careful, an a reach detector's update phase was skipped at the
        // beginning of propagate, then if the reach detector is called during propagate it can push a change onto the list,
        // which can cause the wrong item to be removed here.
        assert(sz == changed.size());

        changed.pop();
    }
    assert(changed.size() == 0);
    return true;
}


void FSMGeneratesDetector::buildReason(Lit p, vec<Lit>& reason, CRef marker){
    if(marker == underprop_marker){
        reason.push(p);
        Var v = var(p);

        int str = getString(v);
        buildGeneratesReason(str, reason);
    }else if(marker == overprop_marker){
        reason.push(p);
        Var v = var(p);
        int str = getString(v);
        buildNonGeneratesReason(str, reason);
    }else{
        assert(false);
    }
}

void FSMGeneratesDetector::buildGeneratesReason(int str, vec<Lit>& conflict){
    //find a path - ideally, the one that traverses the fewest unique transitions - from source to node, learn that one of the transitions on that path must be disabled.

    static vec<NFATransition> path;
    path.clear();
    bool hasPath = underapprox_detector->getPath(str, path);
    assert(hasPath);
    assert(underapprox_detector->generatesString(str));
    for(auto& t:path){
        int edgeID = t.edgeID;

        int output = t.output;
        Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, output);
        assert(outer->value(v) == l_True);
        conflict.push(mkLit(v, true));
    }
    //note: if there are repeated edges in this conflict, they will be cheaply removed by the sat solver anyhow, so that is not a major problem.

}

bool
FSMGeneratesDetector::unique_path_conflict(int s, int string, int str_pos, int emove_count, vec<NFATransition>& path,
                                           vec<Lit>& conflict){
    DynamicFSM& g = g_over;
    if(str_pos == strings[string].size()){
        return true;
    }
    if(emove_count >= g.states()){
        return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
    }

    int l = strings[string][str_pos];
    for(int j = 0; j < g.nIncident(s); j++){
        //now check if the label is active
        int edgeID = g.incident(s, j).id;
        int to = g.incident(s, j).node;
        if(g.transitionEnabled(edgeID, 0, 0)){
            bool set_transition = false;
            if(used_transitions[s].edge == -1){
                used_transitions[s].edge = edgeID;
                used_transitions[s].label = 0;
                set_transition = true;
                //	used_transitions[s].depth = path.size();
            }

            if(used_transitions[s].label == 0 && used_transitions[s].edge == edgeID){

                path.push({edgeID, 0});
                if(unique_path_conflict(to, string, str_pos, emove_count + 1, path,
                                        conflict)){//str_pos is NOT incremented!

                    return true;
                }else{
                    if(set_transition){
                        used_transitions[s].edge = -1;
                        used_transitions[s].label = -1;
                    }
                    path.pop();
                }
            }
        }else{
            Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, 0);
            if(v != var_Undef && used_transitions[s].label == -1){
                assert(outer->value(v) == l_False);
                conflict.push(mkLit(v));
            }
        }

        if(str_pos < strings[string].size()){
            if(g.transitionEnabled(edgeID, 0, l)){
                bool set_transition = false;
                if(used_transitions[s].edge == -1){
                    used_transitions[s].edge = edgeID;
                    used_transitions[s].label = l;
                    set_transition = true;
                    //used_transitions[s].depth = path.size();
                }
                if(used_transitions[s].label == l && used_transitions[s].edge == edgeID){
                    path.push({edgeID, l});
                    if(unique_path_conflict(to, string, str_pos + 1, 0, path, conflict)){//str_pos is incremented

                        return true;
                    }else{
                        if(set_transition){
                            used_transitions[s].edge = -1;
                            used_transitions[s].label = -1;
                        }
                        path.pop();
                    }
                }
            }else{
                Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, l);
                if(v != var_Undef && used_transitions[s].label == -1){
                    assert(outer->value(v) == l_False);
                    conflict.push(mkLit(v));
                }
            }
        }
    }
    return false;
}

void FSMGeneratesDetector::buildNonGeneratesReason(int str, vec<Lit>& conflict){
    vec<int>& string = strings[str];
    used_transitions.clear();
    used_transitions.growTo(g_over.nodes());

    static vec<NFATransition> ignore;
    ignore.clear();
    unique_path_conflict(source, str, 0, 0, ignore, conflict);

}

void FSMGeneratesDetector::printSolution(std::ostream& out){


    static vec<NFATransition> ignore;
    static vec<Lit> ignore_conf;
    for(int str = 0; str < strings.size(); str++){
        ignore.clear();
        ignore_conf.clear();
        used_transitions.clear();
        used_transitions.growTo(g_over.nodes());
        unique_path_conflict(source, str, 0, 0, ignore, ignore_conf);

        printf("digraph{\n");
        if(source >= 0){
            printf("start->%d\n", source);
        }
        for(int n = 0; n < g_under.nodes(); n++){

            int edgeID = used_transitions[n].edge;
            int l = used_transitions[n].label;
            if(l >= 0){
                printf("%d->%d [label=\"", g_under.getEdge(edgeID).from, g_under.getEdge(edgeID).to);


                if(l == 0){
                    printf("{},");
                }else{
                    printf("%c,", 'a' + l - 1);
                }


                printf("\"]\n");
            }
        }


        printf("}\n");
    }
}

bool FSMGeneratesDetector::checkSatisfied(){
    NFAGenerate<> check(g_under, source, strings);
    printSolution(std::cout);
    for(int str = 0; str < generate_lits.size(); str++){
        vec<int>& string = strings[str];


        if(generate_lits[str] != lit_Undef){

            Lit l = generate_lits[str];
            if(outer->value(l) == l_Undef){
                return false;
            }else if(outer->value(l) == l_False && check.generatesString(str)){
                return false;
            }else if(outer->value(l) == l_True && !check.generatesString(str)){
                return false;
            }

        }


    }

    return true;
}


