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


FSMTransducesDetector::FSMTransducesDetector(int detectorID, FSMTheorySolver* outer, DynamicFSM& g_under,
                                             DynamicFSM& g_over, int source, vec<vec<int>>& str, double seed) :
        FSMDetector(detectorID), outer(outer), g_under(g_under), g_over(g_over), source(source), strings(str),
        rnd_seed(seed){

    underReachStatus = new FSMTransducesDetector::TransducesStatus(*this, true);
    overReachStatus = new FSMTransducesDetector::TransducesStatus(*this, false);

    underapprox_detector = new NFATransduce<FSMTransducesDetector::TransducesStatus>(g_under, source, str,
                                                                                     *underReachStatus);
    overapprox_detector = new NFATransduce<FSMTransducesDetector::TransducesStatus>(g_over, source, str,
                                                                                    *overReachStatus);

    underprop_marker = outer->newReasonMarker(getID());
    overprop_marker = outer->newReasonMarker(getID());
}

void FSMTransducesDetector::addTransducesLit(int state, int strID1, int strID2, Var outer_reach_var){


    if(first_destination == -1)
        first_destination = state;

    if(transduce_lits.has(std::tuple<int, int, int>(strID1, strID2, state))){
        Lit r = transduce_lits[std::tuple<int, int, int>(strID1, strID2, state)];
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

    transduce_lits.insert(std::tuple<int, int, int>(strID1, strID2, state), acceptLit);
    all_transduce_lits.push(std::tuple<int, int, int, Lit>(strID1, strID2, state, acceptLit));
    while(accept_lit_map.size() <= accept_var - first_var){
        accept_lit_map.push({-1, -1});
    }
    accept_lit_map[accept_var - first_var] = {strID1, state};

}

void FSMTransducesDetector::TransducesStatus::transduces(int string1, int string2, int state, bool accepts){
    Lit l = detector.transduce_lits[std::tuple<int, int, int>(string1, string2, state)];

    if(l != lit_Undef){
        if(!accepts){
            l = ~l;
        }
        if(polarity == accepts){
            lbool assign = detector.outer->value(l);
            //detector.is_changed[detector.indexOf(var(l))] = true;
            detector.changed.push({l, state, string1, string2});
        }
    }
}


bool FSMTransducesDetector::propagate(vec<Lit>& conflict){

    for(auto& t:all_transduce_lits){

        Lit l = std::get<3>(t);
        int u = std::get<2>(t);
        int str1 = std::get<0>(t);
        int str2 = std::get<1>(t);

        if(underapprox_detector->transducesString(str1, str2, u)){
            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, underprop_marker);
            }else{
                conflict.push(l);
                buildTransducesReason(u, str1, str2, conflict);
                return false;
            }
        }else if(!overapprox_detector->transducesString(str1, str2, u)){
            l = ~l;
            if(outer->value(l) == l_True){
                //do nothing
            }else if(outer->value(l) == l_Undef){
                outer->enqueue(l, overprop_marker);
            }else{
                conflict.push(l);
                buildNonTransducesReason(u, str1, str2, conflict);
                return false;
            }
        }else{

        }
    }

    return true;
}


void FSMTransducesDetector::buildReason(Lit p, vec<Lit>& reason, CRef marker){
    if(marker == underprop_marker){
        reason.push(p);
        Var v = var(p);
        int u = getState(v);
        int str1 = getString1(v);
        int str2 = getString2(v);
        buildTransducesReason(u, str1, str2, reason);
    }else if(marker == overprop_marker){
        reason.push(p);
        Var v = var(p);
        int t = getState(v);
        int str1 = getString1(v);
        int str2 = getString2(v);
        buildNonTransducesReason(t, str1, str2, reason);
    }else{
        assert(false);
    }
}

void FSMTransducesDetector::buildTransducesReason(int node, int str1, int str2, vec<Lit>& conflict){

    static vec<NFATransition> path;
    path.clear();
    bool hasPath = underapprox_detector->getPath(str1, str2, node, path);
    assert(hasPath);
    assert(underapprox_detector->transducesString(str1, str2, node));
    for(auto& t:path){
        int edgeID = t.edgeID;
        int input = t.input;
        int output = t.output;
        Var v = outer->getTransitionVar(g_over.getID(), edgeID, input, output);
        assert(outer->value(v) == l_True);
        conflict.push(mkLit(v, true));
    }
}

bool
FSMTransducesDetector::path_rec(int s, int dest, int string1, int string2, int str1_pos, int str2_pos, int emove_count,
                                vec<NFATransition>& path, vec<Lit>& conflict){
    if(str1_pos == strings[string1].size() && str2_pos == strings[string2].size() && (s == dest || dest < 0)){
        return true;
    }
    DynamicFSM& g = g_over;
    if(emove_count >= g.states()){
        return false;//this is not a great way to solve the problem of avoiding infinite e-move cycles...
    }


    int l_in = 0;
    int l_out = 0;
    if(str1_pos < strings[string1].size()){
        l_in = strings[string1][str1_pos];
    }
    if(str2_pos < strings[string2].size()){
        l_out = strings[string2][str2_pos];
    }

    for(int j = 0; j < g.nIncident(s); j++){
        //now check if the label is active
        int edgeID = g.incident(s, j).id;
        int to = g.incident(s, j).node;

        if(g.transitionEnabled(edgeID, 0, 0)){
            path.push({edgeID, 0, 0});
            if(path_rec(to, dest, string1, string2, str1_pos, str2_pos, emove_count + 1, path,
                        conflict)){//str_pos is NOT incremented!
                return true;
            }else{
                path.pop();
            }
        }else{
            Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, 0);
            if(v != var_Undef){
                assert(outer->value(v) == l_False);
                conflict.push(mkLit(v));
            }
        }

        if(l_in > 0){
            if(g.transitionEnabled(edgeID, l_in, 0)){
                path.push({edgeID, l_in, 0});
                if(path_rec(to, dest, string1, string2, str1_pos + 1, str2_pos, emove_count, path,
                            conflict)){//str_pos is incremented!
                    return true;
                }else{
                    path.pop();
                }
            }else{
                Var v = outer->getTransitionVar(g_over.getID(), edgeID, l_in, 0);
                if(v != var_Undef){
                    assert(outer->value(v) == l_False);
                    conflict.push(mkLit(v));
                }
            }
        }

        if(l_out > 0){
            if(g.transitionEnabled(edgeID, 0, l_out)){
                path.push({edgeID, 0, l_out});
                if(path_rec(to, dest, string1, string2, str1_pos, str2_pos + 1, emove_count, path,
                            conflict)){//str_pos is incremented!
                    return true;
                }else{
                    path.pop();
                }
            }else{
                Var v = outer->getTransitionVar(g_over.getID(), edgeID, 0, l_out);
                if(v != var_Undef){
                    assert(outer->value(v) == l_False);
                    conflict.push(mkLit(v));
                }
            }

        }

        if(l_in >= 0 && l_out >= 0){
            if(g.transitionEnabled(edgeID, l_in, l_out)){
                path.push({edgeID, l_in, l_out});
                if(path_rec(to, dest, string1, string2, str1_pos + 1, str2_pos + 1, 0, path,
                            conflict)){//str_pos is incremented
                    return true;
                }else{
                    path.pop();
                }
            }else{
                Var v = outer->getTransitionVar(g_over.getID(), edgeID, l_in, l_out);
                if(v != var_Undef){
                    assert(outer->value(v) == l_False);
                    conflict.push(mkLit(v));
                }
            }
        }
    }
    return false;
}

void FSMTransducesDetector::buildNonTransducesReason(int node, int str1, int str2, vec<Lit>& conflict){

    static vec<NFATransition> ignore;
    ignore.clear();
    path_rec(source, node, str1, str2, 0, 0, 0, ignore, conflict);

}

void FSMTransducesDetector::printSolution(std::ostream& out){
    g_under.draw(source);
}

bool FSMTransducesDetector::checkSatisfied(){
    NFATransduce<> check(g_under, source, strings);

    g_under.draw(source, first_destination);
    for(int i = 0; i < all_transduce_lits.size(); i++){
        auto& t = all_transduce_lits[i];
        int str1 = std::get<0>(t);
        int str2 = std::get<1>(t);;
        int to = std::get<2>(t);;
        Lit l = std::get<3>(t);;

        vec<int>& string1 = strings[str1];

        vec<int>& string2 = strings[str2];


        if(outer->value(l) == l_Undef){
            return false;
        }else if(outer->value(l) == l_False && check.transducesString(str1, str2, to)){
            return false;
        }else if(outer->value(l) == l_True && !check.transducesString(str1, str2, to)){
            return false;
        }


    }

    return true;
}


