/*
 * DynamicFSM.h
 *
 *  Created on: Dec 15, 2014
 *      Author: sam
 */

#ifndef DYNAMICFSM_H_
#define DYNAMICFSM_H_

#include <vector>
#include "monosat/mtl/Vec.h"
#include "monosat/mtl/Bitset.h"
#include <algorithm>
#include <cassert>
#include "alg/NFATypes.h"
#include "monosat/dgl/DynamicGraph.h"

using namespace dgl;
namespace Monosat {


class DynamicFSM {
    DynamicGraph<bool> g;
    //std::vector<Bitset> edge_status;
    int id;
    bool has_epsilon = true;
    bool is_changed = true;
    bool is_generator = true;
    bool is_acceptor = true;
    //If true, this FSM has non non-constant transitions
    bool is_constant = false;

    bool is_deterministic = false;
    bool must_be_deterministic = false;
    int n_total_input_emoves = 0;
    int n_total_output_emoves = 0;
    int n_enabled_input_emoves = 0;
    int n_enabled_output_emoves = 0;
public:
    vec<Bitset> transitions;

    bool adaptive_history_clear = false;
    int64_t historyClearInterval = 1000;
    int modifications = 0;
    int additions = 0;
    int deletions = 0;
    int in_alphabet = 1;
    int out_alphabet = 1;
    int64_t historyclears = 0;
    struct EdgeChange {
        bool addition;

        int id;
        int input;
        int output;
        int mod;
        int prev_mod;

    };
    std::vector<EdgeChange> history;

public:
    DynamicFSM(int id = -1) : id(id){

    }

    ~DynamicFSM(){

    }

    bool isGenerator() const{
        return is_generator;
    }

    bool isAcceptor() const{
        return is_acceptor;
    }

    /**
     * Return true iff this finite state machine is linear
     * (starting from 'source').
     * A finite state machine is linear if every transition from
     * each state goes only to the next state, and if no state has a self loop or an epsilon transition.
     * @param source
     * @return
     */
    bool isLinear(bool generator=true){
        if (emovesEnabled(!generator)) {
            return false;
        }
        for(int n = 0; n < this->states(); n++){
            int expectedTo = n + 1;
            for(int i = 0; i < nIncident(n); i++){
                int t = incident(n, i).node;
                int edgeID = incident(n, i).id;
                if(transitionEnabled(edgeID, -1, -1)){
                    if(t == n){
                        // no self loops
                        return false;
                    }
                    if(expectedTo != -1 && expectedTo != t){
                        // all edges out of state 'n' go to the same state (which is not equal to n)
                        return false;
                    }
                    assert(t != n);
                    expectedTo = t;
                }
            }

        }
        return true;
    }

    //True if every state _must_ have exactly one outgoing transition in the final solution
    bool mustBeDeterministic() const{
        return must_be_deterministic;
    }

    bool isDeterministic() const{
        return is_deterministic;
    }

    bool isConstant() const{
        return is_constant;
    }

    void setIsConstant(bool isConstant){
        is_constant = isConstant;
    }

    int getID(){
        return id;
    }

    void setEmovesEnabled(bool enabled){
        has_epsilon = enabled;
    }

    bool hasEmoves(bool inputEmoves = true) const{
        return has_epsilon && (inputEmoves ? n_total_input_emoves > 0 : n_total_output_emoves);
    }

    bool emovesEnabled(bool inputEmoves = true) const{
        return has_epsilon && (inputEmoves ? n_enabled_input_emoves > 0 : n_enabled_output_emoves);
    }

/*
	bool emove(int edgeID)const{
		return emovesEnabled() && transitions[edgeID][0];
	}
*/

    int inAlphabet() const{
        return in_alphabet;
    }

    int outAlphabet() const{
        return out_alphabet;
    }

    void addInCharacter(){
        in_alphabet++;
    }

    void addOutCharacter(){
        out_alphabet++;
    }

    // Get an arbitrary enabled input edge from among this set, if any
    // return -1 if no transition among 'usingTransitions' is enabled
    int getAnyEnabledTransition(int edgeID, const Bitset & usingTransitions) const{
        assert(usingTransitions.size() == transitions[edgeID].size());
        return transitions[edgeID].getIndexOfMutualBit(usingTransitions);
    }

    int getAnyEnabledTransition(int edgeID) const{
        return transitions[edgeID].getIndexOfSetBit();
    }


    bool transitionEnabled(int edgeID, int input, int output) const{
        assert(input < inAlphabet());
        assert(output < outAlphabet());
        if(output < 0){
            for(int i = 0; i < out_alphabet; i++){
                if(transitionEnabled(edgeID, input, i))
                    return true;
            }
            return false;
        }else if(input < 0){
            for(int i = 0; i < in_alphabet; i++){
                if(transitionEnabled(edgeID, i, output))
                    return true;
            }
            return false;
        }else{
            int pos = input + output * inAlphabet();
            return transitions[edgeID][pos];
        }
    }

    int addTransition(int from, int to, int edgeID, int input, int output, bool defaultEnabled = true){
        assert(input < inAlphabet());
        assert(output < outAlphabet());
        while(from >= g.nodes() || to >= g.nodes())
            g.addNode();
        if(edgeID == -1){
            edgeID = g.addEdge(from, to, edgeID);
        }
        transitions.growTo(edgeID + 1);
        transitions[edgeID].growTo(inAlphabet() * outAlphabet());
        int pos = input + output * inAlphabet();
        if(defaultEnabled){
            transitions[edgeID].set(pos);
            if(input == 0){
                n_enabled_input_emoves++;
                n_total_input_emoves++;
            }
            if(output == 0){
                n_enabled_output_emoves++;
                n_total_output_emoves++;
            }
        }
        return edgeID;
    }

    void enableTransition(int edgeID, int input, int output){
        assert(edgeID >= 0);
        assert(edgeID < g.edges());
        assert(isEdge(edgeID));
        int pos = input + output * inAlphabet();
        if(!transitions[edgeID][pos]){
            transitions[edgeID].set(pos);
            //edge_status.setStatus(id,true);
            modifications++;
            additions = modifications;
            if(input == 0)
                n_enabled_input_emoves++;
            if(output == 0)
                n_enabled_output_emoves++;
            history.push_back({true, edgeID, input, output, modifications, additions});
        }
    }

    void disableTransition(int edgeID, int input, int output){
        assert(edgeID >= 0);
        assert(edgeID < g.edges());
        assert(isEdge(edgeID));
        int pos = input + output * inAlphabet();
        if(transitions[edgeID][pos]){
            transitions[edgeID].clear(pos);
            modifications++;
            if(input == 0)
                n_enabled_input_emoves--;
            if(output == 0)
                n_enabled_output_emoves--;
            assert(n_enabled_input_emoves >= 0);
            assert(n_enabled_output_emoves >= 0);
            history.push_back({false, edgeID, input, output, modifications, deletions});
            deletions = modifications;
        }
    }

    int addState(){
        return addNode();
    }

    int addNode(){

        g.addNode();
        modifications++;
        additions = modifications;
        deletions = modifications;
        markChanged();
        clearHistory(true);

        return g.nodes() - 1;
    }

    bool edgeEnabled(int edgeID) const{
        return g.edgeEnabled(edgeID);
    }

    bool isEdge(int edgeID) const{
        return g.isEdge(edgeID);
    }

    bool hasEdge(int edgeID) const{
        return isEdge(edgeID);
    }

    //Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
    int addEdge(int from, int to, int nid = -1){ //, int weight=1
        int id = g.addEdge(from, to, nid);

        modifications++;
        additions = modifications;
        markChanged();

        if(transitions.size() <= id){
            transitions.growTo(id + 1);
        }

        return id;
    }

    int nEdgeIDs(){
        return g.nEdgeIDs();
    }

    inline int states() const{
        return g.nodes();
    }

    inline int nodes() const{
        return g.nodes();
    }

    inline int edges() const{
        return g.edges();
    }

    bool hasEdge(int from, int to) const{

        return g.hasEdge(from, to);
    }

    int getEdge(int from, int to) const{
        return g.getEdge(from, to);
    }

    inline int nIncident(int node, bool undirected = false){
        return g.nIncident(node, undirected);
    }

    inline int nDirectedEdges(int node, bool incoming){
        return g.nDirectedEdges(node, incoming);
    }

    inline DynamicGraph<bool>::Edge& directedEdges(int node, int i, bool is_incoming){
        return g.directedEdges(node, i, is_incoming);
    }

    inline int nIncoming(int node, bool undirected = false){
        return g.nIncoming(node, undirected);
    }

    inline DynamicGraph<bool>::Edge& incident(int node, int i, bool undirected = false){
        return g.incident(node, i, undirected);
    }

    inline DynamicGraph<bool>::Edge& incoming(int node, int i, bool undirected = false){
        return g.incoming(node, i, undirected);
    }

    DynamicGraph<bool>::FullEdge getEdge(int id){
        return g.getEdge(id);
    }

    int getCurrentHistory(){
        return modifications;
    }

    int historySize(){
        return history.size();
    }

    void clearHistory(bool forceClear = false){
        //int64_t expect=std::max(1000,historyClearInterval*edges());
        if(history.size()
           && (forceClear
               || (history.size()
                   > (adaptive_history_clear ?
                      std::max((int64_t) 1000, historyClearInterval * edges()) : historyClearInterval)))){//){
            history.clear();
            historyclears++;

        }
        g.clearHistory();
    }

    //force a new modification
    void invalidate(){
        modifications++;
        additions = modifications;
        modifications++;
        deletions = modifications;
        is_changed = true;

    }

    void markChanged(){
        is_changed = true;

    }

    bool changed(){
        return is_changed;
    }

    void clearChanged(){
        is_changed = false;
        g.clearChanged();
    }

    void draw(int source = -1, int dest = -1){

        printf("digraph{\n");
        if(source >= 0){
            printf("start->%d\n", source);
        }
        if(dest >= 0){
            printf("%d [shape=doublecircle]\n", dest);
        }
        for(int i = 0; i < transitions.size(); i++){
            bool any_enabled = false;
            for(int l = 0; l < transitions[i].size(); l++){
                if(transitions[i][l]){
                    any_enabled = true;
                    break;
                }
            }
            if(any_enabled){
                printf("%d->%d [label=\"", g.getEdge(i).from, g.getEdge(i).to);

                for(int in = 0; in < inAlphabet(); in++){
                    for(int out = 0; out < outAlphabet(); out++){
                        int pos = in + inAlphabet() * out;
                        if(transitions[i][pos]){
                            if(out == 0){
                                if(in == 0){
                                    printf("{},");
                                }else{
                                    printf("%d:,", in);
                                }
                            }else{
                                if(in == 0){
                                    printf(":%d,", out);
                                }else{
                                    printf("%d:%d,", in, out);
                                }
                            }

                        }
                    }
                }


                printf("\"]\n");
            }
        }


        printf("}\n");

    }


    bool buildPrefixTable(int startState, int finalState, vec<int>& string, vec<Bitset>& table){
        table.growTo(string.size());
        for(int i = 0; i < table.size(); i++){
            table[i].clear();
            table[i].growTo(this->states());
        }
        //this isn't quite correct, because there may be emoves connecting start to final state...
        if(string.size() == 0){
            return startState == finalState;
        }
        static vec<int> curStates;
        static vec<int> nextStates;
        nextStates.clear();
        curStates.clear();
        curStates.push(startState);
        int pos = 0;
        table[pos].set(startState);

        //initial emove pass:
        if(emovesEnabled()){
            for(int i = 0; i < curStates.size(); i++){
                int s = curStates[i];
                for(int j = 0; j < nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = incident(s, j).id;
                    int to = incident(s, j).node;
                    if(!table[pos][to] && transitionEnabled(edgeID, 0, -1)){
                        table[pos].set(to);
                        curStates.push(to);
                    }

                }
            }
        }

        for(; pos < string.size(); pos++){
            int l = string[pos];
            assert(l > 0);
            for(int i = 0; i < curStates.size(); i++){
                int s = curStates[i];
                for(int j = 0; j < nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = incident(s, j).id;
                    int to = incident(s, j).node;
                    if(!table[pos][to] && transitionEnabled(edgeID, 0, -1)){
                        table[pos].set(to);
                        curStates.push(to);
                    }

                    if(pos + 1 < string.size() && !table[pos + 1][to] && transitionEnabled(edgeID, l, -1)){
                        table[pos + 1].set(to);
                        nextStates.push(to);
                    }
                }
            }

            nextStates.swap(curStates);
            nextStates.clear();
        }

        pos = string.size() - 1;
        //final emove pass:
        if(emovesEnabled()){
            for(int i = 0; i < curStates.size(); i++){
                int s = curStates[i];
                for(int j = 0; j < nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = incident(s, j).id;
                    int to = incident(s, j).node;
                    if(!table[pos][to] && transitionEnabled(edgeID, 0, -1)){
                        table[pos].set(to);
                        curStates.push(to);
                    }
                }
            }
        }
        return table[pos][finalState];
    }

    bool buildSuffixTable(int startState, int finalState, vec<int>& string, vec<Bitset>& table){

        table.growTo(string.size() + 1);
        for(int i = 0; i < table.size(); i++){
            table[i].clear();
            table[i].growTo(this->states() + 1);
        }
        if(string.size() == 0){
            //this isn't quite correct, because there may be emoves connecting start to final state...
            return startState == finalState;
        }
        static vec<int> curStates;
        static vec<int> nextStates;
        nextStates.clear();
        curStates.clear();
        curStates.push(finalState);
        int pos = string.size();
        table[pos].set(finalState);

        //initial emove pass:
        if(emovesEnabled()){
            for(int i = 0; i < curStates.size(); i++){
                int s = curStates[i];
                for(int j = 0; j < nIncoming(s); j++){
                    //now check if the label is active
                    int edgeID = incoming(s, j).id;
                    int to = incoming(s, j).node;
                    if(!table[pos][to] && transitionEnabled(edgeID, 0, -1)){
                        table[pos].set(to);
                        curStates.push(to);
                    }

                }
            }
        }

        for(; pos > 0; pos--){
            int l = string[pos - 1];
            assert(l > 0);
            for(int i = 0; i < curStates.size(); i++){
                int s = curStates[i];
                for(int j = 0; j < nIncoming(s); j++){
                    //now check if the label is active
                    int edgeID = incoming(s, j).id;
                    int to = incoming(s, j).node;
                    if(!table[pos][to] && transitionEnabled(edgeID, 0, -1)){
                        table[pos].set(to);
                        curStates.push(to);
                        //status.reaches(str,to,edgeID,0);
                    }

                    if(pos > 0 && !table[pos - 1][to] && transitionEnabled(edgeID, l, -1)){
                        table[pos - 1].set(to);
                        nextStates.push(to);
                    }
                }
            }

            nextStates.swap(curStates);
            nextStates.clear();

        }
        pos = 0;
        //final emove pass:
        if(emovesEnabled()){
            for(int i = 0; i < curStates.size(); i++){
                int s = curStates[i];
                for(int j = 0; j < nIncoming(s); j++){
                    //now check if the label is active
                    int edgeID = incoming(s, j).id;
                    int to = incoming(s, j).node;
                    if(!table[pos][to] && transitionEnabled(edgeID, 0, -1)){
                        table[pos].set(to);
                        curStates.push(to);
                    }
                }
            }
        }

        return table[0][startState];
    }

private:
    vec<int> next;
    vec<int> cur;

    vec<bool> next_seen;
    vec<bool> cur_seen;

    bool generates_path_rec(int s, int final, int emove_count, vec<NFATransition>& path){
        if(s == final){
            return true;
        }
        if(emove_count >= states()){
            //this is not a great way to solve the problem of avoiding infinite e-move cycles...
            return false;
        }


        for(int j = 0; j < nIncident(s); j++){
            //now check if the label is active
            int edgeID = incident(s, j).id;
            int to = incident(s, j).node;
            if(transitionEnabled(edgeID, 0, 0)){


                path.push({edgeID, 0, 0});
                if(generates_path_rec(to, final, emove_count + 1, path)){
                    //str_pos is NOT incremented!

                    return true;
                }else{

                    path.pop();
                }

            }
            for(int l = 0; l < outAlphabet(); l++){
                if(transitionEnabled(edgeID, 0, l)){
                    bool set_transition = false;


                    path.push({edgeID, 0, l});
                    if(generates_path_rec(to, final, 0, path)){//str_pos is incremented

                        return true;
                    }else{

                        path.pop();
                    }

                }
            }

        }
        return false;
    }

public:
    bool generates(int source, int final, vec<NFATransition>& path){
        return generates_path_rec(source, final, 0, path);
    }


    bool accepts(int source, int final, vec<int>& string){
        return accepts_prefix(source, final, string) == string.size();
    }

    int accepts_prefix(int source, int final, vec<int>& string){
        cur_seen.growTo(states());
        next_seen.growTo(states());
        for(int s:cur){
            assert(cur_seen);
            cur_seen[s] = false;
        }
        cur.clear();
        assert(next.size() == 0);
        cur_seen[source] = true;
        cur.push(source);

        int largest_prefix = 0;

        //initial emove pass:
        if(emovesEnabled()){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = incident(s, j).id;
                    int to = incident(s, j).node;
                    if(!cur_seen[to] && transitionEnabled(edgeID, 0, -1)){
                        cur_seen[to] = true;
                        cur.push(to);

                    }

                }
            }
        }
        if(string.size()){
            for(int k = 0; k < string.size(); k++){
                int l = string[k];
                assert(l > 0);
                for(int i = 0; i < cur.size(); i++){
                    int s = cur[i];
                    for(int j = 0; j < nIncident(s); j++){
                        //now check if the label is active
                        int edgeID = incident(s, j).id;
                        int to = incident(s, j).node;
                        if(!cur_seen[to] && transitionEnabled(edgeID, 0, -1)){
                            cur_seen[to] = true;
                            cur.push(to);
                            //status.reaches(str,to,edgeID,0);
                        }

                        if(!next_seen[to] && transitionEnabled(edgeID, l, -1)){
                            //status.reaches(str,to,edgeID,l);
                            next_seen[to] = true;
                            next.push(to);
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

                if(cur.size()){
                    largest_prefix = k;
                }
            }

            //final emove pass:
            if(emovesEnabled()){
                for(int i = 0; i < cur.size(); i++){
                    int s = cur[i];
                    for(int j = 0; j < nIncident(s); j++){
                        //now check if the label is active
                        int edgeID = incident(s, j).id;
                        int to = incident(s, j).node;
                        if(!cur_seen[to] && transitionEnabled(edgeID, 0, -1)){
                            cur_seen[to] = true;
                            cur.push(to);
                        }

                    }
                }
            }
        }
        if(cur_seen[final]){
            return string.size();
        }else{
            return largest_prefix;
        }

    }

    void clear(){
        g.clear();
        has_epsilon = true;
        is_changed = true;
        transitions.clear();

        in_alphabet = 1;
        out_alphabet = 1;


        next.clear();
        cur.clear();
        next_seen.clear();
        cur_seen.clear();

        invalidate();
        clearHistory(true);
    }

    void copyTo(DynamicFSM& to){
        to.clear();
        g.copyTo(to.g);
        to.has_epsilon = has_epsilon;
        to.in_alphabet = in_alphabet;
        to.out_alphabet = out_alphabet;
        for(Bitset& b:transitions){
            to.transitions.push();
            b.copyTo(to.transitions.last());
        }
    }

};

};


#endif /* DYNAMICFSM_H_ */
