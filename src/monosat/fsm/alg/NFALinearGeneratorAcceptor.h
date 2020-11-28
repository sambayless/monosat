/*
 * NFAGenerateAccept.h
 *
 *  Created on: Jan 3, 2015
 *      Author: sam
 */

#ifndef NFAGENERATEACCEPT_H_
#define NFAGENERATEACCEPT_H_

#include "monosat/fsm/DynamicFSM.h"
#include "monosat/mtl/Bitset.h"
#include "monosat/mtl/Vec.h"
#include "monosat/fsm/alg/NFATypes.h"

using namespace Monosat;
struct ForcedTransition {
    int generator_state;
    int character;
};
struct ChokepointTransition {
    int acceptorEdgeID;
    int character;
};

template<class Status=FSMNullStatus>
class NFALinearGeneratorAcceptor {
    DynamicFSM& gen;
    DynamicFSM& accept;
    Status& status;
    int gen_last_modification = -1;
    int gen_last_addition = -1;
    int gen_last_deletion = -1;
    int gen_history_qhead = 0;
    int gen_last_history_clear = 0;

    int accept_last_modification = -1;
    int accept_last_addition = -1;
    int accept_last_deletion = -1;
    int accept_history_qhead = 0;
    int accept_last_history_clear = 0;


    int stats_full_updates = 0;
    int stats_fast_updates = 0;
    int stats_fast_failed_updates = 0;
    int stats_skip_deletes = 0;
    int stats_skipped_updates = 0;
    int stats_num_skipable_deletions = 0;
    double mod_percentage = 0;

    double stats_full_update_time = 0;
    double stats_fast_update_time = 0;


    vec<int> next;
    vec<int> cur;

    vec<bool> next_seen;
    vec<bool> cur_seen;
    vec<int> cur_from;
    vec<bool> pre_accept_state_seen;

    vec<vec<bool>> all_acceptor_seen;
    vec<Bitset> all_seen_chars;
    vec<int> gen_used_chars;

    vec<int> gen_cur;
    vec<int> gen_next;
    vec<bool> gen_next_seen;
    vec<bool> gen_cur_seen;
    int gen_source;
    int accept_source;

    vec<Bitset> suffixTable;
    vec<int> chars;
    vec<bool> seen_chars;



    struct Check {
        int gen_final;
        int accept_final;
    };

    vec<Check> toCheck;

public:


    NFALinearGeneratorAcceptor(DynamicFSM& gen, DynamicFSM& accept, int gen_source, int accept_source,
                               Status& status = fsmNullStatus) : gen(gen), accept(accept), status(status),
                                                                 gen_source(gen_source), accept_source(accept_source){

        if(!gen.isLinear()){
            //Generator must be _linear_ (but may be non-deterministic.)
            throw std::runtime_error("Generator must be linear (but may be non-deterministic)");
        }

        cur_seen.growTo(accept.states());
        gen_cur_seen.growTo(gen.states());
        next_seen.growTo(accept.states());
        gen_next_seen.growTo(gen.states());
        seen_chars.growTo(gen.outAlphabet() + 1);

    }

private:

    bool stepGeneratorBackward(int final, vec<int>& store, vec<bool>& store_seen, int& cur_gen_state,
                               vec<NFATransition>* path = nullptr){
        DynamicFSM& g = gen;


        for(int i = 0; i < gen_cur.size(); i++){
            int s = gen_cur[i];
            for(int j = 0; j < g.nIncoming(s); j++){
                //now check if the label is active
                int edgeID = g.incoming(s, j).id;
                int to = g.incoming(s, j).node;

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

    bool stepGenerator(int final, vec<int>& store, vec<bool>& store_seen, int& cur_gen_state,
                       vec<NFATransition>* path = nullptr){
        DynamicFSM& g = gen;

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

    void addCheck(int generatorFinalState, int acceptorFinalState){
        toCheck.push({generatorFinalState, acceptorFinalState});
    }

    bool isAttractor(int acceptorState){
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

public:
    bool buildSuffixTable(int gen_final, int accept_final, vec<Bitset>& suffixTable, bool accepting_state_is_attractor,
                          bool invertAcceptance){
        //run the nfa backwards from the end of the generator, and collect the set of reachable fsa states at each step in the (linear) generator.

        suffixTable.growTo(gen.states());
        for(int i = 0; i < suffixTable.size(); i++){
            suffixTable[i].clear();
            suffixTable[i].growTo(accept.states());
        }

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
            suffixTable[gen_pos].set(accept_final);
        }else{
            for(int i = 0; i < accept.states(); i++){
                if(i != accept_final){
                    cur_seen[i] = true;
                    cur.push(i);
                    suffixTable[gen_pos].set(i);
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
        bool any_non_acceptors = !cur_seen[accept_final];//should this be 'accept_source'?
        bool any_non_acceptors_next = false;
        //initial emove pass:
        if(g.emovesEnabled()){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncoming(s); j++){
                    //now check if the label is active
                    int edgeID = g.incoming(s, j).id;
                    int to = g.incoming(s, j).node;
                    if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                        cur_seen[to] = true;
                        suffixTable[gen_pos].set(to);
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
                for(int j = 0; j < gen.nIncoming(s); j++){
                    //now check if the label is active
                    int edgeID = gen.incoming(s, j).id;
                    int to = gen.incoming(s, j).node;
                    if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID, 0, 0)){
                        gen_cur_seen[to] = true;
                        gen_cur.push(to);

                    }
                }
            }
        }

        bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_source];
        bool accepted = false;


        //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
        while(!accepted){
            int ignore;
            bool accepting = stepGeneratorBackward(gen_final, chars, seen_chars, ignore);//get set of next strings
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
                        if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                            cur_seen[to] = true;
                            cur.push(to);
                            suffixTable[gen_pos].set(to);
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
                        for(int j = 0; j < g.nIncoming(s); j++){
                            //now check if the label is active
                            int edgeID = g.incoming(s, j).id;
                            int to = g.incoming(s, j).node;
                            if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                                cur_seen[to] = true;
                                cur.push(to);
                                suffixTable[gen_pos].set(to);
                                if(to != accept_final)
                                    any_non_acceptors = true;
                                //status.reaches(str,to,edgeID,0);
                            }

                            if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0)){
                                if(gen_pos > 0){
                                    suffixTable[gen_pos - 1].set(to);
                                }
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

    bool buildPrefixTable(int gen_final, int accept_final, vec<Bitset>& suffixTable, bool accepting_state_is_attractor,
                          bool invertAcceptance){
        //run the nfa backwards from the end of the generator, and collect the set of reachable fsa states at each step in the (linear) generator.

        suffixTable.growTo(gen.states());
        for(int i = 0; i < suffixTable.size(); i++){
            suffixTable[i].clear();
            suffixTable[i].growTo(accept.states());
        }

        for(int s:cur){
            assert(cur_seen);
            cur_seen[s] = false;
        }
        cur.clear();
        assert(next.size() == 0);
        int gen_pos = gen.states() - 1;
        if(!invertAcceptance){
            cur_seen[accept_source] = true;
            cur.push(accept_source);
            suffixTable[gen_pos].set(accept_source);
        }else{
            for(int i = 0; i < accept.states(); i++){
                if(i != accept_source){
                    cur_seen[i] = true;
                    cur.push(i);
                    suffixTable[gen_pos].set(i);
                }
            }
        }

        for(int s:gen_cur){
            assert(gen_cur_seen);
            gen_cur_seen[s] = false;
        }
        gen_cur.clear();
        assert(next.size() == 0);
        gen_cur_seen[gen_source] = true;
        gen_cur.push(gen_source);
        chars.clear();
        DynamicFSM& g = accept;
        bool any_non_acceptors = !cur_seen[accept_final];//check this!
        bool any_non_acceptors_next = false;
        //initial emove pass:
        if(g.emovesEnabled()){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
                    int to = g.incident(s, j).node;
                    if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                        cur_seen[to] = true;
                        suffixTable[gen_pos].set(to);
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
                    if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID, 0, 0)){
                        gen_cur_seen[to] = true;
                        gen_cur.push(to);

                    }
                }
            }
        }

        bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_source];
        bool accepted = false;


        //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
        while(!accepted){
            int ignore;
            bool accepting = stepGenerator(gen_final, chars, seen_chars, ignore);//get set of next strings
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
                        if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                            cur_seen[to] = true;
                            cur.push(to);
                            suffixTable[gen_pos].set(to);
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
                            if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                                cur_seen[to] = true;
                                cur.push(to);
                                suffixTable[gen_pos].set(to);
                                if(to != accept_final)
                                    any_non_acceptors = true;
                                //status.reaches(str,to,edgeID,0);
                            }

                            if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0)){
                                if(gen_pos > 0){
                                    suffixTable[gen_pos - 1].set(to);
                                }
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
            gen_pos++;
        }

        return accepted;

    }

private:
    bool find_accepts(int gen_final, int accept_final, bool invertAcceptance,
                      vec<vec<ForcedTransition>>* forced_edges = nullptr,
                      vec<ChokepointTransition>* chokepoint_edges = nullptr, vec<int>* pre_accepting_states = nullptr){
        bool accepting_state_is_attractor = !invertAcceptance && isAttractor(accept_final);
        bool hasSuffix = false;


        cur_seen.growTo(accept.states());
        gen_cur_seen.growTo(gen.states());
        next_seen.growTo(accept.states());
        gen_next_seen.growTo(gen.states());
        seen_chars.growTo(gen.outAlphabet() + 1);

        int gen_pos = 0;
        for(int s:cur){
            assert(cur_seen);
            cur_seen[s] = false;
        }
        cur.clear();
        assert(next.size() == 0);
        cur_seen[accept_source] = true;
        cur.push(accept_source);

        if(pre_accepting_states){
            pre_accept_state_seen.clear();
            pre_accept_state_seen.growTo(cur_seen.size());
            pre_accepting_states->clear();
        }
        for(int s:gen_cur){
            assert(gen_cur_seen);
            gen_cur_seen[s] = false;
        }
        gen_cur.clear();
        assert(next.size() == 0);
        gen_cur_seen[gen_source] = true;
        gen_cur.push(gen_source);
        chars.clear();
        DynamicFSM& g = accept;
        bool any_non_acceptors = accept_source != accept_final;
        bool any_non_acceptors_next = false;
        //initial emove pass:
        if(g.emovesEnabled()){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
                    int to = g.incident(s, j).node;
                    if(hasSuffix && !suffixTable[gen_pos][to])
                        continue;
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
                    if(!gen_cur_seen[to] && gen.transitionEnabled(edgeID, 0, 0)){
                        gen_cur_seen[to] = true;
                        gen_cur.push(to);

                    }
                }
            }
        }

        bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_final];
        bool accepted = false;
        if(forced_edges)
            forced_edges->clear();
        //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
        while(!accepted){
            bool has_chokepoint = true;
            int chokepoint_edge = -1;
            int chokepoint_char = -1;
            int cur_gen_state = 0;
            bool any_forced = false;


            bool accepting = stepGenerator(gen_final, chars, seen_chars, cur_gen_state,
                                           nullptr);//get set of next strings



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
                        if(hasSuffix && !suffixTable[gen_pos][to])
                            continue;
                        if(g.transitionEnabled(edgeID, 0, 0)){
                            if(has_chokepoint){
                                if(chokepoint_edge < 0){
                                    chokepoint_edge = edgeID;
                                    chokepoint_char = 0;
                                }else{
                                    has_chokepoint = false;
                                }
                            }

                            if(!cur_seen[to]){
                                cur_seen[to] = true;
                                cur.push(to);
                            }
                            if(to != accept_final){
                                any_non_acceptors = true;
                            }else{
                                if(pre_accepting_states && accepting){
                                    if(!pre_accept_state_seen[s]){
                                        pre_accept_state_seen[s] = true;
                                        pre_accepting_states->push(s);
                                    }
                                    //s is an acceptor fsm state that can lead to the accepting state of the fsm
                                }
                            }
                        }
                    }
                }
            }else{
                if(pre_accepting_states){
                    pre_accept_state_seen.clear();
                    pre_accept_state_seen.growTo(cur_seen.size());
                    pre_accepting_states->clear();
                }
                for(int l:chars){
                    bool character_cannot_lead_to_accepting_state = true;
                    assert(l > 0);
                    for(int i = 0; i < cur.size(); i++){
                        int s = cur[i];
                        for(int j = 0; j < g.nIncident(s); j++){
                            //now check if the label is active
                            int edgeID = g.incident(s, j).id;
                            int to = g.incident(s, j).node;

                            if(g.transitionEnabled(edgeID, 0, 0)){
                                if(!hasSuffix || suffixTable[gen_pos][to]){
                                    if(has_chokepoint){
                                        if(chokepoint_edge < 0){
                                            chokepoint_edge = edgeID;
                                            chokepoint_char = 0;
                                        }else{
                                            has_chokepoint = false;
                                        }
                                    }
                                    if(!cur_seen[to]){
                                        cur_seen[to] = true;
                                        cur.push(to);
                                    }
                                    if(to != accept_final){
                                        any_non_acceptors = true;
                                    }else{
                                        if(pre_accepting_states && accepting){
                                            if(!pre_accept_state_seen[s]){
                                                pre_accept_state_seen[s] = true;
                                                pre_accepting_states->push(s);
                                                //s is an acceptor fsm state that can lead to the accepting state of the fsm
                                            }
                                        }
                                    }
                                }

                                //status.reaches(str,to,edgeID,0);
                            }

                            if(g.transitionEnabled(edgeID, l, 0)){
                                //status.reaches(str,to,edgeID,l);
                                if(!hasSuffix || (gen_pos + 1 < suffixTable.size() && suffixTable[gen_pos + 1][to])){
                                    if(has_chokepoint){
                                        if(chokepoint_edge < 0){
                                            chokepoint_edge = edgeID;
                                            chokepoint_char = l;
                                        }else{
                                            has_chokepoint = false;
                                        }
                                    }
                                    if(!next_seen[to]){
                                        next_seen[to] = true;
                                        next.push(to);
                                    }
                                    character_cannot_lead_to_accepting_state = false;
                                    if(to != accept_final){
                                        any_non_acceptors_next = true;
                                    }else{
                                        if(pre_accepting_states && accepting){
                                            if(!pre_accept_state_seen[s]){
                                                pre_accept_state_seen[s] = true;
                                                pre_accepting_states->push(s);
                                                //s is an acceptor fsm state that can lead to the accepting state of the fsm
                                            }
                                        }
                                    }
                                }


                            }
                        }
                    }

                    if(forced_edges && character_cannot_lead_to_accepting_state){
                        //this is an edge that _must_ be disabled, because it leads to a state in the nfa that cannot reach the acceptor.
                        //forced_edges->push(NFATransition{edgeID,l,0});
                        if(!any_forced){
                            forced_edges->push();
                            any_forced = true;
                        }

                        forced_edges->last().push({cur_gen_state, l});

                    }

                }
            }
            if(!invertAcceptance){
                if(prev_accepting && cur_seen[accept_final]){
                    accepted = true;
                }
            }else{
                if(prev_accepting && any_non_acceptors){
                    accepted = true;
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

            if(has_chokepoint && chokepoint_edge >= 0 && chokepoint_edges){
                assert(chokepoint_char >= 0);
                (*chokepoint_edges).push({chokepoint_edge, chokepoint_char});
            }

            prev_accepting = accepting;
            any_non_acceptors = any_non_acceptors_next;
            any_non_acceptors_next = false;
            gen_pos++;
        }

        return accepted;
    }

    bool find_gen_path(int gen_final, int accept_final, vec<NFATransition>& generator_path, bool invertAcceptance = false,
                       bool all_paths = false){
        bool accepting_state_is_attractor = !invertAcceptance && isAttractor(accept_final);
        generator_path.clear();
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
        DynamicFSM& g = accept;

        bool any_non_acceptors = accept_source != accept_final;
        bool any_non_acceptors_next = false;
        //initial emove pass:
        if(g.emovesEnabled()){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
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
        if(gen.emovesEnabled(false)){
            for(int i = 0; i < gen_cur.size(); i++){
                int s = gen_cur[i];
                for(int j = 0; j < gen.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = gen.incident(s, j).id;
                    int to = gen.incident(s, j).node;
                    if(gen.transitionEnabled(edgeID, 0, 0)){
                        if(!gen_cur_seen[to]){
                            gen_cur_seen[to] = true;
                            gen_cur.push(to);
                        }
                        generator_path.push({edgeID, 0, 0});
                    }
                }
            }
        }

        bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_final];
        bool accepted = false;
        //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
        while(!accepted){
            int prev_path_size = generator_path.size();
            int cur_gen_state = 0;
            bool accepting = stepGenerator(gen_final, chars, seen_chars, cur_gen_state, &generator_path);//get set of next strings
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
                        if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
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
                            if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                                cur_seen[to] = true;
                                cur.push(to);
                                if(to != accept_final)
                                    any_non_acceptors = true;

                            }

                            if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0)){

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
                    generator_path.shrink(generator_path.size() - prev_path_size);
                }
            }else if(invertAcceptance && !all_paths){
                if(prev_accepting && any_non_acceptors){
                    accepted = true;
                    generator_path.shrink(generator_path.size() - prev_path_size);
                }
            }else if(all_paths){
                if(prev_accepting && !any_non_acceptors){
                    accepted = true;
                    generator_path.shrink(generator_path.size() - prev_path_size);
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

    // find an accepting path through both the generator and acceptor (if one exists)
    // else, return false
    bool find_generator_acceptor_path(int gen_final, int accept_final, vec<NFATransition>& generator_path,vec<NFATransition>& acceptor_path){
        acceptor_path.clear();
        generator_path.clear();
        bool accepting_state_is_attractor = isAttractor(accept_final);
        generator_path.clear();
        for(int s:cur){
            assert(cur_seen);
            cur_seen[s] = false;
        }
        cur.clear();
        assert(next.size() == 0);
        cur_seen[accept_source] = true;
        cur.push(accept_source);

        for(vec<bool> & s :all_acceptor_seen ){
            s.clear();
        }
        for(Bitset & seenChars:all_seen_chars){
            seenChars.zero();
        }
        int acceptLength = 0;

        for(int s:gen_cur){
            assert(gen_cur_seen);
            gen_cur_seen[s] = false;
        }

        gen_cur.clear();
        assert(next.size() == 0);
        gen_cur_seen[gen_source] = true;
        gen_cur.push(gen_source);
        chars.clear();
        DynamicFSM& g = accept;

        bool any_non_acceptors = accept_source != accept_final;
        bool any_non_acceptors_next = false;
        //initial emove pass:
        if(g.emovesEnabled()){
            for(int i = 0; i < cur.size(); i++){
                int s = cur[i];
                for(int j = 0; j < g.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = g.incident(s, j).id;
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
        if(gen.emovesEnabled(false)){
            for(int i = 0; i < gen_cur.size(); i++){
                int s = gen_cur[i];
                for(int j = 0; j < gen.nIncident(s); j++){
                    //now check if the label is active
                    int edgeID = gen.incident(s, j).id;
                    int to = gen.incident(s, j).node;
                    if(gen.transitionEnabled(edgeID, 0, 0)){
                        if(!gen_cur_seen[to]){
                            gen_cur_seen[to] = true;
                            gen_cur.push(to);
                        }
                        generator_path.push({edgeID, 0, 0});
                    }
                }
            }
        }

        bool prev_accepting = accepting_state_is_attractor ? true : gen_cur_seen[gen_final];
        bool accepted = false;
        //use the linear generator to produce a (set) of strings. Because the generator is linear, it is only ever in one state, which greatly simplifies the reasoning here...
        while(!accepted){
            int prev_path_size = generator_path.size();
            int cur_gen_state = -1;
            bool accepting = stepGenerator(gen_final, chars, seen_chars, cur_gen_state, &generator_path);//get set of next strings
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
                        if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
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
                            if(!cur_seen[to] && g.transitionEnabled(edgeID, 0, 0)){
                                cur_seen[to] = true;
                                cur.push(to);
                                if(to != accept_final)
                                    any_non_acceptors = true;

                            }

                            if(!next_seen[to] && g.transitionEnabled(edgeID, l, 0)){

                                next_seen[to] = true;
                                next.push(to);
                                if(to != accept_final)
                                    any_non_acceptors_next = true;
                            }
                        }
                    }
                }
            }

            if(prev_accepting && cur_seen[accept_final]){
                accepted = true;
                generator_path.shrink(generator_path.size() - prev_path_size);
            }
            acceptLength++;

            all_acceptor_seen.growTo(acceptLength);
            all_acceptor_seen[acceptLength-1].clear();
            all_acceptor_seen[acceptLength-1].growTo(g.states());
            all_acceptor_seen[acceptLength-1].swap(cur_seen);

            cur.clear();
            next.swap(cur);
            next_seen.swap(cur_seen);

            all_seen_chars.growTo(acceptLength);
            all_seen_chars[acceptLength-1].clear();
            all_seen_chars[acceptLength-1].growTo(g.inAlphabet());

            if(chars.size() == 0){
                //must eventually happen because the generator is linear.
                break;
            }


            for(int l :chars){
                assert(seen_chars[l]);
                all_seen_chars[acceptLength-1].set(l);
                seen_chars[l] = false;
            }
            chars.clear();
            prev_accepting = accepting;
            any_non_acceptors = any_non_acceptors_next;
            any_non_acceptors_next = false;

        }
        for(int s:cur){
            assert(cur_seen);
            cur_seen[s] = false;
        }
        cur.clear();
        if (accepted){
            // walk back from the acceptor's accepting state, narrowing the seen states at each stage to only
            // ones that actually reach the accepting state
            assert(all_seen_chars.size()== all_acceptor_seen.size());
            acceptor_path.clear();
            gen_used_chars.clear();
            int cur_generator_state = gen_final;

            // there is room to improve this to find a guaranteed shortest path
            assert(cur.size()==0);

            cur.push(accept_final);
            cur_seen[accept_final]=true;

            cur_from.clear();
            cur_from.growTo(accept.states(),-1);
            bool foundAcceptingState = false;
            for (int i = acceptLength-1;i>=0;i--){
                // for states reached through emoves, records which state they came from

                assert(i>=0);
                bool foundTransitionToPrevious = false;
                int previousAcceptorState = -1;
                for(int c = 0;c<cur.size();c++){
                    int acceptor_state = cur[c];
                    // first see if there is a non-epsilon transition to an accepting previous state
                    if (i>0){
                        for(int j = 0; j < g.nIncoming(acceptor_state); j++){
                            //now check if the label is active
                            int edgeID = g.incoming(acceptor_state, j).id;
                            int from = g.incoming(acceptor_state, j).node;
                            if(i > 0 && all_acceptor_seen[i - 1][from]){
                                assert(!all_seen_chars[i - 1][0]);
                                int l = g.getAnyEnabledTransition(edgeID, all_seen_chars[i - 1]);
                                if(l > 0){
                                    foundTransitionToPrevious = true;
                                    previousAcceptorState = from;
                                    // this is a transition from a previous state
                                    // add to the path
                                    acceptor_path.push({edgeID, l, 0});
                                    // also include any emoves that were traversed to get to this state
                                    int fromEdge = cur_from[acceptor_state];
                                    int prevAcceptorState = acceptor_state;
                                    while(fromEdge != -1){
                                        assert(prevAcceptorState != g.getEdge(fromEdge).to);
                                        prevAcceptorState = g.getEdge(fromEdge).to;
                                        acceptor_path.push({fromEdge, 0, 0});
                                        fromEdge = cur_from[prevAcceptorState];
                                    }
                                    assert(l > 0);
                                    // find corresponding edge in generator (this code is simplified by the assumption
                                    // that the generator is linear
                                    gen_used_chars.push(l);
                                    break;
                                }
                            }
                        }
                    }

                    if (i==0 && acceptor_state==accept_source){
                        foundAcceptingState = true;
                        // initial emove pass
                        int fromEdge = cur_from[acceptor_state];
                        int prevAcceptorState = acceptor_state;
                        while(fromEdge != -1){
                            assert(prevAcceptorState != g.getEdge(fromEdge).to);
                            prevAcceptorState = g.getEdge(fromEdge).to;
                            acceptor_path.push({fromEdge, 0, 0});
                            fromEdge = cur_from[prevAcceptorState];
                        }
                        break;
                    }

                    if(foundTransitionToPrevious){
                        break;
                    }

                    assert(cur_seen[acceptor_state]);
                    // then there must be an epsilon transition to another state in the acceptor
                    for(int j = 0; j < g.nIncoming(acceptor_state); j++){
                        int edgeID = g.incoming(acceptor_state, j).id;
                        int from = g.incoming(acceptor_state, j).node;
                        if(all_acceptor_seen[i][from] && !cur_seen[from] && g.transitionEnabled(edgeID, 0, 0)){
                            cur.push(from);
                            cur_seen[from] = true;
                            assert(from!=acceptor_state);
                            cur_from[from] = edgeID;
                        }
                    }
                }
                assert(foundTransitionToPrevious || (i==0 && foundAcceptingState));
                assert(previousAcceptorState!=-1 || (i==0 && foundAcceptingState));


                for(int s:cur){
                    assert(cur_seen[s]);
                    cur_seen[s] = false;
                }
                cur.clear();

                cur_from.clear();
                cur_from.growTo(accept.states(),-1);
                cur.push(previousAcceptorState);
                cur_seen[previousAcceptorState]=true;
            }

            for(int s:cur){
                assert(cur_seen[s]);
                cur_seen[s] = false;
            }
            cur.clear();
            assert(foundAcceptingState);
            // narrow generator path to only include transitions that are used by this acceptor path
            // note: this logic is simplified by the assumption that the generator is linear
            int j = 0;
            int i = 0;
            // generator used characters are recorded in reverse order above
            int used_char_pos = gen_used_chars.size()-1;

            for (;used_char_pos>0 && i<generator_path.size();i++){
                assert(used_char_pos>=0);
                int l = gen_used_chars[used_char_pos];
                assert(l>0);
                if (generator_path[i].output==l){
                    generator_path[j++] = generator_path[i];
                    used_char_pos--;
                }
            }
            assert(used_char_pos<=0);

            generator_path.shrink(i-j);

        }
        assert(cur.size()==0);
        assert(next.size()==0);
        return accepted;
    }




public:

    void update(){

        if(gen_last_modification > 0 && gen.modifications == gen_last_modification && accept_last_modification > 0 &&
           accept.modifications == accept_last_modification){
            stats_skipped_updates++;
            return;
        }
        stats_full_updates++;

        if(gen_last_modification <= 0 || gen.changed() || gen_last_history_clear != gen.historyclears ||
           accept_last_modification <= 0 || accept.changed() || accept_last_history_clear != accept.historyclears){
            gen_next_seen.clear();
            gen_next_seen.growTo(gen.states());
            gen_cur_seen.clear();
            gen_cur_seen.growTo(gen.states());

            next_seen.clear();
            next_seen.growTo(accept.states());
            cur_seen.clear();
            cur_seen.growTo(accept.states());
        }

        accept_last_modification = accept.modifications;
        accept_last_deletion = accept.deletions;
        accept_last_addition = accept.additions;
        accept_history_qhead = accept.history.size();
        accept_last_history_clear = accept.historyclears;

        gen_last_modification = gen.modifications;
        gen_last_deletion = gen.deletions;
        gen_last_addition = gen.additions;
        gen_history_qhead = gen.history.size();
        gen_last_history_clear = gen.historyclears;
    }


public:


    bool accepts(int genFinal, int acceptFinal, bool invertAcceptor = false,
                 vec<vec<ForcedTransition>>* forced_edges = nullptr,
                 vec<ChokepointTransition>* chokepoint_edges = nullptr, vec<int>* pre_accepting_states = nullptr){
        return find_accepts(genFinal, acceptFinal, invertAcceptor, forced_edges, chokepoint_edges,
                            pre_accepting_states);

    }

    bool getGeneratorPath(int genFinal, int acceptFinal, vec<NFATransition>& generator_path,
                        bool invert_acceptor = false, bool all_paths = false){
        return find_gen_path(genFinal, acceptFinal, generator_path, invert_acceptor, all_paths);
    }


    bool getGeneratorAceptorPath(int genFinal, int acceptFinal, vec<NFATransition>& generator_path,
                          vec<NFATransition>& acceptor_path){
        return find_generator_acceptor_path(genFinal, acceptFinal, generator_path, acceptor_path);
    }
};


#endif
