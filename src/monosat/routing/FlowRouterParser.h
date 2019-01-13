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

#ifndef FLOW_ROUTER_PARSER_H_
#define FLOW_ROUTER_PARSER_H_

#include <stdio.h>

#include "monosat/utils/ParseUtils.h"
#include "monosat/core/SolverTypes.h"
#include "monosat/routing/FlowRouter.h"
#include "monosat/core/Config.h"
#include "monosat/core/Dimacs.h"
#include "monosat/graph/GraphParser.h"
#include <set>
#include <string>
#include <sstream>

namespace Monosat {


//=================================================================================================

template<class B, class Solver>
class FlowRouterParser : public Parser<B, Solver> {
    using Parser<B, Solver>::mapVar;
    GraphParser<B, Solver>* graph_parser;
    vec<FlowRouter<int64_t>*> flow_routers;
    vec<Var> vars;
    struct RouterStruct {
        int graphID;
        int routerID;
        int sourceNode;
        int destNode;
        Lit maxflowLit;
    };
    vec<RouterStruct> routers;

    struct Nets {
        int graphID;
        int routerID;
        Lit disabledEdgeLit;
        vec<Lit> dest_edge_lits;
        vec<Lit> net_reach_lits;
    };
    vec<Nets> nets;

public:
    FlowRouterParser(GraphParser<B, Solver>* graph_parser) : Parser<B, Solver>("Flow-Router"),
                                                             graph_parser(graph_parser){
        assert(graph_parser);
    }

    bool parseLine(B& in, Solver& S){

        skipWhitespace(in);
        if(*in == EOF)
            return false;
        if(match(in, "f_router_net")){
            int graphID = parseInt(in);
            int routerID = parseInt(in);
            Lit disabledEdgeLit = lit_Undef;
            {
                int parsed_lit = parseInt(in);
                Var var = abs(parsed_lit) - 1;
                var = mapVar(S, var);
                disabledEdgeLit = ((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
            }
            int n_members = parseInt(in);
            assert(n_members > 0);

            vec<Lit> dest_edge_lits;
            vec<Lit> net_reach_lits;

            for(int i = 0; i < n_members; i++){
                {
                    int parsed_lit1 = parseInt(in);
                    Var var1 = abs(parsed_lit1) - 1;
                    var1 = mapVar(S, var1);
                    Lit edge_lit = ((parsed_lit1 > 0) ? mkLit(var1) : ~mkLit(var1));
                    dest_edge_lits.push(edge_lit);
                }
                {
                    int parsed_lit2 = parseInt(in);
                    Var var2 = abs(parsed_lit2) - 1;
                    var2 = mapVar(S, var2);
                    Lit reach_lit = ((parsed_lit2 > 0) ? mkLit(var2) : ~mkLit(var2));
                    net_reach_lits.push(reach_lit);
                }
            }
            nets.push();
            nets.last().graphID = graphID;
            nets.last().routerID = routerID;
            nets.last().disabledEdgeLit = disabledEdgeLit;
            dest_edge_lits.copyTo(nets.last().dest_edge_lits);
            net_reach_lits.copyTo(nets.last().net_reach_lits);

            return true;
        }else if(match(in, "f_router")){
            int graphID = parseInt(in);
            int routerID = parseInt(in);
            int sourceNode = parseInt(in);
            int destNode = parseInt(in);
            assert(graphID >= 0);
            assert(routerID >= 0);

            int parsed_lit = parseInt(in);
            assert (parsed_lit != 0);

            Var var = abs(parsed_lit) - 1;
            var = mapVar(S, var);
            Lit maxflowLit = ((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
            routers.push({graphID, routerID, sourceNode, destNode, maxflowLit});
            return true;
        }

        return false;
    }

    void implementConstraints(Solver& S){
        for(auto& r:routers){
            GraphTheorySolver<int64_t>* G = graph_parser->getGraphTheory(r.graphID);
            assert(G);
            FlowRouter<int64_t>* router = new FlowRouter<int64_t>(&S, G, r.sourceNode, r.destNode, r.maxflowLit);
            flow_routers.growTo(r.routerID + 1);
            flow_routers[r.routerID] = router;
        }
        routers.clear();
        for(auto& n:nets){
            int routerID = n.routerID;
            auto* router = flow_routers[routerID];
            assert(router);
            router->addNet(n.disabledEdgeLit, n.dest_edge_lits, n.net_reach_lits);
        }
        nets.clear();
    }


};

//=================================================================================================
};

#endif /* GRAPH_PARSER_H_ */
