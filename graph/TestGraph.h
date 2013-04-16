/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef TESTGRAPH_H_
#define TESTGRAPH_H_


#include "Graph.h"
namespace Minisat{
class TestGraph:Graph{
	Lit False;
	Lit True;
	TestGraph(Solver * S_):S(S_){
		True = mkLit(S->newVar(),false);
			False=~True;
			S->addClause(True);
	}
    virtual ~TestGraph();
private:
	vec<vec<Edge> > edges;
	virtual Lit newEdge(int from,int to)
    {
    		assert(isNode(from));assert(isNode(to));
    		Lit l = mkLit( S->newVar(),false);
    		Edge e{l,from,to};
    		edges[from].push(e);
    		return l;
    	}
    vec<Lit> c;
	virtual void reachesAny(int from, vec<Lit> reaches,int within_steps){

		//could put these in a separate solver to better test this...

        //bellman-ford:
        //Bellman-ford: For each vertex
		reaches.clear();
		for(int i = 0;i<nNodes();i++){
			reaches.push(False);
		}

		reaches[from]=True;

        for (int i = 0;i<within_steps;i++){
            //For each edge:
        	for(int j = 0;j<edges.size();j++){
        		Edge e = edges[j];
        		if(reaches[e.to]==True){
        			//do nothing
        		}else if (reaches[e.from]==False){
        			//do nothing
        		}else{
        			Lit r = mkLit( S->newVar(), false);
        			c.clear();
        			c.push(~r);c.push(reaches[e.to]);c.push(e.l); //r -> (e.l or reaches[e.to])
        			S->addClause(c);
        			c.clear();
        			c.push(~r);c.push(reaches[e.to]);c.push(reaches[e.from]); //r -> (reaches[e.from]) or reaches[e.to])
					S->addClause(c);
        			c.clear();
        			c.push(r);c.push(~reaches[e.to]); //~r -> ~reaches[e.to]
					S->addClause(c);
					c.clear();
					c.push(r);c.push(~reaches[e.from]);c.push(~e.l); //~r -> (~reaches[e.from] or ~e.l)
					S->addClause(c);
					reaches[e.to]=r   ;//reaches[e.to] == (var & reaches[e.from])| reaches[e.to];
        		}
        	}
        }
    }

};

};

#endif /* TESTGRAPH_H_ */
