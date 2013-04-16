/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef SIMPLEGRAPH_H_
#define SIMPLEGRAPH_H_


#include "Graph.h"
namespace Minisat{
class SimpleGraph:Graph{
	Solver solver;
	Lit False;
	Lit True;
	int offset;
	int super_qhead;
	int local_q;
	SimpleGraph(Solver * S_):S(S_){
		super_qhead=0;
		local_q=0;
		offset=-1;
		True = mkLit(solver.newVar(),false);
				False=~True;
				S->addClause(True);
	}
    virtual ~SimpleGraph();
private:
	struct Edge{
		Lit out_l;
		Lit local_l;
		int from;
		int to;
	};
	vec<vec<Edge> > edges;
	CRef graph_cause;
	virtual void cancelUntil(int level){
		solver.cancelUntil(level);
		super_qhead = S->trail.size();
		local_q = solver.trail.size();
	}

	virtual void newDecisionLevel(){
		solver.newDecisionLevel();
	}

	virtual int newNode(){
		edges.push();
		return nodes++;
	}

	virtual Lit newEdge(int from,int to)
    {
    		assert(isNode(from));assert(isNode(to));
    		Lit out_l = mkLit( S->newVar(),false);
    		Lit local_l;
    		if(offset<0){
    			local_l=  mkLit( solver.newVar(),false);
    			offset = var(out_l)-var(local_l);
    		}else{
    			//make sure the variables are offset apart.
    			local_l = mkLit(var(out_l)-offset,false);
    			assert(var(local_l)>=solver.nVars());
    			while(var(local_l)>=solver.nVars()){
    				solver.newVar();
    			}

    		}
    		Edge e{out_l,local_l,from,to};
    		edges[from].push(e);
    		return out_l;
    	}

    vec<Lit> c;
    vec<Lit> reaches;
	virtual void reachesAny(int from, vec<Lit> properties_out,int within_steps){
		properties_out.clear();
		//could put these in a separate solver to better test this...
		for(int i = 0;i<nNodes();i++){
			properties_out.push(mkLit(S->newVar(),false));
		}
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
        			Lit r = mkLit( solver.newVar(), false);
        			c.clear();
        			c.push(~r);c.push(reaches[e.to]);c.push(e.local_l); //r -> (e.l or reaches[e.to])
        			S->addClause(c);
        			c.clear();
        			c.push(~r);c.push(reaches[e.to]);c.push(reaches[e.from]); //r -> (reaches[e.from]) or reaches[e.to])
					S->addClause(c);
        			c.clear();
        			c.push(r);c.push(~reaches[e.to]); //~r -> ~reaches[e.to]
					S->addClause(c);
					c.clear();
					c.push(r);c.push(~reaches[e.from]);c.push(~e.local_l); //~r -> (~reaches[e.from] or ~e.l)
					solver.addClause(c);
					reaches[e.to]=r   ;//reaches[e.to] == (var & reaches[e.from])| reaches[e.to];
        		}
        	}
        }
    }
	void analyzeFinal(CRef confl, Lit skip_lit, vec<Lit>& out_conflict)
	{
		out_conflict.clear();
	    if (solver.decisionLevel() == 0) return;

	    Clause & c = solver.ca[confl];

	    for (int i = 0; i < c.size(); i++){
	        Var     x = var(c[i]);
	    	if(x==var(skip_lit))
	    		continue;

	        assert(x>=0);
	               assert(x<solver.nVars());
	        if (solver.level(x) > 0)
	        	solver.seen[x] = 1;
	        assert(solver.value(x)!=l_Undef);
	    }

	    int     start = solver.trail.size()-1;
	    int i;
	    for ( i = start; i >= solver.trail_lim[0]; i--){
	        Var     x = var(solver.trail[i]);
	        assert(solver.value(x)!=l_Undef);
	        assert(x>=0);
	        assert(x<solver.nVars());

	        if (solver.seen[x]){
	        	  assert(x!=var(skip_lit));
	            CRef r = solver.reason(x);
	            if (r == CRef_Undef){
	            	Var v = var(solver.trail[i]);
	            	int lev = solver.level(v);
	                out_conflict.push(~solver.trail[i]);
	            }else{
					Clause& c = solver.ca[r];
					for (int j = 0; j < c.size(); j++){

						if (solver.level(var(c[j])) > 0){
							solver.seen[var(c[j])] = 1;
						}
					}
	            }
	            solver.seen[x] = 0;
	        }
	    }
	}
	bool propagate(vec<Lit> & conflict){
		while(super_qhead<S->qhead){
			Lit out_l = S->trail[super_qhead++];
			Lit local_l = mkLit(var(out_l)-offset, sign(out_l));

			if(! solver.enqueue(local_l, CRef_Undef)){
				//this is a conflict
				conflict.clear();
				solver.analyzeFinal(local_l,conflict);
				for(int i = 0;i<conflict.size();conflict++){
					conflict[i]=mkLit(var(conflict[i])+offset,sign(conflict[i]));
				}
				return false;
			}
		}
		CRef confl = solver.propagate();
		if(confl!=CRef_Undef){
			//then we have a conflict which we need to instantiate in S
			conflict.clear();
			analyzeFinal(confl,lit_Undef,conflict);
			for(int i = 0;i<conflict.size();conflict++){
				conflict[i]=mkLit(var(conflict[i])+offset,sign(conflict[i]));
			}
			return false;
		}else{
			//find lits to prop
			while(local_q<solver.qhead){
				Lit local_l =solver.trail[local_q++];
				Lit out_l =mkLit(var(out_l)+offset, sign(out_l));

				if(! S->enqueue(out_l, graph_cause)){
					//this is a conflict
					conflict.clear();
					solver.analyzeFinal(local_l,conflict);
					for(int i = 0;i<conflict.size();conflict++){
						conflict[i]=mkLit(var(conflict[i])+offset,sign(conflict[i]));
					}
					return false;
				}
			}

			return true;
		}

	}

};

};

#endif /* SIMPLEGRAPH_H_ */
