/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef TESTGRAPH_H_
#define TESTGRAPH_H_

#include "core/Theory.h"

#include "Graph.h"
namespace Minisat{
class TestGraph:public GraphTheory{
private:
	int nodes;

	struct Edge{

		Lit l;
		int from;
		int to;
		int weight;
	};
	vec<vec<Edge> > edges;
	vec<int> edge_weights;
	Lit False;
	Lit True;
	Solver * S;
	int id;
    int theory_index;
public:
	TestGraph(Solver * S_, int _id=0):S(S_),id(_id){
		True = mkLit(S->newVar(),false);
			False=~True;
			S->addClause(True);
			nodes=0;


	}
    int getTheoryIndex(){
    	return theory_index;
    }
    void setTheoryIndex(int id){
    	theory_index=id;
    }
	void printStats(int detailLevel){

	}
     ~TestGraph(){};
     int getGraphID(){
    	 return id;
     }
	 int newNode(){
			edges.push();
		return nodes++;
	}
	 void newNodes(int n){
		 for(int i = 0;i<n;i++)
			 newNode();
	 }
	int nNodes(){

		return nodes;
	}
	bool isNode(int n){
		return n>=0 && n<nodes;
	}
	void backtrackUntil(int level){};
	void newDecisionLevel(){};
	bool propagateTheory(vec<Lit> & conflict){return true;};
	bool solveTheory(vec<Lit> & conflict){return true;};
	void enqueueTheory(Lit l){};
	void buildReason(Lit p, vec<Lit> & reason){};
	void implementConstraints(){

	}
	int getWeight(int edgeID){
		return edge_weights[edgeID];
	}
	Lit newEdge(int from,int to)
    {
    		assert(isNode(from));assert(isNode(to));

    		Lit l = mkLit(S->newVar(),false);
    		Edge e{l,from,to,1};
    		edges[from].push(e);
    		edge_weights.push(1);
    		return l;
    	}
	Lit newEdge(int from,int to, Var v, int weight=1 )
    {
    		assert(isNode(from));assert(isNode(to));
    		if(v==var_Undef)
    			v=S->newVar();
    		while(S->nVars()<=v)
    			S->newVar();
    		Lit l = mkLit(v,false);
    		Edge e{l,from,to,weight};
    		edges[from].push(e);
    		edge_weights.push(weight);
    		return l;
    	}
	void connects(int from,int to, Var reach_var,int within_steps=-1){

	}
    vec<Lit> c;
	void reachesAny(int from, vec<Lit> & reaches,int within_steps=-1){
		if(within_steps<0){
			within_steps=nodes;
		}
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
        		for(int k = 0;k<edges.size();k++){
					Edge e = edges[j][k];
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
    }

	void reachesAny(int from, Var firstVar,int within_steps=-1){
			if(within_steps<0){
				within_steps=nodes;
			}
			//could put these in a separate solver to better test this...

	        //bellman-ford:
	        //Bellman-ford: For each vertex
			vec<Lit>  reaches;
			reaches.clear();
			for(int i = 0;i<nNodes();i++){
				reaches.push(False);
			}

			reaches[from]=True;

			while(S->nVars()<=firstVar+nNodes())
				S->newVar();

	        for (int i = 0;i<within_steps;i++){

	            //For each edge:
	        	for(int j = 0;j<edges.size();j++){
	        		for(int k = 0;k<edges[j].size();k++){
						Edge e = edges[j][k];
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


	        for(int i = 0;i<reaches.size();i++){
	        	Var v = firstVar+i;
	        	Lit l = reaches[i];
	        	c.clear();
	        	c.push(mkLit(v,false));
	        	c.push(~l);
	        	S->addClause(c);
	        	c.clear();
	        	c.push(mkLit(v,true));
				c.push(l);
				S->addClause(c);
				c.clear();
	        }

	    }

	vec<bool> has_reaches;
	vec<vec<Lit> >  reaches_vecs;

	void reaches(int from, int to, Var reach_var,int within_steps=-1){

				if(within_steps<0){
					within_steps=nodes;
				}
				//could put these in a separate solver to better test this...

		        //bellman-ford:
		        //Bellman-ford: For each vertex
				has_reaches.growTo(nodes);
				reaches_vecs.growTo(nodes);
				vec<Lit> & reaches = reaches_vecs[from];

				if(!has_reaches[from]){
					has_reaches[from]=true;
					for(int i = 0;i<nNodes();i++){
						reaches.push(False);
					}


					reaches[from]=True;

					while(S->nVars()<=reach_var)
						S->newVar();

					//this is not a good fix
					for(int i = 0;i<1000;i++)
						S->newVar();

					for (int i = 0;i<within_steps;i++){
						//For each edge:
						for(int j = 0;j<edges.size();j++){
							for(int k = 0;k<edges[j].size();k++){
								Edge e = edges[j][k];
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

				}


				Var v = reach_var;
				Lit l = reaches[to];
				c.clear();
				c.push(mkLit(v,false));
				c.push(~l);
				S->addClause(c);
				c.clear();
				c.push(mkLit(v,true));
				c.push(l);
				S->addClause(c);
				c.clear();


		    }

	//v will be true if the minimum weight is <= the specified value
	void minimumSpanningTree(Var v, int minimum_weight){
		while(S->nVars()<=v){
			S->newVar();
		}
			//not implemented, so just leave this unconstrained for now
		printf("Warning: minimum spanning tree constraints are not implemented, and will be left unconstrained\n");


	}


	void edgeInMinimumSpanningTree(Var edgeVar, Var var){
		while(S->nVars()<=var){
			S->newVar();
		}
			//not implemented, so just leave this unconstrained for now
		printf("Warning: minimum spanning tree constraints are not implemented, and will be left unconstrained\n");


	}
	void maxFlow(int from, int to, int max_flow, Var v){
		printf("Warning: max flow/mincut constraints are not implemented, and will be left unconstrained\n");


	}
	void minConnectedComponents(int min_components, Var v){
		printf("Warning: max flow/mincut constraints are not implemented, and will be left unconstrained\n");


	}
};

};

#endif /* TESTGRAPH_H_ */
