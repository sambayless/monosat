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
#include "DynDijkstra.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "MaxFlow.h"
#include "utils/System.h"
namespace Minisat{
class DijGraph:public Theory{
private:


	Lit False;
	Lit True;
	int local_q;

	Solver * S;
	DynamicGraph g;
	DynamicGraph antig;
	DynamicGraph cutGraph;
	vec<Dijkstra*> reach_detectors;
	vec<Dijkstra*> non_reach_detectors;
	vec<vec<Lit> > reach_lits;

	vec<CRef> reach_markers;
	vec<CRef> non_reach_markers;

	vec<int> marker_map;
	vec<int> within;
	struct Edge{
		Var v;
		int from;
		int to;
	};
	vec<EdmondsKarp::Edge> cut;

	//Full matrix
	vec<vec<Edge> > edges;

	//Just a list of the edges
	vec<Edge> edge_list;
	Var min_edge_var;
	int num_edges;
	struct AssignedEdge{
		bool assign:1;
		int from:31;
		int to;
	};
	vec<AssignedEdge> trail;
	vec<int> trail_lim;
	EdmondsKarp mc;
public:

	double mctime;
	double reachtime;
	double unreachtime;
	double pathtime;
	double propagationtime;
	double reachupdatetime;
	double unreachupdatetime;

	DijGraph(Solver * S_):S(S_),mc(cutGraph){
		True = mkLit(S->newVar(),false);
			False=~True;
			S->addClause(True);
			num_edges=0;
			local_q=0;
			mctime=0;
			reachtime=0;
			unreachtime=0;
			pathtime=0;
			propagationtime=0;
			reachupdatetime=0;
			unreachupdatetime=0;
	}

	void printStats(){
		printf("Graph stats:\n");
		printf("Prop Time: %f\n", propagationtime);
		printf("Reach Time: %f (update Time: %f)\n", reachtime,reachupdatetime);
		printf("Unreach Time: %f (Update Time: %f)\n", unreachtime,unreachupdatetime);
		printf("Path Time: %f\n", pathtime);
		printf("Min-cut Time: %f\n", mctime);

	}

     ~DijGraph(){};
	 int newNode(){

		 edges.push();
		 for(int i = 0;i<edges.size();i++)
			 edges[i].growTo(edges.size());

		 antig.addNode();
		 cutGraph.addNode();
		return g.addNode();
	}
	 void newNodes(int n){
		 for(int i = 0;i<n;i++)
			 newNode();
	 }
	int nNodes(){
		return g.nodes;
	}
	bool isNode(int n){
		return n>=0 && n<nNodes();
	}

	void backtrackUntil(int level){
		static int it = 0;
		if(++it==28){
			int a =1;
		}
		//need to remove and add edges in the two graphs accordingly.
		if(trail_lim.size()>level){
			int stop = trail_lim[level];
			for(int i = trail.size()-1;i>=trail_lim[level];i--){
				AssignedEdge e = trail[i];
				if(i==5){
					int a =1;
				}
				if(e.from==0 && e.to==8){
					int  a=1;
				}
				if(e.assign){
					g.removeEdge(e.from,e.to);
				}else{
					antig.addEdge(e.from,e.to);
				}
			}
			trail.shrink(trail.size()-stop);
			trail_lim.shrink(trail_lim.size()-level);
			assert(trail_lim.size()==level);
			int lim = trail_lim.last();

			if(trail.size()<=5){
				int a =1;
			}
		}

		if(local_q>S->qhead)
			local_q=S->qhead;
		assert(dbg_graphsUpToDate());
		for(int i = 0;i<reach_detectors.size();i++){
			reach_detectors[i]->update();
		}

		for(int i = 0;i<non_reach_detectors.size();i++){
			non_reach_detectors[i]->update();
		}
	};
	void newDecisionLevel(){
		trail_lim.push(trail.size());
	};

	void buildReason(Lit p, vec<Lit> & reason){
		CRef marker = S->reason(var(p));
		assert(marker != CRef_Undef);
		int d = marker_map[CRef_Undef- marker];
		reason.push(p);
		assert(d!=0);
		if(d>0){
			double startpathtime = cpuTime();
			d--;
			Dijkstra & detector = *reach_detectors[d];
			//the reason is a path from s to p(provided by d)
			//p is the var for a reachability detector in dijkstra, and corresponds to a node
			Var v = var(p);
			int u =  v- var(reach_lits[d][0]);
			assert(detector.connected(u));

			while(int w= detector.previous(u) > -1){
				Lit l = mkLit( edges[w][u].v,true );
				assert(S->value(l)==l_False);
				reason.push(l);
				u=w;
			}
			double elapsed = cpuTime()-startpathtime;
			pathtime+=elapsed;
		}else{
			d=-d-1;
			Dijkstra & detector = *non_reach_detectors[d];

			//the reason is a cut separating p from s;
			//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

			//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
			//assign the mincut edge weights if they aren't already assigned.

			//set weights

			//compute the mincut
			/*Var v = var(p);
				int t =  v- var(reach_lits[d][0]);
				assert(!detector.connected(t));

			cut.clear();
			mc.minCut (detector.getSource(), t,cut);
			for(int i = 0;i<cut.size();i++){
				reason.push(mkLit( edges[cut[i].u][cut[i].v].v,false ));
			}*/
			Var v = var(p);
							int t =  v- var(reach_lits[d][0]);
			buildNonReachReason(t,d,reason);
			/*for(int i = 0;i<g.nodes;i++){
				for(int j = 0;j<g.adjacency[i];j++){
					int v = g.adjacency[i][j];
					if(mc.parity[i]!=mc.parity[v]){
						//this edge is in the min-cut
						reason.push(mkLit( edges[i][v].v,true ));
					}
				}
			}*/
		}

	}
	void buildReachReason(int node,Dijkstra & d, bool negate,vec<Lit> & conflict){
		//drawFull();
		assert(dbg_reachable(d.source,node));
		double starttime = cpuTime();
		int u = node;
		while(int p = d.previous(u) != -1){
			Var e = min_edge_var+edges[u][p].v;
			conflict.push(mkLit(e, true));
			u = p;
		}
		double elapsed = cpuTime()-starttime;
		pathtime+=elapsed;

	}


	bool dbg_reachable(int from, int to){
		DynamicGraph g;
		for(int i = 0;i<nNodes();i++){
			g.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			Edge e  = edge_list[i];
			if(S->assigns[e.v]==l_True){
				g.addEdge(e.from,e.to);
			}
		}

		Dijkstra d(from,g);

		return d.connected(to);

	}

	bool dbg_notreachable(int from, int to){
		//drawFull(from,to);
		DynamicGraph g;
		for(int i = 0;i<nNodes();i++){
			g.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			Edge e  = edge_list[i];
			if(S->assigns[e.v]!=l_False){
				g.addEdge(e.from,e.to);
			}
		}

		Dijkstra d(from,g);

		return !d.connected(to);

	}

	bool dbg_graphsUpToDate(){
		for(int i = 0;i<edge_list.size();i++){
			Edge e = edge_list[i];
			lbool val = S->value(e.v);
			if(val==l_True || val==l_Undef){
				assert(antig.hasEdge(e.from,e.to));
			}else{
				assert(!antig.hasEdge(e.from,e.to));
			}
			if(val==l_True){
				assert(g.hasEdge(e.from,e.to));
			}else{
				assert(!g.hasEdge(e.from,e.to));
			}
		}
		return true;
	}

	void buildNonReachReason(int node, int detector ,vec<Lit> & conflict){
		int u = node;
		//drawFull( non_reach_detectors[detector]->getSource(),u);
		assert(dbg_notreachable( non_reach_detectors[detector]->getSource(),u));
		double starttime = cpuTime();
		cutGraph.clearChangeSets();
		//ok, set the weights for each edge in the cut graph.
		//We edges to infinite weight if they are undef or true, and weight 1 otherwise.
		for(int i = 0;i<cutGraph.adjacency.size();i++){
			for(int j = 0;j<cutGraph.adjacency[i].size();j++){
				int v = cutGraph.adjacency[i][j];
				Var var = edges[i][v].v;
				if(S->value(var)==l_False){
					mc.setCapacity(i,v,1);
				}else{
					mc.setCapacity(i,v,0xF0F0F0);
				}
			}
		}
		cut.clear();
		int f =mc.minCut(non_reach_detectors[detector]->getSource(),node,cut);
		assert(f<0xF0F0F0); assert(f==cut.size());//because edges are only ever infinity or 1
		for(int i = 0;i<cut.size();i++){
			EdmondsKarp::Edge e = cut[i];

			Lit l = mkLit( edges[e.u][e.v].v,false);
			assert(S->value(l)==l_False);
			conflict.push(l);
		}
		double elapsed = cpuTime()-starttime;
			mctime+=elapsed;
	}
	bool propagateTheory(vec<Lit> & conflict){
		static int itp = 0;
			++itp;
		bool any_change = false;
		double startproptime = cpuTime();
		static vec<int> detectors_to_check;
		detectors_to_check.clear();
		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.
		while(local_q<S->qhead){
			Lit l = S->trail[local_q++];
			Var v = var(l);
			if(v==64){
				int a=1;
			}
			int lev = S->level(v);
			while(lev>trail_lim.size()){
				newDecisionLevel();
			}
			if(v>= min_edge_var && v<min_edge_var+num_edges){
				//this is an edge assignment
				int edge_num = v-min_edge_var;
				int from = edge_list[edge_num].from;
				int to = edge_list[edge_num].to;
				trail.push({!sign(l), from,to});
				AssignedEdge e = trail.last();
				assert(e.from==from);
				assert(e.to==to);
				if (!sign(l)){

					g.addEdge(from,to);
					for(int i = 0;i<reach_detectors.size();i++){
						if(!reach_detectors[i]->marked && (reach_detectors[i]->connected_unsafe(from)||reach_detectors[i]->connected_unsafe(to))){
								detectors_to_check.push((i+1));
								reach_detectors[i]->marked =true;
							}
					}
				}else{
					antig.removeEdge(from,to);
					for(int i = 0;i<non_reach_detectors.size();i++){

						if(!non_reach_detectors[i]->marked && ( non_reach_detectors[i]->connected_unsafe(from)||non_reach_detectors[i]->connected_unsafe(to))){
							detectors_to_check.push(-(i+1));
							non_reach_detectors[i]->marked =true;
						}
					}
				}
			}
		}
			assert(dbg_graphsUpToDate());
			for(int i = 0;i<detectors_to_check.size();i++){
				int d = detectors_to_check[i];
				assert(d!=0);
				if(d>0){
					d--;
					double startdreachtime = cpuTime();
					reach_detectors[d]->marked =false;
					reach_detectors[d]->update();
					double reachUpdateElapsed = cpuTime()-startdreachtime;
					reachupdatetime+=reachUpdateElapsed;
				/*	for(int j =0;j<reach_lits[d].size();j++){
						Lit l = reach_lits[d][j];
						if( reach_detectors[d]->connected(j)){*/
					for(int j = 0;j<reach_detectors[d]->getChanged().size();j++){
							int u = reach_detectors[d]->getChanged()[j];
						
							Lit l = reach_lits[d][u];
							if( reach_detectors[d]->connected(u)){
								assert(dbg_reachable( reach_detectors[d]->getSource(),u));
							if(S->value(l)==l_Undef){
								S->uncheckedEnqueue(l,reach_markers[d]) ;
							}else if (S->value(l)==l_False){
								double elapsed = cpuTime()-startdreachtime;
															reachtime+=elapsed;
								//conflict
								//The reason is all the literals in the shortest path in to s in d
								conflict.push(l);
								buildReachReason(u,*reach_detectors[d],false,conflict);
								//add it to s
								//return it as a conflict

								return false;
							}
						}

					}
					double elapsed = cpuTime()-startdreachtime;
								reachtime+=elapsed;
				}else{
					d=-d-1;
					double startunreachtime = cpuTime();
					non_reach_detectors[d]->marked =false;
					non_reach_detectors[d]->update();
					double unreachUpdateElapsed = cpuTime()-startunreachtime;
					unreachupdatetime+=unreachUpdateElapsed;
					assert(dbg_graphsUpToDate());
					for(int j = 0;j<non_reach_detectors[d]->getChanged().size();j++){
							int u = non_reach_detectors[d]->getChanged()[j];

							Lit l = ~reach_lits[d][u];
							if(! non_reach_detectors[d]->connected(u)){
								assert(dbg_notreachable( non_reach_detectors[d]->getSource(),u));
	/*				for(int j =0;j<reach_lits[d].size();j++){
						Lit l = ~reach_lits[d][j];
						if(! non_reach_detectors[d]->connected(j)){*/
							if(S->value(l)==l_Undef){
								S->uncheckedEnqueue(l,non_reach_markers[d]) ;
							}else if (S->value(l)==l_False){
								double elapsed = cpuTime()-startunreachtime;
												unreachtime+=elapsed;
								//conflict
								//The reason is a cut separating s from t
								conflict.push(l);
								buildNonReachReason(u,d,conflict);
								//add it to s
								//return it as a conflict
								return false;
							}
						}

					}
					double elapsed = cpuTime()-startunreachtime;
					unreachtime+=elapsed;
				}
			}




		g.clearChangeSets();
		antig.clearChangeSets();

		detectors_to_check.clear();
		//really simple dynamic optimization for online checking connectivity:
		//only need to do an update if an added edge has non-infinite endpoints in g, or if a removed edge has non-infinite endpoints in antig
		/*while(detectors_to_check.size()){

		}*/
		double elapsed = cpuTime()-startproptime;
					propagationtime+=elapsed;
		return true;
	};
	bool solve(vec<Lit> & conflict){return true;};
	void drawFull(int from, int to){
			printf("digraph{\n");
			for(int i = 0;i<nNodes();i++){
				if(i==from){
					printf("n%d [label=\"From\", style=filled, fillcolor=blue]\n", i);
				}else if (i==to){
					printf("n%d [label=\"To\", style=filled, fillcolor=red]\n", i);
				}else
					printf("n%d\n", i);
			}

			for(int i = 0;i<edge_list.size();i++){
				Edge & e = edge_list[i];
				char * s = "black";
				if(S->value(e.v)==l_True)
					s="blue";
				else if (S->value(e.v)==l_False)
					s="red";
				printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from,e.to, e.v, s);
			}

			printf("}\n");
		}
	void drawFull(){
		printf("digraph{\n");
		for(int i = 0;i<nNodes();i++){
			printf("n%d\n", i);
		}

		for(int i = 0;i<edge_list.size();i++){
			Edge & e = edge_list[i];
			char * s = "black";
			if(S->value(e.v)==l_True)
				s="blue";
			else if (S->value(e.v)==l_False)
				s="red";
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from,e.to, e.v, s);
		}

		printf("}\n");
	}

	void drawCurrent(){

	}

	Lit newEdge(int from,int to, Var v = var_Undef)
    {
		if(v==var_Undef)
			v = S->newVar();

		if(num_edges>0)
			assert(v==min_edge_var+num_edges);
		else
			min_edge_var=v;
		num_edges++;
		edge_list.push({v,from,to});
		//edges.push({from,to});
		edges[from][to]= {v,from,to};
		antig.addEdge(from,to);
		cutGraph.addEdge(from,to);
    	return mkLit(v,false);
    }
	void reachesAny(int from, Var firstVar,int within_steps=-1){

			assert(from<g.nodes);
			reach_markers.push(S->newReasonMarker());
			int mnum = CRef_Undef- reach_markers.last();
			marker_map.growTo(mnum+1);
			marker_map[mnum] = reach_markers.size();
			//marker_map.insert(reach_markers.last(),reach_markers.size());

			non_reach_markers.push(S->newReasonMarker());
			//marker_map[non_reach_markers.last()]=-non_reach_markers.size();
			//marker_map.insert(non_reach_markers.last(),non_reach_markers.size());

			mnum = CRef_Undef- non_reach_markers.last();
			marker_map.growTo(mnum+1);
			marker_map[mnum] = non_reach_markers.size();

			reach_detectors.push(new Dijkstra(from,g));

			non_reach_detectors.push(new Dijkstra(from,antig));
			reach_lits.push();
			within.push(within_steps);
			for(int i = 0;i<g.nodes;i++){
				Var reachVar = firstVar+i;// S->newVar();
				Lit reachLit=mkLit(reachVar,false);
				//reaches.push(reachLit);
				reach_lits.last().push(reachLit);
			}
			reach_detectors.last()->update();
			non_reach_detectors.last()->update();

	    }
	void reachesAny(int from, vec<Lit> & reaches,int within_steps=-1){
		assert(from<g.nodes);
		reach_markers.push(S->newReasonMarker());
		int mnum = CRef_Undef- reach_markers.last();
		marker_map.growTo(mnum+1);
		marker_map[mnum] = reach_markers.size();
		//marker_map.insert(reach_markers.last(),reach_markers.size());

		non_reach_markers.push(S->newReasonMarker());
		//marker_map[non_reach_markers.last()]=-non_reach_markers.size();
		//marker_map.insert(non_reach_markers.last(),non_reach_markers.size());

		mnum = CRef_Undef- non_reach_markers.last();
		marker_map.growTo(mnum+1);
		marker_map[mnum] = non_reach_markers.size();

		reach_detectors.push(new Dijkstra(from,g));

		non_reach_detectors.push(new Dijkstra(from,antig));
		reach_lits.push();
		within.push(within_steps);
		for(int i = 0;i<g.nodes;i++){
			Var reachVar = S->newVar();
			Lit reachLit=mkLit(reachVar,false);
			reaches.push(reachLit);
			reach_lits.last().push(reachLit);
		}
		reach_detectors.last()->update();
		non_reach_detectors.last()->update();

    }

};

};

#endif /* TESTGRAPH_H_ */
