/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef DGRAPH_H_
#define DGRAPH_H_

#include "core/Theory.h"
#include "Graph.h"
#include "Dijkstra.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "MaxFlow.h"
#include "utils/System.h"
#ifdef DEBUG_GRAPH
#include "TestGraph.h"
#endif
namespace Minisat{



class DijGraph:public GraphTheory{
private:


	Lit False;
	Lit True;
	int local_q;

	Solver * S;
#ifdef DEBUG_GRAPH
	Solver * dbg;
	TestGraph *dbg_graph;
#endif
	DynamicGraph g;
	DynamicGraph antig;
	DynamicGraph cutGraph;

	struct ReachDetector{
		int within;
		int source;
		Dijkstra* positive_reach_detector;
		Dijkstra* negative_reach_detector;
		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;

		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

		ReachDetector():within(-1),source(-1),positive_reach_detector(NULL),negative_reach_detector(NULL){}
	};

	struct ReachInfo{
		int source;
		ReachDetector * detector;
		ReachInfo():source(-1),detector(NULL){}
	};

	vec<ReachInfo> reach_info;
	vec<ReachDetector*> reach_detectors;

	vec<CRef> reach_markers;
	vec<CRef> non_reach_markers;

	vec<int> marker_map;

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
		Var var;
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
#ifdef DEBUG_GRAPH
		dbg=new Solver();
		dbg_graph = new TestGraph(dbg);
#endif
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
#ifdef DEBUG_GRAPH
		 dbg_graph->newNode();
#endif

		 edges.push();
		 for(int i = 0;i<edges.size();i++)
			 edges[i].growTo(edges.size());

		 reach_info.push();
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

#ifdef DEBUG_GRAPH
		bool dbg_clause(const vec<Lit> &conflict){

			static vec<Lit> c;
			c.clear();
			for(int i = 0;i<conflict.size();i++)
				c.push(~conflict[i]);
		//	assert(dbg->solve());
			bool res = dbg->solve(c);
			assert(~res);

			return true;
		}
		bool dbg_propgation(Lit l){

			static vec<Lit> c;
			c.clear();
			for(int i = 0;i<S->trail.size() ;i++){
				Lit l = S->trail[i];
				Var v = var(l);
				if(v>=min_edge_var && v<min_edge_var+num_edges){
					Edge e = edge_list[v-min_edge_var];
					c.push(l);
				}
			}
			c.push(~l);
			bool res = dbg->solve(c);
			assert(~res);

			return true;
		}
#endif

	void dbg_sync(){
#ifdef DEBUG_GRAPH
		static vec<lbool> assigned;
		assigned.clear();
		for(int i = 0;i<num_edges;i++)
			assigned.push(l_Undef);
		int lev = 0;
		int j =0;
		for(int i = 0;i<trail.size();i++){
			while(j<trail_lim.size() && i>=trail_lim[j]){
				lev++;
				j++;
			}
			AssignedEdge e = trail[i];
			Lit l = mkLit( e.var,!e.assign);
			assert(S->value(l)==l_True);
			int expected_level = S->level(var(l));
			assert(S->level(var(l))==lev);
			int edge_num = e.var-min_edge_var;
			assert(assigned[edge_num]==l_Undef);
			assigned[edge_num] = sign(l)?l_False:l_True;
		}

		for(int i = 0;i<S->trail.size();i++){

			Lit l = S->trail[i];
			Var v = var(l);

			int lev = S->level(v);

			if(v>= min_edge_var && v<min_edge_var+num_edges){
				int edge_num = v-min_edge_var;
				lbool assigned_val=assigned[edge_num];
				assert(assigned_val== (sign(l)?l_False:l_True));
			}
		}

#endif
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
				if(e.from==47 && e.to==55){
									int a =1;
								}
				if(e.var==368){
					int  a=1;
				}
				assert(S->value(e.var)==l_Undef);
				if(e.assign){
					g.removeEdge(e.from,e.to);
				}else{
					antig.addEdge(e.from,e.to);
					assert(antig.hasEdge(e.from,e.to));
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
			if(reach_detectors[i]->positive_reach_detector)
				reach_detectors[i]->positive_reach_detector->update();
			if(reach_detectors[i]->negative_reach_detector)
				reach_detectors[i]->negative_reach_detector->update();
		}

		dbg_sync();


	};
	void newDecisionLevel(){
		trail_lim.push(trail.size());
	};

	void buildReason(Lit p, vec<Lit> & reason){
		CRef marker = S->reason(var(p));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef- marker;
		int d = marker_map[pos];
		reason.push(p);
		assert(d!=0);
		if(d>0){
			double startpathtime = cpuTime();
			d--;
			Dijkstra & detector = *reach_detectors[d]->positive_reach_detector;
			//the reason is a path from s to p(provided by d)
			//p is the var for a reachability detector in dijkstra, and corresponds to a node
			Var v = var(p);
			int u =reach_detectors[d]->getNode(v); //reach_detectors[d]->reach_lit_map[v];
			assert(detector.connected(u));
			int w;
			while(( w= detector.previous(u)) > -1){
				Lit l = mkLit( edges[w][u].v,true );
				assert(S->value(l)==l_False);
				reason.push(l);
				u=w;
			}
#ifdef DEBUG_GRAPH
		 assert(dbg_clause(reason));

#endif
			double elapsed = cpuTime()-startpathtime;
			pathtime+=elapsed;
		}else{
			d=-d-1;
			Dijkstra & detector = *reach_detectors[d]->negative_reach_detector;

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
			int t = reach_detectors[d]->getNode(v); // v- var(reach_lits[d][0]);
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
	void buildReachReason(int node,Dijkstra & d,vec<Lit> & conflict){
		drawFull();
		assert(dbg_reachable(d.source,node));
		double starttime = cpuTime();
		int u = node;
		int p;
		while(( p = d.previous(u)) != -1){
			Edge edg = edges[p][u];
			Var e =edges[p][u].v;
			lbool val = S->value(e);
			assert(S->value(e)==l_True);
			conflict.push(mkLit(e, true));
			u = p;
		}
		double elapsed = cpuTime()-starttime;
		pathtime+=elapsed;
#ifdef DEBUG_GRAPH
		 assert(dbg_clause(conflict));

#endif
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
		assert(dbg_notreachable( reach_detectors[detector]->source,u));
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
		int f =mc.minCut(reach_detectors[detector]->source,node,cut);
		assert(f<0xF0F0F0); assert(f==cut.size());//because edges are only ever infinity or 1
		for(int i = 0;i<cut.size();i++){
			EdmondsKarp::Edge e = cut[i];

			Lit l = mkLit( edges[e.u][e.v].v,false);
			assert(S->value(l)==l_False);
			conflict.push(l);
		}
		double elapsed = cpuTime()-starttime;
			mctime+=elapsed;

#ifdef DEBUG_GRAPH
		 assert(dbg_clause(conflict));
#endif
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

		//At level 0, need to propagate constant reaches/source nodes/edges...
		if(local_q==0){
			assert(S->decisionLevel()==0);
			for(int i = 0;i<reach_detectors.size();i++){
				detectors_to_check.push((i+1));
				reach_detectors[i]->positive_reach_detector->marked =true;
			}
			for(int i = 0;i<reach_detectors.size();i++){
				detectors_to_check.push(-(i+1));
				reach_detectors[i]->negative_reach_detector->marked =true;
			}
		}

		while(local_q<S->qhead){
			Lit l = S->trail[local_q++];
			Var v = var(l);

			int lev = S->level(v);
			while(lev>trail_lim.size()){
				newDecisionLevel();
			}
			if(v>= min_edge_var && v<min_edge_var+num_edges){
				if(v==368){
							int a=1;
						}
				//this is an edge assignment
				int edge_num = v-min_edge_var;
				int from = edge_list[edge_num].from;
				int to = edge_list[edge_num].to;
				trail.push({!sign(l), from,to,v});
				AssignedEdge e = trail.last();
				assert(e.from==from);
				assert(e.to==to);
				if(from==47 && to==55){
					int a =1;
				}
				if (!sign(l)){

					g.addEdge(from,to);
					for(int i = 0;i<reach_detectors.size();i++){
						if(!reach_detectors[i]->positive_reach_detector->marked && (reach_detectors[i]->positive_reach_detector->connected_unsafe(from)||reach_detectors[i]->positive_reach_detector->connected_unsafe(to))){
								detectors_to_check.push((i+1));
								reach_detectors[i]->positive_reach_detector->marked =true;
							}
					}
				}else{
					antig.removeEdge(from,to);
					for(int i = 0;i<reach_detectors.size();i++){

						if(!reach_detectors[i]->negative_reach_detector->marked && ( reach_detectors[i]->negative_reach_detector->connected_unsafe(from)||reach_detectors[i]->negative_reach_detector->connected_unsafe(to))){
							detectors_to_check.push(-(i+1));
							reach_detectors[i]->negative_reach_detector->marked =true;
						}
					}
				}
			}
		}
		dbg_sync();
			assert(dbg_graphsUpToDate());
			for(int i = 0;i<detectors_to_check.size();i++){
				int d = detectors_to_check[i];
				assert(d!=0);
				if(d>0){
					d--;
					double startdreachtime = cpuTime();
					reach_detectors[d]->positive_reach_detector->marked =false;
					reach_detectors[d]->positive_reach_detector->update();
					double reachUpdateElapsed = cpuTime()-startdreachtime;
					reachupdatetime+=reachUpdateElapsed;
				/*	for(int j =0;j<reach_lits[d].size();j++){
						Lit l = reach_lits[d][j];
						if( reach_detectors[d]->connected(j)){*/
					for(int j = 0;j<reach_detectors[d]->positive_reach_detector->getChanged().size();j++){
							int u = reach_detectors[d]->positive_reach_detector->getChanged()[j];
							assert(reach_detectors[d]->reach_lits[u]!=lit_Undef);
							Lit l = (reach_detectors[d]->reach_lits[u]); // reach_lits[d][u];

							if( reach_detectors[d]->positive_reach_detector->connected(u)){
								assert(dbg_reachable( reach_detectors[d]->source,u));
							if(S->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
								assert(dbg_propgation(l));
#endif
								S->uncheckedEnqueue(l,reach_markers[d]) ;
							}else if (S->value(l)==l_False){
								double elapsed = cpuTime()-startdreachtime;
															reachtime+=elapsed;
								//conflict
								//The reason is all the literals in the shortest path in to s in d
								conflict.push(l);
								buildReachReason(u,*reach_detectors[d]->positive_reach_detector,conflict);
								//add it to s
								//return it as a conflict
#ifdef DEBUG_GRAPH
								for(int i = 0;i<conflict.size();i++)
											 assert(S->value(conflict[i])==l_False);
#endif

								return false;
							}else{
								int  a=1;
							}
						}else{
							assert(!dbg_reachable( reach_detectors[d]->source,u));
						}

					}
					reach_detectors[d]->positive_reach_detector->clearChanged();
					double elapsed = cpuTime()-startdreachtime;
								reachtime+=elapsed;
				}else{
					d=-d-1;
					double startunreachtime = cpuTime();
					reach_detectors[d]->negative_reach_detector->marked =false;
					reach_detectors[d]->negative_reach_detector->update();
					double unreachUpdateElapsed = cpuTime()-startunreachtime;
					unreachupdatetime+=unreachUpdateElapsed;
					assert(dbg_graphsUpToDate());
					for(int j = 0;j<reach_detectors[d]->negative_reach_detector->getChanged().size();j++){
							int u = reach_detectors[d]->negative_reach_detector->getChanged()[j];
							assert(reach_detectors[d]->reach_lits[u]!=lit_Undef);
							Lit l = ~ (reach_detectors[d]->reach_lits[u]); // ~reach_lits[d][u];
							if(! reach_detectors[d]->negative_reach_detector->connected(u)){
								assert(dbg_notreachable( reach_detectors[d]->source,u));
	/*				for(int j =0;j<reach_lits[d].size();j++){
						Lit l = ~reach_lits[d][j];
						if(! non_reach_detectors[d]->connected(j)){*/
							if(S->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
								assert(dbg_propgation(l));
#endif
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
#ifdef DEBUG_GRAPH
								for(int i = 0;i<conflict.size();i++)
											 assert(S->value(conflict[i])==l_False);
#endif
								return false;
							}
						}

					}
					reach_detectors[d]->negative_reach_detector->clearChanged();
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
					dbg_sync();
		return true;
	};
	bool solveTheory(vec<Lit> & conflict){return true;};
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
		if(from==55 && to==1){
			int a =1;
		}
		if(v==var_Undef)
			v = S->newVar();
#ifdef DEBUG_GRAPH
		 dbg_graph->newEdge(from,to,v);
#endif
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

	void reaches(int from, int to, Var reach_var,int within_steps=-1){

#ifdef DEBUG_GRAPH
		 dbg_graph->reaches(from,  to,reach_var,within_steps);
#endif

			assert(from<g.nodes);

			if (reach_info[from].source<0){

				reach_markers.push(S->newReasonMarker(this));
				int mnum = CRef_Undef- reach_markers.last();
				marker_map.growTo(mnum+1);
				marker_map[mnum] = reach_markers.size();
				//marker_map.insert(reach_markers.last(),reach_markers.size());

				non_reach_markers.push(S->newReasonMarker(this));
				//marker_map[non_reach_markers.last()]=-non_reach_markers.size();
				//marker_map.insert(non_reach_markers.last(),non_reach_markers.size());

				mnum = CRef_Undef- non_reach_markers.last();
				marker_map.growTo(mnum+1);
				marker_map[mnum] = -non_reach_markers.size();

				reach_detectors.push(new ReachDetector());

				reach_detectors.last()->positive_reach_detector = new Dijkstra(from,g);
				reach_detectors.last()->negative_reach_detector = new Dijkstra(from,antig);
				reach_detectors.last()->source=from;

				reach_info[from].source=from;
				reach_info[from].detector=reach_detectors.last();

				reach_detectors.last()->within=within_steps;

				reach_detectors.last()->positive_reach_detector->update();
				reach_detectors.last()->negative_reach_detector->update();

				reach_detectors.last()->first_reach_var = reach_var;

			}

			ReachDetector * d = reach_info[from].detector;
			assert(d);



			while( d->reach_lits.size()<=to)
				d->reach_lits.push(lit_Undef);

			while(S->nVars()<=reach_var)
				S->newVar();

			Lit reachLit=mkLit(reach_var,false);

			if(d->reach_lits[to]==lit_Undef){
				d->reach_lits[to] = reachLit;

				while(d->reach_lit_map.size()<= reach_var- d->first_reach_var ){
					d->reach_lit_map.push(-1);
				}

				d->reach_lit_map[reach_var-d->first_reach_var]=to;

			}else{
				Lit r = d->reach_lits[to];
				//force equality between the new lit and the old reach lit, in the SAT solver
				S->addClause(~r, reachLit);
				S->addClause(r, ~reachLit);
			}




	    }

	void reachesAny(int from, Var firstVar,int within_steps=-1){

		for(int i = 0;i<g.nodes;i++){
			reaches(from,i,firstVar+i,within_steps);
		}

	  }
	void reachesAny(int from, vec<Lit> & reachlits_out,int within_steps=-1){
		for(int i = 0;i<g.nodes;i++){
			Var reachVar = S->newVar();
			reaches(from,i,reachVar,within_steps);
			reachlits_out.push(mkLit(reachVar,false));
		}
    }

};

};

#endif /* DGRAPH_H_ */
