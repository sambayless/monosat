/*
 * SimpleGraph.h
 *
 *  Created on: 2013-04-14
 *      Author: sam
 */

#ifndef DGRAPH_H_
#define DGRAPH_H_

#include "utils/System.h"
#include "core/Theory.h"
#include "Graph.h"
#include "Reach.h"
#include "Dijkstra.h"
#include "Connectivity.h"
#include "Distance.h"
#include "core/SolverTypes.h"
#include "mtl/Map.h"
#include "MaxFlow.h"
#include "IBFS.h"
#include "EdmondsKarp.h"
#include "EdmondsKarpAdj.h"
#include "Chokepoint.h"

#include "utils/System.h"
#ifdef DEBUG_GRAPH
#include "TestGraph.h"
#endif

#ifdef DEBUG_SOLVER
#include "TestGraph.h"
#endif
namespace Minisat{



class GraphTheorySolver:public GraphTheory{
private:


	Lit False;
	Lit True;
	int local_q;

	Solver * S;
#ifdef DEBUG_GRAPH
	Solver * dbg;
	TestGraph *dbg_graph;
#endif
#ifdef DEBUG_SOLVER
	TestGraph * shadow_dbg;
#endif


	vec<lbool> edge_assignments;

	typedef DefaultEdgeStatus PositiveEdgeStatus;
	typedef DefaultEdgeStatus NegativeEdgeStatus;
	typedef DefaultEdgeStatus CutEdgeStatus;


	PositiveEdgeStatus g_status;
	NegativeEdgeStatus antig_status;
	CutEdgeStatus cutGraph_status;

	DynamicGraph<PositiveEdgeStatus> g;
	DynamicGraph<NegativeEdgeStatus> antig;
	DynamicGraph<CutEdgeStatus> cutGraph;

	Var min_edge_var;
	int num_edges;
	struct Assignment{
		bool isEdge:1;
		bool assign:1;
		int from:30;
		int to;
		Var var;
	};

	vec<Assignment> trail;
	vec<int> trail_lim;

	struct ReachDetector{
		GraphTheorySolver & outer;
		int within;
		int source;

		CRef reach_marker;
		CRef non_reach_marker;
		CRef forced_reach_marker;

		Reach * positive_reach_detector;
		Reach * negative_reach_detector;
		Reach *  positive_path_detector;

		vec<Lit>  reach_lits;
		Var first_reach_var;
		vec<int> reach_lit_map;
		vec<int> force_reason;

		struct DistLit{
			Lit l;
			int min_distance;

		};

		vec<vec<DistLit> > dist_lits;


		struct Change{
			Lit l;
			int u;
		};

		vec<Change> changed;

		vec<Change> & getChanged(){
			return changed;
		}

		struct ReachStatus{
			ReachDetector & outer;
			bool polarity;
			void setReachable(int u, bool reachable){
				if(polarity==reachable && u<outer.reach_lits.size()){
					Lit l = outer.reach_lits[u];
					if(l!=lit_Undef){
						lbool assign = outer.outer.S->value(l);
						if(assign!= (reachable? l_True:l_False )){
							 outer.changed.push({reachable? l:~l,u});
						}
					}
				}
			}
			bool isReachable(int u) const{
				return false;
			}

			void setMininumDistance(int u, bool reachable, int distance){
				assert(reachable ==(distance<outer.outer.g.nodes));
				setReachable(u,reachable);

				if(u<outer.dist_lits.size()){
					assert(distance>=0);
					assert(distance<outer.outer.g.nodes);
					for(int i = 0;i<outer.dist_lits[u].size();i++){
						int d =  outer.dist_lits[u][i].min_distance;
						Lit l = outer.dist_lits[u][i].l;
						assert(l!=lit_Undef);
						if(d<distance && !polarity){
							lbool assign = outer.outer.S->value(l);
							if( assign!= l_False ){
								 outer.changed.push({~l,u});
							}
						}else if(d>=distance && polarity){
							lbool assign = outer.outer.S->value(l);
							if( assign!= l_True ){
								 outer.changed.push({l,u});
							}
						}
					}
				}
			}


			ReachStatus(ReachDetector & _outer, bool _polarity):outer(_outer), polarity(_polarity){}
		};
		ReachStatus *positiveReachStatus;
		ReachStatus *negativeReachStatus;


		struct ChokepointStatus{
			ReachDetector & outer;
			bool mustReach(int node){
				Lit l =  outer.reach_lits[node];
				if(l!=lit_Undef){
					return outer.outer.S->value(l)==l_True;
				}
				return false;
			}
			bool operator() (int edge_id){
				return outer.outer.edge_assignments[edge_id]==l_Undef;
			}
			ChokepointStatus(ReachDetector & _outer):outer(_outer){

			}
		}chokepoint_status;

		Chokepoint<ChokepointStatus ,NegativeEdgeStatus> chokepoint;

		int getNode(Var reachVar){
			assert(reachVar>=first_reach_var);
			int index = reachVar-first_reach_var;
			assert(index< reach_lit_map.size());
			assert(reach_lit_map[index]>=0);
			return reach_lit_map[index];
		}

	/*	Lit getLit(int node){

			return reach_lits[node];

		}*/

		void buildReachReason(int node,vec<Lit> & conflict){
			//drawFull();
			Reach & d = *positive_path_detector;


			double starttime = cpuTime();
			d.update();


			if(!outer.dbg_reachable(d.getSource(),node)){
				outer.drawFull();

				d.drawFull();

				assert(false);
			}
			assert(d.connected_unchecked(node));
			if(opt_learn_reaches ==0 || opt_learn_reaches==2){
				int u = node;
				int p;
				while(( p = d.previous(u)) != -1){
					Edge edg = outer.edges[p][u];
					Var e =outer.edges[p][u].v;
					lbool val = outer.S->value(e);
					assert(outer.S->value(e)==l_True);
					conflict.push(mkLit(e, true));
					u = p;

				}
			}else{
				//Instead of a complete path, we can learn reach variables, if they exist
				int u = node;
				int p;
				while(( p = d.previous(u)) != -1){
					Edge edg = outer.edges[p][u];
					Var e =outer.edges[p][u].v;
					lbool val = outer.S->value(e);
					assert(outer.S->value(e)==l_True);
					conflict.push(mkLit(e, true));
					u = p;
					if( u< reach_lits.size() && reach_lits[u]!=lit_Undef && outer.S->value(reach_lits[u])==l_True && outer.S->level(var(reach_lits[u]))< outer.S->decisionLevel()){
						//A potential (fixed) problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
						//This is avoided by ensuring that L is lower level than the conflict.
						Lit l =reach_lits[u];
						assert(outer.S->value(l)==l_True);
						conflict.push(~l);
						break;
					}
				}

			}
			outer.num_learnt_paths++;
			outer.learnt_path_clause_length+= (conflict.size()-1);
			double elapsed = cpuTime()-starttime;
			outer.pathtime+=elapsed;
	#ifdef DEBUG_GRAPH
			 assert(outer.dbg_clause(conflict));

	#endif
		}
		void buildNonReachReason(int node,vec<Lit> & conflict){
			static int it = 0;
			++it;
			int u = node;
			//drawFull( non_reach_detectors[detector]->getSource(),u);
			assert(outer.dbg_notreachable( source,u));
			double starttime = cpuTime();
			outer.cutGraph.clearHistory();
			outer.stats_mc_calls++;
			if(opt_conflict_min_cut){
				if(mincutalg!= ALG_EDKARP_ADJ){
					//ok, set the weights for each edge in the cut graph.
					//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
					for(int u = 0;u<outer.cutGraph.adjacency.size();u++){
						for(int j = 0;j<outer.cutGraph.adjacency[u].size();j++){
							int v = outer.cutGraph.adjacency[u][j].node;
							Var var = outer.edges[u][v].v;
							/*if(S->value(var)==l_False){
								mc.setCapacity(u,v,1);
							}else{*/
							outer.mc->setCapacity(u,v,0xF0F0F0);
							//}
						}
					}

					//find any edges assigned to false, and set their capacity to 1
					for(int i =0;i<outer.trail.size();i++){
						if(outer.trail[i].isEdge && !outer.trail[i].assign){
							outer.mc->setCapacity(outer.trail[i].from, outer.trail[i].to,1);
						}
					}
				}
				outer.cut.clear();

				int f =outer.mc->minCut(source,node,outer.cut);
				assert(f<0xF0F0F0); assert(f==outer.cut.size());//because edges are only ever infinity or 1
				for(int i = 0;i<outer.cut.size();i++){
					MaxFlow::Edge e = outer.cut[i];

					Lit l = mkLit( outer.edges[e.u][e.v].v,false);
					assert(outer.S->value(l)==l_False);
					conflict.push(l);
				}
			}else{
				//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
				//or perhaps we can learn the actual 1-uip cut?


					vec<int>& to_visit  = outer.to_visit;
					vec<char>& seen  = outer.seen;

				    to_visit.clear();
				    to_visit.push(node);
				    seen.clear();
					seen.growTo(outer.nNodes());
				    seen[node]=true;

				    do{

				    	assert(to_visit.size());
				    	int u = to_visit.last();
				    	assert(u!=source);
				    	to_visit.pop();
				    	assert(seen[u]);
				    	assert(!negative_reach_detector->connected_unsafe(u));
				    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
				    	for(int i = 0;i<outer.inv_adj[u].size();i++){
				    		int v = outer.inv_adj[u][i].v;
				    		int from = outer.inv_adj[u][i].from;
				    		assert(from!=u);
				    		assert(outer.inv_adj[u][i].to==u);
				    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				    		int edge_num = v-outer.min_edge_var;

				    		if(outer.edge_assignments[edge_num]==l_False){
				    			//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
				    			//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
				    			//if(!seen[from])
				    				conflict.push(mkLit(v,false));
				    		}else{
				    			assert(from!=source);
				    			//even if it is undef? probably...
				    			if(!seen[from]){
				    				seen[from]=true;
				    				if((opt_learn_reaches ==2 || opt_learn_reaches==3) && from< reach_lits.size() && reach_lits[from]!=lit_Undef && outer.S->value(reach_lits[from])==l_False  && outer.S->level(var(reach_lits[from]))< outer.S->decisionLevel())
				    				{
				    				//The problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
				    					Lit r = reach_lits[from];
				    					assert(var(r)<outer.S->nVars());
				    					assert(outer.S->value(r)==l_False);
				    					conflict.push(r);
				    				}else
				    					to_visit.push(from);
				    			}
				    		}
				    	}
				    }  while (to_visit.size());


			}

			 outer.num_learnt_cuts++;
			 outer.learnt_cut_clause_length+= (conflict.size()-1);

			double elapsed = cpuTime()-starttime;
			 outer.mctime+=elapsed;

	#ifdef DEBUG_GRAPH
			 assert(outer.dbg_clause(conflict));
	#endif
		}

		/**
		 * Explain why an edge was forced (to true).
		 * The reason is that _IF_ that edge is false, THEN there is a cut of disabled edges between source and target
		 * So, create the graph that has that edge (temporarily) assigned false, and find a min-cut in it...
		 */
		void buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict){
					static int it = 0;
					++it;

					assert(outer.edge_assignments[forced_edge_id]==l_True);
					Lit edgeLit =mkLit( outer.edge_list[forced_edge_id].v,false);

					conflict.push(edgeLit);

					int forced_edge_from = outer.edge_list[forced_edge_id].from;
					int forced_edge_to= outer.edge_list[forced_edge_id].to;


					int u = reach_node;
					//drawFull( non_reach_detectors[detector]->getSource(),u);
					assert(outer.dbg_notreachable( source,u));
					double starttime = cpuTime();
					outer.cutGraph.clearHistory();
					outer.stats_mc_calls++;



					if(opt_conflict_min_cut){
						if(mincutalg!= ALG_EDKARP_ADJ){
							//ok, set the weights for each edge in the cut graph.
							//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
							for(int u = 0;u<outer.cutGraph.adjacency.size();u++){
								for(int j = 0;j<outer.cutGraph.adjacency[u].size();j++){
									int v = outer.cutGraph.adjacency[u][j].node;
									Var var = outer.edges[u][v].v;
									/*if(S->value(var)==l_False){
										mc.setCapacity(u,v,1);
									}else{*/
									outer.mc->setCapacity(u,v,0xF0F0F0);
									//}
								}
							}

							//find any edges assigned to false, and set their capacity to 1
							for(int i =0;i<outer.trail.size();i++){
								if(outer.trail[i].isEdge && !outer.trail[i].assign){
									outer.mc->setCapacity(outer.trail[i].from, outer.trail[i].to,1);
								}
							}

							outer.mc->setCapacity(forced_edge_from, forced_edge_to,1);
						}
						outer.cut.clear();

						int f =outer.mc->minCut(source,reach_node,outer.cut);
						assert(f<0xF0F0F0); assert(f==outer.cut.size());//because edges are only ever infinity or 1
						for(int i = 0;i<outer.cut.size();i++){
							MaxFlow::Edge e = outer.cut[i];

							Lit l = mkLit( outer.edges[e.u][e.v].v,false);
							assert(outer.S->value(l)==l_False);
							conflict.push(l);
						}
					}else{
						//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
						//or perhaps we can learn the actual 1-uip cut?


							vec<int>& to_visit  = outer.to_visit;
							vec<char>& seen  = outer.seen;

						    to_visit.clear();
						    to_visit.push(reach_node);
						    seen.clear();
							seen.growTo(outer.nNodes());
						    seen[reach_node]=true;

						    do{

						    	assert(to_visit.size());
						    	int u = to_visit.last();
						    	assert(u!=source);
						    	to_visit.pop();
						    	assert(seen[u]);
						    	assert(!negative_reach_detector->connected_unsafe(u));
						    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
						    	for(int i = 0;i<outer.inv_adj[u].size();i++){
						    		int v = outer.inv_adj[u][i].v;
						    		int from = outer.inv_adj[u][i].from;
						    		assert(from!=u);
						    		assert(outer.inv_adj[u][i].to==u);
						    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
						    		int edge_num = v-outer.min_edge_var;

						    		if(edge_num == forced_edge_id || outer.edge_assignments[edge_num]==l_False){
						    			//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
						    			//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
						    			//if(!seen[from])
						    				conflict.push(mkLit(v,false));
						    		}else{
						    			assert(from!=source);
						    			//even if it is undef? probably...
						    			if(!seen[from]){
						    				seen[from]=true;
						    				if((opt_learn_reaches ==2 || opt_learn_reaches==3) && from< reach_lits.size() && reach_lits[from]!=lit_Undef && outer.S->value(reach_lits[from])==l_False  && outer.S->level(var(reach_lits[from]))< outer.S->decisionLevel())
						    				{
						    				//The problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
						    					Lit r = reach_lits[from];
						    					assert(var(r)<outer.S->nVars());
						    					assert(outer.S->value(r)==l_False);
						    					conflict.push(r);
						    				}else
						    					to_visit.push(from);
						    			}
						    		}
						    	}
						    }  while (to_visit.size());


					}

					 outer.num_learnt_cuts++;
					 outer.learnt_cut_clause_length+= (conflict.size()-1);

					double elapsed = cpuTime()-starttime;
					 outer.mctime+=elapsed;

			#ifdef DEBUG_GRAPH
					 assert(outer.dbg_clause(conflict));
			#endif
				}

		void buildReason(Lit p, vec<Lit> & reason, CRef marker){


				if(marker==reach_marker){
					reason.push(p);
				//	double startpathtime = cpuTime();

					/*Dijkstra & detector = *reach_detectors[d]->positive_dist_detector;
					//the reason is a path from s to p(provided by d)
					//p is the var for a reachability detector in dijkstra, and corresponds to a node
					detector.update();
					Var v = var(p);
					int u =reach_detectors[d]->getNode(v); //reach_detectors[d]->reach_lit_map[v];
					assert(detector.connected(u));
					int w;
					while(( w= detector.previous(u)) > -1){
						Lit l = mkLit( edges[w][u].v,true );
						assert(S->value(l)==l_False);
						reason.push(l);
						u=w;
					}*/
					Var v = var(p);
					int u =getNode(v);
					buildReachReason(u,reason);

		#ifdef DEBUG_GRAPH
				 assert(outer.dbg_clause(reason));

		#endif
					//double elapsed = cpuTime()-startpathtime;
				//	pathtime+=elapsed;
				}else if(marker==non_reach_marker){
					reason.push(p);

					//the reason is a cut separating p from s;
					//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

					//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
					//assign the mincut edge weights if they aren't already assigned.


					Var v = var(p);
					int t = getNode(v); // v- var(reach_lits[d][0]);
					buildNonReachReason(t,reason);

				}else if (marker==forced_reach_marker){
					Var v = var(p);
					//The forced variable is an EDGE that was forced.
					int forced_edge_id =v- outer.min_edge_var;
					//The corresponding node that is the reason it was forced
					int reach_node=force_reason[forced_edge_id];
					 buildForcedEdgeReason( reach_node,  forced_edge_id, reason);
				}else{
					assert(false);
				}
		}



		ReachDetector(GraphTheorySolver & _outer, int _source):outer(_outer),within(-1),source(_source),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),chokepoint_status(*this),chokepoint(chokepoint_status, outer.antig,source){}
	};

	struct ReachInfo{
		int source;
		bool distance;
		ReachDetector * detector;

		ReachInfo():source(-1),detector(NULL){}
	};

	vec<ReachInfo> reach_info;
public:
	vec<ReachDetector*> reach_detectors;


	vec<int> marker_map;

	struct Edge{
		Var v;
		int from;
		int to;
	};
	vec<MaxFlow::Edge> cut;

	//Full matrix
	vec<vec<Edge> > edges;

	//Just a list of the edges
	vec<Edge> edge_list;

	vec<vec<Edge> > inv_adj;


	MaxFlow * mc;
	//MaxFlow * reachprop;

    vec<char> seen;
	vec<int> to_visit;


public:

	double mctime;
	double reachtime;
	double unreachtime;
	double pathtime;
	double propagationtime;
	double reachupdatetime;
	double unreachupdatetime;
	double stats_initial_propagation_time;
	int num_learnt_paths;
	int learnt_path_clause_length;
	int num_learnt_cuts;
	int learnt_cut_clause_length;

	int stats_mc_calls;
	vec<Lit> reach_cut;
	vec<ForceReason> forced_edges;
	struct CutStatus{
		GraphTheorySolver & outer;
		int operator () (int id) const {

			if(outer.edge_assignments[id]==l_False){
				return 1;
			}else{
				return 0xF0F0F0;
			}
		}
		CutStatus(GraphTheorySolver & _outer):outer(_outer){}

	} cutStatus;

	struct PropCutStatus{
		GraphTheorySolver & outer;
		int operator () (int id) const {

			if(outer.edge_assignments[id]==l_Undef){
				return 1;
			}else{
				assert(outer.edge_assignments[id]==l_True);
				return 0xF0F0F0;
			}
		}
		PropCutStatus(GraphTheorySolver & _outer):outer(_outer){}

	}propCutStatus;

	GraphTheorySolver(Solver * S_):S(S_),g(g_status),antig(antig_status) ,cutGraph(cutGraph_status),cutStatus(*this),propCutStatus(*this){
		True = mkLit(S->newVar(),false);
			False=~True;
			S->addClause(True);
			num_edges=0;
			local_q=0;
			mctime=0;
			stats_mc_calls=0;
			reachtime=0;
			unreachtime=0;
			pathtime=0;
			propagationtime=0;
			reachupdatetime=0;
			unreachupdatetime=0;
			stats_initial_propagation_time=0;

			 num_learnt_paths=0;
			 learnt_path_clause_length=0;
			 num_learnt_cuts=0;
			 learnt_cut_clause_length=0;

			if(mincutalg==ALG_IBFS){
				mc = new IBFS<CutEdgeStatus>(cutGraph);

			}else if (mincutalg == ALG_EDKARP_ADJ){

				mc = new EdmondsKarpAdj<CutStatus, CutEdgeStatus>(cutGraph,cutStatus);
				//reachprop = new EdmondsKarpAdj<PropCutStatus, NegativeEdgeStatus>(antig,propCutStatus);
			}else{
				mc = new EdmondsKarp<CutEdgeStatus>(cutGraph);
			}

#ifdef DEBUG_GRAPH
		dbg=new Solver();
		dbg_graph = new TestGraph(dbg);
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg = new TestGraph(S->dbg_solver);
#endif
	}

	void printStats(){

		printf("Graph stats:\n");
		printf("Prop Time: %f (initial: %f)\n", propagationtime,stats_initial_propagation_time);
		printf("Reach Time: %f (update Time: %f)\n", reachtime,reachupdatetime);
		printf("Unreach Time: %f (Update Time: %f)\n", unreachtime,unreachupdatetime);
		printf("Path Time: %f (#Paths: %d, AvgLength %f, total: %d)\n", pathtime, num_learnt_paths, (learnt_path_clause_length /  ((float) num_learnt_paths+1)),learnt_path_clause_length);
		printf("Min-cut Time: %f (%d calls, %f average, #Cuts: %d, AvgLength %f, total: %d)\n", mctime, stats_mc_calls,(mctime/(stats_mc_calls ? stats_mc_calls:1)),  num_learnt_cuts, (learnt_cut_clause_length /  ((float) num_learnt_cuts+1)),learnt_cut_clause_length);

		int stats_full_updates=0;
		int stats_fast_updates=0;
		int stats_failed_fast_updates=0;
		int stats_skip_deletes=0;
		int stats_skipped_updates=0;
		int skipable_deletions = 0;
		double stats_full_update_time=0;
		double stats_fast_update_time=0;

		for(int i = 0;i<reach_detectors.size();i++){
			skipable_deletions+=reach_detectors[i]->positive_reach_detector->stats_num_skipable_deletions;
			skipable_deletions+=reach_detectors[i]->negative_reach_detector->stats_num_skipable_deletions;

			stats_failed_fast_updates+=reach_detectors[i]->positive_reach_detector->stats_fast_failed_updates;
			stats_failed_fast_updates+=reach_detectors[i]->negative_reach_detector->stats_fast_failed_updates;

			stats_full_updates+=reach_detectors[i]->positive_reach_detector->stats_full_updates;
			stats_fast_updates+=reach_detectors[i]->positive_reach_detector->stats_fast_updates;
			stats_skip_deletes+=reach_detectors[i]->positive_reach_detector->stats_skip_deletes;
			stats_skipped_updates+=reach_detectors[i]->positive_reach_detector->stats_skipped_updates;

			stats_full_update_time+=reach_detectors[i]->positive_reach_detector->stats_full_update_time;
			stats_fast_update_time+=reach_detectors[i]->positive_reach_detector->stats_fast_update_time;

			stats_full_updates+=reach_detectors[i]->negative_reach_detector->stats_full_updates;
			stats_fast_updates+=reach_detectors[i]->negative_reach_detector->stats_fast_updates;
			stats_skip_deletes+=reach_detectors[i]->negative_reach_detector->stats_skip_deletes;
			stats_skipped_updates+=reach_detectors[i]->negative_reach_detector->stats_skipped_updates;

			stats_full_update_time+=reach_detectors[i]->negative_reach_detector->stats_full_update_time;
			stats_fast_update_time+=reach_detectors[i]->negative_reach_detector->stats_fast_update_time;
		}

		printf("Dijkstra Full Updates: %d (time: %f, average: %f)\n",stats_full_updates, stats_full_update_time,(stats_full_update_time/(stats_full_updates ? stats_full_updates:1)));
		printf("Dijkstra Fast Updates: %d (time: %f, average: %f, failed: %d)\n",stats_fast_updates, stats_fast_update_time,(stats_fast_update_time/(stats_fast_updates ? stats_fast_updates:1)),stats_failed_fast_updates);
		printf("Dijkstra Skipped Updates: %d (deletionSkips: %d, skipable deletes %d)\n",stats_skipped_updates, stats_skip_deletes,skipable_deletions);

	}

     ~GraphTheorySolver(){};
	 int newNode(){
#ifdef DEBUG_GRAPH
		 dbg_graph->newNode();
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg->newNode();
#endif
		 edges.push();
		 for(int i = 0;i<edges.size();i++)
			 edges[i].growTo(edges.size());
		 inv_adj.push();
		 reach_info.push();
		 antig.addNode();
		 cutGraph.addNode();


		 seen.growTo(nNodes());

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
					if(edge_list[v-min_edge_var].v<0)
								continue;
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
	void	dbg_sync_reachability(){
#ifdef DEBUG_GRAPH

		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			for(int j = 0;j< d->reach_lits.size();j++){
				Lit l = d->reach_lits[j];
				if(l!=lit_Undef){
					int node = d->getNode(var(l));

					if(d->positive_reach_detector->connected(node)){
						assert(S->value(l)==l_True);
					}else if(! d->negative_reach_detector->connected(node)){
						assert(S->value(l)==l_False);
					}
				}

			}
		}
#endif
	}
	void dbg_sync(){
#ifdef DEBUG_GRAPH
		static vec<lbool> assigned;
		assigned.clear();
		for(int i = 0;i<edge_list.size();i++)
			assigned.push(l_Undef);
		int lev = 0;
		int j =0;
		for(int i = 0;i<trail.size();i++){
			while(j<trail_lim.size() && i>=trail_lim[j]){
				lev++;
				j++;
			}

			Assignment e = trail[i];
			if(e.isEdge){

				Lit l = mkLit( e.var,!e.assign);
				assert(S->value(l)==l_True);
				int expected_level = S->level(var(l));
				assert(S->level(var(l))==lev);
				int edge_num = e.var-min_edge_var;
				if(edge_list[edge_num].v<0)
							continue;
				assert(assigned[edge_num]==l_Undef);

				assigned[edge_num] = sign(l)?l_False:l_True;
			}else{
				//Reachability assignment


			}
		}

		for(int i = 0;i<assigned.size();i++){
			assert(edge_assignments[i]== assigned[i]);
		}

		for(int i = 0;i<S->trail.size();i++){

			Lit l = S->trail[i];
			Var v = var(l);

			int lev = S->level(v);

			if(v>= min_edge_var && v<min_edge_var+num_edges){
				int edge_num = v-min_edge_var;
				if(edge_list[edge_num].v<0)
					continue;
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
				Assignment e = trail[i];
				if(e.isEdge){
					assert(S->value(e.var)==l_Undef);
					int edge_num = e.var-min_edge_var;
					assert(edge_assignments[edge_num]!=l_Undef);
					edge_assignments[edge_num]=l_Undef;
					if(e.assign){
						g.disableEdge(e.from,e.to, edge_num);
					}else{
						antig.enableEdge(e.from,e.to,edge_num);
						assert(antig.hasEdge(e.from,e.to));
					}
				}
			}
			trail.shrink(trail.size()-stop);
			trail_lim.shrink(trail_lim.size()-level);
			assert(trail_lim.size()==level);


		}

		if(local_q>S->qhead)
			local_q=S->qhead;
		assert(dbg_graphsUpToDate());
		/*for(int i = 0;i<reach_detectors.size();i++){
			if(reach_detectors[i]->positive_reach_detector)
				reach_detectors[i]->positive_reach_detector->update();
			if(reach_detectors[i]->negative_reach_detector)
				reach_detectors[i]->negative_reach_detector->update();
		}*/

		dbg_sync();


	};

	Lit decideTheory(){
		if(!opt_decide_graph)
			return lit_Undef;
		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector * r = reach_detectors[i];

			Distance<ReachDetector::ReachStatus,NegativeEdgeStatus> * over = (Distance<ReachDetector::ReachStatus,NegativeEdgeStatus>*) r->negative_reach_detector;

			Distance<ReachDetector::ReachStatus,PositiveEdgeStatus> * under = (Distance<ReachDetector::ReachStatus,PositiveEdgeStatus>*) r->positive_reach_detector;

			//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

			//this can be obviously more efficient
			for(int j = 0;j<nNodes();j++){
				Lit l = r->reach_lits[j];
				if(S->value(l)==l_True){
					if(!under->connected(j)){
						//then lets try to connect this
						int a = 1;
						assert(over->connected(j));//Else, we would already be in conflict

						//ok, read back the path from the over to find a candidate edge we can decide
						//find the earliest unconnected node on this path
						int p = j;
						int last = j;
						while(!under->connected(p)){
							last=p;
							assert(p!=r->source);
							int prev = over->previous(p);
							p = prev;

						}
						//ok, now pick some edge p->last that will connect p to last;
						assert(!under->connected(last));
						assert(under->connected(p));

						assert(over->connected(last));
						assert(over->connected(p));

						for(int k = 0;k<antig.adjacency[p].size();k++){
							int to = antig.adjacency[p][k].node;
							if (to==last){
								Var v =edge_list[ antig.adjacency[p][k].id].v;
								if(S->value(v)==l_Undef){
									return mkLit(v,false);
								}else{
									assert(S->value(v)!=l_True);
								}
							}
						}
						assert(false);
					}
				}
			}

		}

		return lit_Undef;
	}

	void backtrackUntil(Lit p){
			//need to remove and add edges in the two graphs accordingly.
			int i = trail.size()-1;
			for(;i>=0;i--){
				Assignment e = trail[i];
				if(e.isEdge){
					int edge_num = e.var-min_edge_var;
					assert(edge_assignments[edge_num]!=l_Undef);
					edge_assignments[edge_num]=l_Undef;
					if(e.assign){
						g.disableEdge(e.from,e.to, edge_num);
					}else{
						antig.enableEdge(e.from,e.to,edge_num);
						assert(antig.hasEdge(e.from,e.to));
					}
				}else{
					if(var(p)==e.var){
						assert(sign(p)!=e.assign);
						break;
					}
				}
			}

			trail.shrink(trail.size()-(i+1));
			//while(trail_lim.size() && trail_lim.last()>=trail.size())
			//	trail_lim.pop();

	/*		for(int i = 0;i<reach_detectors.size();i++){
					if(reach_detectors[i]->positive_reach_detector)
						reach_detectors[i]->positive_reach_detector->update();
					if(reach_detectors[i]->negative_reach_detector)
						reach_detectors[i]->negative_reach_detector->update();
				}*/
		};

	void newDecisionLevel(){
		trail_lim.push(trail.size());
	};

	void buildReason(Lit p, vec<Lit> & reason){
		CRef marker = S->reason(var(p));
		assert(marker != CRef_Undef);
		int pos = CRef_Undef- marker;
		int d = marker_map[pos];

		backtrackUntil(p);


		assert(d<reach_detectors.size());
		reach_detectors[d]->buildReason(p,reason,marker);


	}



	bool dbg_reachable(int from, int to){
#ifdef DEBUG_GRAPH
		DefaultEdgeStatus tmp;
		/*DynamicGraph<> gtest(tmp);
		for(int i = 0;i<nNodes();i++){
			gtest.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
				continue;
			Edge e  = edge_list[i];
			if(S->assigns[e.v]==l_True){
				gtest.addEdge(e.from,e.to);
				assert(g.edgeEnabled(e.v));
			}else{
				assert(!g.edgeEnabled(e.v));
			}
		}*/

		Dijkstra<> d(from,g);
		d.update();
		return d.connected(to);
#else
		return true;
#endif
	}

	bool dbg_notreachable(int from, int to){

#ifdef DEBUG_GRAPH
		//drawFull(from,to);
		DefaultEdgeStatus tmp;
		DynamicGraph<> g(tmp);
		for(int i = 0;i<nNodes();i++){
			g.addNode();
		}

		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge e  = edge_list[i];
			if(S->assigns[e.v]!=l_False){
				g.addEdge(e.from,e.to);
			}
		}

		Dijkstra<> d(from,g);

		return !d.connected(to);
#else
		return true;
#endif
	}

	bool dbg_graphsUpToDate(){
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
				continue;
			Edge e = edge_list[i];
			lbool val = S->value(e.v);
			assert(edge_assignments[i]==val);
			if(val==l_True || val==l_Undef){
				assert(antig.edgeEnabled(i));
				assert(antig.hasEdge(e.from,e.to));
			}else{
				assert(!antig.edgeEnabled(i));
				assert(!antig.hasEdge(e.from,e.to));
			}
			if(val==l_True){
				assert(g.edgeEnabled(i));
				assert(g.hasEdge(e.from,e.to));
			}else{
				assert(!g.edgeEnabled(i));
				assert(!g.hasEdge(e.from,e.to));
			}
		}

#endif
		return true;
	}


	bool propagateTheory(vec<Lit> & conflict){
		static int itp = 0;
		if(	++itp==2){
			int a =1;
		}

		bool any_change = false;
		double startproptime = cpuTime();
		static vec<int> detectors_to_check;

		conflict.clear();
		//Can probably speed this up alot by a) constant propagating reaches that I care about at level 0, and b) Removing all detectors for nodes that appear only in the opposite polarity (or not at all) in the cnf.
		//That second one especially.

		//At level 0, need to propagate constant reaches/source nodes/edges...

		while(local_q<S->qhead){
			Lit l = S->trail[local_q++];
			Var v = var(l);

			int lev = S->level(v);
			while(lev>trail_lim.size()){
				newDecisionLevel();
			}
			if(v>= min_edge_var && v<min_edge_var+edge_list.size()){

				//this is an edge assignment
				int edge_num = v-min_edge_var;
				if(edge_list[edge_num].v<0)
					continue;
				int from = edge_list[edge_num].from;
				int to = edge_list[edge_num].to;
				trail.push({true,!sign(l), from,to,v});
				assert(edge_assignments[edge_num]==l_Undef);
				edge_assignments[edge_num]=sign(l) ? l_False:l_True;
				Assignment e = trail.last();
				assert(e.from==from);
				assert(e.to==to);

				if (!sign(l)){

					g.enableEdge(from,to,edge_num);
				}else{
					antig.disableEdge(from,to,edge_num);

				}
			}
		}
		stats_initial_propagation_time += cpuTime() - startproptime;
		dbg_sync();
			assert(dbg_graphsUpToDate());

			for(int d = 0;d<reach_detectors.size();d++){
				{
					double startdreachtime = cpuTime();
					reach_detectors[d]->changed.clear();
					reach_detectors[d]->positive_reach_detector->update();
					double reachUpdateElapsed = cpuTime()-startdreachtime;
					reachupdatetime+=reachUpdateElapsed;

					double startunreachtime = cpuTime();
					reach_detectors[d]->negative_reach_detector->update();
					double unreachUpdateElapsed = cpuTime()-startunreachtime;
					unreachupdatetime+=unreachUpdateElapsed;

					for(int j = 0;j<reach_detectors[d]->getChanged().size();j++){
							Lit l = reach_detectors[d]->getChanged()[j].l;
							int u =  reach_detectors[d]->getChanged()[j].u;
							bool reach = !sign(l);
							if(S->value(l)==l_True){
								//do nothing
							}else if(S->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
								assert(dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
								if(S->dbg_solver)
									S->dbg_check_propagation(l);
#endif
								trail.push({false,reach,d,0,var(l)});
								if(reach)
									S->uncheckedEnqueue(l,reach_detectors[d]->reach_marker) ;
								else
									S->uncheckedEnqueue(l,reach_detectors[d]->non_reach_marker) ;

							}else if (S->value(l)==l_False){
								conflict.push(l);
								propagationtime+= cpuTime()-startproptime;
								if(reach){

								//conflict
								//The reason is a path in g from to s in d
									reach_detectors[d]->buildReachReason(u,conflict);
								//add it to s
								//return it as a conflict

								}else{
									//The reason is a cut separating s from t
									reach_detectors[d]->buildNonReachReason(u,conflict);

								}
#ifdef DEBUG_GRAPH
								for(int i = 0;i<conflict.size();i++)
											 assert(S->value(conflict[i])==l_False);
#endif
#ifdef DEBUG_SOLVER
								if(S->dbg_solver)
									S->dbg_check(conflict);
#endif


								return false;
							}else{
								int  a=1;
							}

							if(opt_reach_prop){
								forced_edges.clear();
								reach_detectors[d]->chokepoint.collectForcedEdges(forced_edges);
								for(int i = 0;i<forced_edges.size();i++){
									int edge_id = forced_edges[i].edge_id;
									int node = forced_edges[i].node;
									Lit l = mkLit( edge_list[edge_id].v,false);
									if(S->value(l)==l_Undef){
										reach_detectors[d]->force_reason.growTo(edge_id+1);
										reach_detectors[d]->force_reason[edge_id]=node;
										S->enqueue(l,reach_detectors[d]->forced_reach_marker);
									}else if(S->value(l)==l_True){
										//do nothing

									}else{
										//conflict.
										//this actually shouldn't be possible (at this point in the code)
										reach_detectors[d]->buildForcedEdgeReason(node,edge_id,conflict);
										return false;
									}
								}

							}

						}

#ifdef DEBUG_GRAPH
		for(int i = 0;i<reach_detectors[d]->reach_lits.size();i++){
			Lit l = reach_detectors[d]->reach_lits[i];
			if(l!=lit_Undef){
				int u = reach_detectors[d]->getNode(var(l));
				if(reach_detectors[d]->positive_reach_detector->connected(u)){
					assert(S->value(l)==l_True);
				}else if (!reach_detectors[d]->negative_reach_detector->connected(u)){
					assert(S->value(l)==l_False);
				}
			}

		}
#endif
					}

				}






		g.clearHistory();
		antig.clearHistory();

		detectors_to_check.clear();

		double elapsed = cpuTime()-startproptime;
					propagationtime+=elapsed;
					dbg_sync();
					dbg_sync_reachability();
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
				if(edge_list[i].v<0)
							continue;
				Edge & e = edge_list[i];
				const char * s = "black";
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
			const char * s = "black";
			if(S->value(e.v)==l_True)
				s="blue";
			else if (S->value(e.v)==l_False)
				s="red";
			else{
				int  a=1;
			}
			printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", e.from,e.to, e.v, s);
		}

		printf("}\n");
	}

	bool check_solved(){
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge & e = edge_list[i];
			lbool val = S->value(e.v);
			if(val==l_Undef){
				return false;
			}

			if(val==l_True){
				if(!g.hasEdge(e.from,e.to)){
					return false;
				}
				if(!antig.hasEdge(e.from,e.to)){
					return false;
				}
			}else{
				if(g.hasEdge(e.from,e.to)){
					return false;
				}
				if(antig.hasEdge(e.from,e.to)){
					return false;
				}

			}


		}
		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			for(int j = 0;j< d->reach_lits.size();j++){
				Lit l = d->reach_lits[j];
				if(l!=lit_Undef){
					int node = d->getNode(var(l));

					if(S->value(l)==l_True){
						if(!d->positive_reach_detector->connected(node)){
							return false;
						}
					}else if (S->value(l)==l_False){
						if( d->negative_reach_detector->connected(node)){
							return false;
						}
					}else{
						if(d->positive_reach_detector->connected(node)){
							return false;
						}
						if(!d->negative_reach_detector->connected(node)){
							return false;
						}
					}
				}
			}
		}
		return true;
	}

	bool dbg_solved(){
#ifdef DEBUG_GRAPH
		for(int i = 0;i<edge_list.size();i++){
			if(edge_list[i].v<0)
						continue;
			Edge & e = edge_list[i];
			lbool val = S->value(e.v);
			assert(val!=l_Undef);

			if(val==l_True){
				assert(g.hasEdge(e.from,e.to));
				assert(antig.hasEdge(e.from,e.to));
			}else{
				assert(!g.hasEdge(e.from,e.to));
				assert(!antig.hasEdge(e.from,e.to));
			}


		}
		for(int i = 0;i<reach_detectors.size();i++){
			ReachDetector* d  = reach_detectors[i];
			for(int j = 0;j< d->reach_lits.size();j++){
				Lit l = d->reach_lits[j];
				if(l!=lit_Undef){
					int node = d->getNode(var(l));

					if(S->value(l)==l_True){
						assert(d->positive_reach_detector->connected(node));
					}else if (S->value(l)==l_False){
						assert(! d->negative_reach_detector->connected(node));
					}else{
						assert(!d->positive_reach_detector->connected(node));
						assert(d->negative_reach_detector->connected(node));
					}
				}
			}
		}
#endif
		return true;
	}

	void drawCurrent(){

	}

	Lit newEdge(int from,int to, Var v = var_Undef)
    {

		if(v==var_Undef)
			v = S->newVar();
#ifdef DEBUG_GRAPH
		 dbg_graph->newEdge(from,to,v);
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg->newEdge(from,to,v);
#endif
		if(num_edges>0){
		}else
			min_edge_var=v;

		int index = v-min_edge_var;

		while(edge_list.size()<=index){
			edge_list.push({-1,-1,-1});
			edge_assignments.push(l_Undef);
		}

		inv_adj[to].push({v,from,to});

		num_edges++;
		edge_list[index].v =v;
		edge_list[index].from=from;
		edge_list[index].to =to;

		edges[from][to]= {v,from,to};
		g.addEdge(from,to,index);
		g.disableEdge(from,to, index);
		antig.addEdge(from,to,index);
		cutGraph.addEdge(from,to,index);
    	return mkLit(v,false);
    }

	void reaches(int from, int to, Var reach_var,int within_steps=-1){

#ifdef DEBUG_GRAPH
		 dbg_graph->reaches(from,  to,reach_var,within_steps);
#endif
#ifdef DEBUG_SOLVER
		if(S->dbg_solver)
			shadow_dbg->reaches(from,  to,reach_var,within_steps);
#endif
			assert(from<g.nodes);
			if(within_steps>g.nodes)
				within_steps=-1;

			if (reach_info[from].source<0){
				reach_detectors.push(new ReachDetector(*this,from));
				reach_detectors.last()->reach_marker=S->newReasonMarker(this);

				int mnum = CRef_Undef- reach_detectors.last()->reach_marker;
				marker_map.growTo(mnum+1);
				marker_map[mnum] = reach_detectors.size()-1;
				//marker_map.insert(reach_markers.last(),reach_markers.size());

				reach_detectors.last()->non_reach_marker=S->newReasonMarker(this);
				//marker_map[non_reach_markers.last()]=-non_reach_markers.size();
				//marker_map.insert(non_reach_markers.last(),non_reach_markers.size());

				mnum = CRef_Undef- reach_detectors.last()->non_reach_marker;
				marker_map.growTo(mnum+1);
				marker_map[mnum] = reach_detectors.size()-1;

				reach_detectors.last()->forced_reach_marker =S->newReasonMarker(this);
				//marker_map[non_reach_markers.last()]=-non_reach_markers.size();
				//marker_map.insert(non_reach_markers.last(),non_reach_markers.size());

				mnum = CRef_Undef- reach_detectors.last()->forced_reach_marker;
				marker_map.growTo(mnum+1);
				marker_map[mnum] = reach_detectors.size()-1;


				if(within_steps<0 ){
					if(reachalg==ALG_CONNECTIVITY){
						reach_detectors.last()->positiveReachStatus = new ReachDetector::ReachStatus(*reach_detectors.last(),true);
						reach_detectors.last()->negativeReachStatus = new ReachDetector::ReachStatus(*reach_detectors.last(),false);
						reach_detectors.last()->positive_reach_detector = new Connectivity<ReachDetector::ReachStatus,PositiveEdgeStatus>(from,g,*(reach_detectors.last()->positiveReachStatus),1);
						reach_detectors.last()->negative_reach_detector = new Connectivity<ReachDetector::ReachStatus,NegativeEdgeStatus>(from,antig,*(reach_detectors.last()->negativeReachStatus),-1);
						if(opt_conflict_shortest_path)
							reach_detectors.last()->positive_path_detector = new Distance<NullEdgeStatus,PositiveEdgeStatus>(from,g,nullEdgeStatus,1);
						else
							reach_detectors.last()->positive_path_detector = reach_detectors.last()->positive_reach_detector;
					}else if(reachalg==ALG_BFS){
						reach_detectors.last()->positiveReachStatus = new ReachDetector::ReachStatus(*reach_detectors.last(),true);
						reach_detectors.last()->negativeReachStatus = new ReachDetector::ReachStatus(*reach_detectors.last(),false);
						reach_detectors.last()->positive_reach_detector = new Distance<ReachDetector::ReachStatus,PositiveEdgeStatus>(from,g,*(reach_detectors.last()->positiveReachStatus),1);
						reach_detectors.last()->negative_reach_detector = new Distance<ReachDetector::ReachStatus,NegativeEdgeStatus>(from,antig,*(reach_detectors.last()->negativeReachStatus),-1);
						reach_detectors.last()->positive_path_detector = reach_detectors.last()->positive_reach_detector;
					}else{

						reach_detectors.last()->positive_reach_detector = new Dijkstra<PositiveEdgeStatus>(from,g);
						reach_detectors.last()->negative_reach_detector = new Dijkstra<NegativeEdgeStatus>(from,antig);
						reach_detectors.last()->positive_path_detector = reach_detectors.last()->positive_reach_detector;
						//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
					}
				}else{

					if(distalg==ALG_BFS){
						reach_detectors.last()->positiveReachStatus = new ReachDetector::ReachStatus(*reach_detectors.last(),true);
						reach_detectors.last()->negativeReachStatus = new ReachDetector::ReachStatus(*reach_detectors.last(),false);
						reach_detectors.last()->positive_reach_detector = new Distance<ReachDetector::ReachStatus,PositiveEdgeStatus>(from,g,*(reach_detectors.last()->positiveReachStatus),1);
						reach_detectors.last()->negative_reach_detector = new Distance<ReachDetector::ReachStatus,NegativeEdgeStatus>(from,antig,*(reach_detectors.last()->negativeReachStatus),-1);
						reach_detectors.last()->positive_path_detector = reach_detectors.last()->positive_reach_detector;
						/*	if(opt_conflict_shortest_path)
							reach_detectors.last()->positive_dist_detector = new Dijkstra<PositiveEdgeStatus>(from,g);*/
					}else{

						reach_detectors.last()->positive_reach_detector = new Dijkstra<PositiveEdgeStatus>(from,g);
						reach_detectors.last()->negative_reach_detector = new Dijkstra<NegativeEdgeStatus>(from,antig);
						reach_detectors.last()->positive_path_detector = reach_detectors.last()->positive_reach_detector;
						//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
					}


				}
				//reach_detectors.last()->negative_dist_detector = new Dijkstra(from,antig);
				reach_detectors.last()->source=from;

				reach_info[from].source=from;
				reach_info[from].detector=reach_detectors.last();

				reach_detectors.last()->within=within_steps;



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
