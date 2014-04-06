/*
 * ReachDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "ReachDetector.h"
#include "GraphTheory.h"
#include "core/Config.h"
#include "DynamicConnectivity.h"
ReachDetector::ReachDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph<PositiveEdgeStatus> &_g, DynamicGraph<NegativeEdgeStatus> &_antig, int from,double seed):Detector(_detectorID),outer(_outer),g(_g),antig(_antig),within(-1),source(from),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),opt_weight(*this),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){
	check_positive=true;
	check_negative=true;
	rnd_path=NULL;
	opt_path=NULL;
	 if(opt_use_random_path_for_decisions){
		 rnd_weight.clear();
		 rnd_path = new WeightedDijkstra<NegativeEdgeStatus, vec<double> >(from,_antig,rnd_weight);
		 for(int i=0;i<outer->edge_list.size();i++){
			 double w = drand(rnd_seed);

			 rnd_weight.push(w);
		 }

	 }

	 if(opt_use_optimal_path_for_decisions){
		 opt_path = new WeightedDijkstra<NegativeEdgeStatus, OptimalWeightEdgeStatus >(from,_antig,opt_weight);
	 }

	if(reachalg==ReachAlg::ALG_BFS){
							positiveReachStatus = new ReachDetector::ReachStatus(*this,true);
							negativeReachStatus = new ReachDetector::ReachStatus(*this,false);
							positive_reach_detector = new Connectivity<ReachDetector::ReachStatus,PositiveEdgeStatus>(from,_g,*(positiveReachStatus),1);
							negative_reach_detector = new Connectivity<ReachDetector::ReachStatus,NegativeEdgeStatus>(from,_antig,*(negativeReachStatus),-1);
							if(opt_conflict_shortest_path)
								positive_path_detector = new Distance<NullEdgeStatus,PositiveEdgeStatus>(from,_g,nullEdgeStatus,1);
							else
								positive_path_detector =positive_reach_detector;
						}else if(reachalg==ReachAlg::ALG_DFS){
							positiveReachStatus = new ReachDetector::ReachStatus(*this,true);
							negativeReachStatus = new ReachDetector::ReachStatus(*this,false);
							positive_reach_detector = new DFSReachability<ReachDetector::ReachStatus,PositiveEdgeStatus>(from,_g,*(positiveReachStatus),1);
							negative_reach_detector = new DFSReachability<ReachDetector::ReachStatus,NegativeEdgeStatus>(from,_antig,*(negativeReachStatus),-1);
							if(opt_conflict_shortest_path)
								positive_path_detector = new Distance<NullEdgeStatus,PositiveEdgeStatus>(from,_g,nullEdgeStatus,1);
							else
								positive_path_detector =positive_reach_detector;
						}else if(reachalg==ReachAlg::ALG_DISTANCE){
							positiveReachStatus = new ReachDetector::ReachStatus(*this,true);
							negativeReachStatus = new ReachDetector::ReachStatus(*this,false);
							positive_reach_detector = new Distance<ReachDetector::ReachStatus,PositiveEdgeStatus>(from,_g,*(positiveReachStatus),1);
							negative_reach_detector = new Distance<ReachDetector::ReachStatus,NegativeEdgeStatus>(from,_antig,*(negativeReachStatus),-1);
							positive_path_detector = positive_reach_detector;
						}/*else if (reachalg==ReachAlg::ALG_THORUP){
							positiveReachStatus = new ReachDetector::ReachStatus(*this,true);
							negativeReachStatus = new ReachDetector::ReachStatus(*this,false);
							positive_reach_detector = new DynamicConnectivity<ReachDetector::ReachStatus,PositiveEdgeStatus>(_g,*(positiveReachStatus),1);
							negative_reach_detector = new DynamicConnectivity<ReachDetector::ReachStatus,NegativeEdgeStatus>(_antig,*(negativeReachStatus),-1);
							positive_path_detector = positive_reach_detector;
							if(opt_conflict_shortest_path)
								positive_path_detector = new Distance<NullEdgeStatus,PositiveEdgeStatus>(from,_g,nullEdgeStatus,1);
							else
								positive_path_detector =positive_reach_detector;
						}*/else{
							positive_reach_detector = new Dijkstra<PositiveEdgeStatus>(from,_g);
							negative_reach_detector = new Dijkstra<NegativeEdgeStatus>(from,_antig);
							positive_path_detector = positive_reach_detector;
							//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
						}
	positive_reach_detector->setSource(source);
	negative_reach_detector->setSource(source);
	first_reach_var = var_Undef;
	reach_marker=outer->newReasonMarker(getID());
	non_reach_marker=outer->newReasonMarker(getID());
	forced_reach_marker=outer->newReasonMarker(getID());
}

void ReachDetector::addLit(int from, int to, Var outer_reach_var){
	g.invalidate();
	antig.invalidate();
	Var reach_var = outer->newVar(outer_reach_var,getID());
	if(first_reach_var==var_Undef){
		first_reach_var=reach_var;
	}else{
		assert(reach_var>first_reach_var);
	}
	assert(from==source);
	while( reach_lits.size()<=to)
			reach_lits.push(lit_Undef);
/*
	while(outer->S->nVars()<=reach_var)
		outer->S->newVar();*/

	Lit reachLit=mkLit(reach_var,false);

	if(reach_lits[to]==lit_Undef){
		reach_lits[to] = reachLit;

		while(reach_lit_map.size()<= reach_var- first_reach_var ){
			reach_lit_map.push(-1);
		}

		reach_lit_map[reach_var-first_reach_var]=to;
	}else{
		Lit r = reach_lits[to];
		//force equality between the new lit and the old reach lit, in the SAT solver
		outer->makeEqual(r,reachLit);
	/*	outer->S->addClause(~r, reachLit);
		outer->S->addClause(r, ~reachLit);*/
	}
}

void ReachDetector::ReachStatus::setReachable(int u, bool reachable){
				if(reachable){
					assert(!detector.outer->dbg_notreachable( detector.source,u));
				}else{
					assert(detector.outer->dbg_notreachable( detector.source,u));
				}
				if(polarity==reachable && u<detector.reach_lits.size()){
					Lit l = detector.reach_lits[u];
					if(l!=lit_Undef){
						lbool assign = detector.outer->value(l);
						if(assign!= (reachable? l_True:l_False )){
							detector.changed.push({reachable? l:~l,u});
						}
					}
				}
			}

void ReachDetector::ReachStatus::setMininumDistance(int u, bool reachable, int distance){
	assert(reachable ==(distance<detector.outer->g.nodes));
	setReachable(u,reachable);

	if(u<detector.dist_lits.size()){
		assert(distance>=0);
		assert(distance<detector.outer->g.nodes);
		for(int i = 0;i<detector.dist_lits[u].size();i++){
			int d =  detector.dist_lits[u][i].min_distance;
			Lit l = detector.dist_lits[u][i].l;
			assert(l!=lit_Undef);
			if(d<distance && !polarity){
				lbool assign = detector.outer->value(l);
				if( assign!= l_False ){
					detector.changed.push({~l,u});
				}
			}else if(d>=distance && polarity){
				lbool assign = detector.outer->value(l);
				if( assign!= l_True ){
					detector.changed.push({l,u});
				}
			}
		}
	}
}

bool ReachDetector::ChokepointStatus::mustReach(int node){
	Lit l =  detector.reach_lits[node];
	if(l!=lit_Undef){
		return detector.outer->value(l)==l_True;
	}
	return false;
}
bool ReachDetector::ChokepointStatus::operator() (int edge_id){
	return detector.outer->edge_assignments[edge_id]==l_Undef;
}

void ReachDetector::preprocess(){
	//vec<bool> pure;
	//pure.growTo(reach_lits.size());
	//can check if all reach lits appear in only one polarity in the solver constraints; if so, then we can disable either check_positive or check_negative
}

void ReachDetector::buildReachReason(int node,vec<Lit> & conflict){
			//drawFull();
			Reach & d = *positive_path_detector;


			double starttime = cpuTime();
			d.update();

			assert(outer->dbg_reachable(d.getSource(),node));
		/*	if(!outer->dbg_reachable(d.getSource(),node)){
				outer->drawFull();

				//d.drawFull();

				assert(false);
			}*/
			assert(d.connected_unchecked(node));
			if(opt_learn_reaches ==0 || opt_learn_reaches==2){
				int u = node;
				int p;
				while(( p = d.previous(u)) != -1){
					Edge edg = outer->edges[p][u];
					Var e =outer->edges[p][u].v;
					lbool val = outer->value(e);
					assert(outer->value(e)==l_True);
					conflict.push(mkLit(e, true));
					u = p;

				}
			}else{
				//Instead of a complete path, we can learn reach variables, if they exist
				int u = node;
				int p;
				while(( p = d.previous(u)) != -1){
					Edge edg = outer->edges[p][u];
					Var e =outer->edges[p][u].v;
					lbool val = outer->value(e);
					assert(outer->value(e)==l_True);
					conflict.push(mkLit(e, true));
					u = p;
					if( u< reach_lits.size() && reach_lits[u]!=lit_Undef && outer->value(reach_lits[u])==l_True && outer->level(var(reach_lits[u]))< outer->decisionLevel()){
						//A potential (fixed) problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
						//This is avoided by ensuring that L is lower level than the conflict.
						Lit l =reach_lits[u];
						assert(outer->value(l)==l_True);
						conflict.push(~l);
						break;
					}
				}

			}
			outer->num_learnt_paths++;
			outer->learnt_path_clause_length+= (conflict.size()-1);
			double elapsed = cpuTime()-starttime;
			outer->pathtime+=elapsed;
	#ifdef DEBUG_GRAPH
			 assert(outer->dbg_clause(conflict));

	#endif
		}
		void ReachDetector::buildNonReachReason(int node,vec<Lit> & conflict){
			static int it = 0;
			++it;
			int u = node;
			//drawFull( non_reach_detectors[detector]->getSource(),u);
			assert(outer->dbg_notreachable( source,u));
			double starttime = cpuTime();
			outer->cutGraph.clearHistory();
			outer->stats_mc_calls++;
			if(opt_conflict_min_cut){
				if(mincutalg!= MinCutAlg::ALG_EDKARP_ADJ){
					//ok, set the weights for each edge in the cut graph.
					//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
					for(int u = 0;u<outer->cutGraph.adjacency.size();u++){
						for(int j = 0;j<outer->cutGraph.adjacency[u].size();j++){
							int v = outer->cutGraph.adjacency[u][j].node;
							Var var = outer->edges[u][v].v;
							/*if(S->value(var)==l_False){
								mc.setCapacity(u,v,1);
							}else{*/
							outer->mc->setCapacity(u,v,0xF0F0F0);
							//}
						}
					}

					//find any edges assigned to false, and set their capacity to 1
					for(int i =0;i<outer->trail.size();i++){
						if(outer->trail[i].isEdge && !outer->trail[i].assign){
							outer->mc->setCapacity(outer->trail[i].from, outer->trail[i].to,1);
						}
					}
				}
				outer->cut.clear();

				int f =outer->mc->minCut(source,node,outer->cut);
				assert(f<0xF0F0F0); assert(f==outer->cut.size());//because edges are only ever infinity or 1
				for(int i = 0;i<outer->cut.size();i++){
					MaxFlow::Edge e = outer->cut[i];

					Lit l = mkLit( outer->edges[e.u][e.v].v,false);
					assert(outer->value(l)==l_False);
					conflict.push(l);
				}
			}else{
				//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
				//or perhaps we can learn the actual 1-uip cut?


					vec<int>& to_visit  = outer->to_visit;
					vec<char>& seen  = outer->seen;

				    to_visit.clear();
				    to_visit.push(node);
				    seen.clear();
					seen.growTo(outer->nNodes());
				    seen[node]=true;

				    do{

				    	assert(to_visit.size());
				    	int u = to_visit.last();
				    	assert(u!=source);
				    	to_visit.pop();
				    	assert(seen[u]);
				    	assert(!negative_reach_detector->connected_unsafe(u));
				    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
				    	for(int i = 0;i<outer->inv_adj[u].size();i++){
				    		int v = outer->inv_adj[u][i].v;
				    		int from = outer->inv_adj[u][i].from;

				    		assert(from!=u);
				    		assert(outer->inv_adj[u][i].to==u);
				    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				    		int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;

				    		if(outer->edge_assignments[edge_num]==l_False){
				    			//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
				    			//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
				    			//if(!seen[from])
				    				conflict.push(mkLit(v,false));
				    		}else{
				    			assert(from!=source);
				    			//even if it is undef? probably...
				    			if(!seen[from]){
				    				seen[from]=true;
				    				if((opt_learn_reaches ==2 || opt_learn_reaches==3) && from< reach_lits.size() && reach_lits[from]!=lit_Undef && outer->value(reach_lits[from])==l_False  && outer->level(var(reach_lits[from]))< outer->decisionLevel())
				    				{
				    				//The problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
				    					Lit r = reach_lits[from];
				    					assert(var(r)<outer->nVars());
				    					assert(outer->value(r)==l_False);
				    					conflict.push(r);
				    				}else
				    					to_visit.push(from);
				    			}
				    		}
				    	}
				    }  while (to_visit.size());


			}

			 outer->num_learnt_cuts++;
			 outer->learnt_cut_clause_length+= (conflict.size()-1);

			double elapsed = cpuTime()-starttime;
			 outer->mctime+=elapsed;

	#ifdef DEBUG_GRAPH
			 assert(outer->dbg_clause(conflict));
	#endif
		}

		/**
		 * Explain why an edge was forced (to true).
		 * The reason is that _IF_ that edge is false, THEN there is a cut of disabled edges between source and target
		 * So, create the graph that has that edge (temporarily) assigned false, and find a min-cut in it...
		 */
		void ReachDetector::buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict){
					static int it = 0;
					++it;

					assert(outer->edge_assignments[forced_edge_id]==l_True);
					Lit edgeLit =mkLit( outer->edge_list[forced_edge_id].v,false);

					conflict.push(edgeLit);

					int forced_edge_from = outer->edge_list[forced_edge_id].from;
					int forced_edge_to= outer->edge_list[forced_edge_id].to;


					int u = reach_node;
					//drawFull( non_reach_detectors[detector]->getSource(),u);
					assert(outer->dbg_notreachable( source,u));
					double starttime = cpuTime();
					outer->cutGraph.clearHistory();
					outer->stats_mc_calls++;



					if(opt_conflict_min_cut){
						if(mincutalg!= MinCutAlg::ALG_EDKARP_ADJ){
							//ok, set the weights for each edge in the cut graph.
							//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
							for(int u = 0;u<outer->cutGraph.adjacency.size();u++){
								for(int j = 0;j<outer->cutGraph.adjacency[u].size();j++){
									int v = outer->cutGraph.adjacency[u][j].node;
									Var var = outer->edges[u][v].v;
									/*if(S->value(var)==l_False){
										mc.setCapacity(u,v,1);
									}else{*/
									outer->mc->setCapacity(u,v,0xF0F0F0);
									//}
								}
							}

							//find any edges assigned to false, and set their capacity to 1
							for(int i =0;i<outer->trail.size();i++){
								if(outer->trail[i].isEdge && !outer->trail[i].assign){
									outer->mc->setCapacity(outer->trail[i].from, outer->trail[i].to,1);
								}
							}

							outer->mc->setCapacity(forced_edge_from, forced_edge_to,1);
						}
						outer->cut.clear();

						int f =outer->mc->minCut(source,reach_node,outer->cut);
						assert(f<0xF0F0F0); assert(f==outer->cut.size());//because edges are only ever infinity or 1
						for(int i = 0;i<outer->cut.size();i++){
							MaxFlow::Edge e = outer->cut[i];

							Lit l = mkLit( outer->edges[e.u][e.v].v,false);
							assert(outer->value(l)==l_False);
							conflict.push(l);
						}
					}else{
						//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
						//or perhaps we can learn the actual 1-uip cut?


							vec<int>& to_visit  = outer->to_visit;
							vec<char>& seen  = outer->seen;

						    to_visit.clear();
						    to_visit.push(reach_node);
						    seen.clear();
							seen.growTo(outer->nNodes());
						    seen[reach_node]=true;

						    do{

						    	assert(to_visit.size());
						    	int u = to_visit.last();
						    	assert(u!=source);
						    	to_visit.pop();
						    	assert(seen[u]);
						    	assert(!negative_reach_detector->connected_unsafe(u));
						    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
						    	for(int i = 0;i<outer->inv_adj[u].size();i++){
						    		int v = outer->inv_adj[u][i].v;
						    		int from = outer->inv_adj[u][i].from;
						    		assert(from!=u);
						    		assert(outer->inv_adj[u][i].to==u);
						    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
						    		int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;

						    		if(edge_num == forced_edge_id || outer->edge_assignments[edge_num]==l_False){
						    			//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
						    			//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
						    			//if(!seen[from])
						    				conflict.push(mkLit(v,false));
						    		}else{
						    			assert(from!=source);
						    			//even if it is undef? probably...
						    			if(!seen[from]){
						    				seen[from]=true;
						    				if((opt_learn_reaches ==2 || opt_learn_reaches==3) && from< reach_lits.size() && reach_lits[from]!=lit_Undef && outer->value(reach_lits[from])==l_False  && outer->level(var(reach_lits[from]))< outer->decisionLevel())
						    				{
						    				//The problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
						    					Lit r = reach_lits[from];
						    					assert(var(r)<outer->nVars());
						    					assert(outer->value(r)==l_False);
						    					conflict.push(r);
						    				}else
						    					to_visit.push(from);
						    			}
						    		}
						    	}
						    }  while (to_visit.size());


					}

					 outer->num_learnt_cuts++;
					 outer->learnt_cut_clause_length+= (conflict.size()-1);

					double elapsed = cpuTime()-starttime;
					 outer->mctime+=elapsed;

			#ifdef DEBUG_GRAPH
					 assert(outer->dbg_clause(conflict));
			#endif
				}

		void ReachDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){


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
				 assert(outer->dbg_clause(reason));

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
					int forced_edge_id =outer->getEdgeID(v);//v- outer->min_edge_var;
					//The corresponding node that is the reason it was forced
					int reach_node=force_reason[forced_edge_id];
					 buildForcedEdgeReason( reach_node,  forced_edge_id, reason);
				}else{
					assert(false);
				}
		}

		bool ReachDetector::propagate(vec<Assignment> & trail,vec<Lit> & conflict){

			if(check_positive){
				double startdreachtime = cpuTime();
				getChanged().clear();
				positive_reach_detector->update();
				double reachUpdateElapsed = cpuTime()-startdreachtime;
				outer->reachupdatetime+=reachUpdateElapsed;
			}
			if(check_negative){
				double startunreachtime = cpuTime();
				negative_reach_detector->update();
				double unreachUpdateElapsed = cpuTime()-startunreachtime;
				outer->unreachupdatetime+=unreachUpdateElapsed;
			}
			for(int j = 0;j<getChanged().size();j++){
					Lit l = getChanged()[j].l;
					int u =  getChanged()[j].u;
					bool reach = !sign(l);
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
						assert(outer->dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
						if(S->dbg_solver)
							S->dbg_check_propagation(l);
#endif
						trail.push(Assignment(false,reach,detectorID,0,var(l)));
						if(reach)
							outer->enqueue(l,reach_marker) ;
						else
							outer->enqueue(l,non_reach_marker) ;

					}else if (outer->value(l)==l_False){
						conflict.push(l);

						if(reach){

						//conflict
						//The reason is a path in g from to s in d
							buildReachReason(u,conflict);
						//add it to s
						//return it as a conflict

						}else{
							//The reason is a cut separating s from t
							buildNonReachReason(u,conflict);

						}
#ifdef DEBUG_GRAPH
						for(int i = 0;i<conflict.size();i++)
									 assert(outer->value(conflict[i])==l_False);
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
						chokepoint.collectForcedEdges(forced_edges);
						for(int i = 0;i<forced_edges.size();i++){
							int edge_id = forced_edges[i].edge_id;
							int node = forced_edges[i].node;
							Lit l = mkLit( outer->edge_list[edge_id].v,false);
							if(outer->value(l)==l_Undef){
								force_reason.growTo(edge_id+1);
								force_reason[edge_id]=node;
								outer->enqueue(l,forced_reach_marker);
							}else if(outer->value(l)==l_True){
								//do nothing

							}else{
								//conflict.
								//this actually shouldn't be possible (at this point in the code)
								buildForcedEdgeReason(node,edge_id,conflict);
								return false;
							}
						}

					}

				}

			#ifdef DEBUG_GRAPH
					for(int i = 0;i<reach_lits.size();i++){
						Lit l = reach_lits[i];
						if(l!=lit_Undef){
							int u = getNode(var(l));
							if(positive_reach_detector->connected(u)){
								assert(outer->value(l)==l_True);
							}else if (!negative_reach_detector->connected(u)){
								assert(outer->value(l)==l_False);
							}
						}

					}
			#endif
			return true;
		}

bool ReachDetector::checkSatisfied(){

				for(int j = 0;j< reach_lits.size();j++){
					Lit l = reach_lits[j];
					if(l!=lit_Undef){
						int node =getNode(var(l));

						if(outer->value(l)==l_True){
							if(!positive_reach_detector->connected(node)){
								return false;
							}
						}else if (outer->value(l)==l_False){
							if( negative_reach_detector->connected(node)){
								return false;
							}
						}else{
							if(positive_reach_detector->connected(node)){
								return false;
							}
							if(!negative_reach_detector->connected(node)){
								return false;
							}
						}
					}
				}
	return true;
}

void ReachDetector::dbg_sync_reachability(){
		for(int j = 0;j< reach_lits.size();j++){
						Lit l =reach_lits[j];
						if(l!=lit_Undef){
							int node = getNode(var(l));

							if(positive_reach_detector->connected(node)){
								assert(outer->value(l)==l_True);
							}else if(! negative_reach_detector->connected(node)){
								assert(outer->value(l)==l_False);
							}
						}

					}
	}


int ReachDetector::OptimalWeightEdgeStatus::operator [] (int edge) const {
	Var v = detector.outer->edge_list[edge].v;
	lbool val = detector.outer->value(v);
	if(val==l_False){
		assert(false);
		return detector.outer->edge_list.size()*2;
	}else if (val==l_True){
		return 0;
	}else if (val==l_Undef){
		return 1;
	}
}
int ReachDetector::OptimalWeightEdgeStatus::size()const{
	return detector.outer->edge_list.size();
}


Lit ReachDetector::decide(){

	Distance<ReachDetector::ReachStatus,NegativeEdgeStatus> * over = (Distance<ReachDetector::ReachStatus,NegativeEdgeStatus>*)negative_reach_detector;

	Distance<ReachDetector::ReachStatus,PositiveEdgeStatus> * under = (Distance<ReachDetector::ReachStatus,PositiveEdgeStatus>*)positive_reach_detector;

	//we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

	//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

	//this can be obviously more efficient
	//for(int j = 0;j<nNodes();j++){
	for(int k = 0;k<reach_lits.size();k++){
		Lit l =reach_lits[k];
		if(l==lit_Undef)
			continue;
		int j =getNode(var(l));
		if(outer->value(l)==l_True && opt_decide_graph_pos){
			//if(S->level(var(l))>0)
			//	continue;
			assert(over->connected(j));
			if(over->connected(j) && !under->connected(j)){
				//then lets try to connect this
				static vec<bool> print_path;

				assert(over->connected(j));//Else, we would already be in conflict
				int p =j;
				int last=j;
				if(!opt_use_random_path_for_decisions){
					if(opt_use_optimal_path_for_decisions){
						//ok, read back the path from the over to find a candidate edge we can decide
						//find the earliest unconnected node on this path
						opt_path->update();
						 p = j;
						 last = j;
						while(!under->connected(p)){

							last=p;
							assert(p!=source);
							int prev = opt_path->previous(p);
							p = prev;

						}
					}else{
						//ok, read back the path from the over to find a candidate edge we can decide
						//find the earliest unconnected node on this path
						over->update();
						 p = j;
						 last = j;
						while(!under->connected(p)){

							last=p;
							assert(p!=source);
							int prev = over->previous(p);
							p = prev;

						}
					}
				}else{
					//Randomly re-weight the graph sometimes
					if(drand(rnd_seed)<opt_decide_graph_re_rnd){

						for(int i=0;i<outer->edge_list.size() ;i++){
								 double w = drand(rnd_seed);
								/* w-=0.5;
								 w*=w;*/
								 //printf("%f (%f),",w,rnd_seed);
								//rnd_path->setWeight(i,w);
								 rnd_weight[i]=w;
							 }
					}
					rnd_path->update();
					//derive a random path in the graph
					 p = j;
					 last = j;
					 assert(rnd_path->connected(p));
					while(!under->connected(p)){

						last=p;
						assert(p!=source);
						int prev =rnd_path->previous(p);
						p = prev;
						assert(p>=0);
					}

				}

				if(outer->level(var(l))==0 && opt_print_decision_path){
					print_path.clear();
					print_path.growTo(outer->nNodes());

					if(!opt_use_random_path_for_decisions){
											//ok, read back the path from the over to find a candidate edge we can decide
											//find the earliest unconnected node on this path
						over->update();
						int p = j;

						while(p!=source){
							if(opt_print_decision_path)
								print_path[p]=true;


							assert(p!=source);
							int prev = over->previous(p);
							p = prev;

						}
					}else{


						int p = j;

						while(p!=source){
							if(opt_print_decision_path)
								print_path[p]=true;


							assert(p!=source);
							int prev =rnd_path->previous(p);
							p = prev;

						}
					}

					if(opt_print_decision_path){
						printf("From %d to %d:\n",source,j);
						int width = sqrt(outer->nNodes());


						int v = 0;
						//for (int i = 0;i<w;i++){
						//	for(int j = 0;j<w;j++){
						int lasty= 0;
						for(int n = 0;n<outer->nNodes();n++){
							int x = n%width;
							int y = n/width;
							if(y > lasty)
								printf("\n");

							if(n==j){
								printf("* ");
							}else if (n==source){
								printf("#");
							}else{

								if(print_path[n]){
									printf("+ ");
								}else{
									printf("- ");
								}
							}

							lasty=y;
						}
						printf("\n\n");
					}
				}

				//ok, now pick some edge p->last that will connect p to last;
			/*	assert(!under->connected(last));
				assert(under->connected(p));

				assert(over->connected(last));
				assert(over->connected(p));*/
				Var v = outer->edges[p][last].v;
				if(outer->value(v)==l_Undef){
					return mkLit(v,false);
				}else{
					assert(outer->value(v)!=l_True);
				}
			/*	for(int k = 0;k<outer->antig.adjacency[p].size();k++){
					int to = outer->antig.adjacency[p][k].node;
					if (to==last){
						Var v =outer->edge_list[ outer->antig.adjacency[p][k].id].v;
						if(outer->value(v)==l_Undef){
							return mkLit(v,false);
						}else{
							assert(outer->value(v)!=l_True);
						}
					}
				}*/

			}
		}else if(outer->value(l)==l_False && opt_decide_graph_neg){


				//for each negated reachability constraint, we can find a cut through the unassigned edges in the over-approx and disable one of those edges.
				assert(!under->connected(j));
				over->update();
				if(over->connected(j) && ! under->connected(j)){
					//then lets try to disconnect this node from source by walking back along the path in the over approx, and disabling the first unassigned edge we see.
					//(there must be at least one such edge, else the variable would be connected in the under approximation as well - in which case it would already have been propagated.
					int p = j;
					int last = j;
					while(!under->connected(p)){
						last=p;
						assert(p!=source);
						int prev = over->previous(p);
						Var v = outer->edges[prev][p].v;
						if(outer->value(v)==l_Undef){
							return mkLit(v,true);
						}else{
							assert(outer->value(v)!=l_False);
						}
						p = prev;

					}


				}

		}

	}
	return lit_Undef;
};


