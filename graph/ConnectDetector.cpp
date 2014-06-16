/*
 * ConnectDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "ConnectDetector.h"
#include "GraphTheory.h"
#include "core/Config.h"
#include "dgl/DynamicConnectivity.h"
ConnectDetector::ConnectDetector(int _detectorID, GraphTheorySolver * _outer, DynamicGraph &_g, DynamicGraph &_antig, int from,double seed):Detector(_detectorID),outer(_outer),g(_g),antig(_antig),within(-1),source(from),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL),opt_weight(*this),chokepoint_status(*this),chokepoint(chokepoint_status, _antig,source){
	check_positive=true;
	check_negative=true;
	constraintsBuilt=-1;
	first_reach_var = var_Undef;

	rnd_path=nullptr;
	opt_path=nullptr;


	 if(undirectedalg ==ConnectivityAlg::ALG_SAT){
		 positiveReachStatus=nullptr;
		 negativeReachStatus=nullptr;
		 positive_reach_detector=nullptr;
		 negative_reach_detector=nullptr;
		 positive_path_detector=nullptr;
		 reach_marker=CRef_Undef;
		 non_reach_marker=CRef_Undef;
		 forced_reach_marker=CRef_Undef;

		 //we are just going to directly enforce these constraints in the original SAT solver, so _nothing_ will end up happening in this detector (except for creating the clauses needed to enforce these constraints).

		 return;
	 }


	 if(opt_use_random_path_for_decisions){
		 rnd_weight.clear();
		 rnd_path = new WeightedDijkstra< vec<double> >(from,_antig,rnd_weight);
		 for(int i=0;i<outer->edge_list.size();i++){
			 double w = drand(rnd_seed);

			 rnd_weight.push(w);
		 }

	 }

	 if(opt_use_optimal_path_for_decisions){
		 opt_path = new WeightedDijkstra< OptimalWeightEdgeStatus >(from,_antig,opt_weight);
	 }
		positiveReachStatus = new ConnectDetector::ReachStatus(*this,true);
		negativeReachStatus = new ConnectDetector::ReachStatus(*this,false);
	 if(undirectedalg ==ConnectivityAlg::ALG_BFS){

		positive_reach_detector = new BFSReachability<ConnectDetector::ReachStatus,true>(from,_g,*(positiveReachStatus),1);
		negative_reach_detector = new BFSReachability<ConnectDetector::ReachStatus,true>(from,_antig,*(negativeReachStatus),-1);
		if(opt_conflict_shortest_path)
			positive_path_detector = new Distance<Reach::NullStatus,true>(from,_g,Reach::nullStatus,1);
		else
			positive_path_detector =positive_reach_detector;
	}else if(undirectedalg==ConnectivityAlg::ALG_DFS){

		positive_reach_detector = new DFSReachability<ConnectDetector::ReachStatus,true>(from,_g,*(positiveReachStatus),1);
		negative_reach_detector = new DFSReachability<ConnectDetector::ReachStatus,true>(from,_antig,*(negativeReachStatus),-1);
		if(opt_conflict_shortest_path)
			positive_path_detector = new Distance<Reach::NullStatus,true>(from,_g,Reach::nullStatus,1);
		else
			positive_path_detector =positive_reach_detector;
	}else if(undirectedalg==ConnectivityAlg::ALG_DISTANCE){

		positive_reach_detector = new Distance<ConnectDetector::ReachStatus,true>(from,_g,*(positiveReachStatus),1);
		negative_reach_detector = new Distance<ConnectDetector::ReachStatus,true>(from,_antig,*(negativeReachStatus),-1);
		positive_path_detector = positive_reach_detector;
	}else if (undirectedalg==ConnectivityAlg::ALG_THORUP){

		positive_reach_detector = new DynamicConnectivity<ConnectDetector::ReachStatus>(_g,*(positiveReachStatus),1);
		negative_reach_detector = new DynamicConnectivity<ConnectDetector::ReachStatus>(_antig,*(negativeReachStatus),-1);
		positive_path_detector = positive_reach_detector;
		if(opt_conflict_shortest_path)
			positive_path_detector = new Distance<Reach::NullStatus,true>(from,_g,Reach::nullStatus,1);
		else
			positive_path_detector =positive_reach_detector;
	}else{
		positive_reach_detector = new Dijkstra<ConnectDetector::ReachStatus, true>(from,_g,*positiveReachStatus);
		negative_reach_detector = new Dijkstra<ConnectDetector::ReachStatus,true>(from,_antig,*negativeReachStatus);
		positive_path_detector = positive_reach_detector;
		//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
	}
	positive_reach_detector->setSource(source);
	negative_reach_detector->setSource(source);

	reach_marker=outer->newReasonMarker(getID());
	non_reach_marker=outer->newReasonMarker(getID());
	forced_reach_marker=outer->newReasonMarker(getID());
}

void ConnectDetector::buildSATConstraints(int within_steps){
	if(within_steps<0)
		within_steps=g.nodes();
	if(within_steps>g.nodes())
		within_steps=g.nodes();
	if(within_steps>g.edges())
		within_steps=g.edges();
	if(constraintsBuilt>=within_steps)
		return;



	assert(outer->decisionLevel()==0);
	vec<Lit> c;

	if(constraintsBuilt<=0){
		constraintsBuilt=0;
		dist_lits.push();
		Lit True = mkLit(outer->newVar());
		outer->addClause(True);
		assert(outer->value(True)==l_True);
		Lit False = ~True;
		for(int i = 0;i<g.nodes();i++){
			dist_lits[0].push(False);
		}
		dist_lits[0][source]=True;
	}

	vec<Lit>  reaches;


	//bellman-ford, unrolled:
	for (int i = constraintsBuilt;i<within_steps;i++){
		dist_lits.last().copyTo(reaches);
		dist_lits.push();
		reaches.copyTo(dist_lits.last());
		assert(outer->value( reaches[source])==l_True);

		for(int j = 0;j<outer->nNodes();j++){
			Lit r_cur = reaches[j];
			for(Edge & e: outer->undirected_adj[j]){
				//Edge e = outer->edges[j][k];

				int from = e.from;
				int to =e.to;
				if(from==j){
					std::swap(from,to);
				}
				if(outer->value(dist_lits.last()[to])==l_True){
					//do nothing
				}else if (outer->value(reaches[from])==l_False){
					//do nothing
				}else{
					Lit l = mkLit(e.v,false);

					Lit r = mkLit( outer->newVar(), false);

					c.clear();
					c.push(~r);c.push(r_cur);c.push(l); //r -> (e.l or reaches[e.to])
					outer->addClause(c);
					c.clear();
					c.push(~r);c.push(r_cur);c.push(reaches[from]); //r -> (reaches[e.from]) or reaches[e.to])
					outer->addClause(c);
					c.clear();
					c.push(r);c.push(~r_cur); //~r -> ~reaches[e.to]
					outer->addClause(c);
					c.clear();
					c.push(r);c.push(~reaches[from]);c.push(~l); //~r -> (~reaches[e.from] or ~e.l)
					outer->addClause(c);
					r_cur=r;

				}
			}
			dist_lits.last()[j]=r_cur;
		}


	}

	if(within_steps==g.nodes() || within_steps==g.edges()){
		assert(reach_lits.size()==0);
		for(Lit d: dist_lits.last()){
			reach_lits.push(d);
		}
	}
	constraintsBuilt=within_steps;
}

void ConnectDetector::addLit(int from, int to, Var outer_reach_var){
	g.invalidate();
	antig.invalidate();
	Var reach_var = outer->newVar(outer_reach_var,getID());

	if(first_reach_var==var_Undef){
		first_reach_var=reach_var;
	}else{
		assert(reach_var>first_reach_var);
	}
	assert(from==source);


	 if(undirectedalg ==ConnectivityAlg::ALG_SAT){
		 buildSATConstraints();
	 }
		while( reach_lits.size()<=to)
				reach_lits.push(lit_Undef);
	Lit reachLit=mkLit(reach_var,false);

	if(reach_lits[to]==lit_Undef){
		assert(undirectedalg !=ConnectivityAlg::ALG_SAT);
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

void ConnectDetector::ReachStatus::setReachable(int u, bool reachable){
	if(reachable){
		assert(!detector.outer->dbg_notreachable( detector.source,u,true));
	}else{
		assert(detector.outer->dbg_notreachable( detector.source,u,true));
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

void ConnectDetector::ReachStatus::setMininumDistance(int u, bool reachable, int distance){
	assert(reachable ==(distance<detector.outer->g.nodes()));
	setReachable(u,reachable);

/*	if(u<detector.dist_lits.size()){
		assert(distance>=0);
		assert(distance<detector.outer->g.nodes());
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
	}*/
}

bool ConnectDetector::ChokepointStatus::mustReach(int node){
	Lit l =  detector.reach_lits[node];
	if(l!=lit_Undef){
		return detector.outer->value(l)==l_True;
	}
	return false;
}
bool ConnectDetector::ChokepointStatus::operator() (int edge_id){
	return detector.outer->value(detector.outer->edge_list[ edge_id].v)==l_Undef;
}

void ConnectDetector::preprocess(){
	//vec<bool> pure;
	//pure.growTo(reach_lits.size());
	//can check if all reach lits appear in only one polarity in the solver constraints; if so, then we can disable either check_positive or check_negative
}

void ConnectDetector::buildReachReason(int node,vec<Lit> & conflict){
			//drawFull();
			Reach & d = *positive_path_detector;


			double starttime = cpuTime();
			d.update();

			assert(outer->dbg_reachable(d.getSource(),node,true));

			assert(d.connected_unchecked(node));
			if(opt_learn_reaches ==0 || opt_learn_reaches==2){
				int u = node;
				int p;
				while(( p = d.previous(u)) != -1){
					Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
					Var e =edg.v;
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
					Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
					Var e =edg.v;
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
		void ConnectDetector::buildNonReachReason(int node,vec<Lit> & conflict){
			static int it = 0;
			++it;
			int u = node;
			//drawFull( non_reach_detectors[detector]->getSource(),u);
			assert(outer->dbg_notreachable( source,u,true));
			double starttime = cpuTime();
			outer->cutGraph.clearHistory();
			outer->stats_mc_calls++;
			if(opt_conflict_min_cut){
				if(mincutalg!= MinCutAlg::ALG_EDKARP_ADJ){
					//ok, set the weights for each edge in the cut graph.
					//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
					for(int u = 0;u<outer->cutGraph.nodes();u++){
						for(int j = 0;j<outer->cutGraph.nIncident(u,true);j++){
							int v = outer->cutGraph.incident(u,j,true).node;
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
							outer->mc->setCapacity(outer->trail[i].to, outer->trail[i].from,1);
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
				    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
				    	for(int i = 0;i<outer->undirected_adj[u].size();i++){
				    		int v = outer->undirected_adj[u][i].v;
				    		int from = outer->undirected_adj[u][i].from;
				    		int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;
				    		assert(outer->undirected_adj[u][i].to==u);
				    		if(from==u){
								assert(outer->edge_list[edge_num].to == u);
								assert(outer->edge_list[edge_num].from == u);
								continue;//Self loops are allowed, but just make sure nothing got flipped around...
							}
							assert(from!=u);


				    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...

				    		if(outer->value(v)==l_False){
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
		void ConnectDetector::buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict){
					static int it = 0;
					++it;

					assert(outer->value(outer->edge_list[forced_edge_id].v)==l_True);
					Lit edgeLit =mkLit( outer->edge_list[forced_edge_id].v,false);

					conflict.push(edgeLit);

					int forced_edge_from = outer->edge_list[forced_edge_id].from;
					int forced_edge_to= outer->edge_list[forced_edge_id].to;


					int u = reach_node;
					//drawFull( non_reach_detectors[detector]->getSource(),u);
					assert(outer->dbg_notreachable( source,u,true));
					double starttime = cpuTime();
					outer->cutGraph.clearHistory();
					outer->stats_mc_calls++;



					if(opt_conflict_min_cut){
						if(mincutalg!= MinCutAlg::ALG_EDKARP_ADJ){
							//ok, set the weights for each edge in the cut graph.
							//Set edges to infinite weight if they are undef or true, and weight 1 otherwise.
							for(int u = 0;u<outer->cutGraph.nodes();u++){
								for(int j = 0;j<outer->cutGraph.nIncident(u);j++){
									int v = outer->cutGraph.incident(u,j).node;
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
						    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
						    	for(int i = 0;i<outer->undirected_adj[u].size();i++){
						    		int v = outer->undirected_adj[u][i].v;
						    		int from = outer->undirected_adj[u][i].from;
						    		assert(from!=u);
						    		assert(outer->undirected_adj[u][i].to==u);
						    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
						    		int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;

						    		if(edge_num == forced_edge_id || outer->value(v)==l_False){
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

		void ConnectDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){


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
		static int iter = 0;
		bool ConnectDetector::propagate(vec<Lit> & conflict){
			if(!positive_reach_detector)
				return true;

			if(++iter==68305){
				int a=1;
			}
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
			static bool b = false;
			if(!b){
				printf("Warning: Debug code enabled!!\n");
				b=true;
			}
					for(int i = 0;i<reach_lits.size();i++){
						Lit l = reach_lits[i];
						if(l!=lit_Undef){
							int u = getNode(var(l));
							if(positive_reach_detector->connected(u)){
								assert(outer->value(l)==l_True);
								if(outer->value(l)!=l_True){
									exit(3);
								}

							}else if (!negative_reach_detector->connected(u)){
								assert(outer->value(l)==l_False);
								if(outer->value(l)!=l_False){
									exit(3);
								}
							}
						}

					}
			#endif



			return true;
		}

bool ConnectDetector::checkSatisfied(){
	if(positive_reach_detector){
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
	}else{
		Dijkstra<Reach::NullStatus,true>under(source,g) ;
		Dijkstra<Reach::NullStatus,true>over(source,antig) ;
		under.update();
		over.update();

		for(int j = 0;j< reach_lits.size();j++){
			Lit l = reach_lits[j];
			if(l!=lit_Undef){
				int node =j;

				if(outer->value(l)==l_True){
					if(!under.connected(node)){
						return false;
					}
				}else if (outer->value(l)==l_False){
					if( over.connected(node)){
						return false;
					}
				}else{
					if(over.connected(node)){
						return false;
					}
					if(!under.connected(node)){
						return false;
					}
				}
			}
		}
	}
	return true;
}

void ConnectDetector::dbg_sync_reachability(){
#ifndef NDEBUG
	if(!positive_reach_detector)
			return;
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
#endif
	}


int ConnectDetector::OptimalWeightEdgeStatus::operator [] (int edge) const {
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
int ConnectDetector::OptimalWeightEdgeStatus::size()const{
	return detector.outer->edge_list.size();
}


Lit ConnectDetector::decide(){
	if(!negative_reach_detector)
		return lit_Undef;
	auto * over = negative_reach_detector;

	auto * under = positive_reach_detector;

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
				int last_edge=-1;
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
							last_edge= opt_path->incomingEdge(p);
							int prev = opt_path->previous(p);
							p = prev;

						}
					}else{
						//ok, read back the path from the over to find a candidate edge we can decide
						//find the earliest unconnected node on this path
						over->update();
						 p = j;
						 last = j;
						 static int local_it=0;
						 if(++local_it==240){
							 int a=1;
						 }
						while(!under->connected(p)){


							last=p;
							assert(p!=source);
							assert(p!=-1);
							last_edge= over->incomingEdge(p);
							int prev = over->previous(p);
							assert(prev!=-1);
							p = prev;
							if(prev<0){
								printf("%d\n",local_it);
								exit(3);
							}
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
					 last_edge=-1;
					 assert(rnd_path->connected(p));
					while(!under->connected(p)){

						last=p;
						assert(p!=source);
						last_edge=rnd_path->incomingEdge(p);
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

				assert(last_edge>=0);
				assert(outer->edge_list[last_edge].edgeID==last_edge);
				Var v = outer->edge_list[last_edge].v;
				if(outer->value(v)==l_Undef){
					return mkLit(v,false);
				}else{
					assert(outer->value(v)!=l_True);
					if(outer->value(v)!=l_True){
						exit(3);
					}
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
						int incoming_edge = over->incomingEdge(p);
						Var v = outer->edge_list[incoming_edge].v;
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


