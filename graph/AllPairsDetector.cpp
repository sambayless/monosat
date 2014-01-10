/*
 * AllPairsDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "AllPairsDetector.h"
#include "GraphTheory.h"

AllPairsDetector::AllPairsDetector(int _detectorID, GraphTheorySolver * _outer,  DynamicGraph<PositiveEdgeStatus> &_g,DynamicGraph<NegativeEdgeStatus> &_antig, int within_steps ,double seed):
Detector(_detectorID),outer(_outer),within(within_steps),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){


		positiveReachStatus = new AllPairsDetector::ReachStatus(*this,true);
		negativeReachStatus = new AllPairsDetector::ReachStatus(*this,false);
		if(allpairsalg==ALG_FLOYDWARSHALL){
			positive_reach_detector = new FloydWarshall<AllPairsDetector::ReachStatus,PositiveEdgeStatus>(_g,*(positiveReachStatus),1);
			negative_reach_detector = new FloydWarshall<AllPairsDetector::ReachStatus,NegativeEdgeStatus>(_antig,*(negativeReachStatus),-1);
			positive_path_detector = positive_reach_detector;
		}else{
			positive_reach_detector = new DijkstraAllPairs<AllPairsDetector::ReachStatus,PositiveEdgeStatus>(_g,*(positiveReachStatus),1);
			negative_reach_detector = new DijkstraAllPairs<AllPairsDetector::ReachStatus,NegativeEdgeStatus>(_antig,*(negativeReachStatus),-1);
			positive_path_detector = positive_reach_detector;
		}
		/*	if(opt_conflict_shortest_path)
			reach_detectors.last()->positive_dist_detector = new Dijkstra<PositiveEdgeStatus>(from,g);*/


	first_reach_var = var_Undef;

}



void AllPairsDetector::addLit(int from, int to, Var reach_var,int within_steps){
	if(first_reach_var==var_Undef){
		first_reach_var=reach_var;
	}else{
		assert(reach_var>first_reach_var);
	}
	if(within_steps<0)
		within_steps=outer->nNodes();
	if(within_steps>outer->nNodes())
		within_steps=outer->nNodes();

	dist_lits.growTo(outer->nNodes());

	dist_lits[from].growTo(outer->nNodes());
	installed_sources.growTo(outer->nNodes());

	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});

	while(outer->S->nVars()<=reach_var)
		outer->S->newVar();

	if(!installed_sources[from]){
		installed_sources[from]=true;
		sources.push(from);
		positive_reach_detector->addSource(from);
		negative_reach_detector->addSource(from);

	}

	Lit reachLit=mkLit(reach_var,false);
	bool found=false;
	for(int i = 0;i<dist_lits[from][to].size();i++){
		if(dist_lits[from][to][i].min_distance==within_steps){
			found=true;
			Lit r = dist_lits[from][to][within_steps].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->S->addClause(~r, reachLit);
			outer->S->addClause(r, ~reachLit);
		}
	}
	if(!found){
		dist_lits[from][to].push();
		dist_lits[from][to].last().l = reachLit;
		dist_lits[from][to].last().min_distance=within_steps;
		dist_lits[from][to].last().source=from;
		while(reach_lit_map.size()<= reach_var- first_reach_var ){
			reach_lit_map.push({lit_Undef,-1,-1,-1});
		}

		reach_lit_map[reach_var-first_reach_var]={reachLit,within_steps,from,to};
	}



}

void AllPairsDetector::ReachStatus::setReachable(int from,int u, bool reachable){
			/*	if(polarity==reachable && u<detector.reach_lits.size()){
					Lit l = detector.reach_lits[u];
					if(l!=lit_Undef){
						lbool assign = detector.outer->S->value(l);
						if(assign!= (reachable? l_True:l_False )){
							detector.changed.push({reachable? l:~l,u});
						}
					}
				}*/
			}

void AllPairsDetector::ReachStatus::setMininumDistance(int from, int u, bool reachable, int distance){
	assert(reachable ==(distance<detector.outer->g.nodes));
	if(distance<=detector.outer->g.nodes){
		setReachable(from,u,reachable);
	}

		if(u<detector.dist_lits.size()){
			assert(distance>=0);

			for(int i = 0;i<detector.dist_lits[from][u].size();i++){
				int d =  detector.dist_lits[from][u][i].min_distance;

				Lit l = detector.dist_lits[from][u][i].l;
				if(l!=lit_Undef){

					assert(l!=lit_Undef);
					if(d<distance && !polarity){
						lbool assign = detector.outer->S->value(l);
						if( assign!= l_False ){
							detector.changed.push({~l,u,from});
						}
					}else if(d>=distance && polarity){
						lbool assign = detector.outer->S->value(l);
						if( assign!= l_True ){
							detector.changed.push({l,u,from});
						}
					}
				}
			}
		}

}


void AllPairsDetector::buildReachReason(int source, int to,vec<Lit> & conflict){
			//drawFull();
			AllPairs & d = *positive_path_detector;


			double starttime = cpuTime();
			d.update();

			tmp_path.clear();

			assert(d.connected_unchecked(source,to));
			d.getPath(source,to,tmp_path);
			//if(opt_learn_reaches ==0 || opt_learn_reaches==2)
			{
				int u = to;

				//while(( p = d.previous(u)) != -1){
				for(int i = tmp_path.size()-2;i>=0;i--){
					int p = tmp_path[i];
					assert(p!=-1);
					Edge edg = outer->edges[p][u];
					Var e =outer->edges[p][u].v;
					lbool val = outer->S->value(e);
					assert(outer->S->value(e)==l_True);
					conflict.push(mkLit(e, true));
					u = p;

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
		void AllPairsDetector::buildNonReachReason(int source, int node,vec<Lit> & conflict){
			static int it = 0;
			++it;
			int u = node;
			//drawFull( non_reach_detectors[detector]->getSource(),u);
			//assert(outer->dbg_distance( source,u));
			double starttime = cpuTime();
			outer->cutGraph.clearHistory();
			outer->stats_mc_calls++;
			{

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
				    	//assert(negative_reach_detector->distance_unsafe(u)>d);
				    	//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
				    	for(int i = 0;i<outer->inv_adj[u].size();i++){
				    		int v = outer->inv_adj[u][i].v;
				    		int from = outer->inv_adj[u][i].from;
				    		assert(from!=u);
				    		assert(outer->inv_adj[u][i].to==u);
				    		//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				    		int edge_num = v-outer->min_edge_var;

				    		if(outer->edge_assignments[edge_num]==l_False){
				    			//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
				    			//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
				    			//if(!seen[from])
				    				conflict.push(mkLit(v,false));
				    		}else if (from!=source){
				    			//for distance analysis, we _can_ end up reaching source.

				    			//even if it is undef? probably...
				    			if(!seen[from]){
				    				seen[from]=true;
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
		void AllPairsDetector::buildForcedEdgeReason(int source, int to, int forced_edge_id,vec<Lit> & conflict){
					static int it = 0;
					++it;

					assert(outer->edge_assignments[forced_edge_id]==l_True);
					Lit edgeLit =mkLit( outer->edge_list[forced_edge_id].v,false);

					conflict.push(edgeLit);

					int forced_edge_from = outer->edge_list[forced_edge_id].from;
					int forced_edge_to= outer->edge_list[forced_edge_id].to;


					int u = to;
					//drawFull( non_reach_detectors[detector]->getSource(),u);
					assert(outer->dbg_notreachable( source,u));
					double starttime = cpuTime();
					outer->cutGraph.clearHistory();
					outer->stats_mc_calls++;



					//We could learn an arbitrary (non-infinite) cut here, or just the whole set of false edges
					//or perhaps we can learn the actual 1-uip cut?


					vec<int>& to_visit  = outer->to_visit;
					vec<char>& seen  = outer->seen;

					to_visit.clear();
					to_visit.push(to);
					seen.clear();
					seen.growTo(outer->nNodes());
					seen[to]=true;

					do{

						assert(to_visit.size());
						int u = to_visit.last();
						assert(u!=source);
						to_visit.pop();
						assert(seen[u]);
						assert(!negative_reach_detector->connected_unsafe(source,u));
						//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming edges
						for(int i = 0;i<outer->inv_adj[u].size();i++){
							int v = outer->inv_adj[u][i].v;
							int from = outer->inv_adj[u][i].from;
							assert(from!=u);
							assert(outer->inv_adj[u][i].to==u);
							//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
							int edge_num = v-outer->min_edge_var;

							if(edge_num == forced_edge_id || outer->edge_assignments[edge_num]==l_False){
								//note: we know we haven't seen this edge variable before, because we know we haven't visited this node before
								//if we are already planning on visiting the from node, then we don't need to include it in the conflict (is this correct?)
								//if(!seen[from])
									conflict.push(mkLit(v,false));
							}else if (from!=source){
								//assert(from!=source);
								//even if it is undef? probably...
								if(!seen[from]){
									seen[from]=true;
									to_visit.push(from);
								}
							}
						}
					}  while (to_visit.size());




					 outer->num_learnt_cuts++;
					 outer->learnt_cut_clause_length+= (conflict.size()-1);

					double elapsed = cpuTime()-starttime;
					 outer->mctime+=elapsed;

			#ifdef DEBUG_GRAPH
					 assert(outer->dbg_clause(conflict));
			#endif
				}

		void AllPairsDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){


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
					int source = getSource(v);
					buildReachReason(source,u,reason);

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
					int source = getSource(v);
					buildNonReachReason(source,t,reason);

				}else if (marker==forced_reach_marker){
					Var v = var(p);
					//The forced variable is an EDGE that was forced.
					int forced_edge_id =v- outer->min_edge_var;
					//The corresponding node that is the reason it was forced
					int reach_node=force_reason[forced_edge_id];
					int source = getSource(v);
					 buildForcedEdgeReason(source, reach_node,  forced_edge_id, reason);
				}else{
					assert(false);
				}
		}

		bool AllPairsDetector::propagate(vec<Assignment> & trail,vec<Lit> & conflict){


		double startdreachtime = cpuTime();
		getChanged().clear();
		positive_reach_detector->update();
		double reachUpdateElapsed = cpuTime()-startdreachtime;
		outer->reachupdatetime+=reachUpdateElapsed;

		double startunreachtime = cpuTime();
		negative_reach_detector->update();
		double unreachUpdateElapsed = cpuTime()-startunreachtime;
		outer->unreachupdatetime+=unreachUpdateElapsed;

		for(int j = 0;j<getChanged().size();j++){
				Lit l = getChanged()[j].l;
				int u =  getChanged()[j].u;
				int source = getChanged()[j].source;
				bool reach = !sign(l);

				if(outer->S->value(l)==l_True){
					//do nothing
				}else if(outer->S->value(l)==l_Undef){
#ifdef DEBUG_GRAPH
					assert(outer->dbg_propgation(l));
#endif
#ifdef DEBUG_SOLVER
					if(S->dbg_solver)
						S->dbg_check_propagation(l);
#endif
					trail.push(Assignment(false,reach,detectorID,0,var(l)));
					if(reach)
						outer->S->uncheckedEnqueue(l,reach_marker) ;
					else
						outer->S->uncheckedEnqueue(l,non_reach_marker) ;

				}else if (outer->S->value(l)==l_False){
					conflict.push(l);

					if(reach){

					//conflict
					//The reason is a path in g from to s in d
					buildReachReason(source,u,conflict);
					//add it to s
					//return it as a conflict

					}else{
						//The reason is a cut separating s from t
						buildNonReachReason(source,u,conflict);

					}
#ifdef DEBUG_GRAPH
					for(int i = 0;i<conflict.size();i++)
								 assert(outer->S->value(conflict[i])==l_False);
#endif
#ifdef DEBUG_SOLVER
					if(S->dbg_solver)
						S->dbg_check(conflict);
#endif


					return false;
				}else{
					int  a=1;
				}

			}

			#ifdef DEBUG_GRAPH
				for(int m = 0;m<sources.size();m++){
					int s = sources[m];
					for(int i = 0;i<dist_lits[s].size();i++){
						for(int j = 0;j<dist_lits[s][i].size();j++){
							Lit l = dist_lits[s][i][j].l;
							int dist =  dist_lits[s][i][j].min_distance;
							if(l!=lit_Undef){
								int u = getNode(var(l));
								if(positive_reach_detector->distance_unsafe(s,u)<=dist){
									assert(outer->S->value(l)==l_True);
								}else if (negative_reach_detector->distance_unsafe(s,u)>dist){
									int d =negative_reach_detector->distance_unsafe(s,u);
									assert(outer->S->value(l)==l_False);
								}
							}
						}
					}
				}
			#endif
			return true;
		}

bool AllPairsDetector::checkSatisfied(){
	for(int m = 0;m<sources.size();m++){
		int s = sources[m];
				for(int j = 0;j< dist_lits[s].size();j++){
					for(int k = 0;k<dist_lits[s][j].size();k++){
						Lit l = dist_lits[s][j][k].l;
						int dist = dist_lits[s][j][k].min_distance;
						DistLit dl = dist_lits[s][j][k];
						assert(s==dist_lits[s][j][k].source);
						if(l!=lit_Undef){
							int node =getNode(var(l));

							if(outer->S->value(l)==l_True){
								if(positive_reach_detector->distance(s,node)>dist){
									return false;
								}
							}else if (outer->S->value(l)==l_False){
								if( negative_reach_detector->distance(s,node)<=dist){
									return false;
								}
							}else{
								if(positive_reach_detector->distance(s,node)<=dist){
									return false;
								}
								if(!negative_reach_detector->distance(s,node)>dist){
									return false;
								}
							}
						}
					}
				}
	}
	return true;
}
Lit AllPairsDetector::decide(){

	//AllPairs * over = (FloydWarshall<AllPairsDetector::ReachStatus,NegativeEdgeStatus>*) negative_reach_detector;

	//FloydWarshall<AllPairsDetector::ReachStatus,PositiveEdgeStatus> * under = (FloydWarshall<AllPairsDetector::ReachStatus,PositiveEdgeStatus>*)positive_reach_detector;

	//we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

	//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

	//this can be obviously more efficient
	//for(int j = 0;j<nNodes();j++){
	for(int m = 0;m<sources.size();m++){
		int s = sources[m];
		for(int k = 0;k<dist_lits[s].size();k++){
			for(int n = 0;n<dist_lits[s][k].size();n++){
				Lit l = dist_lits[s][k][n].l;
				int min_dist = dist_lits[s][k][n].min_distance;
				if(l==lit_Undef)
					continue;
				int j = getNode(var(l));
				if(outer->S->value(l)==l_True){
					//if(S->level(var(l))>0)
					//	continue;

					assert(negative_reach_detector->distance(s,j)<=min_dist);//else we would already be in conflict before this decision was attempted!
					if(positive_reach_detector->distance(s,j)>min_dist){
						//then lets try to connect this
						static vec<bool> print_path;

						assert(negative_reach_detector->connected(s,j));//Else, we would already be in conflict
						int p =j;
						int last=j;
						//if(!opt_use_random_path_for_decisions)
						{
							//ok, read back the path from the over to find a candidate edge we can decide
							//find the earliest unconnected node on this path
							negative_reach_detector->update();
							negative_reach_detector->getPath(s,j,tmp_path);
							 p = j;
							 last = j;
							 if(tmp_path.size()){
								 tmp_path.pop();
								while(!positive_reach_detector->connected(s,p)){

									last=p;
									assert(p!=s);
									assert(tmp_path.size());
									int prev = tmp_path.last(); //negative_reach_detector->previous(s,p);
									tmp_path.pop();
									p = prev;

								}
							 }
						}

						for(int k = 0;k<outer->antig.adjacency[p].size();k++){
							int to = outer->antig.adjacency[p][k].node;
							if (to==last){
								Var v =outer->edge_list[ outer->antig.adjacency[p][k].id].v;
								if(outer->S->value(v)==l_Undef){
									return mkLit(v,false);
								}else{
									assert(outer->S->value(v)!=l_True);
								}
							}
						}
						assert(false);
					}
				}
			}
		}
	}
	return lit_Undef;
};


