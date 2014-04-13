/*
 * ConnectedComponentsDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "ConnectedComponentsDetector.h"
#include "GraphTheory.h"
#include <limits>
ConnectedComponentsDetector::ConnectedComponentsDetector(int _detectorID, GraphTheorySolver * _outer,  DynamicGraph<PositiveEdgeStatus> &_g,DynamicGraph<NegativeEdgeStatus> &_antig,double seed):
Detector(_detectorID),outer(_outer),g(_g),antig(_antig),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){



		positiveReachStatus = new ConnectedComponentsDetector::ConnectedComponentsStatus(*this,true);
		negativeReachStatus = new ConnectedComponentsDetector::ConnectedComponentsStatus(*this,false);
		//Note: these are _intentionalyl_ swapped
		negative_reach_detector = new DisjointSetsConnectedComponents<ConnectedComponentsDetector::ConnectedComponentsStatus,PositiveEdgeStatus>(_g,*(negativeReachStatus),1);
		positive_reach_detector = new DisjointSetsConnectedComponents<ConnectedComponentsDetector::ConnectedComponentsStatus,NegativeEdgeStatus>(_antig,*(positiveReachStatus),1);
		reach_marker=outer->newReasonMarker(getID());
		non_reach_marker=outer->newReasonMarker(getID());
		//forced_reach_marker=outer->newReasonMarker(getID());
}
void ConnectedComponentsDetector::addConnectedComponentsLit(Var outer_weight_var,int min_components){

	g.invalidate();
	antig.invalidate();
	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});

	Var weight_var = outer->newVar(outer_weight_var,getID());
/*	while(outer->S->nVars()<=weight_var)
		outer->S->newVar();*/

	Lit reachLit=mkLit(weight_var,false);
	bool found=false;
	for(int i = 0;i<connected_components_lits.size();i++){
		if(connected_components_lits[i].min_components==min_components){
			found=true;
			Lit r = connected_components_lits[i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			//outer->S->addClause(~r, reachLit);
			//outer->S->addClause(r, ~reachLit);
			outer->makeEqual(r,reachLit);
		}
	}
	if(!found){
		connected_components_lits.push();
		connected_components_lits.last().l = reachLit;
		connected_components_lits.last().min_components=min_components;
		//weight_lit_map.insert(min_weight,weight_lits.size()-1);
	}



}




/*void ConnectedComponentsDetector::ConnectedComponentsStatus::setConnected(int edgeid, bool in_tree){
	if(edgeid<detector.tree_edge_lits.size()){
		Lit l = detector.tree_edge_lits[edgeid].l;
		if(l!=lit_Undef)
			detector.changed_edges.push({in_tree? l:~l,edgeid});
	}
}*/
void ConnectedComponentsDetector::ConnectedComponentsStatus::setConnected(int u, int v, bool connected){
	/*if(u>v){
		std::swap(u,v);
	}
	if (detector.reachLits[u][v]!=lit_Undef){
		Lit l = detector.reachLits[u][v];
		lbool assign = detector.outer->value(l);
		if(assign==l_True && connected){
			//do nothing
		}else if(assign==l_False && !connected){
			//do nothing
		}else
			detector.changed.push({connected? l:~l,u,v});
	}*/
}

void ConnectedComponentsDetector::ConnectedComponentsStatus::setComponents(int components){

			for(int i = 0;i<detector.connected_components_lits.size();i++){
				int min_components =  detector.connected_components_lits[i].min_components;

				Lit l = detector.connected_components_lits[i].l;
				if(l!=lit_Undef){

					assert(l!=lit_Undef);
					if(components<=min_components && !polarity){
						lbool assign = detector.outer->value(l);
						if( assign!= l_False ){
							detector.changed_weights.push({~l,min_components});
						}
					}else if(components>min_components && polarity){
						lbool assign = detector.outer->value(l);
						if( assign!= l_True ){
							detector.changed_weights.push({l,min_components});
						}
					}
				}
			}


}




		void ConnectedComponentsDetector::buildMinComponentsTooLowReason(int min_components,vec<Lit> & conflict){

			double starttime = rtime(2);


			seen.clear();
			seen.growTo(g.nodes);

			visit.clear();
			//ok,construct spanning forest from among the connected elements using dfs (we don't need to use kruskal's, as the graph is unweighted), and learn that at least one edge in that each must be disabled.
			for(int k = 0;k<negative_reach_detector->numComponents();k++){
			//learn that at least one edge in the tree must be disabled (else, number of connected components cannot be increased)
				int root =negative_reach_detector->getComponent(k);
				seen[root]=true;
				visit.push(root);

				while(visit.size()){
					int u = visit.last();
					visit.pop();

					for(int i = 0;i<g.adjacency_undirected[u].size();i++){
						int v = g.adjacency_undirected[u][i].node;
						int edgeid =  g.adjacency_undirected[u][i].id;
						if( g.edgeEnabled(edgeid) && ! seen[v]){
							seen[v]=true;
							visit.push(v);
							//this edge is in the spanning tree
							Var v = outer->edge_list[edgeid].v;
							assert(outer->value(v)==l_True);
							conflict.push(mkLit(v,true));
						}
					}
				}
			}
			outer->num_learnt_paths++;
			outer->learnt_path_clause_length+= (conflict.size()-1);
			double elapsed = rtime(2)-starttime;
			outer->pathtime+=elapsed;
		}
		void ConnectedComponentsDetector::buildMinComponentsTooHighReason(int min_components,vec<Lit> & conflict){
			static int it = 0;
						++it;

				//drawFull( non_reach_detectors[detector]->getSource(),u);
				//assert(outer->dbg_distance( source,u));

				double starttime = rtime(2);
				int INF=std::numeric_limits<int>::max();
				positive_reach_detector->update();
				int numComponents = positive_reach_detector->numComponents();
				assert(numComponents>=min_components);

					//IF the ConnectedComponents is disconnected, the reason is a separating cut between any two disconnected components.
					//we can find one of these by identifying any two roots

					int min_conflict = INF;

					//walk back down from the each root to find a separating cut of disabled edge.
					//return the smallest such cut.

					edge_in_clause.clear();
					edge_in_clause.growTo(g.nEdgeIDs());
					assert(conflict.size()==1);
					for(int i = 0;i<positive_reach_detector->numComponents();i++){
						int root = positive_reach_detector->getComponent(i);

						//now explore that component.

						visit.clear();
						//ok, now traverse the connected component in the tree below this root.
						seen.clear();
						seen.growTo(g.nodes);

						visit.push(root);
						seen[root]=true;
						while(visit.size()){
							int u = visit.last();
							visit.pop();

							for(int i = 0;i<antig.adjacency_undirected[u].size();i++){
								int edgeid = antig.adjacency_undirected[u][i].id;
								if(antig.edgeEnabled(edgeid)){
									int v = antig.adjacency_undirected[u][i].node;

									if( ! seen[v]){
										//u is a child of node in the minimum spanning tree
										seen[v]=true;
										visit.push(v);
									}
								}else if (!antig.edgeEnabled(edgeid)){
									int v = antig.adjacency_undirected[u][i].node;
									if( ! edge_in_clause[edgeid]){
										edge_in_clause[edgeid]=true;
										Var e =outer->edge_list[edgeid].v;

										assert(outer->value(e)==l_False);
										conflict.push(mkLit(e,false));
									}
								}
							}
						}



					}

					return;

		}
		void ConnectedComponentsDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){

				if(marker==reach_marker){
					reason.push(p);

					Var v = var(p);
					int weight=-1;
					//could swap this out for a map if there are lots of lits..
					for(int i = 0;i<connected_components_lits.size();i++){
						if(var(connected_components_lits[i].l)==v){
							weight=connected_components_lits[i].min_components;
							break;
						}
					}
					assert(weight>=0);
					buildMinComponentsTooHighReason(weight,reason);

					//double elapsed = rtime(2)-startpathtime;
				//	pathtime+=elapsed;
				}else if(marker==non_reach_marker){
					reason.push(p);

					//the reason is a cut separating p from s;
					//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

					//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
					//assign the mincut edge weights if they aren't already assigned.


					Var v = var(p);

					int weight=-1;
					//could swap this out for a map if there are lots of lits..
					for(int i = 0;i<connected_components_lits.size();i++){
						if(var(connected_components_lits[i].l)==v){
							weight=connected_components_lits[i].min_components;
							break;
						}
					}
					assert(weight>=0);

					buildMinComponentsTooLowReason(weight,reason);

				}else{
					assert(false);
				}
		}

		bool ConnectedComponentsDetector::propagate(vec<Assignment> & trail,vec<Lit> & conflict){


		double startdreachtime = rtime(2);
		changed.clear();
		changed_weights.clear();
		positive_reach_detector->update();
		double reachUpdateElapsed = rtime(2)-startdreachtime;
		outer->reachupdatetime+=reachUpdateElapsed;

		double startunreachtime = rtime(2);
		negative_reach_detector->update();
		double unreachUpdateElapsed = rtime(2)-startunreachtime;
		outer->unreachupdatetime+=unreachUpdateElapsed;

		for(int j = 0;j<changed_weights.size();j++){
			Lit l = changed_weights[j].l;
			int components = changed_weights[j].min_components;
			bool reach = !sign(l);
			if(outer->value(l)==l_True){
				//do nothing
			}else if(outer->value(l)==l_Undef){
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
					buildMinComponentsTooHighReason(components,conflict);
				//add it to s
				//return it as a conflict

				}else{
					buildMinComponentsTooLowReason(components,conflict);


				}

				return false;
			}else{
				int  a=1;
			}

		}/*
		for(int j = 0;j<changed.size();j++){
					Lit l = changed[j].l;
					int components = changed_weights[j].min_components;
					bool reach = !sign(l);
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){
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
							buildMinComponentsTooHighReason(components,conflict);
						//add it to s
						//return it as a conflict

						}else{
							buildMinComponentsTooLowReason(components,conflict);


						}

						return false;
					}else{
						int  a=1;
					}

				}*/
			return true;
		}

bool ConnectedComponentsDetector::checkSatisfied(){
	int numConnected = positive_reach_detector->numComponents();
	int numConnectedOver = negative_reach_detector->numComponents();
					for(int k = 0;k<connected_components_lits.size();k++){
						Lit l = connected_components_lits[k].l;
						int moreThanThisManyComponents = connected_components_lits[k].min_components;

						if(l!=lit_Undef){


							if(outer->value(l)==l_True){
								if(negative_reach_detector->numComponents()<=moreThanThisManyComponents){
									return false;
								}
							}else if (outer->value(l)==l_False){
								if( positive_reach_detector->numComponents()>moreThanThisManyComponents){
									return false;
								}
							}else{
								if(negative_reach_detector->numComponents()>moreThanThisManyComponents){
									return false;
								}
								if(!positive_reach_detector->numComponents()<=moreThanThisManyComponents){
									return false;
								}
							}
						}
					}

	return true;
}
Lit ConnectedComponentsDetector::decide(){
	/*ConnectedComponentsDetector *r =this;
	MinimumSpanningTree<ConnectedComponentsDetector::ConnectedComponentsStatus,NegativeEdgeStatus> * over = (MinimumSpanningTree<ConnectedComponentsDetector::ConnectedComponentsStatus,NegativeEdgeStatus>*) r->negative_reach_detector;

	MinimumSpanningTree<ConnectedComponentsDetector::ConnectedComponentsStatus,PositiveEdgeStatus> * under = (MinimumSpanningTree<ConnectedComponentsDetector::ConnectedComponentsStatus,PositiveEdgeStatus>*) r->positive_reach_detector;

	//we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

	//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

	//this can be obviously more efficient
	//for(int j = 0;j<nNodes();j++){
	if(opt_decide_graph_neg){

	}

	for(int k = 0;k<dist_lits.size();k++){
		for(int n = 0;n<dist_lits[k].size();n++){
			Lit l = dist_lits[k][n].l;
			int min_weight = dist_lits[k][n].min_weight;
			if(l==lit_Undef)
				continue;
			int j = r->getNode(var(l));
			if(outer->value(l)==l_True){
				if(opt_decide_graph_pos){
				//if(S->level(var(l))>0)
				//	continue;

				assert(over->distance(j)<=min_dist);//else we would already be in conflict before this decision was attempted!
				if(under->distance(j)>min_dist){
					//then lets try to connect this
					//static vec<bool> print_path;

					assert(over->connected(j));//Else, we would already be in conflict


					int p =j;
					int last=j;
					if(!opt_use_random_path_for_decisions)
					{
						//ok, read back the path from the over to find a candidate edge we can decide
						//find the earliest unconnected node on this path
						over->update();
						 p = j;
						 last = j;
						 int dist = 0;
						while(under->distance(p)>=min_dist-dist){

							last=p;
							assert(p!=r->source);
							int prev = over->previous(p);
							assert(over->distance(p)<=min_dist-dist);
							dist+=1;//should really be weighted
							p = prev;

						}
					}else{
					//This won't work (without modification) because we need to constrain these paths to ones of maximum real distance < min_dist.
						//Randomly re-weight the graph sometimes
						if(drand(rnd_seed)<opt_decide_graph_re_rnd){

							for(int i=0;i<outer->g.nodes;i++){
									 double w = drand(rnd_seed);
									 w-=0.5;
									 w*=w;
									 //printf("%f (%f),",w,rnd_seed);
									 rnd_path->setWeight(i,w);
								 }
						}
						 rnd_path->update();
						//derive a random path in the graph
						 p = j;
						 last = j;
						 assert( rnd_path->connected(p));
						while(!under->connected(p)){

							last=p;
							assert(p!=source);
							int prev = rnd_path->previous(p);
							p = prev;
							assert(p>=0);
						}

					}



					//ok, now pick some edge p->last that will connect p to last;
					assert(!under->connected(last));
					assert(under->connected(p));

					assert(over->connected(last));
					assert(over->connected(p));
					assert(p>-1);
					if(p>-1){
					Var v = outer->edges[p][last].v;
					if(outer->value(v)==l_Undef){
						return mkLit(v,false);
					}else{
						assert(outer->value(v)!=l_True);
					}
					}
					for(int k = 0;k<outer->antig.adjacency[p].size();k++){
						int to = outer->antig.adjacency[p][k].node;
						if (to==last){
							Var v =outer->edge_list[ outer->antig.adjacency[p][k].id].v;
							if(outer->value(v)==l_Undef){
								return mkLit(v,false);
							}else{
								assert(outer->value(v)!=l_True);
							}
						}
					}

				}
				}
			}else if(outer->value(l)==l_False){
				if(opt_decide_graph_neg){

					//assert(over->distance(j)<=min_dist);//else we would already be in conflict before this decision was attempted!


					if(over->distance(j)<=min_dist && under->distance(j)>min_dist){
						//then lets try to disconnect this node from source by walking back along the path in the over approx, and disabling the first unassigned edge we see.
						//(there must be at least one such edge, else the variable would be connected in the under approximation as well - in which case it would already have been propagated.
						int p = j;
						int last = j;
						 int dist = 0;
						// tmp_nodes.clear();
						while(under->distance(p)>min_dist-dist){
							//int d = under->distance(p);
							//int d_over = over->distance(p);
							last=p;
							assert(p!=source);
							int prev = over->previous(p);
							Var v = outer->edges[prev][p].v;
							if(outer->value(v)==l_Undef){
								//if(opt_use_random_path_for_decisions)
								//	tmp_nodes.push(v);
								//else
								return mkLit(v,true);
							}else{
								assert(outer->value(v)!=l_False);
							}
							assert(over->distance(p)<=min_dist-dist);
							dist+=1;//should really be weighted
							p = prev;

						}
						assert(opt_use_random_path_for_decisions);
						assert(tmp_nodes.size()>0);
						int i = irand(rnd_seed,tmp_nodes.size());
						Var v= tmp_nodes[i];
						return mkLit(v,true);

					}
				}
			}
		}
	}*/
	return lit_Undef;
};


