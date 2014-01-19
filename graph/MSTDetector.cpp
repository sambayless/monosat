/*
 * MSTDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "MSTDetector.h"
#include "GraphTheory.h"
#include <limits>
MSTDetector::MSTDetector(int _detectorID, GraphTheorySolver * _outer,  DynamicGraph<PositiveEdgeStatus> &_g,DynamicGraph<NegativeEdgeStatus> &_antig ,vec<int> & _edge_weights,double seed):
Detector(_detectorID),outer(_outer),g(_g),antig(_antig),rnd_seed(seed),edge_weights(_edge_weights),positive_reach_detector(NULL),negative_reach_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){



		positiveReachStatus = new MSTDetector::MSTStatus(*this,true);
		negativeReachStatus = new MSTDetector::MSTStatus(*this,false);
		positive_reach_detector = new Kruskal<MSTDetector::MSTStatus,PositiveEdgeStatus>(_g,*(positiveReachStatus),edge_weights,1);
		negative_reach_detector = new Kruskal<MSTDetector::MSTStatus,NegativeEdgeStatus>(_antig,*(negativeReachStatus),edge_weights,-1);



}

void MSTDetector::addWeightLit(Var weight_var,int min_weight){


	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});

	while(outer->S->nVars()<=weight_var)
		outer->S->newVar();

	Lit reachLit=mkLit(weight_var,false);
	bool found=false;
	for(int i = 0;i<weight_lits.size();i++){
		if(weight_lits[i].min_weight==min_weight){
			found=true;
			Lit r = weight_lits[i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->S->addClause(~r, reachLit);
			outer->S->addClause(r, ~reachLit);
		}
	}
	if(!found){
		weight_lits.push();
		weight_lits.last().l = reachLit;
		weight_lits.last().min_weight=min_weight;
		//weight_lit_map.insert(min_weight,weight_lits.size()-1);
	}



}


void MSTDetector::addTreeEdgeLit(int edge_id, Var reach_var){


	tree_edge_lits.growTo(outer->num_edges);

	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});

	while(outer->S->nVars()<=reach_var)
		outer->S->newVar();

	Lit reachLit=mkLit(reach_var,false);
	bool found=false;

	if(tree_edge_lits[edge_id].l!=lit_Undef){
		found=true;
		Lit r = tree_edge_lits[edge_id].l;
		//force equality between the new lit and the old reach lit, in the SAT solver
		outer->S->addClause(~r, reachLit);
		outer->S->addClause(r, ~reachLit);
	}else{
		tree_edge_lits.push();
		tree_edge_lits.last().l = reachLit;
		tree_edge_lits.last().edgeID=edge_id;

	}


}

void MSTDetector::MSTStatus::inMinimumSpanningTree(int edgeid, bool in_tree){
	if(edgeid<detector.tree_edge_lits.size()){
		Lit l = detector.tree_edge_lits[edgeid].l;
		if(l!=lit_Undef)
			detector.changed_edges.push({in_tree? l:~l,edgeid});
	}
}

void MSTDetector::MSTStatus::setMinimumSpanningTree(int weight){

			for(int i = 0;i<detector.weight_lits.size();i++){
				int min_weight =  detector.weight_lits[i].min_weight;

				Lit l = detector.weight_lits[i].l;
				if(l!=lit_Undef){

					assert(l!=lit_Undef);
					if(min_weight<weight && !polarity){
						lbool assign = detector.outer->S->value(l);
						if( assign!= l_False ){
							detector.changed_weights.push({~l,min_weight});
						}
					}else if(min_weight>=weight && polarity){
						lbool assign = detector.outer->S->value(l);
						if( assign!= l_True ){
							detector.changed_weights.push({l,min_weight});
						}
					}
				}
			}


}


void MSTDetector::buildMinWeightTooSmallReason(int weight,vec<Lit> & conflict){

			MinimumSpanningTree & d = *positive_reach_detector;


			double starttime = cpuTime();
			d.update();
			int minweight = d.weight();
			assert(d.weight()<=weight);
			//learn that at least one edge in the tree must be disabled (else, the minimum weight cannot be higher than the current weight)
			vec<int> & mst = d.getSpanningTree();
			for(int i = 0;i<mst.size();i++){
				int edgeid = mst[i];
				Var v = outer->edge_list[edgeid].v;
				assert(outer->S->value(v)==l_True);
				conflict.push(mkLit(v,true));
			}

			outer->num_learnt_paths++;
			outer->learnt_path_clause_length+= (conflict.size()-1);
			double elapsed = cpuTime()-starttime;
			outer->pathtime+=elapsed;

		}


		bool MSTDetector::walkback(int min_weight, int from, int to){
			int u = from;
			while(u!=to){
				int p = negative_reach_detector->getParent(u);
				int edgeid = negative_reach_detector->getParentEdge(u);
				int weight = edge_weights[edgeid];
				if(weight>min_weight){
					return true;
				}
				u=p;
			}
			return false;
		}

		void MSTDetector::TarjanOLCA(int node,vec<Lit> & conflict){
			ancestors[node]=node;
			for(int i = 0;i<antig.adjacency_undirected[node].size();i++){
				int edgeid = antig.adjacency_undirected[node][i].id;
				if(negative_reach_detector->edgeInTree(edgeid)){
					assert(antig.edgeEnabled(edgeid));
					int v = antig.adjacency_undirected[node][i].node;
					if(negative_reach_detector->getParent(v)==node){
						//u is a child of node in the minimum spanning tree
						TarjanOLCA(v,conflict);
						sets.Union(node,v);
						int set = sets.FindSet(node);
						ancestors[set]=node;

					}
				}
			}
			black[node]=true;
			//now visit all _disabled_ edges of node
			for(int i = 0;i<antig.adjacency_undirected[node].size();i++){
				int edgeid = antig.adjacency_undirected[node][i].id;
				if(!antig.edgeEnabled(edgeid)){
					//this is a disabled edge
					int v = g.adjacency_undirected[node][i].node;
					if(black[v]){
						int set = sets.FindSet(v);
						int lowest_common_ancestor = ancestors[set];

						//ok, now walk back from u and v in the minimum spanning tree until either the lowest common ancestor is seen, or an edge larger than the weight of this disabled edge is found
						int weight =edge_weights[edgeid];
						bool any_larger_weights = walkback(weight,node,lowest_common_ancestor) ||  walkback(weight,v,lowest_common_ancestor) ;
						//if any larger edge was found in either path from u or v to their common ancestor in the minimum spanning tree, then enabling this edge
						//would have replaced that larger edge in the minimum spanning tree, resulting in a smaller minimum spanning tree.
						if(any_larger_weights){
							Var e =outer->edge_list[edgeid].v;
							assert(outer->S->value(e)==l_False);
							conflict.push(mkLit(e,false));
						}
					}
				}
			}
		}

		void MSTDetector::buildMinWeightTooLargeReason(int weight,vec<Lit> & conflict){



			static int it = 0;
			++it;

			//drawFull( non_reach_detectors[detector]->getSource(),u);
			//assert(outer->dbg_distance( source,u));
			double starttime = cpuTime();
			int INF=std::numeric_limits<int>::max();
			negative_reach_detector->update();
			int mstweight = negative_reach_detector->weight();

			if(negative_reach_detector->weight()==INF){
				//IF the mst is disconnected, then we define it's weight to be infinite. In this case, the reason is a separating cut between any two disconnected components.
				//we can find one of these by identifying any two roots

				int min_conflict = INF;

				//walk back down from the each root to find a separating cut of disabled edge.
				//return the smallest such cut.



				assert(conflict.size()==1);
				for(int i = 0;i<negative_reach_detector->numComponents();i++){
					int root = negative_reach_detector->getRoot(i);
					tmp_conflict.clear();
					visit.clear();
					//ok, now traverse the connected component in the tree below this root.
					visit.push(root);
					seen.clear();
					seen.growTo(g.nodes);
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

								Var e =outer->edge_list[edgeid].v;
								assert(outer->S->value(e)==l_False);
								tmp_conflict.push(mkLit(e,false));
							}
						}
					}

						if(tmp_conflict.size()<min_conflict){
							min_conflict = tmp_conflict.size();
							conflict.shrink(conflict.size()-1);//keep only the first conflict element, which is the constraint lit
							assert(conflict.size()==1);
							for(int i = 0;i<tmp_conflict.size();i++){
								conflict.push(tmp_conflict[i]);
							}
						}
					if(!opt_mst_min_cut){
						break;//stop looking as soon as we find a cut, instead of trying to find the smallest one.
					}
				}



				return;
			}

			ancestors.clear();
			ancestors.growTo(g.nodes,-1);
			black.clear();
			black.growTo(g.nodes);
			sets.Reset();
			sets.AddElements(g.nodes);

			//Noah's algorithm for finding a minimal neccesary set of edges to add to decrease the weight of the tree
			//all we are doing here is visiting each disabled edge, then examining the cycle that would be created in the minimum spanning tree if we were to
			//add that edge to the tree. If any edge in that cycle is larger than the disabled edge, then enabling that edge would result in a smaller minimum spanning tree.
			//so we have to add that edge to the conflict as part of the reason that the minimum spanning tree is larger than it should be.

			//conceptually, thats really simple. The problem is that visiting that cycle requires computing the lowest common ancestor of the two nodes on either side of the disabled edge,
			//and to make that efficient we are going to use Tarjan's offline lca algorithm. This results in the convoluted code below.

			int root = 0;
			while(negative_reach_detector->getParent(root)!=-1){
				root = negative_reach_detector->getParent(root);
			}
			TarjanOLCA(root,conflict);
#ifndef NDEBUG
			for(int i = 0;i<black.size();i++)
				assert(black[i]);
#endif


			 outer->num_learnt_cuts++;
			 outer->learnt_cut_clause_length+= (conflict.size()-1);

			double elapsed = cpuTime()-starttime;
			 outer->mctime+=elapsed;


		}


		int MSTDetector::walkback_edge(int min_weight, int edge_id, int from, int to, bool & found){
			int u = from;
			int maxweight = 0;
			while(u!=to){
				int p = negative_reach_detector->getParent(u);
				int edgeid = negative_reach_detector->getParentEdge(u);
				int weight = edge_weights[edgeid];
				if(edgeid == edge_id){
					assert(!found);//can't enounter the edge twice while traversing the cycle
					found=true;
				}
				if(weight>maxweight){
					maxweight=weight;
				}
				u=p;
			}
			return maxweight;
		}


		void MSTDetector::TarjanOLCA_edge(int node,int check_edgeid,int lowest_endpoint,vec<Lit> & conflict){
					ancestors[node]=node;
					for(int i = 0;i<antig.adjacency_undirected[node].size();i++){
						int edgeid = antig.adjacency_undirected[node][i].id;
						if(negative_reach_detector->edgeInTree(edgeid)){
							int v = antig.adjacency_undirected[node][i].node;
							if(negative_reach_detector->getParent(v)==node){
								//u is a child of node in the minimum spanning tree
								TarjanOLCA_edge(v,edgeid,lowest_endpoint,conflict);
								sets.Union(node,v);
								int set = sets.FindSet(node);
								ancestors[set]=node;

							}
						}
					}
					black[node]=true;
					//now visit all _disabled_ edges of node
					for(int i = 0;i<g.adjacency_undirected[node].size();i++){
						int edgeid = g.adjacency_undirected[node][i].id;
						if(!antig.edgeEnabled(edgeid)){
							//this is a disabled edge
							int v = g.adjacency_undirected[node][i].node;
							if(black[v]){
								int set = sets.FindSet(v);
								int lowest_common_ancestor = ancestors[set];

								//ok, now walk back from u and v in the minimum spanning tree until either the lowest common ancestor is seen, or an edge larger than the weight of this disabled edge is found
								int weight =edge_weights[edgeid];
								bool found=false;

								int largest_weight = walkback_edge(weight,check_edgeid,node,lowest_common_ancestor, found);

								int largest_weight_other = walkback_edge(weight,check_edgeid,v,lowest_common_ancestor,found) ;
								if(found   && (largest_weight>weight || largest_weight_other >weight )){
									//if any edge larger than the disabled edge's weight was found in either path from u or v to their common ancestor in the minimum spanning tree, then enabling this edge
									//would lead  to replacing that larger edge in the minimum spanning tree, resulting in a smaller minimum spanning tree, possibly resulting in the edge we really care about being removed from the tree.

									Var e =outer->edge_list[edgeid].v;
									assert(outer->S->value(e)==l_False);
									conflict.push(mkLit(e,false));

								}
							}
						}
					}
				}

				void MSTDetector::buildNonEdgeReason(int edgeid,vec<Lit> & conflict){
					static int it = 0;
					++it;
					ancestors.clear();
					ancestors.growTo(g.nodes,-1);

					black.clear();
					black.growTo(g.nodes);
					//drawFull( non_reach_detectors[detector]->getSource(),u);
					//assert(outer->dbg_distance( source,u));
					double starttime = cpuTime();
					int INF=std::numeric_limits<int>::max();
					negative_reach_detector->update();

					sets.Reset();
					sets.AddElements(g.nodes);


					//Noah's algorithm for finding a minimal neccesary set of edges to add to decrease the weight of the tree
					//all we are doing here is visiting each disabled edge, then examining the cycle that would be created in the minimum spanning tree if we were to
					//add that edge to the tree. If any edge in that cycle is larger than the disabled edge, then enabling that edge would result in a smaller minimum spanning tree.
					//so we have to add that edge to the conflict as part of the reason that the minimum spanning tree is larger than it should be.

					//conceptually, thats really simple. The problem is that visiting that cycle requires computing the lowest common ancestor of the two nodes on either side of the disabled edge,
					//and to make that efficient we are going to use Tarjan's offline lca algorithm. This results in the convoluted code below.
					int u = g.all_edges[edgeid].from;
					int v = g.all_edges[edgeid].to;
					int lower_endpoint = u;
					if(negative_reach_detector->getParent(v)==u){
						lower_endpoint=v;
					}else{
						assert(negative_reach_detector->getParent(u)==v);
					}

					int root = 0;
					while(negative_reach_detector->getParent(root)!=-1){
						root = negative_reach_detector->getParent(root);
					}
					TarjanOLCA_edge(root,edgeid, lower_endpoint,conflict);//run tarjan's off-line lowest common ancestor query from node 0, arbitrarily.



					 outer->num_learnt_cuts++;
					 outer->learnt_cut_clause_length+= (conflict.size()-1);

					double elapsed = cpuTime()-starttime;
					 outer->mctime+=elapsed;


				}
				void MSTDetector::buildEdgeReason(int edgeid,vec<Lit> & conflict){
					//the reason that an edge is NOT in the minimum spanning tree is the paths to lca from either edge. So long as all those edges are in the tree, each of which is <= weight to this edge,
					//this edge cannot be in the mst.
					seen.clear();
					seen.growTo(g.nodes);
					int u = g.all_edges[edgeid].from;
					int v = g.all_edges[edgeid].to;
					positive_reach_detector->update();
					assert(!g.edgeEnabled(edgeid));
					int r = u;
					while(r != -1 ){
						seen[r]=true;
						r= positive_reach_detector->getParent(r);

					}
					r=v;
					while(! seen[r]){
						r = positive_reach_detector->getParent(r);
						assert(r!=-1);
					}
					int lca = r;
					r=u;
					while(r!=lca){
						assert(seen[u]);
						seen[r]=false;
						int p =  positive_reach_detector->getParent(r);
						int edge = positive_reach_detector->getParentEdge(r);
						assert(g.edgeEnabled(edge));
						Var v = outer->edge_list[edge].v;
						assert(outer->S->value(v)==l_True);
						conflict.push(mkLit(v,true));
						r= p;
					}

					r=v;
					while(r!=lca){
						assert(!seen[u]);
						int p =  positive_reach_detector->getParent(r);
						int edge = positive_reach_detector->getParentEdge(r);
						assert(g.edgeEnabled(edge));
						Var v = outer->edge_list[edge].v;
						assert(outer->S->value(v)==l_True);
						conflict.push(mkLit(v,true));
						r= p;
					}

					r = lca;
					while(r!=-1){
						assert(seen[r] || r==lca);
						seen[r]=false;
						r = positive_reach_detector->getParent(r);
					}
				}
		/**
		 * Explain why an edge was forced (to true).
		 * The reason is that _IF_ that edge is false, THEN there is a cut of disabled edges between source and target
		 * So, create the graph that has that edge (temporarily) assigned false, and find a min-cut in it...
		 */
	/*	void MSTDetector::buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict){
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




					 outer->num_learnt_cuts++;
					 outer->learnt_cut_clause_length+= (conflict.size()-1);

					double elapsed = cpuTime()-starttime;
					 outer->mctime+=elapsed;

			#ifdef DEBUG_GRAPH
					 assert(outer->dbg_clause(conflict));
			#endif
				}
*/
		void MSTDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){

				if(marker==reach_marker){
					reason.push(p);

					Var v = var(p);
					int weight=-1;
					//could swap this out for a map if there are lots of lits..
					for(int i = 0;i<weight_lits.size();i++){
						if(var(weight_lits[i].l)==v){
							weight=weight_lits[i].min_weight;
							break;
						}
					}
					assert(weight>=0);
					buildMinWeightTooSmallReason(weight,reason);

					//double elapsed = cpuTime()-startpathtime;
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
					for(int i = 0;i<weight_lits.size();i++){
						if(var(weight_lits[i].l)==v){
							weight=weight_lits[i].min_weight;
							break;
						}
					}
					assert(weight>=0);

					buildMinWeightTooLargeReason(weight,reason);

				}else if (marker==forced_reach_marker){
					//not implemented yet
					/*Var v = var(p);
					//The forced variable is an EDGE that was forced.
					int forced_edge_id =v- outer->min_edge_var;
					//The corresponding node that is the reason it was forced
					int reach_node=force_reason[forced_edge_id];
					 buildForcedEdgeReason( reach_node,  forced_edge_id, reason);*/
				}else{
					assert(false);
				}
		}

		bool MSTDetector::propagate(vec<Assignment> & trail,vec<Lit> & conflict){


		double startdreachtime = cpuTime();
		changed_edges.clear();
		changed_weights.clear();
		positive_reach_detector->update();
		double reachUpdateElapsed = cpuTime()-startdreachtime;
		outer->reachupdatetime+=reachUpdateElapsed;

		double startunreachtime = cpuTime();
		negative_reach_detector->update();
		double unreachUpdateElapsed = cpuTime()-startunreachtime;
		outer->unreachupdatetime+=unreachUpdateElapsed;

		for(int j = 0;j<changed_weights.size();j++){
			Lit l = changed_weights[j].l;
			int weight = changed_weights[j].weight;
			bool reach = !sign(l);
			if(outer->S->value(l)==l_True){
				//do nothing
			}else if(outer->S->value(l)==l_Undef){
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
					buildMinWeightTooSmallReason(weight,conflict);
				//add it to s
				//return it as a conflict

				}else{

					buildMinWeightTooLargeReason(weight,conflict);

				}

				return false;
			}else{
				int  a=1;
			}

		}

		for(int j = 0;j<changed_edges.size();j++){
					Lit l = changed_edges[j].l;
					int edge = changed_edges[j].edgeID;
					bool reach = !sign(l);
					if(outer->S->value(l)==l_True){
						//do nothing
					}else if(outer->S->value(l)==l_Undef){
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
							buildEdgeReason(edge,conflict);
						//add it to s
						//return it as a conflict

						}else{
							//The reason is a cut separating s from t
							buildNonEdgeReason(edge,conflict);

						}

						return false;
					}else{
						int  a=1;
					}

				}
			return true;
		}

bool MSTDetector::checkSatisfied(){


					for(int k = 0;k<weight_lits.size();k++){
						Lit l = weight_lits[k].l;
						int dist = weight_lits[k].min_weight;

						if(l!=lit_Undef){


							if(outer->S->value(l)==l_True){
								if(positive_reach_detector->weight()>dist){
									return false;
								}
							}else if (outer->S->value(l)==l_False){
								if( negative_reach_detector->weight()<=dist){
									return false;
								}
							}else{
								if(positive_reach_detector->weight()<=dist){
									return false;
								}
								if(!negative_reach_detector->weight()>dist){
									return false;
								}
							}
						}
					}

	return true;
}
Lit MSTDetector::decide(){
	/*MSTDetector *r =this;
	MinimumSpanningTree<MSTDetector::MSTStatus,NegativeEdgeStatus> * over = (MinimumSpanningTree<MSTDetector::MSTStatus,NegativeEdgeStatus>*) r->negative_reach_detector;

	MinimumSpanningTree<MSTDetector::MSTStatus,PositiveEdgeStatus> * under = (MinimumSpanningTree<MSTDetector::MSTStatus,PositiveEdgeStatus>*) r->positive_reach_detector;

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
			if(outer->S->value(l)==l_True){
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
					if(outer->S->value(v)==l_Undef){
						return mkLit(v,false);
					}else{
						assert(outer->S->value(v)!=l_True);
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

				}
				}
			}else if(outer->S->value(l)==l_False){
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
							if(outer->S->value(v)==l_Undef){
								//if(opt_use_random_path_for_decisions)
								//	tmp_nodes.push(v);
								//else
								return mkLit(v,true);
							}else{
								assert(outer->S->value(v)!=l_False);
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


