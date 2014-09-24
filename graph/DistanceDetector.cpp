/*
 * DistanceDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */




#include <core/Config.h>
#include <dgl/RamalReps.h>
#include <dgl/UnweightedRamalReps.h>
#include <graph/DistanceDetector.h>
#include <graph/GraphTheory.h>
#include <mtl/Rnd.h>
#include <mtl/Vec.h>
#include <utils/Options.h>
//#include "dgl/UnweightedDistance.h"
#include <cassert>
#include <cstdlib>
#include <iostream>


using namespace Minisat;
template<typename Weight>
DistanceDetector<Weight>::DistanceDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,std::vector<Weight> & weights,  DynamicGraph &_g,DynamicGraph &_antig, int from, double seed):
Detector(_detectorID),outer(_outer),weights(weights),g(_g),antig(_antig),source(from),rnd_seed(seed),positive_reach_detector(NULL),negative_reach_detector(NULL),positive_path_detector(NULL),positiveReachStatus(NULL),negativeReachStatus(NULL){
	max_unweighted_distance=0;
	rnd_path=NULL;

	constraintsBuilt=-1;
	first_reach_var = var_Undef;
	stats_pure_skipped=0;
	 if(distalg ==DistAlg::ALG_SAT){
		 positiveReachStatus=nullptr;
		 negativeReachStatus=nullptr;
		 positive_reach_detector=nullptr;
		 negative_reach_detector=nullptr;
		 positive_path_detector=nullptr;
		 reach_marker=CRef_Undef;
		 non_reach_marker=CRef_Undef;
		 forced_reach_marker=CRef_Undef;
		 weighted_reach_marker=CRef_Undef;
		 weighted_non_reach_marker=CRef_Undef;
		 //we are just going to directly enforce these constraints in the original SAT solver, so _nothing_ will end up happening in this detector (except for creating the clauses needed to enforce these constraints).

		 return;
	 }

	 if(opt_use_random_path_for_decisions){
		 rnd_weight.clear();
		 rnd_path = new WeightedDijkstra< double>(from,_antig,rnd_weight);
		 for(int i=0;i<outer->edge_list.size();i++){
			 double w = drand(rnd_seed);

			 rnd_weight.push_back(w);
		 }

	 }
	/* if(opt_use_optimal_path_for_decisions){
		 opt_path = new WeightedDijkstra< OptimalWeightEdgeStatus >(from,_antig,opt_weight);
	 }*/
	 positiveReachStatus = new DistanceDetector<Weight>::ReachStatus(*this,true);
	 negativeReachStatus = new DistanceDetector<Weight>::ReachStatus(*this,false);

	//select the unweighted distance detectors
	if(distalg==DistAlg::ALG_DISTANCE){
			positive_reach_detector = new UnweightedBFS<DistanceDetector<Weight>::ReachStatus>(from,_g,*(positiveReachStatus),0);
			negative_reach_detector = new UnweightedBFS<DistanceDetector<Weight>::ReachStatus>(from,_antig,*(negativeReachStatus),0);
			positive_path_detector = positive_reach_detector;
		/*	if(opt_conflict_shortest_path)
			reach_detectors.last()->positive_dist_detector = new Dijkstra<PositiveEdgeStatus>(from,g);*/
	}else if (distalg==DistAlg::ALG_RAMAL_REPS){

		positive_reach_detector = new UnweightedRamalReps<typename DistanceDetector<Weight>::ReachStatus>(from,_g,*(positiveReachStatus),0);
		negative_reach_detector = new UnweightedRamalReps<typename DistanceDetector<Weight>::ReachStatus>(from,_antig,*(negativeReachStatus),0);
		positive_path_detector =  new UnweightedBFS<Reach::NullStatus>(from,_g,Reach::nullStatus,0);

	}else{
		positive_reach_detector = new UnweightedDijkstra<typename DistanceDetector<Weight>::ReachStatus>(from,_g,*positiveReachStatus,0);
		negative_reach_detector = new UnweightedDijkstra<typename DistanceDetector<Weight>::ReachStatus>(from,_antig,*negativeReachStatus,0);
		positive_path_detector = positive_reach_detector;
		//reach_detectors.last()->positive_dist_detector = new Dijkstra(from,g);
	}

	//select the _weighted_ distance detectors
	 positiveDistanceStatus = new DistanceDetector<Weight>::DistanceStatus(*this,true);
	 negativeDistanceStatus = new DistanceDetector<Weight>::DistanceStatus(*this,false);


	 if (distalg==DistAlg::ALG_RAMAL_REPS){
		positive_weighted_distance_detector = new RamalReps<Weight,typename DistanceDetector<Weight>::DistanceStatus>(from,_g,weights,*(positiveDistanceStatus),0);
		negative_weighted_distance_detector = new RamalReps<Weight,typename DistanceDetector<Weight>::DistanceStatus>(from,_antig,weights,*(negativeDistanceStatus),0);
		 positive_weighted_path_detector =  new Dijkstra<Weight>(from,_g,weights);
	 }else{
		 positive_weighted_distance_detector = new Dijkstra<Weight,typename DistanceDetector<Weight>::DistanceStatus>(from,_g,weights,*positiveDistanceStatus,0);
		 negative_weighted_distance_detector = new Dijkstra<Weight,typename DistanceDetector<Weight>::DistanceStatus>(from,_antig,weights,*negativeDistanceStatus,0);
		 positive_weighted_path_detector = positive_weighted_distance_detector;
	 }


	reach_marker=outer->newReasonMarker(getID());
	non_reach_marker=outer->newReasonMarker(getID());
	forced_reach_marker=outer->newReasonMarker(getID());
	weighted_reach_marker=outer->newReasonMarker(getID());
	weighted_non_reach_marker=outer->newReasonMarker(getID());
}

template<typename Weight>
void DistanceDetector<Weight>::buildSATConstraints(int within_steps){
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
		unweighted_dist_lits.push();
		Lit True = mkLit(outer->newVar());
		outer->addClause(True);
		assert(outer->value(True)==l_True);
		Lit False = ~True;
		for(int i = 0;i<g.nodes();i++){
			unweighted_sat_lits[0].push(False);
		}
		unweighted_sat_lits[0][source]=True;
	}

	vec<Lit> reaches;

	//bellman-ford:
	for (int i = constraintsBuilt;i<within_steps;i++){
		unweighted_sat_lits.last().copyTo(reaches);
		unweighted_sat_lits.push();
		reaches.copyTo(unweighted_sat_lits.last());
		assert(outer->value( reaches[source])==l_True);
		//For each edge:
		for(int j = 0;j<g.nodes();j++){
			Lit r_cur = reaches[j];

			for(Edge & e: outer->inv_adj[j]){
					if(outer->value(unweighted_sat_lits.last()[e.to])==l_True){
						//do nothing
					}else if (outer->value(reaches[e.from])==l_False){
						//do nothing
					}else{
						Lit l = mkLit(e.v,false);
						Lit r = mkLit( outer->newVar(), false);

						c.clear();
						c.push(~r);c.push(reaches[e.to]);c.push(l); //r -> (e.l or reaches[e.to])
						outer->addClause(c);
						c.clear();
						c.push(~r);c.push(reaches[e.to]);c.push(reaches[e.from]); //r -> (reaches[e.from]) or reaches[e.to])
						outer->addClause(c);
						c.clear();
						c.push(r);c.push(~reaches[e.to]); //~r -> ~reaches[e.to]
						outer->addClause(c);
						c.clear();
						c.push(r);c.push(~reaches[e.from]);c.push(~l); //~r -> (~reaches[e.from] or ~e.l)
						outer->addClause(c);
						r_cur=r;

					}
			}
			unweighted_sat_lits.last()[j]=r_cur   ;//reaches[e.to] == (var & reaches[e.from])| reaches[e.to];
		}

	}

	constraintsBuilt=within_steps;
}


template<typename Weight>
void DistanceDetector<Weight>::addWeightedShortestPathLit(int from, int to, Var outer_reach_var,Weight within_distance){
	g.invalidate();
	antig.invalidate();
	Var reach_var= outer->newVar(outer_reach_var,getID());
	assert(from==source);
	weighted_dist_lits.push(WeightedDistLit{mkLit(reach_var),to,within_distance});
	//sort(weighted_dist_lits);
}

template<typename Weight>
void DistanceDetector<Weight>::addUnweightedShortestPathLit(int from, int to, Var outer_reach_var,int within_steps){
	g.invalidate();
	antig.invalidate();
	Var reach_var= outer->newVar(outer_reach_var,getID());

	if(first_reach_var==var_Undef){
		first_reach_var=reach_var;
	}else{
		assert(reach_var>first_reach_var);
	}
	assert(from==source);
	assert(within_steps>=0);
	if(within_steps>g.nodes())
		within_steps=g.nodes();
	if(within_steps<0)
		within_steps=g.nodes();

	if(within_steps>=max_unweighted_distance){
		max_unweighted_distance=within_steps;
		if(opt_compute_max_distance){
			positive_reach_detector->setMaxDistance(max_unweighted_distance);
			negative_reach_detector->setMaxDistance(max_unweighted_distance);
			//positive_path_detector->setMaxDistance(max_distance);
		}
	}

	Lit reachLit=mkLit(reach_var,false);

	 if(distalg ==DistAlg::ALG_SAT){
		 buildSATConstraints(within_steps);

		if(unweighted_sat_lits.size() >=within_steps){
			//then we are using the sat encoding and can directly look up its corresponding distance lit

			Lit r = unweighted_sat_lits[within_steps][to];
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->makeEqual(r,reachLit);
		}
		return;
	 }

		while( unweighted_dist_lits.size()<=to)
			unweighted_dist_lits.push();

	bool found=false;
	for(int i = 0;i<unweighted_dist_lits[to].size();i++){
		if(unweighted_dist_lits[to][i].min_unweighted_distance==within_steps){
			found=true;
			Lit r = unweighted_dist_lits[to][i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->makeEqual(r,reachLit);
			/*outer->S->addClause(~r, reachLit);
			outer->S->addClause(r, ~reachLit);*/
		}
	}



	if(!found){
		unweighted_dist_lits[to].push();
		unweighted_dist_lits[to].last().l = reachLit;
		unweighted_dist_lits[to].last().min_unweighted_distance=within_steps;
		while(reach_lit_map.size()<= reach_var- first_reach_var ){
			reach_lit_map.push(-1);
		}

		reach_lit_map[reach_var-first_reach_var]=to;
	}



}
template<typename Weight>
void DistanceDetector<Weight>::ReachStatus::setReachable(int u, bool reachable){
/*	if(polarity==reachable && u<detector.reach_lits.size()){
		Lit l = detector.reach_lits[u];
		if(l!=lit_Undef){
			lbool assign = detector.outer->value(l);
			if(assign!= (reachable? l_True:l_False )){
				detector.changed.push({reachable? l:~l,u});
			}
		}
	}*/
}
template<typename Weight>
void DistanceDetector<Weight>::ReachStatus::setMininumDistance(int u, bool reachable, int  distance){
	assert(reachable ==(distance<detector.outer->g.nodes()));
	if(distance<=detector.outer->g.nodes()){
		setReachable(u,reachable);
	}

		if(u<detector.unweighted_dist_lits.size()){
			assert(distance>=0);

			for(int i = 0;i<detector.unweighted_dist_lits[u].size();i++){
				int& min_distance =  detector.unweighted_dist_lits[u][i].min_unweighted_distance;

				Lit l = detector.unweighted_dist_lits[u][i].l;
				if(l!=lit_Undef){

					assert(l!=lit_Undef);
					if(!polarity && (!reachable || (distance>min_distance))){
						lbool assign = detector.outer->value(l);
						if( assign!= l_False ){
							detector.changed.push({~l,u});
						}
					}else if(polarity && reachable && (distance<=min_distance)){
						lbool assign = detector.outer->value(l);
						if( assign!= l_True ){
							detector.changed.push({l,u});
						}
					}
				}
			}
		}

}

template<typename Weight>
void DistanceDetector<Weight>::DistanceStatus::setReachable(int u, bool reachable){

}
template<typename Weight>
void DistanceDetector<Weight>::DistanceStatus::setMininumDistance(int u, bool reachable, Weight & distance){




}

template<typename Weight>
void DistanceDetector<Weight>::buildReachReason(int node,vec<Lit> & conflict){
			//drawFull();
			Reach & d = *positive_path_detector;
			stats_unweighted_leq_reasons++;

			double starttime = rtime(2);
			d.update();
			assert(outer->dbg_reachable(d.getSource(),node));
/*
			if(!outer->dbg_reachable(d.getSource(),node)){
				outer->drawFull();

				d.drawFull();

				assert(false);
			}*/
			assert(positive_reach_detector->connected(node));
			/*if(!d.connected_unchecked(node)){
				exit(3);
			}*/
			assert(d.connected_unchecked(node));
			//if(opt_learn_reaches ==0 || opt_learn_reaches==2)
			{
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
			}/*else{
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
					if( u< reach_lits.size() && reach_lits[u]!=lit_Undef && outer->value(reach_lits[u])==l_True && outer->S->level(var(reach_lits[u]))< outer->S->decisionLevel()){
						//A potential (fixed) problem with the above: reach lit can be false, but have been assigned after r in the trail, messing up clause learning if this is a reason clause...
						//This is avoided by ensuring that L is lower level than the conflict.
						Lit l =reach_lits[u];
						assert(outer->value(l)==l_True);
						conflict.push(~l);
						break;
					}
				}

			}*/
			outer->num_learnt_paths++;
			outer->learnt_path_clause_length+= (conflict.size()-1);
			double elapsed = rtime(2)-starttime;
			outer->pathtime+=elapsed;

		}
template<typename Weight>
void DistanceDetector<Weight>::buildNonReachReason(int node,vec<Lit> & conflict){
	static int it = 0;
	stats_unweighted_gt_reasons++;
	++it;
	int u = node;
	//drawFull( non_reach_detectors[detector]->getSource(),u);
	//assert(outer->dbg_distance( source,u));
	double starttime = rtime(2);
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
				//Ok, then add all its incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
				for(int i = 0;i<outer->inv_adj[u].size();i++){
					int v = outer->inv_adj[u][i].v;
					int from = outer->inv_adj[u][i].from;
					assert(outer->inv_adj[u][i].to==u);
					//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
					int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;

					if(from==u){
						assert(outer->edge_list[edge_num].to == u);
						assert(outer->edge_list[edge_num].from == u);
						continue;//Self loops are allowed, but just make sure nothing got flipped around...
					}
					assert(from!=u);

					if(outer->value(v)==l_False){
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

	double elapsed = rtime(2)-starttime;
	 outer->mctime+=elapsed;

}

template<typename Weight>
void DistanceDetector<Weight>::printSolution(){

		 vec<bool> to_show;
		 to_show.growTo(g.nodes());
		 for(auto & w:weighted_dist_lits){
			 to_show[w.u]=true;
		 }
		 positive_weighted_path_detector->update();
		 for(int to = 0;to<g.nodes();to++){
			 if(!to_show[to])
				 continue;

			Distance<Weight> & d = *positive_weighted_path_detector;
			d.update();
			Weight & actual_dist = d.distance(to);
			std::cout<<"Shortest Weighted Path " << source <<"->"<<to<<" is " << actual_dist<<": ";
			 //std::cout<< "Weighted Distance Constraint " <<dimacs(w.l) << " (" << source <<"->" << w.u << ") <=" << w.min_distance ;


			 if(d.connected(to) ){
				 vec<int> path;
				int u = to;
				path.push(u);
				int p;
				while(( p = d.previous(u)) != -1){
					Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
					path.push(p);
					u = p;
				}


				for(int i = path.size()-1;i>=0;i--){
					std::cout<<path[i] <<",";
				}
				std::cout<<'\n';
			 }else{
				 std::cout<<": FALSE\n";
			 }
		 }



}

template<typename Weight>
void DistanceDetector<Weight>::buildDistanceLEQReason(int to,Weight & min_distance,vec<Lit> & conflict){
	stats_distance_leq_reasons++;
	Distance<Weight> & d = *positive_weighted_path_detector;
	positive_weighted_distance_detector->update();
	Weight & actual_dist = positive_weighted_distance_detector->distance(to);
	double starttime = rtime(2);
	d.update();
	assert(outer->dbg_reachable(d.getSource(),to));
	assert(positive_reach_detector->connected(to));
	assert(positive_weighted_distance_detector->distance(to)<=min_distance && positive_weighted_distance_detector->distance(to)!=positive_weighted_distance_detector->unreachable());
	assert(d.connected_unchecked(to));
	if(!d.connected(to) || d.distance(to)>min_distance){
		fprintf(stderr,"Error in shortest path detector, aborting\n");
		exit(4);
	}
	//the reason that the distance is less than or equal to min_distance is because the shortest path is less than this weight
	{
		int u = to;
		int p;
		while(( p = d.previous(u)) != -1){
			Edge & edg = outer->edge_list[d.incomingEdge(u)]; //outer->edges[p][u];
			Var e =edg.v;
			lbool val = outer->value(e);
			assert(outer->value(e)==l_True);
			conflict.push(mkLit(e, true));
			u = p;
		}
	}
	outer->num_learnt_paths++;
	outer->learnt_path_clause_length+= (conflict.size()-1);
	double elapsed = rtime(2)-starttime;
	outer->pathtime+=elapsed;

}

template<typename Weight>
void DistanceDetector<Weight>::buildDistanceGTReason(int to,Weight & min_distance,vec<Lit> & conflict){
	static int it = 0;
	stats_distance_gt_reasons++;
	++it;
	int u = to;

	Weight & actual_dist = negative_weighted_distance_detector->distance(to);
	/*for(auto & e:antig.all_edges){
		if(antig.edgeEnabled(e.id)){
			printf("nxg.add_edge(%d,%d,weight=",e.from,e.to);
			std::cout<<outer->getWeight(e.id)<<")\n";
		}
	}*/
	bool connected = negative_weighted_distance_detector->connected(to);
#ifndef NDEBUG
	Dijkstra<Weight> d(source,antig,weights);
	Weight & expected = d.distance(to);
	assert(expected==actual_dist);
	assert(!negative_weighted_distance_detector->connected(to) || actual_dist>min_distance);
#endif
	double starttime = rtime(2);

	{

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
			//add all of this node's incoming disabled edges to the cut, and visit any unseen, non-disabled incoming.edges()
			for(int i = 0;i<outer->inv_adj[u].size();i++){
				int v = outer->inv_adj[u][i].v;
				int from = outer->inv_adj[u][i].from;
				assert(outer->inv_adj[u][i].to==u);
				//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
				int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;

				if(from==u){
					assert(outer->edge_list[edge_num].to == u);
					assert(outer->edge_list[edge_num].from == u);
					continue;//Self loops are allowed, but just make sure nothing got flipped around...
				}
				assert(from!=u);

				if(outer->value(v)==l_False){
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

	double elapsed = rtime(2)-starttime;
	 outer->mctime+=elapsed;


}

		/**
		 * Explain why an edge was forced (to true).
		 * The reason is that _IF_ that edge is false, THEN there is a cut of disabled edges between source and target
		 * So, create the graph that has that edge (temporarily) assigned false, and find a min-cut in it...
		 */
template<typename Weight>
		void DistanceDetector<Weight>::buildForcedEdgeReason(int reach_node, int forced_edge_id,vec<Lit> & conflict){
					static int it = 0;
					++it;

					assert(outer->value(outer->edge_list[forced_edge_id].v)==l_True);
					Lit edgeLit =mkLit( outer->edge_list[forced_edge_id].v,false);

					conflict.push(edgeLit);

					int forced_edge_from = outer->edge_list[forced_edge_id].from;
					int forced_edge_to= outer->edge_list[forced_edge_id].to;


					int u = reach_node;
					//drawFull( non_reach_detectors[detector]->getSource(),u);
					assert(outer->dbg_notreachable( source,u));
					double starttime = rtime(2);
					outer->cutGraph.clearHistory();
					outer->stats_mc_calls++;




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
					for(int i = 0;i<outer->inv_adj[u].size();i++){
						int v = outer->inv_adj[u][i].v;
						int from = outer->inv_adj[u][i].from;
						assert(from!=u);
						assert(outer->inv_adj[u][i].to==u);
						//Note: the variable has to not only be assigned false, but assigned false earlier in the trail than the reach variable...
						int edge_num =outer->getEdgeID(v);// v-outer->min_edge_var;

						if(edge_num == forced_edge_id || outer->value(v)==l_False){
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

					double elapsed = rtime(2)-starttime;
					 outer->mctime+=elapsed;


				}
template<typename Weight>
		void DistanceDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker){


				if(marker==reach_marker){
					reason.push(p);
				//	double startpathtime = rtime(2);

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


					//double elapsed = rtime(2)-startpathtime;
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
template<typename Weight>
		bool DistanceDetector<Weight>::propagate(vec<Lit> & conflict){
			if(!positive_reach_detector)
				return true;

			static int iter = 0;
			if(++iter==1624){//18303
				int a=1;
			}

		//printf("iter %d\n",iter);

		getChanged().clear();
		if(!opt_detect_pure_theory_lits || unassigned_positives>0){
			double startdreachtime = rtime(2);
			stats_under_updates++;
			positive_reach_detector->update();
			double reachUpdateElapsed = rtime(2)-startdreachtime;
			stats_under_update_time+=reachUpdateElapsed;
		}else
			stats_skipped_under_updates++;
			//stats_pure_skipped++;


		//  outer->reachupdatetime+=reachUpdateElapsed;


		if(!opt_detect_pure_theory_lits || unassigned_negatives>0){
			double startunreachtime = rtime(2);
			stats_over_updates++;
  		    negative_reach_detector->update();
  			double unreachUpdateElapsed = rtime(2)-startunreachtime;
  			stats_over_update_time+=unreachUpdateElapsed;
		}else
			stats_skipped_over_updates++;

		//outer->unreachupdatetime+=unreachUpdateElapsed;

		if(opt_rnd_shuffle){
			randomShuffle(rnd_seed, changed);
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

			}

		if(opt_rnd_shuffle && weighted_dist_lits.size()){
			randomShuffle(rnd_seed, weighted_dist_lits);
		}
		//now, check for weighted distance lits
		for(auto & dist_lit:weighted_dist_lits){
			Lit l = dist_lit.l;
			int to = dist_lit.u;
			Weight & min_dist =  dist_lit.min_distance;
			Weight & over_dist = positive_weighted_distance_detector->distance(to);
			Weight & under_dist = negative_weighted_distance_detector->distance(to);
			if(positive_weighted_distance_detector->connected(to) && positive_weighted_distance_detector->distance(to)<=min_dist){
				if(outer->value(l)==l_True){
					//do nothing
				}else if (outer->value(l)==l_Undef){
					outer->enqueue(l,weighted_reach_marker) ;
				}else if (outer->value(l)==l_False){
					//conflict
					conflict.push(l);
					buildDistanceLEQReason(to,min_dist,conflict);
					return false;
				}
			}
			if(! negative_weighted_distance_detector->connected(to) || negative_weighted_distance_detector->distance(to)>min_dist){
				if(outer->value(~l)==l_True){
					//do nothing
				}else if (outer->value(~l)==l_Undef){
					outer->enqueue(~l,weighted_non_reach_marker) ;
				}else if (outer->value(~l)==l_False){
					//conflict
					conflict.push(~l);
					buildDistanceGTReason(to,min_dist,conflict);
					return false;
				}
			}

		}



			#ifdef DEBUG_DIJKSTRA
					for(int i = 0;i<unweighted_dist_lits.size();i++){
						for(int j = 0;j<unweighted_dist_lits[i].size();j++){
							Lit l = unweighted_dist_lits[i][j].l;
							int dist =  unweighted_dist_lits[i][j].min_unweighted_distance;
							if(l!=lit_Undef){
								int u = getNode(var(l));
								if((!opt_detect_pure_theory_lits || unassigned_positives>0) && positive_reach_detector->distance_unsafe(u)<=dist){
									if(outer->dbg_value(l)!=l_True){
										assert(false);
										exit(3);
									}
								}else if ((!opt_detect_pure_theory_lits || unassigned_negatives>0) && negative_reach_detector->distance_unsafe(u)>dist){
									int d =negative_reach_detector->distance_unsafe(u);
									if(outer->dbg_value(l)!=l_False){
										assert(false);
										exit(3);
									}
								}
							}
						}
					}
			#endif
			return true;
		}
template<typename Weight>
bool DistanceDetector<Weight>::checkSatisfied(){
			UnweightedDijkstra<>under(source,g) ;
			UnweightedDijkstra<>over(source,antig) ;
			under.update();
			over.update();
			for(int j = 0;j< unweighted_dist_lits.size();j++){
				for(int k = 0;k<unweighted_dist_lits[j].size();k++){
					Lit l = unweighted_dist_lits[j][k].l;
					int dist = unweighted_dist_lits[j][k].min_unweighted_distance;

					if(l!=lit_Undef){
						int node =getNode(var(l));

						if(outer->value(l)==l_True){
							if(!under.connected(node) || under.distance(node)>dist){
								return false;
							}
						}else if (outer->value(l)==l_False){
							if(under.connected(node)  &&  over.distance(node)<=dist){
								return false;
							}
						}else{
							if(under.connected(node) && under.distance(node)<=dist){
								return false;
							}
							if(!over.connected(node) || over.distance(node)>dist){
								return false;
							}
						}
					}
				}
		}

	/*else{/
		/*{
			UnweightedDijkstra<>under(source,g) ;
			UnweightedDijkstra<>over(source,antig) ;
			under.update();
			over.update();
			for(int j = 0;j< unweighted_sat_lits.size();j++){
					for(int k = 0;k<unweighted_sat_lits[j].size();k++){
						Lit l = unweighted_sat_lits[j][k];
						int dist =j;

						if(l!=lit_Undef){
							int node =k;

							if(outer->value(l)==l_True){
								if(under.distance(node)>dist){
									return false;
								}
							}else if (outer->value(l)==l_False){
								if( over.distance(node)<=dist){
									return false;
								}
							}else{
								if(under.distance(node)<=dist){
									return false;
								}
								if(!over.distance(node)>dist){
									return false;
								}
							}
						}
					}
			}
		}*/
		{
			Dijkstra<Weight>under(source,g,weights) ;
			Dijkstra<Weight>over(source,antig,weights) ;
			under.update();
			over.update();
			//now, check for weighted distance lits
			for(auto & dist_lit:weighted_dist_lits){
				Lit l = dist_lit.l;
				int to = dist_lit.u;
				Weight & min_dist =  dist_lit.min_distance;
				Weight & under_dist = under.distance(to);
				Weight & over_dist = over.distance(to);
				if(under.connected(to) && under_dist<=min_dist){
					if(outer->value(l)==l_True){
						//do nothing
					}else if (outer->value(l)==l_Undef){
						outer->enqueue(l,weighted_reach_marker) ;
					}else if (outer->value(l)==l_False){
						return false;
					}
				}
				if(!over.connected(to) || over_dist>min_dist){
					if(outer->value(~l)==l_True){
						//do nothing
					}else if (outer->value(~l)==l_Undef){
						outer->enqueue(~l,weighted_non_reach_marker) ;
					}else if (outer->value(~l)==l_False){

						return false;
					}
				}
			}

		}

		//	}
	return true;
}

/*
int DistanceDetector<Weight>::OptimalWeightEdgeStatus::operator [] (int edge) const {
	Var v = detector.outer->edge_list[edge].v;
	lbool val = detector.outer->value(v);
	if(val==l_False){
		assert(false);
		return detector.outer->edge_list.size()*2;
	}else if (val==l_True){
		return 0;
	}else{
		assert(val==l_Undef);
		return 1;
	}
}

int DistanceDetector<Weight>::OptimalWeightEdgeStatus::size()const{
	return detector.outer->edge_list.size();
}
*/

template<typename Weight>
Lit DistanceDetector<Weight>::decide(){
	if(!opt_decide_graph_distance || !negative_reach_detector)
		return lit_Undef;
	DistanceDetector *r =this;
	auto * over =  r->negative_reach_detector;

	auto *  under =  r->positive_reach_detector;

	//we can probably also do something similar, but with cuts, for nodes that are decided to be unreachable.

	//ok, for each node that is assigned reachable, but that is not actually reachable in the under approx, decide an edge on a feasible path

	//this can be obviously more efficient
	//for(int j = 0;j<nNodes();j++){
/*	if(opt_decide_graph_neg){

	}*/

	for(int k = 0;k<unweighted_dist_lits.size();k++){
		for(int n = 0;n<unweighted_dist_lits[k].size();n++){
			Lit l = unweighted_dist_lits[k][n].l;
			int min_dist = unweighted_dist_lits[k][n].min_unweighted_distance;
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

							for(int i=0;i<outer->edge_list.size();i++){
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
				/*	assert(!under->connected(last));
					assert(under->connected(p));

					assert(over->connected(last));
					assert(over->connected(p));*/
					assert(p>-1);
					if(p>-1){
					Var v = outer->edges[p][last].v;
					if(outer->value(v)==l_Undef){
						return mkLit(v,false);
					}else{
						assert(outer->value(v)!=l_True);
					}
					}
			/*		for(int k = 0;k<outer->antig.adjacency[p].size();k++){
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
				/*		assert(opt_use_random_path_for_decisions);
						assert(tmp_nodes.size()>0);
						int i = irand(rnd_seed,tmp_nodes.size());
						Var v= tmp_nodes[i];
						return mkLit(v,true);*/

					}
				}
			}
		}
	}
	return lit_Undef;
};

template class DistanceDetector<int>;
template class DistanceDetector<double>;
#include <gmpxx.h>
template class DistanceDetector<mpq_class>;
