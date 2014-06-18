/*
 * SteinerDetector.c
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */



#include "SteinerDetector.h"
#include "GraphTheory.h"
#include "dgl/SteinerApprox.h"
#include "dgl/graph/DynamicNodes.h"
#include <limits>
#include <set>
#include "core/Config.h"
using namespace Minisat;
SteinerDetector::SteinerDetector(int detectorID, GraphTheorySolver * outer,  DynamicGraph &g,DynamicGraph &antig ,double seed):
Detector(detectorID),outer(outer),g(g),antig(antig),rnd_seed(seed),edge_weights(edge_weights),positive_reach_detector(NULL),negative_reach_detector(NULL){
	checked_unique=false;
	all_unique=true;
	positiveStatus = new SteinerDetector::SteinerStatus(*this,true);
	negativeStatus = new SteinerDetector::SteinerStatus(*this,false);
	positive_reach_detector =  new SteinerApprox<DynamicNodes,SteinerDetector::SteinerStatus>(g,underTerminalSet,*positiveStatus,1);//new SpiraPan<SteinerDetector::MSTStatus>(_g,*(positiveReachStatus),1);
	negative_reach_detector = new SteinerApprox<DynamicNodes,SteinerDetector::SteinerStatus>(antig,overTerminalSet,*negativeStatus,-1);


	reach_marker=outer->newReasonMarker(getID());
	non_reach_marker=outer->newReasonMarker(getID());

	reach_edge_marker=outer->newReasonMarker(getID());
	non_reach_edge_marker=outer->newReasonMarker(getID());
	first_reach_var=var_Undef;
}

void SteinerDetector::addTerminalNode(int node, Var outer_Var){
	Var var = outer->newVar(outer_Var,getID());
	underTerminalSet.addNode(node);
	underTerminalSet.setNodeEnabled(node,false);
	overTerminalSet.addNode(node);
	overTerminalSet.setNodeEnabled(node,true);
	terminal_map.growTo(g.nodes());
	terminal_map[node]=var;
	//vec<int> terminal_var_map;
	terminal_var_map.growTo(var+1,var_Undef);
	terminal_var_map[var]=node;
}

void SteinerDetector::addWeightLit(Var outer_weight_var,int min_weight){
	g.invalidate();
	antig.invalidate();

	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});
	Var weight_var = outer->newVar(outer_weight_var,getID());

	Lit reachLit=mkLit(weight_var,false);
	bool found=false;
	for(int i = 0;i<weight_lits.size();i++){
		if(weight_lits[i].min_weight==min_weight){
			found=true;
			Lit r = weight_lits[i].l;
			//force equality between the new lit and the old reach lit, in the SAT solver
			outer->makeEqual(r,reachLit);
			/*outer->S->addClause(~r, reachLit);
			outer->S->addClause(r, ~reachLit);*/
		}
	}
	if(!found){
		weight_lits.push();
		weight_lits.last().l = reachLit;
		weight_lits.last().min_weight=min_weight;

		//weight_lit_map.insert(min_weight,weight_lits.size()-1);
	}



}

void SteinerDetector::SteinerStatus::setMinimumSteinerTree(int weight){

	for(int i = 0;i<detector.weight_lits.size();i++){
		int min_weight =  detector.weight_lits[i].min_weight;
		Lit l = detector.weight_lits[i].l;
		if(l!=lit_Undef){
			assert(l!=lit_Undef);
			if(min_weight<weight && !polarity){
				lbool assign = detector.outer->value(l);
				if( assign!= l_False ){
					detector.changed_weights.push({~l,min_weight});
				}
			}else if(min_weight>=weight && polarity){
				lbool assign = detector.outer->value(l);
				if( assign!= l_True ){
					detector.changed_weights.push({l,min_weight});
				}
			}
		}
	}

}
void SteinerDetector::buildMinWeightTooSmallReason(int weight,vec<Lit> & conflict){


		}


		void SteinerDetector::buildMinWeightTooLargeReason(int weight,vec<Lit> & conflict){



		}




		void SteinerDetector::buildReason(Lit p, vec<Lit> & reason, CRef marker){

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
					for(int i = 0;i<weight_lits.size();i++){
						if(var(weight_lits[i].l)==v){
							weight=weight_lits[i].min_weight;
							break;
						}
					}
					assert(weight>=0);

					buildMinWeightTooLargeReason(weight,reason);

				}else{
					assert(false);
				}
		}

		bool SteinerDetector::propagate(vec<Lit> & conflict){
			static int it = 0;
			if(++it==7){
				int a = 1;
			}

			double startdreachtime = rtime(2);
			stats_under_updates++;
			positive_reach_detector->update();
			double reachUpdateElapsed = rtime(2)-startdreachtime;
			stats_under_update_time+=reachUpdateElapsed;
		//}else
		//	stats_skipped_under_updates++;

		//if(negative_reach_detector && (!opt_detect_pure_theory_lits || unassigned_negatives>0)){
			double startunreachtime = rtime(2);
			stats_over_updates++;
			negative_reach_detector->update();
			double unreachUpdateElapsed = rtime(2)-startunreachtime;
			stats_over_update_time+=unreachUpdateElapsed;
		//}else
		//	stats_skipped_over_updates++;



		for(int j = 0;j<changed_weights.size();j++){
			Lit l = changed_weights[j].l;
			//printf("mst: %d\n",dimacs(l));
			int weight = changed_weights[j].weight;

			bool reach = !sign(l);
			if(outer->value(l)==l_True){
				//do nothing
			}else if(outer->value(l)==l_Undef){
				//trail.push(Assignment(false,reach,detectorID,0,var(l)));
				if(reach)
					outer->enqueue(l,reach_marker) ;
				else
					outer->enqueue(l,non_reach_marker) ;

			}else if (outer->value(l)==l_False){
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

			return true;
		}

bool SteinerDetector::checkSatisfied(){
	Kruskal<MinimumSpanningTree::NullStatus> positive_checker(g,MinimumSpanningTree::nullStatus,0);
	Kruskal<MinimumSpanningTree::NullStatus> negative_checker(antig,MinimumSpanningTree::nullStatus,0);
	positive_checker.update();
	negative_checker.update();
	for(int k = 0;k<weight_lits.size();k++){
		Lit l = weight_lits[k].l;
		int dist = weight_lits[k].min_weight;

		if(l!=lit_Undef){


			if(outer->value(l)==l_True){
				if(positive_checker.weight()>dist){
					return false;
				}
			}else if (outer->value(l)==l_False){
				if( negative_checker.weight()<=dist){
					return false;
				}
			}else{
				if(positive_checker.weight()<=dist){
					return false;
				}
				if(!negative_checker.weight()>dist){
					return false;
				}
			}
		}
	}

	return true;
}
Lit SteinerDetector::decide(){

	return lit_Undef;
};


