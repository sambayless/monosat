/****************************************************************************************[Solver.h]
The MIT License (MIT)

Copyright (c) 2014, Sam Bayless

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


#include "MaxflowDetector.h"
#include "GraphTheory.h"
#include "dgl/EdmondsKarpAdj.h"
#include "dgl/KohliTorr.h"
#include "dgl/EdmondsKarpDynamic.h"
#include "dgl/Dinics.h"
#include "dgl/DinicsLinkCut.h"
using namespace Monosat;

template<>
void MaxflowDetector<int>::buildDinitzLinkCut(){
	positive_detector =  new DinitzLinkCut<std::vector<int>>(g,capacities);
	negative_detector = new DinitzLinkCut<std::vector<int>>(antig,capacities);
}

template<typename Weight>
void MaxflowDetector<Weight>::buildDinitzLinkCut(){
	positive_detector =  nullptr;
	negative_detector = nullptr;
	assert(false);
	fprintf(stderr,"Warning: Dinitz Link/Cut Tree implementation only supports unweighted or integer weight graphs, aborting!\n");
	exit(1);
}


template<typename Weight>
MaxflowDetector<Weight>::MaxflowDetector(int _detectorID, GraphTheorySolver<Weight> * _outer,std::vector<Weight> & capacities,  DynamicGraph &_g,DynamicGraph &_antig, int from, int _target,double seed):
Detector(_detectorID),outer(_outer),capacities(capacities),over_graph(_g),g(_g),antig(_antig),source(from),target(_target),rnd_seed(seed),positive_detector(NULL),negative_detector(NULL){
	if(mincutalg==MinCutAlg::ALG_EDKARP_DYN){
		positive_detector = new EdmondsKarpDynamic<std::vector<Weight>,Weight>(_g,capacities);
		negative_detector = new EdmondsKarpDynamic<std::vector<Weight>,Weight>(_antig,capacities);
		positive_conflict_detector =positive_detector;//new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_conflict_detector =negative_detector;// new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if(opt_conflict_min_cut)
				learn_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps);
	}else if (mincutalg==MinCutAlg::ALG_EDKARP_ADJ){
		positive_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		positive_conflict_detector = positive_detector;
		negative_conflict_detector = negative_detector;
		if(opt_conflict_min_cut)
				learn_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps);
	}/*else if (mincutalg==MinCutAlg::ALG_IBFS){
		positive_detector = new IBFS(_g);
		negative_detector = new IBFS(_antig);
		positive_conflict_detector =new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
	}*/else if (mincutalg==MinCutAlg::ALG_DINITZ){
		positive_detector = new Dinitz<std::vector<Weight>,Weight>(_g,capacities);
		negative_detector = new Dinitz<std::vector<Weight>,Weight>(_antig,capacities);
		positive_conflict_detector = positive_detector;// new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_conflict_detector = negative_detector;//new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if(opt_conflict_min_cut)
				learn_cut = new Dinitz<std::vector<int>,int>(learn_graph,learn_caps);
	}else if (mincutalg==MinCutAlg::ALG_DINITZ_LINKCUT){
		//link-cut tree currently only supports ints (enforcing this using tempalte specialization...).
		buildDinitzLinkCut();
		positive_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if(opt_conflict_min_cut)
				learn_cut = new DinitzLinkCut<std::vector<int>>(learn_graph,learn_caps);
	}else if (mincutalg==MinCutAlg::ALG_KOHLI_TORR){
		positive_detector = new KohliTorr<std::vector<Weight>,Weight>(_g,capacities);
		negative_detector = new KohliTorr<std::vector<Weight>,Weight>(_antig,capacities);
		positive_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_conflict_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		if(opt_conflict_min_cut)
				learn_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps);
	}else{
		positive_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_g,capacities);
		negative_detector = new EdmondsKarpAdj<std::vector<Weight>,Weight>(_antig,capacities);
		positive_conflict_detector = positive_detector;
		negative_conflict_detector = negative_detector;
		if(opt_conflict_min_cut)
				learn_cut = new EdmondsKarpAdj<std::vector<int>,int>(learn_graph,learn_caps);
	}






	first_reach_var = var_Undef;
	reach_marker=outer->newReasonMarker(getID());
	non_reach_marker=outer->newReasonMarker(getID());
	forced_reach_marker=outer->newReasonMarker(getID());
}


template<typename Weight>
void MaxflowDetector<Weight>::addFlowLit(Weight maxflow, Var outer_reach_var){
	g.invalidate();
	antig.invalidate();
	Var reach_var = outer->newVar(outer_reach_var,getID());
	if(first_reach_var==var_Undef){
		first_reach_var=reach_var;
	}else{
		assert(reach_var>first_reach_var);
	}

	assert(maxflow>=0);



	//while( dist_lits[to].size()<=within_steps)
	//	dist_lits[to].push({lit_Undef,-1});

/*	while(outer->S->nVars()<=reach_var)
		outer->S->newVar();*/

	Lit reachLit=mkLit(reach_var,false);
	bool found=false;

	flow_lits.push();
	flow_lits.last().l = reachLit;
	flow_lits.last().max_flow=maxflow;
	while(reach_lit_map.size()<= reach_var- first_reach_var ){
		reach_lit_map.push(-1);
	}

	reach_lit_map[reach_var-first_reach_var]=flow_lits.size()-1;

}





template<typename Weight>
void MaxflowDetector<Weight>::buildMaxFlowTooHighReason(Weight flow,vec<Lit> & conflict){
			//drawFull();
			double starttime = rtime(2);
			tmp_cut.clear();
			Weight actual_flow = positive_conflict_detector->maxFlow(source, target);

			//just collect the set of edges which have non-zero flow, and return them
			//(Or, I could return a cut, probably)
			for(int i = 0;i<g.edges();i++){
				if(g.edgeEnabled(i)){
					if(positive_conflict_detector->getEdgeFlow(i)>0){
						Var v = outer->edge_list[i].v;
						assert(outer->value(v)==l_True);
						conflict.push(mkLit(v,true));
					}
				}
			}

#ifdef RECORD
			std::sort(conflict.begin(),conflict.end());
		if(g.outfile){
			fprintf(g.outfile,"toohigh ");
			for(int i = 0;i<conflict.size();i++){
				fprintf(g.outfile,"%d,",dimacs(conflict[i]));
			}
			fprintf(g.outfile,"\n");
			fflush(g.outfile);
		}
		if(antig.outfile){
			fprintf(antig.outfile,"toohigh ");
			for(int i = 0;i<conflict.size();i++){
				fprintf(antig.outfile,"%d,",dimacs(conflict[i]));
			}
			fprintf(antig.outfile,"\n");
				fflush(antig.outfile);
			}
#endif
			stats_under_conflicts++;
			outer->num_learnt_paths++;
			outer->learnt_path_clause_length+= (conflict.size()-1);
			double elapsed = rtime(2)-starttime;
			outer->pathtime+=elapsed;
		}
template<typename Weight>
		void MaxflowDetector<Weight>::buildMaxFlowTooLowReason(Weight maxflow,vec<Lit> & conflict){
			static int it = 0;
			++it;
			if(it==22){
				int a=1;
			}
			double starttime = rtime(2);
			if(opt_conflict_min_cut){
				Weight foundflow = negative_conflict_detector->maxFlow(source,target);
				//find a minimal s-t cut in the residual graph.
				//to do so, first construct the residual graph
				int INF=0x0FF0F0;
				//for each edge in the original graph, need to add a forward and backward edge here.
				if(learn_graph.nodes()<g.nodes()){
					while(learn_graph.nodes()<g.nodes())
						learn_graph.addNode();

					for (auto & e:g.all_edges){
						learn_graph.addEdge(e.from,e.to);
					}
					back_edges.growTo(g.edges());
					for (auto & e:g.all_edges){
						back_edges[e.id] = learn_graph.addEdge(e.to,e.from);
					}
					learn_caps.resize(learn_graph.edges());
				}

				//now, set learn_graph to the residual graph
				for (auto & e:g.all_edges){
					int from = e.from;
					int to = e.to;
					int v = outer->getEdgeVar(e.id);
					lbool val =outer->value(v);
					if(val!=l_False){
						Weight flow = negative_conflict_detector->getEdgeFlow(e.id);
						Weight capacity =  negative_conflict_detector->getEdgeCapacity(e.id);
						if(capacity==0){
							learn_graph.disableEdge(e.id);
							learn_graph.disableEdge(back_edges[e.id]);
						}else{
							if(flow>0){
								//then there is capacity in the backward edge in the residual graph
								int back_edge = back_edges[e.id];
								learn_graph.enableEdge(back_edges[e.id]);
								learn_caps[back_edge] = INF;
							}else{
								learn_graph.disableEdge(back_edges[e.id]);
							}
							if (flow<capacity){
								//then there is capacity in the forward edge in the residual graph
								learn_caps[e.id] = INF;
								learn_graph.enableEdge(e.id);
							}else{
								learn_graph.disableEdge(e.id);
							}
						}
					}else{

						learn_graph.enableEdge(e.id);
						learn_graph.disableEdge(back_edges[e.id]);
						learn_caps[e.id] = 1;
					}
				}
				learn_graph.invalidate();
				learn_graph.drawFull(true);
				outer->cut.clear();
				antig.drawFull(true);
				int f =learn_cut->minCut(source,target,outer->cut);
				if(f<0x0FF0F0){
					assert(f<0xF0F0F0); assert(f==outer->cut.size());//because edges are only ever infinity or 1
					for(int i = 0;i<outer->cut.size();i++){
						MaxFlowEdge e = outer->cut[i];

						Lit l = mkLit( outer->getEdgeVar(e.id),false);
						assert(outer->value(l)==l_False);
						conflict.push(l);
					}
				}else{
					int a=1;
					//there is no way to increase the max flow.
				}

				 outer->num_learnt_cuts++;
				 outer->learnt_cut_clause_length+= (conflict.size()-1);
				 stats_over_conflicts++;
				return;
			}

			//drawFull( non_reach_detectors[detector]->getSource(),u);
			//assert(outer->dbg_distance( source,u));

			//The reason why we can't reach this assignment is a cut through the disabled edges in the residual graph from the overapprox.
			//we could search for a min cut, but instead we will just step back in the RESIDUAL graph, from u to s, collecting disabled edges.
			seen.clear();
			seen.growTo(outer->nNodes());
			visit.clear();
			Weight foundflow = negative_conflict_detector->maxFlow(source,target);
			visit.push(target);
			for(int k = 0;k<visit.size();k++){
				int u = visit[k];
				for(int i = 0;i<g.nIncoming(u);i++){
					int p = g.incoming(u,i).node;
					if(!seen[p]){


						int edgeid = g.incoming(u,i).id;
						int v = outer->getEdgeVar(edgeid);

						//assert( g.incoming(u,i).to==u);
						if(outer->value(v)!=l_False){
							//this is an enabled edge in the overapprox

							Weight residual_capacity = negative_conflict_detector->getEdgeResidualCapacity(edgeid);
							if(residual_capacity>0){
								seen[p]=true;
								visit.push(p);
							}
						}else{
							//this is a disabled edge, and we can add it to the cut.
							//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
							conflict.push(mkLit(v,false));
						}
					}
					//pred
				}
				//we are walking back in the RESIDUAL graph, which can mean backwards traversal of edges with flow!
				for(int i = 0;i<g.nIncident(u);i++){
					int p = g.incident(u,i).node;
					if(!seen[p]){
						int edgeid = g.incident(u,i).id;
						int v = outer->getEdgeVar(edgeid);

						//assert( g.incoming(u,i).to==u);
						if(outer->value(v)!=l_False){
							//this is the residual capacity of the backwards edge in the residual graph - which is equal to the forwards flow on this edge!
							Weight residual_capacity = negative_conflict_detector->getEdgeFlow(edgeid);
							if(residual_capacity ){
								seen[p]=true;
								visit.push(p);
							}
						}else{
							//this is a disabled edge, and we can add it to the cut.
							//we're going to assume the edge has non-zero capacity here, otherwise we could exclude it (but it shouldn't really even be in this graph in that case, anyways).
							conflict.push(mkLit(v,false));
						}
					}
					//pred
				}
			}

#ifdef RECORD
			if(negative_detector!=negative_conflict_detector){
				Weight foundflow2 = negative_detector->maxFlow(source,target);
				assert(foundflow==foundflow2);
			}
			std::sort(conflict.begin(),conflict.end());
		if(g.outfile){
			fprintf(g.outfile,"toolow ");
			for(int i = 0;i<conflict.size();i++){
				fprintf(g.outfile,"%d,",dimacs(conflict[i]));
			}
			fprintf(g.outfile,"\n");
			fflush(g.outfile);
		}
		if(antig.outfile){
			fprintf(antig.outfile,"toolow ");
			for(int i = 0;i<conflict.size();i++){
				fprintf(antig.outfile,"%d,",dimacs(conflict[i]));
			}
			fprintf(antig.outfile,"\n");
				fflush(antig.outfile);
			}
#endif
			 outer->num_learnt_cuts++;
			 outer->learnt_cut_clause_length+= (conflict.size()-1);
			 stats_over_conflicts++;
			double elapsed = rtime(2)-starttime;
			 outer->mctime+=elapsed;

		}

template<typename Weight>
		void MaxflowDetector<Weight>::buildReason(Lit p, vec<Lit> & reason, CRef marker){


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
					Weight flow = flow_lits[reach_lit_map[v-first_reach_var]].max_flow;
					//int u =getNode(v);
					buildMaxFlowTooHighReason(flow,reason);


					//double elapsed = rtime(2)-startpathtime;
				//	pathtime+=elapsed;
				}else if(marker==non_reach_marker){
					reason.push(p);

					//the reason is a cut separating p from s;
					//We want to find a min-cut in the full graph separating, where activated edges (ie, those still in antig) are weighted infinity, and all others are weighted 1.

					//This is a cut that describes a minimal set of edges which are disabled in the current graph, at least one of which would need to be activated in order for s to reach p
					//assign the mincut edge weights if they aren't already assigned.


					Var v = var(p);
					Weight flow = flow_lits[reach_lit_map[v-first_reach_var]].max_flow;
					//int t = getNode(v); // v- var(reach_lits[d][0]);
					buildMaxFlowTooLowReason(flow,reason);

				}else{
					assert(false);
				}
		}
		template<typename Weight>
		bool MaxflowDetector<Weight>::propagate(vec<Lit> & conflict){
			if(flow_lits.size()==0){
				return true;
			}
			static int iter1=0;
			if(++iter1==135){
				int a=1;
			}
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"iter %d\n",iter1);
			fflush(g.outfile);
		}
		if(antig.outfile){
					fprintf(antig.outfile,"iter %d\n",iter1);
					fflush(antig.outfile);
				}
#endif

			Weight over_maxflow=-1;
			Weight under_maxflow=-1;

			if(positive_detector && (!opt_detect_pure_theory_lits || unassigned_positives>0)){
				double startdreachtime = rtime(2);
				stats_under_updates++;
				under_maxflow = positive_detector->maxFlow(source,target);
				assert(under_maxflow==positive_conflict_detector->maxFlow(source,target));
				double reachUpdateElapsed = rtime(2)-startdreachtime;
				stats_under_update_time+=reachUpdateElapsed;
			}else
				stats_skipped_under_updates++;

			if(negative_detector && (!opt_detect_pure_theory_lits || unassigned_negatives>0)){
				double startunreachtime = rtime(2);
				stats_over_updates++;
				over_maxflow = negative_detector->maxFlow(source,target);
				assert(over_maxflow==negative_conflict_detector->maxFlow(source,target));
				double unreachUpdateElapsed = rtime(2)-startunreachtime;
				stats_over_update_time+=unreachUpdateElapsed;
			}else
				stats_skipped_over_updates++;

			for(int j = 0;j<flow_lits.size();j++){

				DistLit f = flow_lits[j];
				Lit l = f.l;
				Weight& maxflow = f.max_flow;
				//int u = getNode(var(l));


				if(under_maxflow>=maxflow){
					if(outer->value(l)==l_True){
						//do nothing
					}else if(outer->value(l)==l_Undef){
						//trail.push(Assignment(false,true,detectorID,0,var(l)));
						outer->enqueue(l,reach_marker) ;

					}else if(outer->value(l)==l_False){
						conflict.push(l);
						buildMaxFlowTooHighReason(maxflow,conflict);
						return false;
					}

				}else if(over_maxflow<maxflow){
					if(outer->value(l)==l_False){
						//do nothing
					}else if(outer->value(l)==l_Undef){
						//trail.push(Assignment(false,false,detectorID,0,var(l)));
						outer->enqueue(~l,reach_marker) ;

					}else if(outer->value(l)==l_True){
						conflict.push(~l);
						buildMaxFlowTooLowReason(maxflow,conflict);
						return false;
					}

				}



			}

			return true;
		}
template<typename Weight>
bool MaxflowDetector<Weight>::checkSatisfied(){
	EdmondsKarpAdj<std::vector<Weight>,Weight> positiveCheck(g,capacities);
	EdmondsKarpAdj<std::vector<Weight>,Weight> negativeCheck(antig,capacities);
		for(int j = 0;j< flow_lits.size();j++){

				Lit l = flow_lits[j].l;
				Weight dist = flow_lits[j].max_flow;

				if(l!=lit_Undef){
					//int node =getNode(var(l));

					if(outer->value(l)==l_True){
						if(positiveCheck.maxFlow(source,target)<dist){
							return false;
						}
					}else if (outer->value(l)==l_False){
						if( negativeCheck.maxFlow(source,target)>=dist){
							return false;
						}
					}else{
						if(positiveCheck.maxFlow(source,target)>=dist){
							return false;
						}
						if(!negativeCheck.maxFlow(source,target)<dist){
							return false;
						}
					}
				}

		}
	return true;
}
template<typename Weight>
void MaxflowDetector<Weight>::printSolution(){
	Weight f = positive_detector->maxFlow(source,target);
	printf("Maximum %d->%d flow in graph %d is ", this->source,this->target,this->outer->getGraphID());
	std::cout<<f <<"\n";
	if(opt_verb>0){

		int maxw = log10(outer->nNodes() )+1;
		int width = sqrt(outer->nNodes());
		if(opt_width>0){
				width=opt_width;
		}
			int height =width;
			if(opt_height>0){
				height = opt_height;
			}

			int lasty= 0;
			int extra =  outer->nNodes() % width ? (width- outer->nNodes() % width ):0;
			for(int n = 0;n<outer->nNodes();n++){
				int x = n%width;

				int y = (n )/width;
				if(y > lasty)
					printf("\n");
				Weight total_flow = 0;
				for(int e = 0;e<g.edges();e++){
					if(g.getEdge(e).to==n && g.edgeEnabled(e)){
						total_flow+=positive_detector->getEdgeFlow(e);

					}
				}

				//printf("%*d ",maxw,total_flow);
				std::cout<<total_flow<<" ";
					lasty=y;
				}
				printf("\n");

				for(int n = 0;n<outer->nNodes();n++){

					int total_flow = 0;
					for(int e = 0;e<g.edges();e++){
						if(g.getEdge(e).to==n && g.edgeEnabled(e)){
							Weight flow = positive_detector->getEdgeFlow(e);
							if(flow>0){
								printf("flow  %d to %d = ", g.getEdge(e).from,g.getEdge(e).to);
								std::cout<<flow<<"\n";
							}
						}
					}

				}
				printf("\n");
	}

}

template<typename Weight>
Lit MaxflowDetector<Weight>::decide(){
	static int it =0;
	++it;
	double startdecidetime = rtime(2);
	auto * over =negative_conflict_detector;
		auto * under = positive_conflict_detector;

		Weight under_flow = under->maxFlow(source,target) ;
		Weight over_flow = over->maxFlow(source,target) ;

		if(to_decide.size() && last_decision_status== over->numUpdates()){
			while(to_decide.size()){
				Lit l = to_decide.last();
				to_decide.pop();
				if(outer->value(l)==l_Undef){
					stats_decide_time+=  rtime(2)-startdecidetime;
					return l;
				}
			}
		}


		for(int k = 0;k<flow_lits.size();k++){
			Lit l =flow_lits[k].l;
			if(l==lit_Undef)
				continue;

			Weight & required_flow = flow_lits[k].max_flow;

			if(outer->value(l)==l_True && opt_decide_graph_pos){
				assert(over_flow>=required_flow);

#ifndef NDEBUG
				static vec<bool> dbg_expect;
				int dbg_count = 0;
				dbg_expect.clear();
				dbg_expect.growTo(g.edges());
				for(int edgeID = 0;edgeID<g.edges();edgeID++){
					lbool val = outer->value(outer->getEdgeVar(edgeID));
					if(val==l_Undef){
						if(over->getEdgeFlow(edgeID)>0){
							dbg_expect[edgeID]=true;
							dbg_count++;
						}
					}else if (val == l_False){
						assert(over->getEdgeFlow(edgeID)==0);
					}
				}
#endif

				if(under_flow <required_flow){

					//then decide an unassigned edge of the currently selected flow
					to_decide.clear();
					last_decision_status= over->numUpdates();

					seen.clear();
					seen.growTo(antig.nodes(),false);
					q.clear();

					if(opt_conflict_dfs){
						//do a dfs to find that edge. Could start from either the source or the target.




						if(opt_conflict_from_source){

							q.push_back(source);

							while(q.size()){
								int u =q.back();
								q.pop_back();
								for(int i = 0;i<antig.nIncident(u) ;i++){
									if(!antig.edgeEnabled(antig.incident(u,i).id))
										continue;
									int v =antig.incident(u,i).node;
									int edgeID= antig.incident(u,i).id;
									if (over->getEdgeFlow(edgeID)>0){
										Var var = outer->getEdgeVar(edgeID);
										if(outer->value(var)==l_Undef){
											to_decide.push(mkLit(var,false));
										}
										if(!seen[v]){
											seen[v]=true;
											q.push_back(v);
										}
									}
								}
							}
						}else{
							q.push_back(target);
							while(q.size()){
								int u =q.back();
								q.pop_back();
								for(int i = 0;i<antig.nIncoming(u) ;i++){
									if(!antig.edgeEnabled(antig.incoming(u,i).id))
										continue;
									int v =antig.incoming(u,i).node;
									int edgeID= antig.incoming(u,i).id;
									if (over->getEdgeFlow(edgeID)>0){
										Var var = outer->getEdgeVar(edgeID);
										if(outer->value(var)==l_Undef){
											to_decide.push(mkLit(var,false));
										}
										if(!seen[v]){
											seen[v]=true;
											q.push_back(v);
										}
									}
								}
							}
						}
					}else{
						if(opt_conflict_from_source){
							q.push_back(source);
							//do a bfs to find that edge. Could start from either the source or the target.
							for (int i = 0;i<q.size();i++){
								int u = q[i];
								for(int i = 0;i<antig.nIncident(u);i++){
									if(!antig.edgeEnabled(antig.incident(u,i).id))
										continue;
									int edgeID = antig.incident(u,i).id;
									int v = antig.incident(u,i).node;
									if (over->getEdgeFlow(edgeID)>0){
										Var var = outer->getEdgeVar(edgeID);
										assert(outer->value(var)!=l_False);
										if(outer->value(var)==l_Undef){
											to_decide.push(mkLit(var,false));
										}
										if(!seen[v]){
											seen[v]=true;
											q.push_back(v);
										}
									}
								}
							}
						}else{
							q.push_back(target);
							//do a bfs to find that edge. Could start from either the source or the target.
							for (int i = 0;i<q.size();i++){
								int u = q[i];
								for(int i = 0;i<antig.nIncoming(u);i++){
									if(!antig.edgeEnabled(antig.incoming(u,i).id))
										continue;
									int edgeID = antig.incoming(u,i).id;
									int v = antig.incoming(u,i).node;
									if (over->getEdgeFlow(edgeID)>0){
										Var var = outer->getEdgeVar(edgeID);
										assert(outer->value(var)!=l_False);
										if(outer->value(var)==l_Undef){
											to_decide.push(mkLit(var,false));
										}
										if(!seen[v]){
											seen[v]=true;
											q.push_back(v);
										}
									}
								}
							}
						}
					}
					assert(to_decide.size()==dbg_count);


				}

			}else if(outer->value(l)==l_False && opt_decide_graph_neg){



			}
			if(to_decide.size() && last_decision_status== over->numUpdates()){
				while(to_decide.size()){
					Lit l = to_decide.last();
					to_decide.pop();
					if(outer->value(l)==l_Undef){
						stats_decide_time+= rtime(2)-startdecidetime;
						return l;
					}
				}
			}
		}
		stats_decide_time+= rtime(2)-startdecidetime;
		return lit_Undef;
};


template class MaxflowDetector<int>;
template class MaxflowDetector<double>;
#include <gmpxx.h>
template class MaxflowDetector<mpq_class>;
