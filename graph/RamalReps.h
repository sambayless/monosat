/*
 * Dijkstra.h
 *
 *  Created on: 2013-05-28
 *      Author: sam
 */

#ifndef RAMAL_REPS_H_
#define RAMAL_REPS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "Reach.h"
#include "Dijkstra.h"
#include "core/Config.h"

namespace Minisat{
template<class Status,class EdgeStatus=DefaultEdgeStatus>
class RamalReps:public Reach{
public:
	DynamicGraph<EdgeStatus> & g;
	Status &  status;
	int reportPolarity;
	bool reportDistance;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	int INF;





	vec<int> old_dist;
	vec<int> changed;
	vec<bool> node_changed;
	vec<int> dist;
	vec<int> prev;
	struct DistCmp{
		vec<int> & _dist;
		 bool operator()(int a, int b)const{
			return _dist[a]<_dist[b];
		}
		 DistCmp(vec<int> & d):_dist(d){};
	};
	Heap<DistCmp> q;

	vec<int> edgeInShortestPathGraph;
	vec<int> delta;
	vec<int> changeset;
public:

	int stats_full_updates;
	int stats_fast_updates;
	int stats_fast_failed_updates;
	int stats_skip_deletes;
	int stats_skipped_updates;
	int stats_num_skipable_deletions;
	double mod_percentage;

	double stats_full_update_time;
	double stats_fast_update_time;
	RamalReps(int s,DynamicGraph<EdgeStatus> & graph,	Status &  status, int reportPolarity=0,bool reportDistance=false):g(graph),status(status),reportPolarity(reportPolarity),reportDistance(reportDistance), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),q(DistCmp(dist)){

		mod_percentage=0.2;
		stats_full_updates=0;
		stats_fast_updates=0;
		stats_skip_deletes=0;
		stats_skipped_updates=0;
		stats_full_update_time=0;
		stats_fast_update_time=0;
	}
	//Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};


	void setSource(int s){
		source = s;
		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}
	int getSource(){
		return source;
	}

	vec<int> & getChanged(){
		return changed;
	}
	void clearChanged(){
		changed.clear();
	}

	void drawFull(){

	}

	void dbg_delta(){
#ifndef NDEBUG
		dbg_delta_lite();
		assert(delta.size()==g.nodes);

		for(int i = 0;i<g.nEdgeIDs();i++){
			if(!g.edgeEnabled(i)){
				assert(!edgeInShortestPathGraph[i]);
			}
		}

		vec<int> dbg_delta;
		vec<int> dbg_dist;
		dbg_dist.growTo(g.nodes,INF);
		dbg_delta.growTo(g.nodes);
		dbg_dist[getSource()]=0;
		struct DistCmp{
			vec<int> & _dist;
			 bool operator()(int a, int b)const{
				return _dist[a]<_dist[b];
			}
			 DistCmp(vec<int> & d):_dist(d){};
		};
		Heap<DistCmp> q(dbg_dist);

		q.insert(getSource());

		while(q.size()){
			int u = q.removeMin();
			if(dbg_dist[u]==INF)
				break;
			dbg_delta[u]=0;

			for(int i = 0;i<g.inverted_adjacency[u].size();i++){
					if(!g.edgeEnabled( g.inverted_adjacency[u][i].id))
						continue;

					int edgeID = g.inverted_adjacency[u][i].id;
					int v = g.all_edges[edgeID].from;
					int alt = dbg_dist[v]+ 1;
					assert(alt>=dbg_dist[u]);
				/*	if (alt==dbg_dist[u]){
						dbg_delta[u]++;
					}*/
				}

			for(int i = 0;i<g.adjacency[u].size();i++){
				if(!g.edgeEnabled( g.adjacency[u][i].id))
					continue;

				int edgeID = g.adjacency[u][i].id;
				int v = g.all_edges[edgeID].to;
				int alt = dbg_dist[u]+ 1;
				if(alt<dbg_dist[v]){

					dbg_dist[v]=alt;

					if(!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}/*else if (alt==dbg_dist[v]){
					dbg_delta[v]++;
				}*/
			}
		}

		for(int u = 0;u<g.nodes;u++){
			int d = dist[u];
			if(u==53){
				int a=1;
			}
			assert(dbg_dist[u]==dist[u]);
			for(int i = 0;i<g.inverted_adjacency[u].size();i++){
				if(!g.edgeEnabled( g.inverted_adjacency[u][i].id))
					continue;

				int edgeID = g.inverted_adjacency[u][i].id;
				int v = g.all_edges[edgeID].from;
				int alt = dbg_dist[v]+ 1;
				assert(alt>=dbg_dist[u]);
				if(alt==dbg_dist[u]){
					dbg_delta[u]++;
					assert(edgeInShortestPathGraph[edgeID]);
				}else{
					assert(!edgeInShortestPathGraph[edgeID]);
				}

			}

		}
		for(int u = 0;u<g.nodes;u++){
			int d = dist[u];
			int d_expect = dbg_dist[u];
			assert(d==dbg_dist[u]);

		}
		for(int u = 0;u<g.nodes;u++){
			int d = delta[u];
			int d_expect = dbg_delta[u];
			assert(d==dbg_delta[u]);

		}
		dbg_delta_lite();
#endif
	}

	void GRRInc(int edgeID){
		static int iter = 0;
		++iter;
		dbg_delta_lite();
		assert(g.edgeEnabled(edgeID));
		if(edgeInShortestPathGraph[edgeID])
			return;
		int ru =  g.all_edges[edgeID].from;
		int rv =  g.all_edges[edgeID].to;
		if(rv==63){
				int a=1;
			}
		int rdv = dist[rv];
		int rdu = dist[ru];
		if(rv==1 ){
			int a =1;
		}
		int weight = 1;
		if(dist[rv]<dist[ru]+weight)
			return;
		else if (dist[rv]==dist[ru]+weight ){
			assert(!edgeInShortestPathGraph[edgeID]);
			edgeInShortestPathGraph[edgeID]=true;
			delta[rv]++;//we have found an alternative shortest path to v
			return;
		}
		edgeInShortestPathGraph[edgeID]=true;
		delta[rv]++;
		dist[rv]=dist[ru]+weight;
		q.clear();
		q.insert(rv);

		while(q.size()){

			int u = q.removeMin();
			if(u==63 || u==55){
					int a=1;
				}
			if(!node_changed[u]){
				node_changed[u]=true;
				changed.push(u);
			}
			delta[u]=0;
			for(auto & e:g.inverted_adjacency[u]){
				int adjID = e.id;
				if (g.edgeEnabled(adjID)){

						assert(g.all_edges[adjID].to==u);
						int v =g.all_edges[adjID].from;
						int w = 1;//assume a weight of one for now
						int du = dist[u];
						int dv = dist[v];
						if(dist[u]==dist[v]+w ){
							edgeInShortestPathGraph[adjID]=true;
							delta[u]++;
						}else if (dist[u]<(dist[v]+w)){
							//This doesn't hold for us, because we are allowing multiple edges to be added at once.
							//assert(dist[u]<(dist[v]+w));

							edgeInShortestPathGraph[adjID]=false;
						}else{
							//don't do anything. This will get corrected in a future call to GRRInc.
							//assert(false);

						}
				}else{
					edgeInShortestPathGraph[adjID]=false;//need to add this, because we may have disabled multiple edges at once.
				}
			}

			for(auto & e:g.adjacency[u]){
				int adjID = e.id;
				if (g.edgeEnabled(adjID)){
						assert(g.all_edges[adjID].from==u);
						int s =g.all_edges[adjID].to;
						int w = 1;//assume a weight of one for now
						int du = dist[u];
						int ds = dist[s];
						if(dist[s]>dist[u]+w){
							dist[s]=dist[u]+w;
							q.update(s);
						}else if (dist[s]==dist[u]+w && !edgeInShortestPathGraph[adjID]){
							edgeInShortestPathGraph[adjID]=true;
							delta[s]++;
						}
				}
			}
		}
		dbg_delta_lite();
	}
	void dbg_delta_lite(){
#ifndef NDEBUG
		for(int u = 0;u<g.nodes;u++){
			int del = delta[u];
			int d = dist[u];
			int num_in = 0;
			for(auto & e:g.inverted_adjacency[u]){
				int adjID = e.id;
				int from = g.all_edges[adjID].from;

				int dfrom = dist[from];
				if(edgeInShortestPathGraph[adjID])
					num_in++;
			}
			assert(del==num_in);
		}
#endif

	}

	void GRRDec(int edgeID){
		dbg_delta_lite();
		assert(!g.edgeEnabled(edgeID));
		//First, check if this edge is actually in the shortest path graph
		if(!edgeInShortestPathGraph[edgeID])
			return;
		edgeInShortestPathGraph[edgeID]=false;//remove this edge from the shortest path graph

		int ru =  g.all_edges[edgeID].from;
		int rv =  g.all_edges[edgeID].to;
		if(rv==63){
				int a=1;
			}
		assert(delta[rv]>0);
		delta[rv]--;
		if(delta[rv]>0)
			return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.

		q.clear();
		changeset.clear();
		changeset.push(rv);

		//find all effected nodes whose shortest path lengths may now be increased (or that may have become unreachable)
		for(int i = 0;i<changeset.size();i++){
			int u = changeset[i];
			if(u==63 || u==55){
					int a=1;
				}
			dist[u] = INF;
			for(auto & e:g.adjacency[u]){
				int adjID = e.id;
				if (g.edgeEnabled(adjID)){
					if(edgeInShortestPathGraph[adjID]){
						edgeInShortestPathGraph[adjID]=false;
						assert(g.all_edges[adjID].from==u);
						int s =g.all_edges[adjID].to;
						if(s==63 || s==55){
								int a=1;
							}
						assert(delta[s]>0);
						delta[s]--;
						if(delta[s]==0){
							changeset.push(s);
						}
					}
				}
			}
		}

		for (int i = 0;i<changeset.size();i++){
			int u = changeset[i];
			if(u==63 || u==55){
					int a=1;
				}
			assert(dist[u]==INF);
			for(auto & e:g.inverted_adjacency[u]){
				int adjID = e.id;

				if (g.edgeEnabled(adjID)){
					assert(g.all_edges[adjID].to==u);
					int v =g.all_edges[adjID].from;
					int w = 1;//assume a weight of one for now
					int alt = dist[v]+w;
					assert(!edgeInShortestPathGraph[adjID]);
					if(dist[u]>alt){
						dist[u]=alt;
					}
				}

			}
			if(dist[u]!=INF){
				//q.insert(u);
				//dbg_Q_add(q,u);
				q.insert(u);

				if(!reportDistance && reportPolarity>=0){
					if(!node_changed[u]){
						node_changed[u]=true;
						changed.push(u);
					}
				}
			}else if ( reportPolarity<=0){
				//have to mark this change even if we are reporting distanec, as u has not been added to the queue.
				if(!node_changed[u]){
					node_changed[u]=true;
					changed.push(u);
				}
			}
		}

		while(q.size()>0){
			int u = q.removeMin();
			if(u==53 ){
					int a=1;
				}
			if(reportDistance){
				if(dist[u]!=INF){
					if(reportPolarity>=0){
						if(!node_changed[u]){
							node_changed[u]=true;
							changed.push(u);
						}
					}
				}else if (reportPolarity<=0){
					if(!node_changed[u]){
						node_changed[u]=true;
						changed.push(u);
					}
				}
			}
			for(auto & e:g.adjacency[u]){
				int adjID = e.id;
				if (g.edgeEnabled(adjID)){
					assert(g.all_edges[adjID].from==u);
					int s =g.all_edges[adjID].to;
					if(s==63 || s==55){
							int a=1;
						}
					int w = 1;//assume a weight of one for now
					int alt = dist[u]+w;
					if(dist[s]>alt){
						dist[s]=alt;
						q.update(s);
					}else if (dist[s]==alt && !edgeInShortestPathGraph[adjID] ){
						edgeInShortestPathGraph[adjID]=true;
						delta[s]++;//added by sam... not sure if this is correct or not.
					}
				}
			}

			for(auto & e:g.inverted_adjacency[u]){
				int adjID = e.id;
				if (g.edgeEnabled(adjID)){

						assert(g.all_edges[adjID].to==u);
						int v =g.all_edges[adjID].from;
						int dv = dist[v];
						int du = dist[u];
						bool edgeIn = edgeInShortestPathGraph[adjID];
						int w = 1;//assume a weight of one for now
						if(dist[u]==dist[v]+w && ! edgeInShortestPathGraph[adjID]){
							assert(!edgeInShortestPathGraph[adjID]);
							edgeInShortestPathGraph[adjID]=true;
							delta[u]++;
						}else if (dist[u]<dist[v]+w &&  edgeInShortestPathGraph[adjID]){
							 edgeInShortestPathGraph[adjID]=false;
							 delta[u]--;
							assert(!edgeInShortestPathGraph[adjID]);
						}else if(dist[u]>dist[v]+w) {
							//assert(false);
						}
				}
			}
		}
		dbg_delta_lite();
	}




	void update( ){
#ifdef RECORD
		if(g.outfile){
			fprintf(g.outfile,"r %d\n", getSource());
		}
#endif
		static int iteration = 0;
		int local_it = ++iteration ;
		if(local_it==7668){
			int a =1;
		}
		if(last_modification>0 && g.modifications==last_modification)
				return;
		if(last_modification<=0 || g.changed()){
			INF=g.nodes+1;
			dist.growTo(g.nodes,INF);
			dist[getSource()]=0;
			delta.growTo(g.nodes);
			node_changed.growTo(g.nodes,true);

			for(int i = 0;i<g.nodes;i++){
				node_changed[i]=true;
				changed.push(i);//On the first round, report status of all nodes.
			}
		}
		edgeInShortestPathGraph.growTo(g.nEdgeIDs());
		if(last_history_clear!=g.historyclears){
			history_qhead=0;
			last_history_clear=g.historyclears;

		}
		double startdupdatetime = rtime(2);
		for (int i = history_qhead;i<g.history.size();i++){
			int edgeid = g.history[i].id;
			if(g.history[i].addition && g.edgeEnabled(edgeid)){
				GRRInc(edgeid);
			}else if (!g.history[i].addition &&  !g.edgeEnabled(edgeid)){
				GRRDec(edgeid);
			}
		}

		//for(int i = 0;i<g.nodes;i++){
		//	int u=i;
		for(int u:changed){
			//int u = changed[i];
			node_changed[u]=false;
			//CANNOT clear the change flag here, because if we backtrack and then immediately re-propagate an edge before calling theoryPropagate, the change to this node may be missed.

			if(reportPolarity<=0 && dist[u]>=INF){
				status.setReachable(u,false);
				status.setMininumDistance(u,dist[u]<INF,dist[u]);
			}else if (reportPolarity>=0 && dist[u]<INF){
				status.setReachable(u,true);
				status.setMininumDistance(u,dist[u]<INF,dist[u]);
			}
		}
		changed.clear();
		//}
		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;

		stats_full_update_time+=rtime(2)-startdupdatetime;;
	}


	bool dbg_path(int to){
#ifdef DEBUG_DIJKSTRA
		assert(connected(to));
		if(to == source){
			return true;
		}
		int p = previous(to);

		if(p<0){
			return false;
		}
		if(p==to){
			return false;
		}

		return dbg_path(p);


#endif
		return true;
	}
	bool dbg_uptodate(){
//#ifdef DEBUG_GRAPH
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
			return true;
		dbg_delta();
		Dijkstra<EdgeStatus,false> d(source,g);

		for(int i = 0;i<g.nodes;i++){
			int distance = dist[i];
			int dbgdist = d.distance(i);
			if(distance!=dbgdist){
				exit(4);
			}
		}
//#endif
#endif
		return true;
	}

	bool connected_unsafe(int t){
		dbg_uptodate();
		return t<dist.size() && dist[t]<INF;
	}
	bool connected_unchecked(int t){
		assert(last_modification==g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return dist[t]<INF;
	}
	int distance(int t){
		if(last_modification!=g.modifications)
					update();
		return dist[t];
	}
	int distance_unsafe(int t){
		if(connected_unsafe(t))
			return dist[t];
		else
			return INF;
	}
	int incomingEdge(int t){

		assert(false);//not yet implemented...
		assert(t>=0 && t<prev.size());
		assert(prev[t]>=-1 );
		return prev[t];
	}
	int previous(int t){
		if(prev[t]<0)
			return -1;

		assert(g.all_edges[incomingEdge(t)].to==t);
		return g.all_edges[incomingEdge(t)].from;
	}
};
};
#endif /* DIJKSTRA_H_ */
