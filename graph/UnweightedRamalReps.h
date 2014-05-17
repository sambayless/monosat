/*
 * Dijkstra.h
 *
 *  Created on: 2013-05-28
 *      Author: sam
 */

#ifndef UNWEIGHTED_RAMAL_REPS_H_
#define UNWEIGHTED_RAMAL_REPS_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "Reach.h"
#include "Dijkstra.h"
#include "core/Config.h"
#include "mtl/Sort.h"
namespace Minisat{
template<class Status,class EdgeStatus=DefaultEdgeStatus>
class UnweightedRamalReps:public Reach{
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

	int maxDistance;



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
	//Heap<DistCmp> q;
	vec<bool> in_queue;
	vec<bool> in_queue2;
	vec<int> q;
	vec<int> q2;
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
	UnweightedRamalReps(int s,DynamicGraph<EdgeStatus> & graph,	Status &  status, int reportPolarity=0,bool reportDistance=true):g(graph),status(status),reportPolarity(reportPolarity),reportDistance(reportDistance), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(-1){
		maxDistance=-1;
		mod_percentage=0.2;
		stats_full_updates=0;
		stats_fast_updates=0;
		stats_skip_deletes=0;
		stats_skipped_updates=0;
		stats_full_update_time=0;
		stats_fast_update_time=0;
	}
	//Dijkstra(const Dijkstra& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),q(DistCmp(dist)),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};
	void setMaxDistance(int _maxDistance){
		if(_maxDistance != maxDistance){
			last_modification=-1;//force the next update to recompute from scratch
			if(_maxDistance<0){
				maxDistance=INF;
			}else
				maxDistance=_maxDistance;
		}
	}

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
				if(edgeInShortestPathGraph[i]){
					exit(3);
				}
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
					if(maxDistance>=0 && alt>maxDistance)
						alt=INF;
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
				if(maxDistance>=0 && alt>maxDistance)
					alt=INF;
				if(alt<dbg_dist[v]){

					dbg_dist[v]=alt;

					if(!q.inHeap(v))
						q.insert(v);
					else
						q.decrease(v);
				}
			}
		}

		for(int u = 0;u<g.nodes;u++){
			int d = dist[u];

			int db = dbg_dist[u];
			assert(dbg_dist[u]==dist[u]);

			for(int i = 0;i<g.inverted_adjacency[u].size();i++){
				if(!g.edgeEnabled( g.inverted_adjacency[u][i].id))
					continue;

				int edgeID = g.inverted_adjacency[u][i].id;
				int v = g.all_edges[edgeID].from;
				int alt = dbg_dist[v]+ 1;
				int du = dbg_dist[u];
				if(maxDistance>=0 && alt>maxDistance)
					alt=INF;
				assert(alt>=dbg_dist[u]);
				if(alt==dbg_dist[u] && alt<INF){
					dbg_delta[u]++;
					assert(edgeInShortestPathGraph[edgeID]);
				}else if (alt<INF){
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
			int du = dist[u];
			if(dist[u]<INF){
				int d = delta[u];
				int d_expect = dbg_delta[u];
				assert(d==dbg_delta[u]);
			}
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

		int rdv = dist[rv];
		int rdu = dist[ru];

		int weight = 1;
		int altw = dist[ru]+weight;
		if(altw>maxDistance)
			altw=INF;
		if(dist[rv]<altw)
			return;
		else if (dist[rv]==altw && dist[rv]<INF ){
			assert(!edgeInShortestPathGraph[edgeID]);
			edgeInShortestPathGraph[edgeID]=true;
			delta[rv]++;//we have found an alternative shortest path to v
			return;
		}else if (altw==INF){
			return;//don't do anything
		}
		assert(altw<INF);
		edgeInShortestPathGraph[edgeID]=true;
		delta[rv]++;
		dist[rv]=altw;

		q.clear();
		in_queue.clear();in_queue.growTo(g.nodes);
		in_queue2.growTo(g.nodes);
		q.push(rv);
		in_queue[rv]=true;

		for(int i = 0;i<q.size();i++){
			int u =q[i];
			dbg_Q_order(q);
			assert(dist[u]<INF);
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
						int alt = dist[v]+w;
						if(alt>maxDistance)
							alt=INF;
						if(dist[u]==alt ){
							assert(alt<INF);
							edgeInShortestPathGraph[adjID]=true;
							delta[u]++;
						}else if (dist[u]<alt){
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
						int alt = dist[u]+w;
						if(alt>maxDistance)
							alt=INF;
						if(dist[s]>alt){
							dist[s]=alt;
							if(!in_queue[s]){
							//q.update(s);
								dbg_Q_add(q,s);
								q.push(s);
								in_queue[s]=true;
							}
							dbg_not_seen_q(q,s,i);
						}else if (dist[s]==alt && !edgeInShortestPathGraph[adjID] && alt<INF){
							assert(alt<INF);
							edgeInShortestPathGraph[adjID]=true;
							delta[s]++;
						}
				}
			}
		}
		dbg_delta_lite();
	}
	void dbg_not_seen_q(vec<int> & q, int u, int from){
		bool found = false;
		for(int i = from;i<q.size();i++){
			if(q[i]==u){
				found=true;
				break;
			}
		}
		assert(found);
	}
	void dbg_Q_add(vec<int> & q,int u){
#ifndef NDEBUG
		//assert(!in_queue[u]);
		for(int v:q){
			assert(u!=v);
			assert(dist[v]<=dist[u]);
			//assert(in_queue[v]);
		}

#endif
	}
	void dbg_Q_order(vec<int> & _q){
#ifndef NDEBUG

		for(int i = 1;i<_q.size();i++){
			int v = _q[i];
			int u = _q[i-1];

			if(&_q==&q){
			assert(in_queue[v]);
			assert(in_queue[u]);
			if(!(in_queue2[u] || in_queue2[v])){

				assert(dist[u]<=dist[v]);
			}
			}else{
				assert(in_queue2[u]);
				assert(in_queue2[v]);
				assert(dist[u]<=dist[v]);
			}
		}

#endif
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
				if(edgeInShortestPathGraph[adjID]){

					assert(dist[u]<INF);
					num_in++;
				}
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

		assert(delta[rv]>0);
		delta[rv]--;
		if(delta[rv]>0)
			return; //the shortest path hasn't changed in length, because there was an alternate route of the same length to this node.

		in_queue.clear();
		in_queue.growTo(g.nodes);
		in_queue2.clear();
		in_queue2.growTo(g.nodes);
		q.clear();
		q2.clear();

		changeset.clear();
		changeset.push(rv);

		//find all effected nodes whose shortest path lengths may now be increased (or that may have become unreachable)
		for(int i = 0;i<changeset.size();i++){
			int u = changeset[i];

			dist[u] = INF;
			for(auto & e:g.adjacency[u]){
				int adjID = e.id;
				if (g.edgeEnabled(adjID)){
					if(edgeInShortestPathGraph[adjID]){
						edgeInShortestPathGraph[adjID]=false;
						assert(g.all_edges[adjID].from==u);
						int s =g.all_edges[adjID].to;

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

			assert(dist[u]==INF);
			for(auto & e:g.inverted_adjacency[u]){
				int adjID = e.id;

				if (g.edgeEnabled(adjID)){
					assert(g.all_edges[adjID].to==u);
					int v =g.all_edges[adjID].from;
					int w = 1;//assume a weight of one for now
					int alt = dist[v]+w;
					if(alt>maxDistance)
						alt=INF;
					assert(!edgeInShortestPathGraph[adjID]);
					if(dist[u]>alt){
						dist[u]=alt;
					}
				}

			}

			if(dist[u]<INF){
				//q.insert(u);
				//dbg_Q_add(q,u);
				q.push(u);
				in_queue[u]=true;

				if( reportPolarity>=0){
					if(!node_changed[u]){
						node_changed[u]=true;
						changed.push(u);
					}
				}
			}else if (reportPolarity<=0){
				//call this even if we are reporting distance, because u hasn't been placed in the queue!
				if(!node_changed[u]){
					node_changed[u]=true;
					changed.push(u);
				}
			}

		}
		sort(q,DistCmp(dist));
		int i=0,j=0;
		while(i<q.size()||j<q2.size()){
				int u;
				if(i==q.size()){
					assert(j<q2.size());
					u=q2[j++];
				}else if (j==q2.size()){
					assert(i<q.size());
					u=q[i++];
					if(in_queue2[u]){
						continue;
					}
				}else if(dist[q[i]]<dist[q2[j]]){
					u = q[i++];
					if(in_queue2[u]){
						continue;
					}
				}else{
					assert(dist[q2[j]]<=dist[q[i]]);
					u=q2[j++];
				}
				assert(dist[u]<INF);
				if(reportDistance){
					if(reportPolarity>=0){
						if(!node_changed[u]){
							node_changed[u]=true;
							changed.push(u);
						}
					}
				}
				dbg_Q_order(q);
				dbg_Q_order(q2);

				for(auto & e:g.adjacency[u]){
					int adjID = e.id;
					if (g.edgeEnabled(adjID)){
						assert(g.all_edges[adjID].from==u);
						int s =g.all_edges[adjID].to;
						int w = 1;//assume a weight of one for now
						int alt = dist[u]+w;
						if(alt>maxDistance)
							alt=INF;
						if(dist[s]>alt){
							assert(alt<INF);
							if(reportPolarity>=0 && dist[s]>=0){
								//This check is needed (in addition to the above), because even if we are NOT reporting distances, it is possible for a node that was previously not reachable
								//to become reachable here. This is ONLY possible because we are batching multiple edge incs/decs at once (otherwise it would be impossible for removing an edge to decrease the distance to a node).
								if(!node_changed[s]){
									node_changed[s]=true;
									changed.push(s);
								}
							}
							dist[s]=alt;
							if(!in_queue2[s]){
								dbg_Q_add(q2,s);
								q2.push(s);
								in_queue2[s]=true;
							}

							dbg_Q_order(q2);
							//dbg_not_seen_q(q2,s,j);
						}else if (dist[s]==alt && !edgeInShortestPathGraph[adjID] && dist[s]<INF ){
							assert(dist[s]<INF);
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
							int alt = dist[v]+w;
							if(alt>maxDistance)
								alt=INF;
							if(dist[u]==alt && ! edgeInShortestPathGraph[adjID]){
								assert(!edgeInShortestPathGraph[adjID]);
								assert(alt<INF);
								edgeInShortestPathGraph[adjID]=true;
								delta[u]++;
							}else if (dist[u]<alt &&  edgeInShortestPathGraph[adjID]){
								 edgeInShortestPathGraph[adjID]=false;
								 delta[u]--;
								assert(!edgeInShortestPathGraph[adjID]);
							}else if(dist[u]>alt) {
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
		if(local_it==671){
			int a =1;
		}
		if(last_modification>0 && g.modifications==last_modification)
				return;
		if(last_modification<=0 || g.changed()){//Note for the future: there is probably room to improve this further.
			stats_full_updates++;
			INF=g.nodes+1;
			dist.growTo(g.nodes,INF);
			dist[getSource()]=0;
			delta.growTo(g.nodes);
			edgeInShortestPathGraph.growTo(g.nEdgeIDs());
			node_changed.growTo(g.nodes);
			changed.clear();
			if(maxDistance<0)
				maxDistance=INF;
			for(int i = 0;i<g.nodes;i++){
				if((dist[i]>=INF && reportPolarity<=0) || (dist[i]<INF && reportPolarity>=0)){
				node_changed[i]=true;
				changed.push(i);//On the first round, report status of all nodes.
				}
			}
		}

		if(last_history_clear!=g.historyclears){
			history_qhead=0;
			last_history_clear=g.historyclears;
			for(int edgeid = 0;edgeid<g.edges;edgeid++){
				if(g.edgeEnabled(edgeid)){
					GRRInc(edgeid);
				}else{
					GRRDec(edgeid);
				}
			}
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
		for(int u:changed){
			//int u=i;
			//int u = changed[i];
			node_changed[u]=false;

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
		dbg_delta();

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
		assert(dbg_uptodate());
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
		if(last_modification<0)
			return true;
		dbg_delta();
		Dijkstra<NullReachStatus,EdgeStatus,false> d(source,g);

		for(int i = 0;i<g.nodes;i++){
			int distance = dist[i];
			if(maxDistance>=0 && distance>maxDistance){
				distance=INF;
			}
			int dbgdist = d.distance(i);
			if(maxDistance>=0 && dbgdist>maxDistance){
				dbgdist=INF;
			}
			if(distance!=dbgdist){
				assert(false);
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
