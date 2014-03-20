
#ifndef DYNAMIC_CONNECTIVITY_H_
#define DYNAMIC_CONNECTIVITY_H_

#include "mtl/Vec.h"
#include "mtl/Heap.h"
#include "DynamicGraph.h"
#include "EulerTree.h"
#include "TreapCustom.h"
#include "core/Config.h"
#include "Reach.h"
using namespace Minisat;



template<class Status,class EdgeStatus=DefaultEdgeStatus>
class DynamicConnectivity:public Reach{
public:

	DynamicGraph<EdgeStatus> & g;
	Status &  status;
	int last_modification;
	int last_addition;
	int last_deletion;
	int history_qhead;

	int last_history_clear;

	int source;
	int INF;


	vec<int> q;
	vec<int> check;
	const int reportPolarity;

	//vec<char> old_seen;
	vec<char> seen;;
//	vec<int> changed;


	vec<int> prev;

	struct DefaultReachStatus{
			vec<bool> stat;
				void setReachable(int u, bool reachable){
					stat.growTo(u+1);
					stat[u]=reachable;
				}
				bool isReachable(int u) const{
					return stat[u];
				}
				DefaultReachStatus(){}
			};

	vec<int> trail;
	vec<int> levels;
	vec<int> trail_lim;
public:


	DynamicConnectivity(int s,DynamicGraph<EdgeStatus> & graph, Status & _status, int _reportPolarity=0 ):g(graph), status(_status), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(s),INF(0),reportPolarity(_reportPolarity){
		marked=false;
		mod_percentage=0.2;
		stats_full_updates=0;
		stats_fast_updates=0;
		stats_skip_deletes=0;
		stats_skipped_updates=0;
		stats_full_update_time=0;
		stats_fast_update_time=0;
		stats_num_skipable_deletions=0;
		stats_fast_failed_updates=0;
	}
	//Connectivity(const Connectivity& d):g(d.g), last_modification(-1),last_addition(-1),last_deletion(-1),history_qhead(0),last_history_clear(0),source(d.source),INF(0),mod_percentage(0.2),stats_full_updates(0),stats_fast_updates(0),stats_skip_deletes(0),stats_skipped_updates(0),stats_full_update_time(0),stats_fast_update_time(0){marked=false;};


	void setSource(int s){
		source = s;
		last_modification=-1;
		last_addition=-1;
		last_deletion=-1;
	}
	int getSource(){
		return source;
	}

	/*void updateFast(){
		stats_fast_updates++;
		double start_time = cpuTime();

		assert(last_deletion==g.deletions);
		last_modification=g.modifications;
		last_addition=g.additions;
		INF=g.nodes+1;
		seen.growTo(g.nodes);
		prev.growTo(g.nodes);

		if(lastaddlist!=g.addlistclears){
			addition_qhead=0;
			lastaddlist=g.addlistclears;
		}
		int start = q.size();
		//ok, now check if any of the added edges allow for new connectivity
		for (int i = addition_qhead;i<g.addition_list.size();i++){
			int u=g.addition_list[i].u;
			int v=g.addition_list[i].v;

			if(!seen[v]){
				q.push(v);
				seen[v]=1;
				prev[v]=u;
			}
		}
		addition_qhead=g.addition_list.size();

		for(int i = start;i<q.size();i++){
			int u = q[i];
			assert(seen[u]);
			for(int i = 0;i<g.adjacency[u].size();i++){
				int v = g.adjacency[u][i];

				if(!seen[v]){
					//this was changed
					changed.push(v);
					seen[v]=1;
					prev[v]=u;
					q.push(v);
				}
			}
		}
		stats_fast_update_time+=cpuTime()-start_time;
	}*/
/*	vec<int> & getChanged(){
		return changed;
	}
	void clearChanged(){
		changed.clear();
	}*/


	/*
	 * WARNING: THIS FUNDAMENTALLY WONT WORK if there are any cycles in the graph!
	 * inline void delete_update(int to){
		q.clear();
		q.push(to);
		seen[to]=0;
		//Is this really safe? check it very carefully, it could easily be wrong
		while(q.size()){
			int u = q.last();
			q.pop();
			assert(!seen[u]);
			for(int i = 0;i<g.inverted_adjacency[u].size();i++){
				int v = g.inverted_adjacency[u][i];
				if(seen[v]){
					seen[v]=1;
					//Then since to is still seen, we are up to date
					break;
				}
			}
			if(!seen[u]){
				for(int i = 0;i<g.adjacency[u].size();i++){
					int v = g.adjacency[u][i];
					if(seen[v] && prev[v]==to){
						seen[v]=0;
					}
				}
			}else{
#ifdef GRAPH_DEBUG
				for(int i = 0;i<g.adjacency[u].size();i++){
						int v = g.adjacency[u][i];
						assert(seen[v]);
				}
#endif
			}
		}
	}*/

	void setNodes(int n){
		q.capacity(n);
		check.capacity(n);
		seen.growTo(n);
		prev.growTo(n);
		INF=g.nodes+1;
	}

/*
	int decisionLevel(){
		return trail_lim.size()-1;
	}
	void cancelUntil(int level){
	    if (decisionLevel() > level){
	        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
	            int node = trail[c];
	            seen[node]=false;
	        }

	        trail.shrink(trail.size() - trail_lim[level]);
	        trail_lim.shrink(trail_lim.size() - level);
	    }
	}

*/
	vec<EulerTree::EulerVertex*> nodes;

	bool next(EulerTree::EulerVertex* node) {
	  TreapCustom::Node * n = node->node;
	  if(n) {
	    n = n->next;
	  }
	  while(n) {
	    if(n->value.type == "vertex") {
	      break;
	    }

	    n = n->next;
	  }
	  node->node = n;
	  return n;
	}

	bool hasNext(EulerTree::EulerVertex* node) {
		 TreapCustom::Node * n = node->node;
	  if(n) {
	    n = n.next
	  }
	  while(n) {
	    if(n->value.type == "vertex") {
	      break;
	    }

	    n = n->next;
	  }
	  return n;
	}

/*
	proto.prev = function() {
	  var n = this.node
	  if(n) {
	    n = n.prev
	  }
	  while(n) {
	    if(n.value.type === "vertex") {
	      break
	    }

	    n = n.prev
	  }
	  this.node = n
	  return !!n
	}

	proto.hasPrev = function() {
	  var n = this.node
	  if(n) {
	    n = n.prev
	  }
	  while(n) {
	    if(n.value.type === "vertex") {
	      break
	    }

	    n = n.prev
	  }
	  return !!n
	}
	*/

	int KEY_COUNTER = 0;
	struct DynamicEdge;
	struct DynamicVertex;
	EulerTree<DynamicVertex*> et;
	struct DynamicVertex {
	  int value;
	  vec<DynamicVertex*> adjacent;
	  vec<EulerTree<DynamicVertex*> ::EulerVertex*> euler;


	};

	struct DynamicEdge {
		 int value;
		 int key;
		 DynamicVertex* s;
		 DynamicVertex* t;
		 int edgeID;
		 int level;
		 vec<int> euler_half_edges;
	};
	//Raise the level of an edge, optionally inserting into higher level trees
	void raiseLevel(DynamicEdge* edge) {
	  DynamicVertex* s = edge->s;
	  DynamicVertex* t = edge->t;

	  //Update position in edge lists
	  removeEdge(s, edge);
	  removeEdge(t, edge);
	  edge->level += 1;
	  elist.insert(s.adjacent, edge);
	  elist.insert(t.adjacent, edge);

	  //Update flags for s
	  if(s->euler.size() <= edge->level) {
		s->euler.push(et.createVertex(s));
	  }
	  EulerTree<DynamicVertex*>::EulerVertex * vs = s->euler[edge->level];
	  et.setFlag(vs,true);

	  //Update flags for t
	  if(t->euler.size() <= edge->level) {
		t->euler.push(et.createVertex(t));
	  }
	  EulerTree<DynamicVertex*>::EulerVertex * vt = t->euler[edge->level];
	  et.setFlag(vt, true);

	  //Relink if necessary
	  if(edge->euler) {
		edge->euler.push(et.link(vs,vt,edge->edgeID, edge));
	  }
	}

	//Remove edge from list and update flags
	void removeEdge(DynamicVertex *vertex,DynamicEdge* edge) {
	 // var adj = vertex->adjacent;
	//  var idx = elist->index(adj, edge);
	 // adj.splice(idx, 1);
	  //there has got to be a better way to do this...
		vertex->adjacent.remove(edge);
	  //Check if flag needs to be updated
	  if(!((idx < adj.length && adj[idx].level == edge.level) ||
		   (idx > 0 && adj[idx-1].level == edge.level))) {
		vertex->euler[edge->level]->setFlag(false);
	  }
	}

	//Add an edge to all spanning forests with level <= edge.level
	void link(DynamicEdge* edge) {
	  vec<DynamicVertex*> & vs = edge->s->euler;
	  vec<DynamicVertex*> & vt = edge->t->euler;
	  //var euler = new Array(edge.level+1);
	  for(var i=0; i<edge->level+1; ++i) {
		if(vs.size() <= i) {
		  vs.push(et.createVertex (edge->s));
		}
		if(vt.size() <= i) {
		  vt.push(et.createVertex(edge->t));
		}
		edge->euler.push( link( vs[i], vt[i], edge));
	  }
	  //edge.euler = euler;
	}
	//Search over tv for edge connecting to tw
	  bool visit(TreapCustom::Node *node) {
		if(node->flag) {
		  int v = node->value.value;
		  var adj = v.adjacent;
		  for(var ptr=elist.level(adj, level); ptr<adj.length && adj[ptr].level === level; ++ptr) {
			var e = adj[ptr];
			var es = e.s;
			var et = e.t;
			if(es.euler[level].path(et.euler[level])) {
			  raiseLevel(e);
			  ptr -= 1;
			} else {
			  //Found the edge, relink components
			  link(e);
			  return true;
			}
		  }
		}
		if(node->left && node->left->flagAggregate) {
		  if(visit(node->left)) {
			return true;
		  }
		}
		if(node->right && node->right->flagAggregate) {
		  if(visit(node->right)) {
			return true;
		  }
		}
		return false;
	  }
	void cut(DynamicEdge* edge) {
		  int level;

		  //Don't double cut an edge
		  if(!edge->s) {
		    return;
		  }



		  removeEdge(edge->s, edge);
		  removeEdge(edge->t, edge);
		  if(edge->euler) {
		    //Cut edge from tree
		    for(var i=0; i<edge->euler.size(); ++i) {
		      cut(edge->euler[i]);
		    }

		    //Find replacement, looping over levels
		    for(var i=edge->level; i>=0; --i) {
		      TreapCustom::Node * tv = et.t.findRoot(edge->s.euler[i]->node);
		      TreapCustom::Node * tw =  et.t.findRoot(edge->t.euler[i]->node);
		      level = i;
		      if(tv->count > tw->count) {
		        visit(tw);
		      } else {
		        visit(tv);
		      }
		    }
		  }
		  edge->s = NULL;
		  edge->t =NULL;
		  edge->euler.clear();
		  edge->level=32;
		}

		bool connected(DynamicVertex * a, DynamicVertex* b) {
		  return path(a->euler[0],b->euler[0]);
		}


		DynamicEdge * link(DynamicVertex * node, DynamicVertex * other, int value) {
		  DynamicEdge* e = new DynamicEdge(value, (KEY_COUNTER++), this, other, 0, NULL);
		  if(!path( node->euler[0],other->euler[0])) {
		    link(e);
		  }
		  et.setFlag( node->euler[0],true);
		  et.setFlag( other->euler[0],true);
		  elist->insert(node->adjacent, e);
		  elist->insert(other->adjacent, e);
		  return e;
		}

		//Returns the number of vertices in this connected component
		int componentSize(DynamicVertex * node) {
		  return count(node->euler[0]);
		}

		//Removes the vertex from the graph
			void cut(DynamicVertex * v) {
			  while(v->adjacent.size() > 0) {
			   cut( v->adjacent.last());
			  }
			}

			int component(DynamicVertex * node) {
				  return node->euler[0]->node;//should convert this to an integer...
			}
			DynamicVertex * createVertex(int value) {
				  //var euler = [null]
				DynamicVertex * v = new DynamicVertex();
				v->value= value;
				v->euler.push(et.createVertex (v));
				  return v;
				}
	void update( ){
		static int iteration = 0;
		int local_it = ++iteration ;

		if(last_modification>0 && g.modifications==last_modification){
			stats_skipped_updates++;
			return;
		}

		if(last_deletion==g.deletions){
			stats_num_skipable_deletions++;
		}

		setNodes(g.nodes);

		if(reportPolarity<1){
			for(int u = 0;u<g.nodes;u++){
				if(!seen[u]){
					status.setReachable(u,false);
				}else if(reportPolarity==0){
					status.setReachable(u,true);
				}
			}
		}
		assert(dbg_uptodate());

		last_modification=g.modifications;
		last_deletion = g.deletions;
		last_addition=g.additions;

		history_qhead=g.history.size();
		last_history_clear=g.historyclears;
	}

	bool dbg_path(int to){
#ifdef DEBUG_DIJKSTRA
		assert(connected(to));
		if(to == source){
			return true;
		}
		int p = prev[to];

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
	void drawFull(){
				printf("digraph{\n");
				for(int i = 0;i< g.nodes;i++){

					if(seen[i]){
						printf("n%d [fillcolor=blue style=filled]\n", i);
					}else{
						printf("n%d \n", i);
					}


				}

				for(int i = 0;i< g.adjacency.size();i++){
					for(int j =0;j<g.adjacency[i].size();j++){
					int id  =g.adjacency[i][j].id;
					int u =  g.adjacency[i][j].node;
					const char * s = "black";
					if( g.edgeEnabled(id))
						s="blue";
					else
						s="red";



					printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i,u, id, s);
					}
				}

				printf("}\n");
			}
	bool dbg_uptodate(){
#ifdef DEBUG_DIJKSTRA
		if(last_modification<=0)
			return true;
		Dijkstra<EdgeStatus> d(source,g);
		d.update();
		//drawFull();
		for(int i = 0;i<g.nodes;i++){

			int dbgdist = d.dist[i];
			if(!seen[i])
				assert(dbgdist==d.INF  );
			else{

				if(!(dbgdist<d.INF)){
					drawFull();
				}
				assert(dbgdist<d.INF);
			}
		}
#endif
		return true;
	}

	bool connected_unsafe(int t)const{
		return t<seen.size() && seen[t];
	}
	bool connected_unchecked(int t)const{
		assert(last_modification==g.modifications);
		return connected_unsafe(t);
	}
	bool connected(int t){
		if(last_modification!=g.modifications)
			update();

		assert(dbg_uptodate());

		return seen[t];
	}
	int distance(int t){
		if(connected(t))
			return 1;
		else
			return INF;
	}
	int distance_unsafe(int t){
		if(connected_unsafe(t))
			return 1;
		else
			return INF;
	}
	int previous(int t){
		assert(t>=0 && t<prev.size());
		assert(prev[t]>=-1 && prev[t]<prev.size());
		return prev[t];
	}

};

#endif
