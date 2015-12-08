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

#ifndef SPIRA_PAN_H_
#define SPIRA_PAN_H_

#include <vector>
#include "alg/Heap.h"
#include "mtl/Sort.h"
#include "DynamicGraph.h"
#include "core/Config.h"
#include "MinimumSpanningTree.h"
#include "Kruskal.h"
#include <algorithm>
#include <limits>
#include <iostream>

namespace dgl {
/**
 * This is an implementation of Spira and Pan's (1975) dynamic minimum spanning tree algorithm.
 * Unlike in the original (1975) version, we initialize the MST in new graphs using Prim's, and also use Prim's to connect separated components back
 * together after a string of edge deletions (the original paper only considers a single edge deletion at a time, which
 * is inefficient if multiple edges are deleted at once).
 */
template<class Status, typename Weight = int>
class SpiraPan: public MinimumSpanningTree<Weight>, public DynamicGraphAlgorithm {
public:
	
	DynamicGraph<Weight> & g;

	Status & status;
	int last_modification=-1;
	Weight min_weight=0;
	int last_addition=0;
	int last_deletion=0;
	int history_qhead=0;

	int last_history_clear=0;
	int alg_id=-1;
	Weight INF;

	std::vector<int> mst;
	std::vector<int> q;
	std::vector<int> check;
	const int reportPolarity;

	//std::vector<char> old_seen;
	std::vector<bool> in_tree;
	std::vector<bool> keep_in_tree;
	std::vector<int> parents;
	std::vector<int> parent_edges;
	std::vector<int> components_to_visit;
	std::vector<bool> component_needs_visit;
	std::vector<int> component_member; //pointer to one arbitrary member of each non-empty component
	std::vector<Weight> component_edge_weight;
	std::vector<bool> edge_enabled;
	struct VertLt {
		const std::vector<Weight>& keys;

		bool operator ()(int x, int y) const {
			return keys[x] < keys[y];
		}
		VertLt(const std::vector<Weight>& _key) :
				keys(_key) {
		}
	};

	Heap<VertLt> Q;

	std::vector<bool> seen;

	int num_sets = 0;

	std::vector<int> edge_to_component;
	std::vector<int> empty_components; //list of component ids with no member nodes
	std::vector<int> components;

#ifndef NDEBUG
	Kruskal<typename MinimumSpanningTree<Weight>::NullStatus, Weight> dbg;
#endif
public:
	
	int stats_full_updates=0;
	int stats_fast_updates=0;
	int stats_fast_failed_updates=0;
	int stats_skip_deletes=0;
	int stats_skipped_updates=0;
	int stats_num_skipable_deletions=0;
	double mod_percentage=0;

	double stats_full_update_time=0;
	double stats_fast_update_time=0;

	SpiraPan(DynamicGraph<Weight> & graph,  Status & status, int reportPolarity = 0) :
			g(graph),  status(status), last_modification(-1), last_addition(-1), last_deletion(-1), history_qhead(
					0), last_history_clear(0), INF(0), reportPolarity(reportPolarity), Q(VertLt(component_edge_weight))
#ifndef NDEBUG
					, dbg(g,  MinimumSpanningTree<Weight>::nullStatus, 0)
#endif
	{
		alg_id=g.addDynamicAlgorithm(this);
		mod_percentage = 0.2;

		min_weight = -1;
		
	}
	
	void setNodes(int n) {
		q.reserve(n);
		check.reserve(n);
		in_tree.resize(g.nEdgeIDs());
		seen.resize(n);
		//INF=std::numeric_limits<int>::max();
		while (component_edge_weight.size() <= g.nodes()) {
			component_edge_weight.push_back(INF);
		}
		
		parents.resize(n, -1);
		edge_to_component.resize(n, -1);
		
		parent_edges.resize(n, -1);
	}

	void dbg_printSpanningTree(bool showWeights=true){
#ifndef NDEBUG

#ifndef NDEBUG
		printf("graph{\n");
	/*	for (int i = 0; i < g.nodes(); i++) {
			printf("n%d\n", i);
		}*/

		for (int i = 0; i < g.adjacency_list.size(); i++) {
			for (int j = 0; j < g.adjacency_list[i].size(); j++) {
				int id = g.adjacency_list[i][j].id;

				int u = g.adjacency_list[i][j].node;
				const char * s = "black";
				if(in_tree[id]){
					s="green";
					assert(g.edgeEnabled(id));
				}else if (g.edgeEnabled(id))
					s = "red";
				else
					s = "blue";
				if(showWeights){
					std::stringstream ss;
					ss<<g.getWeight(id);
					printf("n%d -- n%d [label=\"v%d w=%s\",color=\"%s\"]\n", i,u, id,ss.str().c_str(), s);
				}else{
					printf("n%d -- n%d [label=\"v%d\",color=\"%s\"]\n", i, u, id, s);
				}
			}
		}
		printf("}\n");
#endif


#endif

	}

	/*
	 * 	//chin houck insertion
	 void insert(int newNode){
	 if(g.adjacency_undirected[newNode].size()==0)
	 return;
	 marked.clear();
	 marked.resize(g.nodes());
	 incident_edges.clear();
	 incident_edges.resize(g.nodes(),-1);
	 for(auto & edge:g.adjacency_undirected[newNode]){
	 incident_edges[edge.node]=edge.id;
	 }
	 insert(g.adjacency_undirected[newNode][0]);
	 }

	 //Insert a node into an _existing_ minimum spanning tree
	 void insert(int r, int z, int & t){
	 int m = incident_edges[r];
	 marked[r]=true;

	 for(auto edge:g.adjacency_undirected[r]){//IF this is a dense graph, then this loop is highly sub-optimal!
	 //for each incident edge in the old MST:
	 if(in_tree[edge.id] && !marked[ edge.node]){
	 insert(edge.node,z,t);
	 assert(dbg_is_largest_edge_on_path(m,r,z));
	 assert(dbg_is_largest_edge_on_path(t,edge.node,z));
	 int k = t;
	 if(t<0 || g.g.getWeight(edge.id)>g.weights[t])
	 k=edge.id;
	 int h = t;
	 if(t<0 || g.g.getWeight(edge.id)<g.weights[t])
	 h=edge.id;
	 if(h>=0){
	 keep_in_tree[h]=true;
	 }
	 if(m<0 || (k>=0 && g.weights[k]<g.weights[m])){
	 m=k;
	 }
	 }
	 }
	 t=m;
	 }*/

	bool dbg_is_largest_edge_on_path(int edge, int from, int to) {
#ifndef NDEBUG
		
#endif
		return true;
	}
	
	void dbg_parents() {
#ifndef NDEBUG
		//check that the parents don't cycle
		for (int i = 0; i < g.nodes(); i++) {
			int p = i;
			int num_parents = 0;
			while (p != -1) {
				num_parents++;
				if (parents[p] > -1) {
					int e = parent_edges[p];
					int u =g.getEdge(e).from;
					int v = g.getEdge(e).to;
					assert(u == p || v == p);
					assert(u == parents[p] || v == parents[p]);
					assert(in_tree[e]);
					int compi = components[i];
					int compp = components[p];
					assert(components[p] == components[i]);
					int parent = parents[p];
					int pc = components[parents[p]];
					assert(components[parents[p]] == components[i]);
				} else {
					assert(parent_edges[p] == -1);
				}
				
				p = parents[p];
				assert(num_parents <= g.nodes());
				
			}
			
		}
		
		std::vector<int> used_components;
		for (int i = 0; i < g.nodes(); i++) {
			if (!std::count(used_components.begin(), used_components.end(), components[i])) {
				used_components.push_back(components[i]);
				
			}
		}
		
		for (int c : used_components) {
			assert(!std::count(empty_components.begin(), empty_components.end(), c));
			assert(component_member[c] != -1);
			int m = component_member[c];
			assert(components[m] == c);
		}
		
		for (int c : empty_components) {
			//assert(!used_components.contains(c));
			assert(!std::count(used_components.begin(), used_components.end(), c));
			assert(component_member[c] == -1);
		}
		
#endif
	}
	
	void addEdgeToMST(int edgeid) {
		int u = g.getEdge(edgeid).from;
		int v = g.getEdge(edgeid).to;
		
		dbg_parents();
		if (components[u] != components[v]) {
			
			//If u,v are in separate components, then this edge must be in the mst (and we need to fix the component markings)


			int higher_component = v;
			int lower_component = u;
			int new_c = components[u];
			int old_c = components[v];
			if (new_c > old_c) {
				std::swap(higher_component, lower_component);
				std::swap(new_c, old_c);
			}
			
			if(component_needs_visit[old_c] || component_needs_visit[new_c]){
				return;//let prims connect these two components.
			}

			if (component_needs_visit[old_c]) {
#ifndef NDEBUG
				bool found = false;
				for (int i : components_to_visit) {
					if (i == old_c)
						found = true;
				}
				assert(found);
#endif
				if (!component_needs_visit[new_c]) {
					component_needs_visit[new_c] = true;
					components_to_visit.push_back(new_c);
				}else{

				}
				//component_needs_visit[old_c]=false; //dont do this, because the component has not yet been removed from the components to visit list.
			}
			num_sets--;
			in_tree[edgeid] = true;
			//ok, now set every node in the higher component to be in the lower component with a simple dfs.
			//fix the parents at the same time.
			min_weight += g.getWeight(edgeid);
			assert(components[lower_component] == new_c);
			assert(components[higher_component] == old_c);
			components[higher_component] = new_c;
			parents[higher_component] = lower_component;
			parent_edges[higher_component] = edgeid;
			q.clear();
			q.push_back(higher_component);
			while (q.size()) {
				int n = q.back();
				q.pop_back();
				for (int i = 0; i < g.nIncident(n, true); i++) {
					auto & edge = g.incident(n, i, true);
					if (in_tree[edge.id]) {
						
						int t = edge.node;
						if (components[t] == old_c) {
							components[t] = new_c;
							parents[t] = n;
							parent_edges[t] = edge.id;
							q.push_back(t);
						}
						
					}
				}
			}
			component_member[old_c] = -1;
			empty_components.push_back(old_c);
			dbg_parents();
		} else {
			dbg_parents();
			assert(components[u] == components[v]);
			if (parents[u] == v || parents[v] == u) {
				//If there is already another edge (u,v) that is in the tree, then we at most need to swap that edge out for this one.
				//note that this can only be the case if u is the parent of v or vice versa
				int p_edge = parent_edges[u];
				if (parents[v] == u) {
					p_edge = parent_edges[v];
				}
				if ( g.getWeight(p_edge) > g.getWeight(edgeid)) {
					//then swap these edges without changing anything else.
					in_tree[p_edge] = false;
					in_tree[edgeid] = true;
					Weight delta = g.getWeight(p_edge) - g.getWeight(edgeid);
					min_weight -= delta;
					if (parents[v] == u) {
						assert(parent_edges[v] == p_edge);
						parent_edges[v] = edgeid;
					} else {
						assert(parent_edges[u] == p_edge);
						parent_edges[u] = edgeid;
					}
				}
				dbg_parents();
			} else {
				//otherwise, find the cycle induced by adding this edge into the MST (by walking up the tree to find the LCA - if we are doing many insertions, could we swap this out for tarjan's OLCA?).
				int p = u;
				
				while (p > -1) {
					seen[p] = true;
					p = parents[p];
				}
				Weight max_edge_weight = -1;
				int max_edge = -1;
				bool edge_on_left = false;	//records which branch of the tree rooted at p the edge we are replacing is
				p = v;
				while (true) {
					assert(p > -1);	//u and v must share a parent, because u and v are in the same connected component and we have already computed the mst.
					if (seen[p]) {
						break;
					} else {
						if (parents[p] > -1 && g.getWeight(parent_edges[p]) > max_edge_weight) {
							assert(parent_edges[p] > -1);
							max_edge_weight = g.getWeight(parent_edges[p]);
							max_edge = parent_edges[p];
						}
						p = parents[p];
					}
				}
				assert(seen[p]);
				int lca = p;				//this is the lowest common parent of u and v.
				p = u;
				while (p != lca) {
					assert(seen[p]);
					seen[p] = false;
					if (parents[p] > -1 && g.getWeight(parent_edges[p]) > max_edge_weight) {
						assert(parent_edges[p] > -1);
						max_edge_weight = g.getWeight(parent_edges[p]);
						max_edge = parent_edges[p];
						edge_on_left = true;
					}
					p = parents[p];
				}
				
				//reset remaining 'seen' vars
				assert(p == lca);
				while (p > -1) {
					assert(seen[p]);
					seen[p] = false;
					p = parents[p];
				}
				assert(max_edge > -1);
				if (max_edge_weight > g.getWeight(edgeid)) {
					//then swap out that edge with the new edge.
					//this will also require us to repair parent edges in the cycle
					min_weight -= max_edge_weight;
					min_weight += g.getWeight(edgeid);
					in_tree[edgeid] = true;
					in_tree[max_edge] = false;
					
					int last_p;
					if (edge_on_left) {
						p = u;
						last_p = v;
					} else {
						p = v;
						last_p = u;
					}
					
					int last_edge = edgeid;
					while (parent_edges[p] != max_edge) {
						assert(p > -1);
						int next_p = parents[p];
						std::swap(parent_edges[p], last_edge);					//re-orient the parents.
						parents[p] = last_p;
						last_p = p;
						p = next_p;
					}
					assert(parent_edges[p] == max_edge);
					std::swap(parent_edges[p], last_edge);					//re-orient the parents.
					parents[p] = last_p;
				}
				dbg_parents();
			}
			
		}
		dbg_parents();
	}
	
	void removeEdgeFromMST(int edgeid) {
		dbg_parents();
		if (!in_tree[edgeid]) {
			
			//If an edge is disabled that was NOT in the MST, then no update is required.
		} else {
			//this is the 'tricky' case for dynamic mst.
			//following Spira & Pan, each removed edge splits the spanning tree into separate components that are MST's for those components.
			//after all edges that will be removed have been removed, we will run Prim's to stitch those components back together, if possible.
			
			int u = g.getEdge(edgeid).from;
			int v = g.getEdge(edgeid).to;
			
			in_tree[edgeid] = false;
			min_weight -= g.getWeight(edgeid);
			num_sets++;
			
			assert(components[u] == components[v]);
			assert(parents[u] == v || parents[v] == u);
			
			if (parents[u] == v) {
				parents[u] = -1;
				parent_edges[u] = -1;
			} else {
				parents[v] = -1;
				parent_edges[v] = -1;
			}
			int start_node = v;
			//if we want to maintain the guarantee components are always assigned the lowest node number that they contain, we'd need to modify the code below a bit.
			assert(empty_components.size());
			int new_c = empty_components.back();
			empty_components.pop_back();
			int old_c = components[u];
			assert(component_member[new_c] == -1);
			component_member[new_c] = start_node;
			component_member[old_c] = u;
			components_to_visit.push_back(new_c);
			component_needs_visit[new_c] = true;
			components[start_node] = new_c;
			assert(q.size() == 0);
			//relabel the components of the tree that has been split off.
			q.clear();
			q.push_back(start_node);
			while (q.size()) {
				int n = q.back();
				q.pop_back();
				//for(auto & edge:g.adjacency_undirected[n]){
				for (int i = 0; i < g.nIncident(n, true); i++) {
					auto & edge = g.incident(n, i, true);
					if (in_tree[edge.id]) {
						//assert(g.edgeEnabled(edge.id));
						int t = edge.node;
						if (components[t] == old_c) {
							components[t] = new_c;
							q.push_back(t);
						}
					}
				}
			}
			
		}
		dbg_parents();
	}
	
	void prims() {
		dbg_parents();
		//component_weight.clear();
		
		for (int i = 0; i < components_to_visit.size(); i++) {
			int c = components_to_visit[i];
			int start_node = component_member[c];
			
			if (start_node == -1) {
				//then this component has already been merged into another one, no need to visit it.
				continue;
			}
			component_needs_visit[c] = false;			//because we just visited this component
			assert(c >= 0);
#ifndef NDEBUG
			for (auto w : component_edge_weight) {
				assert(w == INF);
			}
#endif
			//ok, try to connect u's component to v
			//ideally, we'd use the smaller of these two components...
			int smallest_edge = -1;
			Weight smallest_weight = INF;
			Q.insert(c);
			
			//do a dfs to find all the edges leading out of this component.
			while (Q.size()) {
				int cur_component = Q.removeMin();
				
				int last_p = -1;
				int last_edge = -1;
				
				if (cur_component != c) {
					//connect these two components together
					int edgeid = edge_to_component[cur_component];
					assert(g.getWeight(edgeid) == component_edge_weight[cur_component]);
					int u = g.getEdge(edgeid).from;
					int v = g.getEdge(edgeid).to;
					assert(components[u] == c || components[v] == c);
					assert(components[u] == cur_component || components[v] == cur_component);
					
					last_edge = edgeid;
					if (components[u] == cur_component) {
						//attach these components together using this edge.
						//this will force us to re-root cur_component at u (this will happen during the dfs below)
						start_node = u;
						last_p = v;
					} else {
						start_node = v;
						last_p = u;
					}
					
					num_sets--;
					in_tree[edgeid] = true;
					min_weight += g.getWeight(edgeid);
					parents[start_node] = last_p;
					parent_edges[start_node] = edgeid;
					assert(components[start_node] == cur_component);
					components[start_node] = c;
					component_member[cur_component] = -1;
					empty_components.push_back(cur_component);
				}
				component_edge_weight[cur_component] = INF;
				q.clear();
				q.push_back(start_node);
				seen[start_node] = true;
				//do a bfs over the component, finding all edges that leave the component. fix the parent edges in the same pass, if needed.
				//could also do a dfs here - would it make a difference?
				for (int i = 0; i < q.size(); i++) {
					int n = q[i];
					assert(components[n] == c);
					for (int i = 0; i < g.nIncident(n, true); i++) {
						auto & edge = g.incident(n, i, true);
						if (edge_enabled[edge.id]) {
							int t = edge.node;
							if (!in_tree[edge.id]) {
								//assert(g.edgeEnabled(edge.id));
								int ncomponent = components[t];
								if (ncomponent != c && ncomponent != cur_component) {
									Weight w = g.getWeight(edge.id);
									if (w < component_edge_weight[ncomponent]) {
										
										edge_to_component[ncomponent] = edge.id;
										component_edge_weight[ncomponent] = w;
										
										Q.update(ncomponent);
									}
								}
							} else if (components[t] == cur_component) {
								if (components[t] != c) {
									components[t] = c;
									parents[t] = n;
									parent_edges[t] = edge.id;
								}
								if (!seen[t]) {				//can we avoid this seen check?
									seen[t] = true;
									q.push_back(t);
								}
							}
						}
					}
				}
				for (int s : q) {
					assert(seen[s]);
					seen[s] = false;
				}
				q.clear();
#ifndef NDEBUG
				for (bool b : q)
					assert(!b);
#endif
				
			}
		}
		components_to_visit.clear();
	}
	long num_updates = 0;
	int numUpdates() const {
		return num_updates;
	}
	void update() {
		static int iteration = 0;
		int local_it = ++iteration;

		if (g.outfile) {
			fprintf(g.outfile, "m\n");
			fflush(g.outfile);
		}

		if (last_modification > 0 && g.modifications == last_modification) {
			assert(min_weight == dbg.forestWeight());
			assert(num_sets == dbg.numComponents());
			return;
		}

		assert(components_to_visit.size() == 0);
		if (last_modification <= 0 || g.changed() || last_history_clear != g.historyclears) {
			INF = 1;				//g.nodes()+1;
			setNodes(g.nodes());

			for (auto & w : g.getWeights())
				INF += w;
			edge_enabled.clear();
			edge_enabled.resize(g.edges());
			for (int i = 0;i<g.edges();i++)
				if (g.hasEdge(i) && g.edgeEnabled(i)){
					edge_enabled[i]=true;//this can be improved (see Kohli-Torr!)
				}
			seen.clear();
			seen.resize(g.nodes());
			min_weight = 0;
			num_sets = g.nodes();
			empty_components.clear();
			component_needs_visit.clear();
			for (int i = 0; i < g.nodes(); i++)
				component_needs_visit.push_back(false);
			components.clear();
			for (int i = 0; i < g.nodes(); i++) {
				components.push_back(i);
				components_to_visit.push_back(i);
				component_needs_visit[i] = true;
			}
			component_edge_weight.clear();
			while (component_edge_weight.size() <= g.nodes()) {
				component_edge_weight.push_back(INF);
			}
			component_member.clear();
			for (int i = 0; i < g.nodes(); i++)
				component_member.push_back(i);
			
			mst.clear();
			parents.clear();
			parents.resize(g.nodes(), -1);
			parent_edges.clear();
			parent_edges.resize(g.edges(), -1);
			for (int i = 0; i < in_tree.size(); i++)
				in_tree[i] = false;
			last_history_clear = g.historyclears;
			history_qhead = g.historySize();//have to skip any additions or deletions that are left here, as otherwise the tree wont be an MST at the beginning of the addEdgeToMST method, which is an error.
					
		}

		//std::cout<<"Weight " << min_weight << " Components " << num_sets << " Dbg Weight: " << dbg.forestWeight() << " Components " << dbg.numComponents() <<"\n";
		for (int i = history_qhead; i < g.historySize(); i++) {
			
			int edgeid = g.getChange(i).id;
			if (g.getChange(i).addition && g.edgeEnabled(edgeid) && !edge_enabled[edgeid]) {
				prims();//to maintain correctness in spirapan, prims apparently must be called before addEdgeToMST.
				//however, the current implementation can likely be improved by only running prims on the components of the endpoints of edgeid...
				edge_enabled[edgeid]=true;
				addEdgeToMST(edgeid);
			} else if (!g.getChange(i).addition && !g.edgeEnabled(edgeid) && edge_enabled[edgeid]) {
				removeEdgeFromMST(edgeid);
				edge_enabled[edgeid]=false;
			}
			//std::cout<<"Weight " << min_weight << " Components " << num_sets << " Dbg Weight: " << dbg.forestWeight() << " Components " << dbg.numComponents() <<"\n";
		}
#ifndef NDEBUG
		for(int i = 0;i<g.edges();i++)
			assert(g.edgeEnabled(i)==edge_enabled[i]);
#endif
		prims();
		//std::cout<<"Weight " << min_weight << " Components " << num_sets << " Dbg Weight: " << dbg.forestWeight() << " Components " << dbg.numComponents() <<"\n";
		//g.drawFull(true);
		dbg_parents();
		//dbg_printSpanningTree();
#ifndef NDEBUG
		Weight expect = dbg.forestWeight();
		//dbg.dbg_printSpanningTree();
		assert(min_weight == dbg.forestWeight());
		assert(num_sets == dbg.numComponents());
#endif
		
		status.setMinimumSpanningTree(num_sets > 1 ? INF : min_weight, num_sets <= 1);
		
		//if(reportPolarity>-1){
		for (int i = 0; i < in_tree.size(); i++) {
			
			//Note: for the tree edge detector, polarity is effectively reversed.
			if (reportPolarity < 1 && (!g.edgeEnabled(i) || in_tree[i])) {
				status.inMinimumSpanningTree(i, true);
				assert(!g.edgeEnabled(i) || dbg.edgeInTree(i));
			} else if (reportPolarity > -1 && (g.edgeEnabled(i) && !in_tree[i])) {
				status.inMinimumSpanningTree(i, false);
				assert(!dbg.edgeInTree(i));
			}
		}
		assert(dbg_uptodate());
		num_updates++;
		last_modification = g.modifications;
		last_deletion = g.deletions;
		last_addition = g.additions;
		
		history_qhead = g.historySize();
		g.updateAlgorithmHistory(this,alg_id,history_qhead);
		last_history_clear = g.historyclears;
		
		;
	}

	void updateHistory(){
		update();
	}

	std::vector<int> & getSpanningTree() {
		update();
		return mst;
	}
	
	int getParent(int node) {
		update();
		
		return parents[node];
	}
	int getParentEdge(int node) {
		if (getParent(node) != -1)
			return parent_edges[node];
		else
			return -1;
	}
	bool edgeInTree(int edgeid) {
		update();
		return in_tree[edgeid];
	}
	bool dbg_mst() {
		
		return true;
	}
	
	Weight & weight() {
		update();
		assert(dbg_uptodate());
		if (num_sets == 1)
			return min_weight;
		else
			return INF;
	}
	
	Weight & forestWeight() {
		update();
		assert(dbg_uptodate());
		return min_weight;
	}
	int numComponents() {
		update();
		return num_sets;
	}
	int getComponent(int node) {
		update();
		return components[node];
	}
	int getRoot(int component = 0) {
		update();
		int c = 0;
		int nc = 0;
		while (component_member[c] == -1 || nc++ < component) {
			c++;
		}
		assert(component_member[c] > -1);
		return component_member[c];
		//return parents[component];
	}
	
	bool dbg_uptodate() {
#ifndef NDEBUG
		
		dbg_parents();
		//int n_components = 0;
		//check that each component has a unique root
		for (int c = 0; c < g.nodes(); c++) {
			if (component_member[c] > -1) {
				//n_components++;
				int root = component_member[c];
				while (parents[root] > -1) {
					root = parents[root];
				}
				
				for (int i = 0; i < g.nodes(); i++) {
					if (components[i] == c) {
						//check that the root of i is root
						int p = i;
						while (parents[p] > -1) {
							p = parents[p];
						}
						assert(p == root);
					}
				}
				
			}
		}
		
		assert(num_sets == g.nodes() - empty_components.size());
		Weight sumweight = 0;
		in_tree.resize(g.nEdgeIDs());
		for (int i = 0; i < g.nEdgeIDs(); i++) {
			if (in_tree[i]) {
				sumweight += g.getWeight(i);
			}
		}
		assert(sumweight == min_weight);
#endif
		return true;
	}
	;
	
	/*
	 #ifndef NDEBUG
	 int rootcount =0;
	 for(int i = 0;i<parents.size();i++){
	 if(parents[i]==-1)
	 rootcount++;
	 }
	 assert(rootcount==1);
	 #endif
	 */

};
}
;
#endif
