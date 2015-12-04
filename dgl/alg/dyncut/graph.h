/* graph.h */
/*
 This software library implements the maxflow algorithm
 described in

 "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision."
 Yuri Boykov and Vladimir Kolmogorov.
 In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 
 September 2004

 This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
 at Siemens Corporate Research. To make it available for public use,
 it was later reimplemented by Vladimir Kolmogorov based on open publications.

 If you use this software for research purposes, you should cite
 the aforementioned paper in any resulting publication.

 ----------------------------------------------------------------------

 REUSING TREES:

 Starting with version 3.0, there is a also an option of reusing search
 trees from one maxflow computation to the next, as described in

 "Efficiently Solving Dynamic Markov Random Fields Using Graph Cuts."
 Pushmeet Kohli and Philip H.S. Torr
 International Conference on Computer Vision (ICCV), 2005

 If you use this option, you should cite
 the aforementioned paper in any resulting publication.
 */

/*
 For description, license, example usage see README.TXT.
 */

/*
 This code is from Reduce, Reuse and Recycle (2008), Karteek Alahari, Pushmeet Kohli and Philip Torr
 (http://research.microsoft.com/en-us/um/people/pkohli/code/rrr.zip)

 Relevant README section:

 Dynamic Graph Cuts, version 2
 Copyright 2005 Pushmeet Kohli (pushmeet.kohli@brookes.ac.uk), Philip HS Torr
 (philiptorr@brookes.ac.uk).
 [For files block.h graph.h]

 This software library implements the dynamic maxflow algorithm described in:

 Efficiently Solving Dynamic Markov Random Fields using Graph Cuts
 Pushmeet Kohli and Philip H. S. Torr
 In the Tenth IEEE International Conference on Computer Vision (ICCV 2005).

 The algorithm uses the maxflow algorithm code described in:

 An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy
 Minimization in Vision, Yuri Boykov and Vladimir Kolmogorov.
 In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI),
 September 2004.

 The source code also comes under the following license:

 Copyright 2001 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk), Yuri Boykov
 (yuri@csd.uwo.ca).

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef __KOHLI_TORR_GRAPH_H__
#define __KOHLI_TORR_GRAPH_H__
#include <vector>
#include <string.h>
#include "block.h"
#include <cstdlib>
#include <assert.h>
#include <new>
#include <algorithm>
#include <stdexcept>
// NOTE: in UNIX you need to use -DNDEBUG preprocessor option to suppress assert's!!!

namespace kohli_torr {

/*
 special constants for node->parent
 */
#define TERMINAL ( (arc *) 1 )		/* to terminal */
#define ORPHAN   ( (arc *) 2 )		/* orphan */
#define INFINITE_D ((int)(((unsigned)-1)/2))		/* infinite distance to the terminal */

// captype: type of edge capacities (excluding t-links)
// tcaptype: type of t-links (edges between nodes and terminals)
// flowtype: type of total flow
//
template<typename captype, typename tcaptype, typename flowtype> class Graph {
	
	struct node;
	struct arc;

public:
	typedef enum {
		SOURCE = 0, SINK = 1
	} termtype; // terminals 
	typedef int node_id;
	typedef arc* arc_id;
	bool preserve_backward_order = false;
	std::vector<node_id> t_edge_nodes;//added by sam; used to find the flows on each edge
	std::vector<node_id> s_edge_nodes;//added by sam; used to find the flows on each edge
	std::vector<int> changed_edges;//added by sam; keeps track of edges whose flows might have changed during the current call
	std::vector<node_id> nodes_partition_changed; //added by sam: keeps track of nodes whose partition might have changed.
											  //See also changed_list, from the original implementation, which is similar but has slightly different behaviour.

	/////////////////////////////////////////////////////////////////////////
	//                     BASIC INTERFACE FUNCTIONS                       //
	//              (should be enough for most applications)               //
	/////////////////////////////////////////////////////////////////////////
	
	// Constructor. 
	// The first argument gives an estimate of the maximum number of nodes that can be added
	// to the graph, and the second argument is an estimate of the maximum number of edges.
	// The last (optional) argument is the pointer to the function which will be called 
	// if an error occurs; an error message is passed to this function. 
	// If this argument is omitted, an exception will be thrown.
	//
	// IMPORTANT: It is possible to add more nodes to the graph than node_num_max 
	// (and node_num_max can be zero). However, if the count is exceeded, then 
	// the internal memory is reallocated (increased by 50%) which is expensive. 
	// Also, temporarily the amount of allocated memory would be more than twice than needed.
	// Similarly for edges.
	// If you wish to avoid this overhead, you can download version 2.2, where nodes and edges are stored in blocks.
	Graph(int node_num_max, int edge_num_max, void (*err_function)(const char *) = NULL);

	// Destructor
	~Graph();

	// Adds node(s) to the graph. By default, one node is added (num=1); then first call returns 0, second call returns 1, and so on. 
	// If num>1, then several nodes are added, and node_id of the first one is returned.
	// IMPORTANT: see note about the constructor 
	node_id add_node(int num = 1);

	// Adds a bidirectional edge between 'i' and 'j' with the weights 'cap' and 'rev_cap'.
	// IMPORTANT: see note about the constructor 
	int add_edge(node_id i, node_id j, captype cap, captype rev_cap);

	// Adds new edges 'SOURCE->i' and 'i->SINK' with corresponding weights.
	// Can be called multiple times for each node.
	// Weights can be negative.
	// NOTE: the number of such edges is not counted in edge_num_max.
	//       No internal memory is allocated by this call.
	void add_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink);

	// Computes the maxflow. Can be called several times.
	// FOR DESCRIPTION OF reuse_trees, SEE mark_node().
	// FOR DESCRIPTION OF changed_list, SEE remove_from_changed_list().
	flowtype maxflow(bool reuse_trees = false, Block<node_id>** changed_list = NULL);

	// After the maxflow is computed, this function returns to which
	// segment the node 'i' belongs (Graph<captype,tcaptype,flowtype>::SOURCE or Graph<captype,tcaptype,flowtype>::SINK).
	//
	// Occasionally there may be several minimum cuts. If a node can be assigned
	// to both the source and the sink, then default_segm is returned.
	termtype what_segment(node_id i, termtype default_segm = SOURCE);

	//////////////////////////////////////////////
	//       ADVANCED INTERFACE FUNCTIONS       //
	//      (provide access to the graph)       //
	//////////////////////////////////////////////
	
	////////////////////////////
	// 1. Reallocating graph. //
	////////////////////////////
	
	// Removes all nodes and edges. 
	// After that functions add_node() and add_edge() must be called again. 
	//
	// Advantage compared to deleting Graph and allocating it again:
	// no calls to delete/new (which could be quite slow).
	//
	// If the graph structure stays the same, then an alternative
	// is to go through all nodes/edges and set new residual capacities
	// (see functions below).
	void reset();

	////////////////////////////////////////////////////////////////////////////////
	// 2. Functions for getting pointers to arcs and for reading graph structure. //
	//    NOTE: adding new arcs may invalidate these pointers (if reallocation    //
	//    happens). So it's best not to add arcs while reading graph structure.   //
	////////////////////////////////////////////////////////////////////////////////
	
	// The following two functions return arcs in the same order that they
	// were added to the graph. NOTE: for each call add_edge(i,j,cap,cap_rev)
	// the first arc returned will be i->j, and the second j->i.
	// If there are no more arcs, then the function can still be called, but
	// the returned arc_id is undetermined.
	
	arc_id get_first_arc();
	arc_id get_arc(int index) {
		return &arcs[index];
	}
	
	arc* get_reverse(arc* from) {
		return from->sister;
	}
	
	arc_id get_next_arc(arc_id a);

	// other functions for reading graph structure
	int get_node_num() {
		return node_num;
	}
	//returns number of arcs
	int get_arc_num() {
		return (int) (arc_last - arcs);
	}
	void get_arc_ends(arc_id a, node_id& i, node_id& j); // returns i,j to that a = i->j
			
	///////////////////////////////////////////////////
	// 3. Functions for reading residual capacities. //
	///////////////////////////////////////////////////
	
	// returns residual capacity of SOURCE->i minus residual capacity of i->SINK
	tcaptype get_trcap(node_id i);
	// returns residual capacity of arc a
	captype get_rcap(arc* a);
	captype get_ecap(arc* a);
	captype get_flow(arc* a);	//added by Sam
			
	/////////////////////////////////////////////////////////////////
	// 4. Functions for setting residual capacities.               //
	//    NOTE: If these functions are used, the value of the flow //
	//    returned by maxflow() will not be valid!                 //
	/////////////////////////////////////////////////////////////////
	
	void set_trcap(node_id i, tcaptype trcap);
	void set_rcap(arc* a, captype rcap);

	tcaptype getTweight(node_id i) {
		return nodes[i].t_cap;
	}
	
	// Edit capacity of t-edge when "using" tree-recycling 
	void edit_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink);

	// Edit capacity of t-edge when "not using" tree-recycling :		
	// If yoy are editing capacities using this function, "maxflow(false)" needs to be called
	void edit_tweights_wt(node_id i, tcaptype cap_source, tcaptype cap_sink);

	//(Added by Sam)
	//Don't use - linear search!
	bool has_edge(node_id from, node_id to) {

		arc *a, *a_rev;
		a = nodes[from].first;
		
		while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
			a = a->next;
		
		if (!a || a->head != &nodes[to]) {
			return false;
		} else {
			return true;
		}
	}
	
	//(Added by Sam)
	//Don't use - linear search!
	captype get_edge_capacity(node_id from, node_id to,arc*ac=nullptr) {
		arc *a, *a_rev;
		if(!ac){
			assert(false);
			a = nodes[from].first;

			while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
				a = a->next;
		}else{
#ifndef NDEBUG
			a = nodes[from].first;
			while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
				a = a->next;
			assert(a==ac);
#endif
			a = ac;
			assert(a->head==&nodes[to]);

		}
		
		if (a->head != &nodes[to]) {
			throw std::invalid_argument("Specified edge doesn't exist");
		} else {
			return a->e_cap;
		}
	}

	// Edit capacity of n-edge when "using" tree-recycling, by adding (or subtracting) from the existing capacity
	//(Added by Sam)
	void edit_edge_inc(node_id from, node_id to, captype added_cap, captype added_rev_cap,arc*a=nullptr);

	// Edit capacity of n-edge when "using" tree-recycling 
	void edit_edge(node_id from, node_id to, captype cap, captype rev_cap,arc*a=nullptr);

	// Edit capacity of n-edge when "not using" tree-recycling :		
	// If you are editing capacities using this function, "maxflow(false)" needs to be called
	void edit_edge_wt(node_id from, node_id to, captype cap, captype rev_cap,arc*a=nullptr);

	tcaptype MIN(tcaptype a, tcaptype b);
	tcaptype MAX(tcaptype a, tcaptype b);

	////////////////////////////////////////////////////////////////////
	// 5. Functions related to reusing trees & list of changed nodes. //
	////////////////////////////////////////////////////////////////////
	
	// If flag reuse_trees is true while calling maxflow(), then search trees
	// are reused from previous maxflow computation (unless it's the first call to maxflow()).
	// In this case BEFORE calling maxflow() the user must
	// specify which parts of the graph have changed by calling mark_node():
	//   add_tweights(i),set_trcap(i)    => call mark_node(i)
	//   add_edge(i,j),set_rcap(a)       => call mark_node(i); mark_node(j)
	//
	// This option makes sense only if a small part of the graph is changed.
	// The initialization procedure goes only through marked nodes then.
	// 
	// mark_node(i) can either be called before or after graph modification.
	// Can be called more than once per node, but calls after the first one
	// do not have any effect.
	// 
	// NOTE: 
	//   1. It is not necessary to call mark_node() if the change is ``not essential'',
	//      i.e. sign(trcap) is preserved for a node and zero/nonzero status is preserved for an arc.
	//   2. To check that you marked all necessary nodes, you can call maxflow(true) after calling maxflow(false).
	//      If everything is correct, the two calls must return the same value of flow. (Useful for debugging).
	void mark_node(node_id i);

	// If changed_list is not NULL while calling maxflow(), then the algorithm
	// keeps a list of nodes which could potentially have changed their segmentation label
	// (unless it's the first call to maxflow).
	// In this case AFTER calling maxflow() the user must call remove_from_changed_list()
	// for every node in the list. (Exception: this is necessary only if the next 
	// maxflow computation uses option reuse_trees).
	//
	// Nodes which are not in the list are guaranteed to keep their old segmentation label (SOURCE or SINK).
	//
	// Pointer to the list is returned in changed_list. (See block.h on how to read it).
	//
	// NOTE: The user should not deallocate the returned list. (This will be done by Graph destructor).
	void remove_from_changed_list(node_id i) {
		assert(i >= 0 && i < node_num && nodes[i].is_in_changed_list);
		nodes[i].is_in_changed_list = 0;
	}
	
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	
private:
	// internal variables and functions
	
	struct node {
		arc *first = nullptr;		// first outcoming arc
		
		arc *parent = nullptr;	// node's parent
		node *next = nullptr;		// pointer to the next active node
		//   (or to itself if it is the last node in the list)
		int TS = 0;		// timestamp showing when DIST was computed
		int DIST = 0;		// distance to the terminal
		int is_sink :1;	// flag showing whether the node is in the source or in the sink tree (if parent!=NULL)
		int is_marked :1;	// set by mark_node()
		int is_in_changed_list :1; // set by maxflow if node has changed partitions since last assignment
		int is_in_partition_list :1; //added by sam: set if a node has changed partitions, ever (including first assignment)
		int in_t_edges_set :1;//added by sam
		int in_s_edges_set :1;//added by sam
		tcaptype tr_cap;		// if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
		// otherwise         -tr_cap is residual capacity of the arc node->SINK 
		
		tcaptype t_cap;
		tcaptype con_flow;
	};

	struct arc {
		node *head = nullptr;		// node the arc points to
		arc *next = nullptr;		// next arc with the same originating node
		arc *sister = nullptr;	// reverse arc
		bool might_have_flow = false; //added by sam
		captype r_cap;		// residual capacity
		captype e_cap;		// original capacity
	};

	struct nodeptr {
		node *ptr = nullptr;
		nodeptr *next = nullptr;
	};
	static const int NODEPTR_BLOCK_SIZE = 128;

	node *nodes, *node_last, *node_max; // node_last = nodes+node_num, node_max = nodes+node_num_max;
	arc *arcs, *arc_last, *arc_max; // arc_last = arcs+2*edge_num, arc_max = arcs+2*edge_num_max;
			
	int node_num;

	DBlock<nodeptr> *nodeptr_block;

	flowtype flow;			// total flow
	
	// reusing trees & list of changed pixels
	int maxflow_iteration=0; // counter
	bool keep_changed_list=false;
	Block<node_id> *changed_list;

	void (*error_function)(const char *);	// this function is called if a error occurs,
	// with a corresponding error message
	// An exception will be thrown if this is NULL.
	/////////////////////////////////////////////////////////////////////////
	
	node *queue_first[2], *queue_last[2];	// list of active nodes
	nodeptr *orphan_first=nullptr, *orphan_last=nullptr;		// list of pointers to orphans
	int TIME=0;					// monotonically increasing global counter
	
	/////////////////////////////////////////////////////////////////////////
	
	void reallocate_nodes(int num); // num is the number of new nodes
	void reallocate_arcs();

	// functions for processing active list
	void set_active(node *i);
	node *next_active();

	// functions for processing orphans list
	void set_orphan_front(node* i); // add to the beginning of the list
	void set_orphan_rear(node* i);  // add to the end of the list
			
	void add_to_changed_list(node* i);

	void maxflow_init();             // called if reuse_trees == false
	void maxflow_reuse_trees_init(); // called if reuse_trees == true
	void augment(arc *middle_arc);
	void process_source_orphan(node *i);
	void process_sink_orphan(node *i);

	std::vector<int> Q;
	std::vector<arc*> prev;
	std::vector<flowtype> M;
	arc fake_arc;
	arc fake_prev;
	int edmonds_karp_bfs(flowtype & store_flow, int source_node, int sink_node, captype & bridge_capacity,
			bool backward) {
		prev.resize(this->get_node_num(), nullptr);
		M.resize(this->get_node_num(), 0);
		for (int i = 0; i < Q.size(); i++) {
			int u = Q[i];
			
			prev[u] = nullptr;
			
		}
#ifndef NDEBUG
		for (int i = 0; i < this->get_node_num(); i++)
			assert(prev[i] == nullptr);
#endif
		Q.clear();
		if (!backward) {
			for (int i = 0; i < t_edge_nodes.size(); i++) {
				int u = t_edge_nodes[i];
				if (source_node == sink_node && u == source_node)
					continue;
				assert(nodes[u].in_t_edges_set);
				tcaptype & t_edge_flow = nodes[u].tr_cap;
				if (t_edge_flow < 0) {
					//there is flow from the source node to this node.
					//which we want to eliminate.
					//so treat this as capacity on this edge
					prev[u] = &fake_prev;
					M[u] = -t_edge_flow;
					Q.push_back(u);
				}
			}
		} else {
			for (int i = 0; i < s_edge_nodes.size(); i++) {
				int u = s_edge_nodes[i];
				if (source_node == sink_node && u == source_node)
					continue;
				assert(nodes[u].in_s_edges_set);
				tcaptype & s_edge_flow = nodes[u].tr_cap;
				if (s_edge_flow > 0) {
					//there is flow from the source node to this node.
					//which we want to eliminate.
					//so treat this as capacity on this edge
					prev[u] = &fake_prev;
					M[u] = s_edge_flow;
					Q.push_back(u);
				}
			}
		}
		
		for (int j = 0; j < Q.size(); j++) {
			int u = Q[j];
			
			//first check if this node has capacity to the sink.
			tcaptype & t_edge_flow = nodes[u].tr_cap;
			if (!backward && ((source_node != sink_node && t_edge_flow > 0))) {
				store_flow = std::min(M[u], t_edge_flow);
				return u;
			} else if (!backward && source_node == sink_node && u == sink_node) {
				// assert(t_edge_flow<0);
				// tcaptype s_edge_flow = -t_edge_flow;
				store_flow = M[u]; //std::min(M[u],s_edge_flow);
				return u;
			} else if (backward && source_node == sink_node && u == sink_node) {
				store_flow = M[u]; // std::min(M[u],t_edge_flow);
				return u;
			}
			
			if (u == sink_node) {
				//connect the sink node to the source node with infinite capacity
				int to = source_node;
				if (prev[to] == nullptr && bridge_capacity > 0) {
					prev[to] = &fake_arc;
					M[to] = std::min(M[u], bridge_capacity);
					Q.push_back(to);
				}
			}
			
			arc * a = nodes[u].first;
			int cpos = Q.size();
			while (a) {
				int to = a->head - nodes;
				if (prev[to] == nullptr) {
					
					flowtype f = a->e_cap - a->r_cap;
					if (!backward && f > 0) {
						prev[to] = a;
						M[to] = std::min(M[u], f);
						Q.push_back(to);
					} else if (backward && f < 0) {
						tcaptype fb = -f;
						prev[to] = a;
						M[to] = std::min(M[u], fb);
						Q.push_back(to);
					}
				}
				a = a->next;
			}
			if (preserve_backward_order && backward && Q.size() - cpos > 1) {
				//reverse the order of the edges added to the queue if walking backward through the graph, so that they are visited in the same order as if we had walked forwards through the graph
				std::reverse(Q.begin() + cpos, Q.end());
			}
		}
		store_flow = 0;
		return -2;
	}
	void unmarkPartitionAssignment(node * n){
		n->is_in_partition_list=false;
	}
public:
	inline void partitionAssignment(node * n){
		if(!n->is_in_partition_list){
			n->is_in_partition_list=true;
			nodes_partition_changed.push_back( (node_id) (n - nodes));
		}
	}

	void unmarkPartitionAssignment(int node_id){
		unmarkPartitionAssignment(nodes+node_id);
	}
	void markFlowEdge(arc * edge) {
		if (!edge->might_have_flow) {
			assert(!edge->sister->might_have_flow);
			edge->might_have_flow = true;
			edge->sister->might_have_flow = true;
			changed_edges.push_back(edge - get_first_arc());
			changed_edges.push_back(edge->sister - get_first_arc());
		}
		assert(edge->might_have_flow);
		assert(edge->sister->might_have_flow);
	}
	
	void unmarkFlowEdge(arc * edge) {
		edge->might_have_flow = false;
		edge->sister->might_have_flow = false;
	}

	
	//Run edmonds-karp to remove any excess flow on t-edges
	void clear_t_edges(int source_node, int sink_node) {
		static int iter = 0;
		if (++iter == 96) {
			int a = 1;
		}
		flowtype total_flow = maxflow(true, nullptr);
		flowtype f = 0;
#ifndef NDEBUG
		for (int i = 0; i < get_node_num(); i++) {
			tcaptype & t_edge_flow = nodes[i].tr_cap;
			if (t_edge_flow < 0) {
				//if(i!=sink_node){//the sink node will be skipped anyhow, so it doesnt need to be in the t_edges_set.
				assert(nodes[i].in_t_edges_set);
				assert(std::find(t_edge_nodes.begin(), t_edge_nodes.end(), i) != t_edge_nodes.end());
				//}
			}
			if (t_edge_flow > 0) {
				//if(i!=sink_node){//the sink node will be skipped anyhow, so it doesnt need to be in the t_edges_set.
				assert(nodes[i].in_s_edges_set);
				assert(std::find(s_edge_nodes.begin(), s_edge_nodes.end(), i) != s_edge_nodes.end());
				//}
			}
			if (!nodes[i].in_t_edges_set) {
				assert(std::find(t_edge_nodes.begin(), t_edge_nodes.end(), i) == t_edge_nodes.end());
			}
			if (!nodes[i].in_s_edges_set) {
				assert(std::find(s_edge_nodes.begin(), s_edge_nodes.end(), i) == s_edge_nodes.end());
			}
		}
		
#endif
		captype bridge_capacity = -(nodes[sink_node].t_cap - nodes[sink_node].tr_cap + total_flow);
		nodes[sink_node].tr_cap = 0; //( nodes[sink_node].t_cap-nodes[sink_node].tr_cap+ total_flow);
		nodes[source_node].tr_cap = 0; // ( nodes[source_node].t_cap-nodes[source_node].tr_cap- total_flow);
		while (true) {
			int node = edmonds_karp_bfs(f, source_node, sink_node, bridge_capacity, false);
			
			if (node < 0)
				break;
			assert(f > 0);
			assert(nodes[node].tr_cap >= f);
			nodes[node].tr_cap -= f;
			int v = node;
			while (prev[v] != &fake_prev) {
				arc* edge = prev[v];
				if (edge == &fake_arc) {
					//then this is the inserted, infinite capacity arc between the sink and source nodes.
					assert(v == source_node);
					assert(bridge_capacity >= f);
					bridge_capacity -= f;
					v = sink_node;
					continue;
				}
				assert(edge->e_cap - edge->r_cap >= f);
				edge->r_cap += f; //remove this flow from this edge by adding it to its remaining capacity;
						
				markFlowEdge(edge->sister);
				assert(edge->sister->r_cap >= f);
				edge->sister->r_cap -= f;
				int u = edge->sister->head - nodes;
				v = u;
			}
			assert(nodes[v].tr_cap <= -f);
			nodes[v].tr_cap += f;
		}
		
		while (true) {
			int node = edmonds_karp_bfs(f, source_node, source_node, bridge_capacity, true);
			
			if (node < 0)
				break;
			assert(f > 0);
			/*        assert( nodes[node].tr_cap>=f);
			 nodes[node].tr_cap-=f;*/
			int v = node;
			while (prev[v] != &fake_prev) {
				arc* edge = prev[v];
				if (edge == &fake_arc) {
					//then this is the inserted, infinite capacity arc between the sink and source nodes.
					assert(v == source_node);
					assert(bridge_capacity >= f);
					bridge_capacity -= f;
					v = sink_node;
					continue;
				}
				assert(edge->r_cap >= f);
				edge->r_cap -= f; //remove this flow from this edge by adding it to its remaining capacity;
				markFlowEdge(edge);
				edge->sister->r_cap += f;
				int u = edge->sister->head - nodes;
				v = u;
			}
			assert(nodes[v].tr_cap >= f);
			nodes[v].tr_cap -= f;
		}
		while (true) {
			int node = edmonds_karp_bfs(f, sink_node, sink_node, bridge_capacity, false);
			
			if (node < 0)
				break;
			assert(f > 0);
			/*	   assert( nodes[node].tr_cap>=f);
			 nodes[node].tr_cap-=f;*/
			int v = node;
			while (prev[v] != &fake_prev) {
				arc* edge = prev[v];
				if (edge == &fake_arc) {
					//then this is the inserted, infinite capacity arc between the sink and source nodes.
					assert(v == source_node);
					assert(bridge_capacity >= f);
					bridge_capacity -= f;
					v = sink_node;
					continue;
				}
				assert(edge->e_cap - edge->r_cap >= f);
				edge->r_cap += f; //remove this flow from this edge by adding it to its remaining capacity;
				markFlowEdge(edge->sister);
				assert(edge->sister->r_cap >= f);
				edge->sister->r_cap -= f;
				int u = edge->sister->head - nodes;
				v = u;
			}
			assert(nodes[v].tr_cap <= -f);
			nodes[v].tr_cap += f;
		}
		nodes[sink_node].tr_cap = (nodes[sink_node].t_cap + total_flow);
		nodes[source_node].tr_cap = (nodes[source_node].t_cap - total_flow);
		
		int i, j = 0;
		for (i = 0; i < t_edge_nodes.size(); i++) {
			int u = t_edge_nodes[i];
			assert(nodes[u].in_t_edges_set);
			if (nodes[u].tr_cap >= 0) {
				nodes[u].in_t_edges_set = 0;
			} else {
				t_edge_nodes[j++] = u;
			}
		}
		t_edge_nodes.resize(j);
		
		j = 0;
		for (i = 0; i < s_edge_nodes.size(); i++) {
			int u = s_edge_nodes[i];
			assert(nodes[u].in_s_edges_set);
			if (nodes[u].tr_cap <= 0) {
				nodes[u].in_s_edges_set = 0;
			} else {
				s_edge_nodes[j++] = u;
			}
		}
		s_edge_nodes.resize(j);
		
	}
	
	void test_consistency(node* current_node = NULL); // debug function
};

///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////

template<typename captype, typename tcaptype, typename flowtype>
inline typename Graph<captype, tcaptype, flowtype>::node_id Graph<captype, tcaptype, flowtype>::add_node(int num) {
	assert(num > 0);
	
	if (node_last + num > node_max)
		reallocate_nodes(num);
	
	if (num == 1) {
		node_last->first = NULL;
		node_last->tr_cap = 0;
		node_last->t_cap = 0;
		node_last->con_flow = 0;
		node_last->is_marked = 0;
		node_last->is_in_changed_list = 0;
		node_last->is_in_partition_list = 0;
		node_last->in_t_edges_set = 0;
		node_last->in_s_edges_set = 0;
		node_last++;
		return node_num++;
	} else {
		memset(node_last, 0, num * sizeof(node));
		
		node_id i = node_num;
		node_num += num;
		node_last += num;
		return i;
	}
}

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::add_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink) {
	assert(i >= 0 && i < node_num);
	
	nodes[i].tr_cap = cap_source - cap_sink;
	nodes[i].t_cap = cap_source - cap_sink;
	nodes[i].con_flow = MIN(cap_source, cap_sink);
	flow += MIN(cap_source, cap_sink);
}

template<typename captype, typename tcaptype, typename flowtype>
inline int Graph<captype, tcaptype, flowtype>::add_edge(node_id _i, node_id _j, captype cap, captype rev_cap) {
	assert(_i >= 0 && _i < node_num);
	assert(_j >= 0 && _j < node_num);
	assert(_i != _j);
	assert(cap >= 0);
	assert(rev_cap >= 0);
	
	if (arc_last == arc_max)
		reallocate_arcs();
	
	arc *a = arc_last++;
	arc *a_rev = arc_last++;
	
	node* i = nodes + _i;
	node* j = nodes + _j;
	
	a->sister = a_rev;
	a_rev->sister = a;
	a->next = i->first;
	i->first = a;
	a_rev->next = j->first;
	j->first = a_rev;
	a->head = j;
	a_rev->head = i;
	a->r_cap = cap;
	a_rev->r_cap = rev_cap;
	a->e_cap = cap;
	a_rev->e_cap = rev_cap;
	return a - arcs;
}

template<typename captype, typename tcaptype, typename flowtype>
inline typename Graph<captype, tcaptype, flowtype>::arc* Graph<captype, tcaptype, flowtype>::get_first_arc() {
	return arcs;
}

template<typename captype, typename tcaptype, typename flowtype>
inline typename Graph<captype, tcaptype, flowtype>::arc* Graph<captype, tcaptype, flowtype>::get_next_arc(arc* a) {
	return a + 1;
}

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::get_arc_ends(arc* a, node_id& i, node_id& j) {
	assert(a >= arcs && a < arc_last);
	i = (node_id) (a->sister->head - nodes);
	j = (node_id) (a->head - nodes);
}

template<typename captype, typename tcaptype, typename flowtype>
inline tcaptype Graph<captype, tcaptype, flowtype>::get_trcap(node_id i) {
	assert(i >= 0 && i < node_num);
	return nodes[i].tr_cap;
}

template<typename captype, typename tcaptype, typename flowtype>
inline captype Graph<captype, tcaptype, flowtype>::get_rcap(arc* a) {
	assert(a >= arcs && a < arc_last);
	return a->r_cap;
}
template<typename captype, typename tcaptype, typename flowtype>
inline captype Graph<captype, tcaptype, flowtype>::get_ecap(arc* a) {
	assert(a >= arcs && a < arc_last);
	return a->e_cap;
}
template<typename captype, typename tcaptype, typename flowtype>
inline captype Graph<captype, tcaptype, flowtype>::get_flow(arc* a) {
	assert(a >= arcs && a < arc_last);
	return a->e_cap - a->r_cap;
}

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::set_trcap(node_id i, tcaptype trcap) {
	assert(i >= 0 && i < node_num);
	nodes[i].tr_cap = trcap;
}

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::set_rcap(arc* a, captype rcap) {
	assert(a >= arcs && a < arc_last);
	a->r_cap = rcap;
}

template<typename captype, typename tcaptype, typename flowtype>
inline typename Graph<captype, tcaptype, flowtype>::termtype Graph<captype, tcaptype, flowtype>::what_segment(node_id i,
		termtype default_segm) {
	if (nodes[i].parent) {
		return (nodes[i].is_sink) ? SINK : SOURCE;
	} else {
		return default_segm;
	}
}

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::mark_node(node_id _i) {
	node* i = nodes + _i;
	if (!i->next) {
		/* it's not in the list yet */
		if (queue_last[1])
			queue_last[1]->next = i;
		else
			queue_first[1] = i;
		queue_last[1] = i;
		i->next = i;
	}
	i->is_marked = 1;
}

template<typename captype, typename tcaptype, typename flowtype>
inline tcaptype Graph<captype, tcaptype, flowtype>::MIN(tcaptype a, tcaptype b) {
	if (a < b)
		return a;
	return b;
}

template<typename captype, typename tcaptype, typename flowtype>
inline tcaptype Graph<captype, tcaptype, flowtype>::MAX(tcaptype a, tcaptype b) {
	if (a > b)
		return a;
	return b;
}

template<typename captype, typename tcaptype, typename flowtype>
Graph<captype, tcaptype, flowtype>::Graph(int node_num_max, int edge_num_max, void (*err_function)(const char *)) :
		node_num(0), nodeptr_block(NULL), changed_list(NULL), error_function(err_function) {
	if (node_num_max < 16)
		node_num_max = 16;
	if (edge_num_max < 16)
		edge_num_max = 16;
	
	nodes = (node*) malloc(node_num_max * sizeof(node));
	arcs = (arc*) malloc(2 * edge_num_max * sizeof(arc));
	if (!nodes || !arcs) {
		if (error_function)
			(*error_function)("Not enough memory!");
		throw  std::bad_alloc();
	}
	
	for (int i = 0; i < node_num_max; i++) {
		new (&nodes[i]) node();
	}
	for (int i = 0; i < 2 * edge_num_max; i++) {
		new (&arcs[i]) arc();
	}
	
	node_last = nodes;
	node_max = nodes + node_num_max;
	arc_last = arcs;
	arc_max = arcs + 2 * edge_num_max;
	
	maxflow_iteration = 0;
	flow = 0;
}

template<typename captype, typename tcaptype, typename flowtype>
Graph<captype, tcaptype, flowtype>::~Graph() {
	if (nodeptr_block) {
		delete nodeptr_block;
		nodeptr_block = NULL;
	}
	if (changed_list) {
		delete changed_list;
		changed_list = NULL;
	}
	free(nodes);
	free(arcs);
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::reset() {
	node_last = nodes;
	arc_last = arcs;
	node_num = 0;
	
	if (nodeptr_block) {
		delete nodeptr_block;
		nodeptr_block = NULL;
	}
	
	maxflow_iteration = 0;
	flow = 0;
	
	t_edge_nodes.clear();
	s_edge_nodes.clear();
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::reallocate_nodes(int num) {
	int node_num_max = (int) (node_max - nodes);
	int num_nodes = (int) (node_max - nodes);
	node* nodes_old = nodes;
	
	node_num_max += node_num_max / 2;
	if (node_num_max < node_num + num)
		node_num_max = node_num + num;
	nodes = (node*) realloc(nodes_old, node_num_max * sizeof(node));
	if (!nodes) {
		if (error_function)
			(*error_function)("Not enough memory!");
		throw  std::bad_alloc();
	}
	
	for (int i = num_nodes; i < node_num_max; i++) {
		new (&nodes[i]) node();
	}
	
	node_last = nodes + node_num;
	node_max = nodes + node_num_max;
	
	if (nodes != nodes_old) {
		arc* a;
		for (a = arcs; a < arc_last; a++) {
			a->head = (node*) ((char*) a->head + (((char*) nodes) - ((char*) nodes_old)));
		}
	}
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::reallocate_arcs() {
	int arc_num_max = (int) (arc_max - arcs);
	int arc_num = (int) (arc_last - arcs);
	arc* arcs_old = arcs;
	
	arc_num_max += arc_num_max / 2;
	if (arc_num_max & 1)
		arc_num_max++;
	arcs = (arc*) realloc(arcs_old, arc_num_max * sizeof(arc));
	if (!arcs) {
		if (error_function)
			(*error_function)("Not enough memory!");
		throw  std::bad_alloc();
	}
	
	for (int i = arc_num; i < arc_num_max; i++) {
		new (&arcs[i]) arc();
	}
	
	arc_last = arcs + arc_num;
	arc_max = arcs + arc_num_max;
	
	if (arcs != arcs_old) {
		node* i;
		arc* a;
		for (i = nodes; i < node_last; i++) {
			if (i->first)
				i->first = (arc*) ((char*) i->first + (((char*) arcs) - ((char*) arcs_old)));
		}
		for (a = arcs; a < arc_last; a++) {
			if (a->next)
				a->next = (arc*) ((char*) a->next + (((char*) arcs) - ((char*) arcs_old)));
			a->sister = (arc*) ((char*) a->sister + (((char*) arcs) - ((char*) arcs_old)));
		}
	}
}

/*********************************************************************************************/
/***** Dynamic MAXFLOW CODE *******/
/*********************************************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::edit_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink) {
	tcaptype oldRes = nodes[i].t_cap;
	assert(cap_source>=0);assert(cap_sink>=0);
	if (nodes[i].t_cap != cap_source - cap_sink) {
		if (nodes[i].t_cap > 0)
			flow -= MIN(nodes[i].t_cap - nodes[i].tr_cap, nodes[i].t_cap);
		else
			flow += MAX(0, nodes[i].tr_cap);
		
		nodes[i].tr_cap = (cap_source - cap_sink) - (nodes[i].t_cap - nodes[i].tr_cap);
		nodes[i].t_cap = cap_source - cap_sink;
		
		if (nodes[i].t_cap > 0)
			flow += MIN(nodes[i].t_cap - nodes[i].tr_cap, nodes[i].t_cap);
		else
			flow -= MAX(0, nodes[i].tr_cap);
		
		if (!(((oldRes > 0) && (nodes[i].tr_cap > 0)) || ((oldRes < 0) && (nodes[i].tr_cap < 0))))
			mark_node(i);
	}
	flow -= nodes[i].con_flow;
	nodes[i].con_flow = MIN(cap_source, cap_sink);
	flow += MIN(cap_source, cap_sink);
	if (!nodes[i].in_t_edges_set && nodes[i].tr_cap < 0) {
		nodes[i].in_t_edges_set = 1;
		t_edge_nodes.push_back(i);
	} else if (!nodes[i].in_s_edges_set && nodes[i].tr_cap > 0) {
		nodes[i].in_s_edges_set = 1;
		s_edge_nodes.push_back(i);
	}
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::edit_tweights_wt(node_id i, tcaptype cap_source, tcaptype cap_sink) {
	tcaptype oldRes = nodes[i].t_cap;
	
	if (nodes[i].t_cap != cap_source - cap_sink) {
		if (nodes[i].t_cap > 0)
			flow -= MIN(nodes[i].t_cap - nodes[i].tr_cap, nodes[i].t_cap);
		else
			flow += MAX(0, nodes[i].tr_cap);
		
		nodes[i].tr_cap = (cap_source - cap_sink) - (nodes[i].t_cap - nodes[i].tr_cap);
		nodes[i].t_cap = cap_source - cap_sink;
		
		if (nodes[i].t_cap > 0)
			flow += MIN(nodes[i].t_cap - nodes[i].tr_cap, nodes[i].t_cap);
		else
			flow -= MAX(0, nodes[i].tr_cap);
	}
	flow -= nodes[i].con_flow;
	nodes[i].con_flow = MIN(cap_source, cap_sink);
	flow += MIN(cap_source, cap_sink);
	if (!nodes[i].in_t_edges_set && nodes[i].tr_cap < 0) {
		nodes[i].in_t_edges_set = 1;
		t_edge_nodes.push_back(i);
	} else if (!nodes[i].in_s_edges_set && nodes[i].tr_cap > 0) {
		nodes[i].in_s_edges_set = 1;
		s_edge_nodes.push_back(i);
	}
}

/***********************************************************************************************/
/***********************************************************************************************/

// Edit capacity of n-edge when "using" tree-recycling, by adding (or subtracting) from the existing capacity
template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::edit_edge_inc(node_id from, node_id to, captype added_cap,
		captype added_rev_cap,arc*ac) {
	arc *a, *a_rev;
	
	if(!ac){
		assert(false);
		a = nodes[from].first;
		while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
			a = a->next;
	}else{
#ifndef NDEBUG
		a = nodes[from].first;
		while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
			a = a->next;
		assert(a==ac);
#endif
		a = ac;
		assert(a->head==&nodes[to]);

	}

	
	if (!a || a->head != &nodes[to]) {
		throw std::invalid_argument("Specified edge doesn't exist");
	} else {
		// Modifying flow value
		
		if (nodes[from].t_cap > 0)
			flow -= MIN(nodes[from].t_cap - nodes[from].tr_cap, nodes[from].t_cap);
		else
			flow += MAX(0, nodes[from].tr_cap);
		
		if (nodes[to].t_cap > 0)
			flow -= MIN(nodes[to].t_cap - nodes[to].tr_cap, nodes[to].t_cap);
		else
			flow += MAX(0, nodes[to].tr_cap);
		
		captype eflow, excess;
		a_rev = a->sister;
		eflow = a->e_cap - a->r_cap;
		
		captype cap = a->e_cap + added_cap;
		captype rev_cap = a_rev->e_cap + added_rev_cap;
		
		if (eflow == 0) {
			if (((a->e_cap == 0) && (cap > a->e_cap)) || ((a_rev->e_cap == 0) && (rev_cap > a_rev->e_cap))) {
				mark_node(from);
				mark_node(to);
			}
			a->r_cap += (cap - a->e_cap);
			a_rev->r_cap += (rev_cap - a_rev->e_cap);
		} else if (eflow > 0) {
			markFlowEdge(a);//IMPORTANT: the capacity of this edge has changed (possibly enough to reduce the total flow, possibly not).
			//REGARDLESS of whether this is a capacity increase or decrease, and regardless of whether there is now enough capacity for the flow on this edge, because this edge previously carried flow,
			//the assignment of which DGL edge is assigned that flow may have changed. So the edge must be marked.

			if (cap >= a->e_cap) {
				if (eflow >= a->e_cap) {
					mark_node(from);
					mark_node(to);
				}
				a->r_cap += (cap - a->e_cap);
				a_rev->r_cap += (rev_cap - a_rev->e_cap);
			} else if (cap < a->e_cap) {
				if (eflow <= cap) {
					a->r_cap -= (a->e_cap - cap);
					a_rev->r_cap -= (a_rev->e_cap - rev_cap);
					
					if (eflow == cap) {
						mark_node(from);
						mark_node(to);
					}
				} else {
					excess = eflow - cap;
					a->r_cap = 0;
					a_rev->r_cap = rev_cap + cap;
					
					nodes[from].tr_cap += excess;
					nodes[to].tr_cap -= excess;
					if (!nodes[to].in_t_edges_set && nodes[to].tr_cap < 0) {
						nodes[to].in_t_edges_set = 1;
						t_edge_nodes.push_back(to);
					}
					if (!nodes[from].in_s_edges_set && nodes[from].tr_cap > 0) {
						nodes[from].in_s_edges_set = 1;
						s_edge_nodes.push_back(from);
					}
					
					if (nodes[from].tr_cap != 0 || nodes[from].parent == a)
						mark_node(from);
					if (nodes[to].tr_cap != 0 || nodes[to].parent == a)
						mark_node(to);
				}
			}
		} else {
			markFlowEdge(a_rev);//IMPORTANT: the capacity of this edge has changed (possibly enough to reduce the total flow, possibly not).
			//REGARDLESS of whether this is a capacity increase or decrease, and regardless of whether there is now enough capacity for the flow on this edge, because this edge previously carried flow,
			//the assignment of which DGL edge is assigned that flow may have changed. So the edge must be marked.
			eflow *= -1;
			if (rev_cap >= a_rev->e_cap) {
				if (eflow == a_rev->e_cap) {
					mark_node(from);
					mark_node(to);
				}
				a->r_cap += (cap - a->e_cap);
				a_rev->r_cap += (rev_cap - a_rev->e_cap);
			} else if (rev_cap < a_rev->e_cap) {
				if (eflow <= rev_cap) {
					a->r_cap -= (a->e_cap - cap);
					a_rev->r_cap -= (a_rev->e_cap - rev_cap);
					
					if (eflow == cap) {
						mark_node(from);
						mark_node(to);
					}
				} else {
					excess = eflow - rev_cap;
					a_rev->r_cap = 0;
					a->r_cap = rev_cap + cap;
					
					nodes[from].tr_cap -= excess;
					nodes[to].tr_cap += excess;
					if (!nodes[to].in_t_edges_set && nodes[to].tr_cap < 0) {
						nodes[to].in_t_edges_set = 1;
						t_edge_nodes.push_back(to);
					}
					if (!nodes[from].in_s_edges_set && nodes[from].tr_cap > 0) {
						nodes[from].in_s_edges_set = 1;
						s_edge_nodes.push_back(from);
					}
					if (nodes[from].tr_cap != 0 || nodes[from].parent == a)
						mark_node(from);
					if (nodes[to].tr_cap != 0 || nodes[to].parent == a)
						mark_node(to);
				}
			}
		}
		a->e_cap = cap;
		a_rev->e_cap = rev_cap;
		
		// Modifying flow value
		
		if (nodes[from].t_cap > 0)
			flow += MIN(nodes[from].t_cap - nodes[from].tr_cap, nodes[from].t_cap);
		else
			flow -= MAX(0, nodes[from].tr_cap);
		
		if (nodes[to].t_cap > 0)
			flow += MIN(nodes[to].t_cap - nodes[to].tr_cap, nodes[to].t_cap);
		else
			flow -= MAX(0, nodes[to].tr_cap);
	}
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::edit_edge(node_id from, node_id to, captype cap, captype rev_cap, arc*ac) {
	arc *a, *a_rev;
	if(!ac){
		assert(false);
		a = nodes[from].first;

		while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
			a = a->next;
	}else{
#ifndef NDEBUG
		a = nodes[from].first;
		while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
			a = a->next;
		assert(a==ac);
#endif
		a = ac;
		assert(a->head==&nodes[to]);
	}
	
	if (a->head != &nodes[to]) {
		throw std::invalid_argument("Specified edge doesn't exist");
	} else {
		// Modifying flow value 
		
		if (nodes[from].t_cap > 0)
			flow -= MIN(nodes[from].t_cap - nodes[from].tr_cap, nodes[from].t_cap);
		else
			flow += MAX(0, nodes[from].tr_cap);
		
		if (nodes[to].t_cap > 0)
			flow -= MIN(nodes[to].t_cap - nodes[to].tr_cap, nodes[to].t_cap);
		else
			flow += MAX(0, nodes[to].tr_cap);
		
		captype eflow, excess;
		a_rev = a->sister;
		eflow = a->e_cap - a->r_cap;
		
		if (eflow == 0) {
			if (((a->e_cap == 0) && (cap > a->e_cap)) || ((a_rev->e_cap == 0) && (rev_cap > a_rev->e_cap))) {
				mark_node(from);
				mark_node(to);
			}
			a->r_cap += (cap - a->e_cap);
			a_rev->r_cap += (rev_cap - a_rev->e_cap);
		} else if (eflow > 0) {
			markFlowEdge(a);//IMPORTANT: the capacity of this edge has changed (possibly enough to reduce the total flow, possibly not).
			//REGARDLESS of whether this is a capacity increase or decrease, and regardless of whether there is now enough capacity for the flow on this edge, because this edge previously carried flow,
			//the assignment of which DGL edge is assigned that flow may have changed. So the edge must be marked.

			if (cap >= a->e_cap) {
				if (eflow >= a->e_cap) {
					mark_node(from);
					mark_node(to);
				}
				a->r_cap += (cap - a->e_cap);
				a_rev->r_cap += (rev_cap - a_rev->e_cap);
			} else if (cap < a->e_cap) {
				if (eflow <= cap) {
					a->r_cap -= (a->e_cap - cap);
					a_rev->r_cap -= (a_rev->e_cap - rev_cap);
					
					if (eflow == cap) {
						mark_node(from);
						mark_node(to);
					}
				} else {
					excess = eflow - cap;
					a->r_cap = 0;
					a_rev->r_cap = rev_cap + cap;
					
					nodes[from].tr_cap += excess;
					nodes[to].tr_cap -= excess;
					if (!nodes[to].in_t_edges_set && nodes[to].tr_cap < 0) {
						nodes[to].in_t_edges_set = 1;
						t_edge_nodes.push_back(to);
					}
					if (!nodes[from].in_s_edges_set && nodes[from].tr_cap > 0) {
						nodes[from].in_s_edges_set = 1;
						s_edge_nodes.push_back(from);
					}
					if (nodes[from].tr_cap != 0 || nodes[from].parent == a)
						mark_node(from);
					if (nodes[to].tr_cap != 0 || nodes[to].parent == a)
						mark_node(to);
				}
			}
		} else {
			markFlowEdge(a_rev);//IMPORTANT: the capacity of this edge has changed (possibly enough to reduce the total flow, possibly not).
			//REGARDLESS of whether this is a capacity increase or decrease, and regardless of whether there is now enough capacity for the flow on this edge, because this edge previously carried flow,
			//the assignment of which DGL edge is assigned that flow may have changed. So the edge must be marked.

			eflow *= -1;
			if (rev_cap >= a_rev->e_cap) {
				if (eflow == a_rev->e_cap) {
					mark_node(from);
					mark_node(to);
				}
				a->r_cap += (cap - a->e_cap);
				a_rev->r_cap += (rev_cap - a_rev->e_cap);
			} else if (rev_cap < a_rev->e_cap) {
				if (eflow <= rev_cap) {
					a->r_cap -= (a->e_cap - cap);
					a_rev->r_cap -= (a_rev->e_cap - rev_cap);
					
					if (eflow == cap) {
						mark_node(from);
						mark_node(to);
					}
				} else {
					excess = eflow - rev_cap;
					a_rev->r_cap = 0;
					a->r_cap = rev_cap + cap;
					
					nodes[from].tr_cap -= excess;
					nodes[to].tr_cap += excess;
					if (!nodes[to].in_t_edges_set && nodes[to].tr_cap < 0) {
						nodes[to].in_t_edges_set = 1;
						t_edge_nodes.push_back(to);
					}
					if (!nodes[from].in_s_edges_set && nodes[from].tr_cap > 0) {
						nodes[from].in_s_edges_set = 1;
						s_edge_nodes.push_back(from);
					}
					if (nodes[from].tr_cap != 0 || nodes[from].parent == a)
						mark_node(from);
					if (nodes[to].tr_cap != 0 || nodes[to].parent == a)
						mark_node(to);
				}
			}
		}
		a->e_cap = cap;
		a_rev->e_cap = rev_cap;
		
		// Modifying flow value 
		
		if (nodes[from].t_cap > 0)
			flow += MIN(nodes[from].t_cap - nodes[from].tr_cap, nodes[from].t_cap);
		else
			flow -= MAX(0, nodes[from].tr_cap);
		
		if (nodes[to].t_cap > 0)
			flow += MIN(nodes[to].t_cap - nodes[to].tr_cap, nodes[to].t_cap);
		else
			flow -= MAX(0, nodes[to].tr_cap);
	}
}

/***********************************************************************************************/
/***********************************************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::edit_edge_wt(node_id from, node_id to, captype cap, captype rev_cap, arc * ac) {
	arc *a, *a_rev;
	if(!ac){
		assert(false);
		a = nodes[from].first;

		while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
			a = a->next;
	}else{
	#ifndef NDEBUG
			a = nodes[from].first;
			while ((a != NULL) && (a != a->next) && (a->head != &nodes[to]))
				a = a->next;
			assert(a==ac);
	#endif
			a = ac;
			assert(a->head==&nodes[to]);
		}
	if (a->head != &nodes[to]) {
		throw std::invalid_argument("Specified edge doesn't exist");
	} else {
		if (nodes[from].t_cap > 0)
			flow -= MIN(nodes[from].t_cap - nodes[from].tr_cap, nodes[from].t_cap);
		else
			flow += MAX(0, nodes[from].tr_cap);
		if (nodes[to].t_cap > 0)
			flow -= MIN(nodes[to].t_cap - nodes[to].tr_cap, nodes[to].t_cap);
		else
			flow += MAX(0, nodes[to].tr_cap);
		
		captype eflow, excess;
		a_rev = a->sister;
		eflow = a->e_cap - a->r_cap;
		if (eflow>0){
			markFlowEdge(a);//IMPORTANT: the capacity of this edge has changed (possibly enough to reduce the total flow, possibly not).
			//REGARDLESS of whether this is a capacity increase or decrease, and regardless of whether there is now enough capacity for the flow on this edge, because this edge previously carried flow,
			//the assignment of which DGL edge is assigned that flow may have changed. So the edge must be marked.
		}else if (eflow<0){
			markFlowEdge(a_rev);//IMPORTANT: the capacity of this edge has changed (possibly enough to reduce the total flow, possibly not).
			//REGARDLESS of whether this is a capacity increase or decrease, and regardless of whether there is now enough capacity for the flow on this edge, because this edge previously carried flow,
			//the assignment of which DGL edge is assigned that flow may have changed. So the edge must be marked.

		}
		if ((eflow > 0 && eflow > cap) || (eflow < 0 && -eflow > rev_cap)) {
			if (eflow > 0) {
				excess = eflow - cap;
				a->r_cap = 0;
				a_rev->r_cap = rev_cap + cap;
			} else {
				excess = eflow + rev_cap;
				a->r_cap = rev_cap + cap;
				a_rev->r_cap = 0;
			}
			
			nodes[from].tr_cap += excess;
			nodes[to].tr_cap -= excess;
			if (!nodes[to].in_t_edges_set && nodes[to].tr_cap < 0) {
				nodes[to].in_t_edges_set = 1;
				t_edge_nodes.push_back(to);
			}
			if (!nodes[from].in_s_edges_set && nodes[from].tr_cap > 0) {
				nodes[from].in_s_edges_set = 1;
				s_edge_nodes.push_back(from);
			}
		} else {
			a->r_cap += (cap - a->e_cap);
			a_rev->r_cap += (rev_cap - a_rev->e_cap);
		}
		
		a->e_cap = cap;
		a_rev->e_cap = rev_cap;
		
		// Modifying flow value 
		
		if (nodes[from].t_cap > 0)
			flow += MIN(nodes[from].t_cap - nodes[from].tr_cap, nodes[from].t_cap);
		else
			flow -= MAX(0, nodes[from].tr_cap);
		
		if (nodes[to].t_cap > 0)
			flow += MIN(nodes[to].t_cap - nodes[to].tr_cap, nodes[to].t_cap);
		else
			flow -= MAX(0, nodes[to].tr_cap);
		
	}
}

/*********************************************************************************************/
/***** MAXFLOW CODE *******/
/*********************************************************************************************/

/*
 Functions for processing active list.
 i->next points to the next node in the list
 (or to i, if i is the last node in the list).
 If i->next is NULL iff i is not in the list.

 There are two queues. Active nodes are added
 to the end of the second queue and read from
 the front of the first queue. If the first queue
 is empty, it is replaced by the second queue
 (and the second queue becomes empty).
 */

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::set_active(node *i) {
	if (!i->next) {
		/* it's not in the list yet */
		if (queue_last[1])
			queue_last[1]->next = i;
		else
			queue_first[1] = i;
		queue_last[1] = i;
		i->next = i;
	}
}

/*
 Returns the next active node.
 If it is connected to the sink, it stays in the list,
 otherwise it is removed from the list
 */
template<typename captype, typename tcaptype, typename flowtype>
inline typename Graph<captype, tcaptype, flowtype>::node* Graph<captype, tcaptype, flowtype>::next_active() {
	node *i;
	
	while (1) {
		if (!(i = queue_first[0])) {
			queue_first[0] = i = queue_first[1];
			queue_last[0] = queue_last[1];
			queue_first[1] = NULL;
			queue_last[1] = NULL;
			if (!i)
				return NULL;
		}
		
		/* remove it from the active list */
		if (i->next == i)
			queue_first[0] = queue_last[0] = NULL;
		else
			queue_first[0] = i->next;
		i->next = NULL;
		
		/* a node in the list is active iff it has a parent */
		if (i->parent)
			return i;
	}
}

/***********************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::set_orphan_front(node *i) {
	nodeptr *np;
	i->parent = ORPHAN;
	np = nodeptr_block->New();
	np->ptr = i;
	np->next = orphan_first;
	orphan_first = np;
}

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::set_orphan_rear(node *i) {
	nodeptr *np;
	i->parent = ORPHAN;
	np = nodeptr_block->New();
	np->ptr = i;
	if (orphan_last)
		orphan_last->next = np;
	else
		orphan_first = np;
	orphan_last = np;
	np->next = NULL;
}

/***********************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
inline void Graph<captype, tcaptype, flowtype>::add_to_changed_list(node *i) {
	if (keep_changed_list && !i->is_in_changed_list) {
		node_id* ptr = changed_list->New();
		*ptr = (node_id) (i - nodes);
		i->is_in_changed_list = true;
	}
	partitionAssignment(i);
}

/***********************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::maxflow_init() {
	node *i;
	
	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = NULL;
	
	TIME = 0;
	
	for (i = nodes; i < node_last; i++) {
		i->next = NULL;
		i->is_marked = 0;
		i->is_in_changed_list = 0;
		i->TS = TIME;
		if (i->tr_cap > 0) {
			/* i is connected to the source */
			i->is_sink = 0;
			i->parent = TERMINAL;
			partitionAssignment(i);
			set_active(i);
			i->DIST = 1;
		} else if (i->tr_cap < 0) {
			/* i is connected to the sink */
			i->is_sink = 1;
			i->parent = TERMINAL;
			partitionAssignment(i);
			set_active(i);
			i->DIST = 1;
		} else {
			i->parent = NULL;
		}
	}
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::maxflow_reuse_trees_init() {
	node* i;
	node* j;
	node* queue = queue_first[1];
	arc* a;
	nodeptr* np;
	
	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = orphan_last = NULL;
	
	TIME++;
	
	while ((i = queue)) {
		queue = i->next;
		if (queue == i)
			queue = NULL;
		i->next = NULL;
		i->is_marked = 0;
		set_active(i);
		
		if (i->tr_cap == 0) {
			if (i->parent)
				set_orphan_rear(i);
			continue;
		}
		
		if (i->tr_cap > 0) {
			if (!i->parent || i->is_sink) {
				i->is_sink = 0;
				for (a = i->first; a; a = a->next) {
					j = a->head;
					if (!j->is_marked) {
						if (j->parent == a->sister)
							set_orphan_rear(j);
						if (j->parent && j->is_sink && a->r_cap > 0)
							set_active(j);
					}
				}
				add_to_changed_list(i);

			}
		} else {
			if (!i->parent || !i->is_sink) {
				i->is_sink = 1;
				for (a = i->first; a; a = a->next) {
					j = a->head;
					if (!j->is_marked) {
						if (j->parent == a->sister)
							set_orphan_rear(j);
						if (j->parent && !j->is_sink && a->sister->r_cap > 0)
							set_active(j);
					}
				}
				add_to_changed_list(i);

			}
		}
		i->parent = TERMINAL;
		i->TS = TIME;
		i->DIST = 1;

	}
	
	//test_consistency();
	
	/* adoption */
	while ((np = orphan_first)) {
		orphan_first = np->next;
		i = np->ptr;
		nodeptr_block->Delete(np);
		if (!orphan_first)
			orphan_last = NULL;
		if (i->is_sink)
			process_sink_orphan(i);
		else
			process_source_orphan(i);
	}
	/* adoption end */

	//test_consistency();
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::augment(arc *middle_arc) {
	node *i;
	arc *a;
	tcaptype bottleneck;
	
	/* 1. Finding bottleneck capacity */
	/* 1a - the source tree */
	bottleneck = middle_arc->r_cap;
	for (i = middle_arc->sister->head;; i = a->head) {
		a = i->parent;
		if (a == TERMINAL)
			break;
		if (bottleneck > a->sister->r_cap)
			bottleneck = a->sister->r_cap;
	}
	if (bottleneck > i->tr_cap)
		bottleneck = i->tr_cap;
	/* 1b - the sink tree */
	for (i = middle_arc->head;; i = a->head) {
		a = i->parent;
		if (a == TERMINAL)
			break;
		if (bottleneck > a->r_cap)
			bottleneck = a->r_cap;
	}
	if (bottleneck > -i->tr_cap)
		bottleneck = -i->tr_cap;
	
	/* 2. Augmenting */
	/* 2a - the source tree */
	middle_arc->sister->r_cap += bottleneck;
	middle_arc->r_cap -= bottleneck;
	markFlowEdge(middle_arc);
	for (i = middle_arc->sister->head;; i = a->head) {
		a = i->parent;
		if (a == TERMINAL)
			break;
		a->r_cap += bottleneck;
		a->sister->r_cap -= bottleneck;
		markFlowEdge(a->sister);
		if (!a->sister->r_cap) {
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i->tr_cap -= bottleneck;
	assert(bottleneck >= 0);
	if (!i->tr_cap) {
		set_orphan_front(i); // add i to the beginning of the adoption list
	} else if (!i->in_t_edges_set && i->tr_cap < 0) {
		
		i->in_t_edges_set = 1;
		t_edge_nodes.push_back(i - nodes);
	}
	
	/* 2b - the sink tree */
	for (i = middle_arc->head;; i = a->head) {
		a = i->parent;
		if (a == TERMINAL)
			break;
		a->sister->r_cap += bottleneck;
		a->r_cap -= bottleneck;
		markFlowEdge(a);
		if (!a->r_cap) {
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	assert(bottleneck >= 0);
	i->tr_cap += bottleneck;
	if (!i->tr_cap) {
		set_orphan_front(i); // add i to the beginning of the adoption list
	} else if (!i->in_s_edges_set && i->tr_cap > 0) {
		i->in_s_edges_set = 1;
		s_edge_nodes.push_back(i - nodes);
	}
	
	flow += bottleneck;
}

/***********************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::process_source_orphan(node *i) {
	node *j;
	arc *a0, *a0_min = NULL, *a;
	int d, d_min = INFINITE_D;
	
	/* trying to find a new parent */
	for (a0 = i->first; a0; a0 = a0->next)
		if (a0->sister->r_cap) {
			j = a0->head;
			if (!j->is_sink && (a = j->parent)) {
				/* checking the origin of j */
				d = 0;
				while (1) {
					if (j->TS == TIME) {
						d += j->DIST;
						break;
					}
					a = j->parent;
					d++;
					if (a == TERMINAL) {
						j->TS = TIME;
						j->DIST = 1;
						break;
					}
					if (a == ORPHAN) {
						d = INFINITE_D;
						break;
					}
					j = a->head;
				}
				if (d < INFINITE_D) /* j originates from the source - done */
				{
					if (d < d_min) {
						a0_min = a0;
						d_min = d;
					}
					/* set marks along the path */
					for (j = a0->head; j->TS != TIME; j = j->parent->head) {
						j->TS = TIME;
						j->DIST = d--;
					}
				}
			}
		}
	
	if ((i->parent = a0_min)) //is this assignment intentional
			{
		i->TS = TIME;
		i->DIST = d_min + 1;
	} else {
		/* no parent is found */
		add_to_changed_list(i);

		/* process neighbors */
		for (a0 = i->first; a0; a0 = a0->next) {
			j = a0->head;
			if (!j->is_sink && (a = j->parent)) {
				if (a0->sister->r_cap)
					set_active(j);
				if (a != TERMINAL && a != ORPHAN && a->head == i) {
					set_orphan_rear(j); // add j to the end of the adoption list
				}
			}
		}
	}
}

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::process_sink_orphan(node *i) {
	node *j;
	arc *a0, *a0_min = NULL, *a;
	int d, d_min = INFINITE_D;
	
	/* trying to find a new parent */
	for (a0 = i->first; a0; a0 = a0->next)
		if (a0->r_cap) {
			j = a0->head;
			if (j->is_sink && (a = j->parent)) {
				/* checking the origin of j */
				d = 0;
				while (1) {
					if (j->TS == TIME) {
						d += j->DIST;
						break;
					}
					a = j->parent;
					d++;
					if (a == TERMINAL) {
						j->TS = TIME;
						j->DIST = 1;
						break;
					}
					if (a == ORPHAN) {
						d = INFINITE_D;
						break;
					}
					j = a->head;
				}
				if (d < INFINITE_D) /* j originates from the sink - done */
				{
					if (d < d_min) {
						a0_min = a0;
						d_min = d;
					}
					/* set marks along the path */
					for (j = a0->head; j->TS != TIME; j = j->parent->head) {
						j->TS = TIME;
						j->DIST = d--;
					}
				}
			}
		}
	
	if ((i->parent = a0_min)) //is this assignment intentional?
			{
		i->TS = TIME;
		i->DIST = d_min + 1;
	} else {
		/* no parent is found */
		add_to_changed_list(i);

		/* process neighbors */
		for (a0 = i->first; a0; a0 = a0->next) {
			j = a0->head;
			if (j->is_sink && (a = j->parent)) {
				if (a0->r_cap)
					set_active(j);
				if (a != TERMINAL && a != ORPHAN && a->head == i) {
					set_orphan_rear(j); // add j to the end of the adoption list
				}
			}
		}
	}
}

/***********************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
flowtype Graph<captype, tcaptype, flowtype>::maxflow(bool reuse_trees, Block<node_id>** _changed_list) {
	node *i, *j, *current_node = NULL;
	arc *a;
	nodeptr *np, *np_next;
	
	if (!nodeptr_block) {
		nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);
	}
	
	if (maxflow_iteration == 0) {
		reuse_trees = false;
		_changed_list = NULL;
	}
	
	if (_changed_list) {
		keep_changed_list = true;
		if (changed_list)
			changed_list->Reset();
		else
			changed_list = new Block<node_id>(NODEPTR_BLOCK_SIZE, error_function);
		*_changed_list = changed_list;
	} else
		keep_changed_list = false;
	
	if (reuse_trees)
		maxflow_reuse_trees_init();
	else
		maxflow_init();
	
	// main loop
	while (1) {
		// test_consistency(current_node);
		
		if ((i = current_node)) {
			i->next = NULL; /* remove active flag */
			if (!i->parent)
				i = NULL;
		}
		if (!i) {
			if (!(i = next_active()))
				break;
		}
		
		/* growth */
		if (!i->is_sink) {
			/* grow source tree */
			for (a = i->first; a; a = a->next)
				if (a->r_cap) {
					j = a->head;
					if (!j->parent) {
						j->is_sink = 0;
						j->parent = a->sister;
						j->TS = i->TS;
						j->DIST = i->DIST + 1;
						set_active(j);
						add_to_changed_list(j);

					} else if (j->is_sink)
						break;
					else if (j->TS <= i->TS && j->DIST > i->DIST) {
						/* heuristic - trying to make the distance from j to the source shorter */
						j->parent = a->sister;
						j->TS = i->TS;
						j->DIST = i->DIST + 1;
					}
				}
		} else {
			/* grow sink tree */
			for (a = i->first; a; a = a->next)
				if (a->sister->r_cap) {
					j = a->head;
					if (!j->parent) {
						j->is_sink = 1;
						j->parent = a->sister;
						j->TS = i->TS;
						j->DIST = i->DIST + 1;
						set_active(j);
						add_to_changed_list(j);

					} else if (!j->is_sink) {
						a = a->sister;
						break;
					} else if (j->TS <= i->TS && j->DIST > i->DIST) {
						/* heuristic - trying to make the distance from j to the sink shorter */
						j->parent = a->sister;
						j->TS = i->TS;
						j->DIST = i->DIST + 1;
					}
				}
		}
		
		TIME++;
		
		if (a) {
			i->next = i; /* set active flag */
			current_node = i;
			
			/* augmentation */
			augment(a);
			/* augmentation end */

			/* adoption */
			while ((np = orphan_first)) {
				np_next = np->next;
				np->next = NULL;
				
				while ((np = orphan_first)) {
					orphan_first = np->next;
					i = np->ptr;
					nodeptr_block->Delete(np);
					if (!orphan_first)
						orphan_last = NULL;
					if (i->is_sink)
						process_sink_orphan(i);
					else
						process_source_orphan(i);
				}
				
				orphan_first = np_next;
			}
			/* adoption end */
		} else
			current_node = NULL;
	}
	// test_consistency();
	
	if (!reuse_trees || (maxflow_iteration % 64) == 0) {
		delete nodeptr_block;
		nodeptr_block = NULL;
	}
	
	maxflow_iteration++;
	return flow;
}

/***********************************************************************/

template<typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::test_consistency(node* current_node) {
	node *i;
	arc *a;
	int r;
	int num1 = 0, num2 = 0;
	
	// test whether all nodes i with i->next!=NULL are indeed in the queue
	for (i = nodes; i < node_last; i++) {
		if (i->next || i == current_node)
			num1++;
	}
	for (r = 0; r < 3; r++) {
		i = (r == 2) ? current_node : queue_first[r];
		if (i)
			for (;; i = i->next) {
				num2++;
				if (i->next == i) {
					if (r < 2)
						assert(i == queue_last[r]);
					else
						assert(i == current_node);
					break;
				}
			}
	}
	assert(num1 == num2);
	
	for (i = nodes; i < node_last; i++) {
		// test whether all edges in seach trees are non-saturated
		if (i->parent == NULL) {
		} else if (i->parent == ORPHAN) {
		} else if (i->parent == TERMINAL) {
			if (!i->is_sink)
				assert(i->tr_cap > 0);
			else
				assert(i->tr_cap < 0);
		} else {
			if (!i->is_sink)
				assert(i->parent->sister->r_cap > 0);
			else
				assert(i->parent->r_cap > 0);
		}
		// test whether passive nodes in search trees have neighbors in
		// a different tree through non-saturated edges
		if (i->parent && !i->next) {
			if (!i->is_sink) {
				assert(i->tr_cap >= 0);
				for (a = i->first; a; a = a->next) {
					if (a->r_cap > 0)
						assert(a->head->parent && !a->head->is_sink);
				}
			} else {
				assert(i->tr_cap <= 0);
				for (a = i->first; a; a = a->next) {
					if (a->sister->r_cap > 0)
						assert(a->head->parent && a->head->is_sink);
				}
			}
		}
		// test marking invariants
		if (i->parent && i->parent != ORPHAN && i->parent != TERMINAL) {
			assert(i->TS <= i->parent->head->TS);
			if (i->TS == i->parent->head->TS)
				assert(i->DIST > i->parent->head->DIST);
		}
	}
}
}
;
#endif
