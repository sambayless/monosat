/*
 * DynamicConnect.h
 *
 *  Created on: Mar 26, 2014
 *      Author: sam
 */

#ifndef THORUPDYNAMICCONNECT_H_
#define THORUPDYNAMICCONNECT_H_

#include "DynamicConnectivityImpl.h"
#include "alg/EulerTree.h"
#include <cmath>
#ifndef NDEBUG
#include "mtl/Sort.h"
#include "NaiveDynamicConnectivity.h"
#endif

class ThorupDynamicConnectivity:public DynamicConnectivityImpl{

#ifndef NDEBUG
	NaiveDynamicConnectivity dbg;
#endif

struct Edge{
	int edgeID;
	int from;
	int to;
	int in_forest:1;
	int enabled:1;
	int incident_from:1;
	int incident_to:1;
	int level:28;
	Edge():edgeID(-1),from(-1),to(-1),in_forest(false),enabled(false),incident_from(false),incident_to(false),level(0){

	}
};
public:
vec<Edge> edges;
vec<bool> seen;
vec<int> t_seen;

int nodes;

vec<EulerTree> et;
//This is a list of _ALL_ incident edges to each node, except the ones that are actually in the forest, regardless of their level.
vec<vec<int> > incident_edges;
private:

int levels;

bool insert(int edgeID){
	Edge & e = edges[edgeID];
	assert(e.edgeID==edgeID);
	e.level=levels-1;
	//fix these later


	if (!et.last().connected(e.from,e.to)){
		et.last().link(e.from,e.to,edgeID);
		e.in_forest=true;
		dbg_levels();
		return true;
	}else{
		e.in_forest=false;
		if(!e.incident_from){
			e.incident_from=true;
			assert(!incident_edges[e.from].contains(edgeID));
			incident_edges[e.from].push(edgeID);
		}
		if(!e.incident_to){
			e.incident_to=true;
			assert(!incident_edges[e.to].contains(edgeID));
			incident_edges[e.to].push(edgeID);
		}
		dbg_levels();
		return false;
	}


}

bool insert_unchecked(int edgeID, bool already_connected){
	Edge & e = edges[edgeID];
	assert(e.edgeID==edgeID);
	e.level=levels-1;
	//fix these later
	//assert(!incident_edges[e.from].contains(edgeID));
	//assert(!incident_edges[e.to].contains(edgeID));
	assert(et.last().connected(e.from,e.to)==already_connected);
	if (!already_connected){
		et.last().link(e.from,e.to,edgeID);
		e.in_forest=true;
		dbg_levels();
		return true;
	}else{
		e.in_forest=false;
		if(!e.incident_from){
			e.incident_from=true;
			assert(!incident_edges[e.from].contains(edgeID));
			incident_edges[e.from].push(edgeID);
		}
		if(!e.incident_to){
			e.incident_to=true;
			assert(!incident_edges[e.to].contains(edgeID));
			incident_edges[e.to].push(edgeID);
		}
		dbg_levels();
		return false;
	}

}

void dbg_tree(int level){
#ifndef NDEBUG
	for(int n = 0;n<nodes;n++){
		int sz = et[level].getFullTreeSize(n);
		assert(sz<=pow(2,level+1));//invariant from paper (note that we are adding one because our levels are offset by one (and reversed!) from the papers)

		//ensure that the forest at this level is contained in the forest at each higher level
		for(auto it = et[level].begin_half_edge_tour(n);it != et[level].end_half_edge_tour();++it){
			int cur_edge = *it;
			Edge & treeEdge = edges[cur_edge];
			for(int l = level+1;l<levels;l++){
				assert(et[l].connected(n,treeEdge.from));
				assert(et[l].connected(n,treeEdge.to));
				assert(et[l].edgeInTree(treeEdge.edgeID));
			}
		}
	}
#endif
}
void dbg_checkGraph(){
#ifndef NDEBUG
	for(int i = 0;i<nodes;i++){
	/*	for(int e:incident_edges[i]){
			assert(!edges[e].in_forest);
		}*/
		for(int j = 0;j<nodes;j++){
			assert(connected(i,j)==dbg.connected(i,j));

		}
	}
	for(int i = 0;i<edges.size();i++){
		if(edges[i].enabled){
			assert(et.last().connected(edges[i].from,edges[i].to));
			//This is invariant 1 from the paper.
			assert(et[edges[i].level].connected(edges[i].from,edges[i].to));

			if(edges[i].in_forest){
				for(int j = edges[i].level;j<levels;j++){
					assert(et[j].connected(edges[i].from,edges[i].to));
					assert(et[j].edgeInTree(edges[i].edgeID));//because the forests are subsets of
				}

			}
		}else{
			assert(!edges[i].enabled);
		}
	}
	dbg_levels();
#endif
}

void dbg_incident(){
#ifndef NDEBUG
	for(int n = 0;n<nodes;n++){
		vec<int> & incident = incident_edges[n];
		vec<int> seen;
		for(int edge:incident){
			//ensure no duplicates
			assert(!seen.contains(edge));
			seen.push(edge);
			Edge & e =edges[edge];
			assert(e.incident_from || e.incident_to);
			if(n==e.from){
				assert(e.incident_from);
			}else{
				assert(n==e.to);
				assert(e.incident_to);
			}
		}

	}

	for(Edge & e:edges){
		if(e.incident_from){
			assert(incident_edges[e.from].contains(e.edgeID));
		}

		if(e.incident_to){
			assert(incident_edges[e.to].contains(e.edgeID));
		}
	}



#endif
}

void dbg_print(){
#ifndef NDEBUG
	vec<bool> seen;
	seen.growTo(nodes);
	for(int n = 0;n<nodes;n++){
		if(!seen[n]){
			seen[n]=true;
			vec<int> e;
			getConnectedComponentEdges(n,e);

			vec<int> node_set;
			getConnectedComponent(n,node_set);
			sort(node_set);
			if(node_set.size()>1){
			for(int n:node_set){
				seen[n]=true;
				printf("%d,",n);
			}
			printf("\n");

			for(int i =0;i<e.size();i++){
				printf("(%d->%d)", edges[e[i]].from,edges[e[i]].to );
			}
			printf("\n");
			}
		}
	}
#endif
}

bool visit(int w,int otherTreeVertex, int removedEdge, int level, int & replacementEdge){
	if(seen[w]){
		return false;
	}
	seen[w]=true;
	int k,j=0;

	//Note that because I am only storing one list of incident edges per vertex - rather than one per vertex per level -
	//this loop may do more work than it needs to. In sparsely connected graphs, I am hoping the improved space/locality still makes this a win.
	for(k = 0; k<incident_edges[w].size();k++){
		Edge & e = edges[incident_edges[w][k]];
		if(incident_edges[w][k]==removedEdge ||e.in_forest || !e.enabled){
			//remove this edge
			if(w==e.from){
				assert(e.incident_from);
				e.incident_from=false;
			}else{
				assert(w==e.to);
				assert(e.incident_to);
				e.incident_to=false;
			}
			continue;
		}
		assert(!edges[incident_edges[w][k]].in_forest);


		assert(e.from==w || e.to==w);
		if(e.level==level){
			int otherNode = e.from==w?e.to:e.from;
			if(et[level].connected(otherTreeVertex,otherNode)){
				assert(!e.in_forest);

				//remove this edge from v's incident edges.
				//wrong. The incident edge list must include the edges in the minimum spanning forest as well.
	/*			for (int h = k+1;h<incident_edges[v].size();h++){
					incident_edges[v][h-1]=incident_edges[v][h];
				}*/

				//remove the current edge
				if(w==e.from){
					assert(e.incident_from);
					e.incident_from=false;
				}else{
					assert(w==e.to);
					assert(e.incident_to);
					e.incident_to=false;
				}

				for(k=k+1;k<incident_edges[w].size();k++){
					Edge & e = edges[ incident_edges[w][k]];
					if(incident_edges[w][k]==removedEdge  ||e.in_forest  || !e.enabled){
						//remove this edge
						if(w==e.from){
							assert(e.incident_from);
							e.incident_from=false;
						}else{
							assert(w==e.to);
							assert(e.incident_to);
							e.incident_to=false;
						}
						continue;
					}
					incident_edges[w][j++]=incident_edges[w][k];
				}
				incident_edges[w].shrink(k-j);
				assert(!incident_edges[w].contains(e.edgeID));
				replacementEdge=e.edgeID;
				return true;
			}else{
				//note: because Tv was a spanning tree connecting its connected component, each incident edge was connected to Tv U Tu before the edge was cut, and hence is still connected to either Tv or Tu (because that edge was part of the component), and hence each incident node is in Tv.

				assert(et[level].connected(w,otherNode));

				e.level--;
				assert(e.level>=0);

				//invariant 1 in the paper... but this only holds once we have finished moving tv down a level, which we are doing while this function is being called, so it doesn't hold quite yet.
				//assert(et[e.level].connected(e.from,e.to));

				//invariant 2 in the paper
				assert(et[e.level].getFullTreeSize(otherNode) <= pow(2,e.level+1));
			}

		}else{
/*			if(e.level>level)
				assert(!et[level].connected(e.from,e.to));
			else
				assert(et[level].connected(e.from,e.to));*/
		}
		incident_edges[w][j++]=incident_edges[w][k];
	}
	incident_edges[w].shrink(k-j);
	return false;

}

bool cut(int edgeID){
	Edge & e = edges[edgeID];
	if(!e.in_forest){
		return false;//this edge is not in the top level forest (and hence also not in _any_ level forest), so we don't need to do anything special to remove it
	}

	int u = e.from;
	int v = e.to;
	assert(e.level>0);
#ifndef NDEBUG
	for(int i = 0;i<e.level;i++){
		assert(!et[i].connected(u,v));
	}
#endif
	assert(et.last().connected(u,v));
	//find a replacement edge to connect u and v, if one exists
	bool foundReplacement=false;
	for(int i = e.level;!foundReplacement && i<levels;i++){
		dbg_tree(i);

		//apparently this doesn't hold.. but shouldn't it??
		assert(et[i].connected(u,v));

		et[i].cut(edgeID);
		//int sumsize=et[i].getFullTreeSize(tu) + et[i].getFullTreeSize(tv);
		assert(et[i].getFullTreeSize(u) + et[i].getFullTreeSize(v) <= pow(2,i+1));//invariant from paper... (note that we are adding one to i, because our levels start at log2(n)-1 and decrease, rather than starting at 0 and increasing)
		assert(!et[i].connected(u,v));//these must be disconnected now that we have cut them

		if(et[i].getFullTreeSize(v)>et[i].getFullTreeSize(u)){
			std::swap(u,v);
		}
		//int sz = et[i].getFullTreeSize(v);
		assert(et[i].getFullTreeSize(v)<=et[i].getFullTreeSize(u));
		assert(et[i].getFullTreeSize(v) <= pow(2,i));

		edges[edgeID].in_forest=false;



		//The paper says we 'can' set all the in_forest edges of Tv at level i to level i-1, but it doesn't say that we must do so, or should do so, or whether we are expected to do so or not.

		//It appears that we _must_ do so, in order to avoid violating condition 1 when we move other incident edges p below
		//the easiest option is probably to just visit the whole et tour.
		int replacementEdge=-1;
		dbg_incident();
		foundReplacement=visit(v,u,edgeID,i, replacementEdge);//start with this vertex, so that we handle singleton nodes correctly.
		dbg_incident();
		//iterate through all the nodes in this tree in the forest at level i, and see if any have an incident edge at or above this level
		for(auto it = et[i].begin_half_edge_tour(v);it != et[i].end_half_edge_tour();++it){
			int cur_edge = *it;
			Edge & treeEdge = edges[cur_edge];
			//ok, we must move the level of this tree edge down. Must do this for all tv's tree edges.
			if(treeEdge.level==i){
				treeEdge.level--;
				assert(treeEdge.level>=0);
				assert(!et[treeEdge.level].connected(treeEdge.from,treeEdge.to));
				/*int oldsz1 = et[treeEdge.level].getFullTreeSize(treeEdge.from);
				int oldsz2 = et[treeEdge.level].getFullTreeSize(treeEdge.to);
				int psize =  et[i].getFullTreeSize(v);*/
				//assert(et[treeEdge.level].getFullTreeSize(treeEdge.from) + et[treeEdge.level].getFullTreeSize(treeEdge.to) <= pow(2,treeEdge.level+1));
				dbg_tree(treeEdge.level);
				dbg_printTree(treeEdge.level,treeEdge.from);
				dbg_printTree(treeEdge.level,treeEdge.to);
				dbg_printTree(i,v);
				assert(et[treeEdge.level].getFullTreeSize(treeEdge.from) + et[treeEdge.level].getFullTreeSize(treeEdge.to) <= et[i].getFullTreeSize(v) );
				et[treeEdge.level].link(treeEdge.from,treeEdge.to,treeEdge.edgeID);
				dbg_tree(treeEdge.level);
				//invariant 1 in the paper
				assert(et[treeEdge.level].connected(treeEdge.from,treeEdge.to));
				//int nsz = et[treeEdge.level].getFullTreeSize(treeEdge.from);
				//invariant 2 in the paper
				assert(et[treeEdge.level].getFullTreeSize(treeEdge.from) <= pow(2,treeEdge.level+1));
			}

			for(int n = 0;!foundReplacement && n<2;n++){
				//note: we only visit this search loop if we have not already found a replacement edge.
				int w = n? treeEdge.to:treeEdge.from;
				dbg_incident();
				static int iter = 0;
				if(++iter ==93){
					int a=1;
					int b =2;
				}
				foundReplacement = visit(w,u,edgeID,i,replacementEdge);
				dbg_incident();
			}
		}

		seen[v]=false;//needed to handle singleton nodes
		//run through the tour again and clear all the seen markers. it would be nice to avoid this...
		for(auto it = et[i].begin_half_edge_tour(v);it != et[i].end_half_edge_tour();++it){
			int cur_edge = *it;
			Edge & treeEdge = edges[cur_edge];
			seen[treeEdge.from]=false;
			seen[treeEdge.to]=false;
		}
		if(foundReplacement){
			assert(replacementEdge!=-1);
			Edge & e = edges[replacementEdge];
			assert(e.level==i);
			//this is a replacement edge to keep the two components connected
			et[e.level].link(e.from,e.to,e.edgeID);
			e.in_forest=true;
			//now link all the higher level trees
			for(int h = e.level+1;h<levels;h++){
				et[h].cut(edgeID);
				assert(!et[h].connected(u,v));
				et[h].link(e.from,e.to,e.edgeID);
			}
		}
		dbg_tree(i);
	}

	e.level=levels;
	dbg_levels();
	return !foundReplacement;

}
void dbg_levels(){
#ifndef NDEBUG
	for(int l = 0;l<levels;l++){
		dbg_tree(l);

	}
#endif

}

void dbg_printTree(int level, int fromnode){
#ifndef NDEBUG
	vec<bool> dbg_seen;
	dbg_seen.growTo(nodes);
	vec<int> treenodes;
	treenodes.push(fromnode);
	dbg_seen[fromnode]=true;
	for(auto it = et[level].begin_half_edge_tour(fromnode);it != et[level].end_half_edge_tour();++it){
		int cur_edge = *it;
		Edge & treeEdge = edges[cur_edge];
		if(!dbg_seen[treeEdge.from]){
			dbg_seen[treeEdge.from]=true;
			treenodes.push(treeEdge.from);
		}
		if(!dbg_seen[treeEdge.to]){
			dbg_seen[treeEdge.to]=true;
			treenodes.push(treeEdge.to);
		}
	}
	sort(treenodes);
	printf("tree level %d: ",level);
	for(int i:treenodes){
		printf("%d,",i);
	}
	printf("\n");
	assert(treenodes.size()==et[level].getFullTreeSize(fromnode));
#endif
}

void setEdgeLevel(int edgeID, int level){
	assert(level==edges[edgeID].level);
	//et[level].setHasIncidentEdges(edges[edgeID].from,true);
	//et[level].setHasIncidentEdges(edges[edgeID].to,true);
}

public:
ThorupDynamicConnectivity():nodes(0),levels(0){

}
bool connected(int u, int v){
#ifndef NDEBUG
/*	bool c = et.last().connected(u,v);
	bool d= dbg.connected(u,v);
	if(c!=d){
	dbg.dbg_print();
	dbg_print();
	}*/
#endif
	assert(et.last().connected(u,v)==dbg.connected(u,v));
	return et.last().connected(u,v);
}

int numComponents(){
	assert(et.last().numComponents()==dbg.numComponents());
	return et.last().numComponents();
}

int findRoot(int node){
	return et.last().findRoot(node);
}

void addNode(){
	nodes++;

	incident_edges.push();
	seen.push();
	levels = (int)(floor(log(nodes)/log(2))+1);
	et.growTo(levels);
	for(EulerTree & t:et){
		while(t.nVertices()<nodes)
			t.createVertex();
		assert(t.nVertices()==nodes);
	}
#ifndef NDEBUG
	dbg.addNode();
#endif
}

void addEdge( int from, int to,int edgeID){
	edges.growTo(edgeID+1);
	if(edges[edgeID].edgeID==-1){
		edges[edgeID].from=from;
		edges[edgeID].to=to;
		edges[edgeID].level=levels-1;
		edges[edgeID].edgeID=edgeID;
		setEdgeLevel(edgeID,edges[edgeID].level);
	#ifndef NDEBUG
		dbg.addEdge(from,to,edgeID);
	#endif
	}
}

bool edgeEnabled(int edgeid)const{
	return edges[edgeid].enabled;
}

int nNodes()const{
	return nodes;
}
int nEdges()const{
	return edges.size();
}

//If from and to are connected, finds an ARBITRARY connecting path.
void getPath(int from, int to, vec<int> & nodes_out){
	nodes_out.clear();
	if(from==to){
		nodes_out.push(from);
		return;
	}else{
		if(connected(from,to)){
			//ok, walk through the tour until a path is found. But we also want to drop from this nodes that are not needed...
			t_seen.clear();
			t_seen.growTo(nodes);
			int curLevel = 1;
			t_seen[from]=curLevel;
			nodes_out.push(from);


			for(auto it = et.last().begin_half_edge_tour(from,true);it != et.last().end_half_edge_tour();++it){
				int cur_edge = *it;
				Edge & treeEdge = edges[cur_edge];
				assert(t_seen[treeEdge.to]==curLevel || t_seen[treeEdge.from]==curLevel);
				if(t_seen[treeEdge.to] && t_seen[treeEdge.from]){
					//this is a step back up the tree
					nodes_out.pop();
					curLevel--;
					assert(t_seen[treeEdge.to]==curLevel || t_seen[treeEdge.from]==curLevel);
				}else if(t_seen[treeEdge.to]==curLevel){
					assert(!t_seen[treeEdge.from]);
					curLevel++;
					t_seen[treeEdge.from]=curLevel;
					nodes_out.push(treeEdge.from);
				}else if(t_seen[treeEdge.from]==curLevel){
					assert(!t_seen[treeEdge.to]);
					curLevel++;
					t_seen[treeEdge.to]=curLevel;
					nodes_out.push(treeEdge.to);
				}
				if(nodes_out.last()==to){
					break;
				}
			}
			t_seen.clear();
			assert(nodes_out.last()==to);
		}
	}
}
void getPathEdges(int from, int to, vec<int> & edges_out){
	edges_out.clear();
	if(from==to){

		return;
	}else{
		if(connected(from,to)){
			//ok, walk through the tour until a path is found. But we also want to drop from this nodes that are not needed...
			t_seen.clear();
			t_seen.growTo(nodes);
			int curLevel = 1;
			t_seen[from]=curLevel;


			for(auto it = et.last().begin_half_edge_tour(from,true);it != et.last().end_half_edge_tour();++it){
				int cur_edge = *it;
				Edge & treeEdge = edges[cur_edge];
				assert(t_seen[treeEdge.to]==curLevel || t_seen[treeEdge.from]==curLevel);
				if(t_seen[treeEdge.to] && t_seen[treeEdge.from]){
					//this is a step back up the tree
					edges_out.pop();
					curLevel--;
					assert(t_seen[treeEdge.to]==curLevel || t_seen[treeEdge.from]==curLevel);
				}else if(t_seen[treeEdge.to]==curLevel){
					assert(!t_seen[treeEdge.from]);
					curLevel++;
					t_seen[treeEdge.from]=curLevel;
					edges_out.push(treeEdge.edgeID);
				}else if(t_seen[treeEdge.from]==curLevel){
					assert(!t_seen[treeEdge.to]);
					curLevel++;
					t_seen[treeEdge.to]=curLevel;
					edges_out.push(treeEdge.edgeID);
				}
				if(treeEdge.from==to || treeEdge.to==to){
						break;
					}
			}
			t_seen.clear();

		}
	}
}
void getConnectedComponent(int from, vec<int> & nodes_out){
	t_seen.clear();
	t_seen.growTo(nodes);
	t_seen[from]=true;
	nodes_out.clear();
	nodes_out.push(from);
	for(auto it = et.last().begin_half_edge_tour(from);it != et.last().end_half_edge_tour();++it){
		int cur_edge = *it;
		Edge & treeEdge = edges[cur_edge];
		if(!t_seen[treeEdge.to]){
			t_seen[treeEdge.to]=true;
			nodes_out.push(treeEdge.to);
		}
		if(!t_seen[treeEdge.from]){
			t_seen[treeEdge.from]=true;
			nodes_out.push(treeEdge.from);
		}
	}
#ifndef NDEBUG
/*	for(int i =0;i<nodes;i++){
		if(t_seen[i]){
			assert(nodes_out.contains(i));
			assert(dbg.connected(from,i));
		}else{
			assert(!dbg.connected(from,i));
		}
	}*/
#endif

	t_seen.clear();
}


/**
 * If strict is true, then the edges of the component are oriented so they start at from
 */
void getConnectedComponentEdges(int from, vec<int>& edges_out, bool strict=false){
	t_seen.clear();
	t_seen.growTo(edges.size());
	for(auto it = et.last().begin_half_edge_tour(from,strict);it != et.last().end_half_edge_tour();++it){
		int cur_edge = *it;
		Edge & treeEdge = edges[cur_edge];
		if(!t_seen[treeEdge.edgeID]){
			t_seen[treeEdge.edgeID]=true;
			edges_out.push(treeEdge.edgeID);
		}
	}
	t_seen.clear();
}

//disable all edges
void clear(){

	for(int i = 0;i<edges.size();i++){
		edges[i].enabled=false;
		edges[i].in_forest=false;
		edges[i].level=levels;
	}
	for(EulerTree & t:et){
		t.clear();
	}

}

/**
 * Returns true if the set of connected components have changed
 */
bool setEdgeEnabled(int from, int to,int edgeID, bool enabled){
	static int iter = 0;
	if(++iter==531){
		int a=1;
	}
	addEdge(from,to,edgeID);
	bool changed = false;
/*	if(enabled)
		printf("Enable Edge (%d->%d)\n", edges[edgeID].from,edges[edgeID].to);
	else
		printf("Disable Edge (%d->%d)\n", edges[edgeID].from,edges[edgeID].to);*/
	dbg_checkGraph();
	assert( edges[edgeID].edgeID==edgeID);assert( edges[edgeID].to==to);assert( edges[edgeID].from==from);
	if(enabled && ! edges[edgeID].enabled){
		//printf("%d: Enable Edge (%d->%d)\n",iter, edges[edgeID].from,edges[edgeID].to);
		edges[edgeID].enabled=true;
		changed = insert(edgeID);
	}else if(!enabled && edges[edgeID].enabled){
		//printf("%d: Disable Edge (%d->%d)\n",iter, edges[edgeID].from,edges[edgeID].to);
		edges[edgeID].enabled=false;
		changed = cut(edgeID);
	}
#ifndef NDEBUG
	dbg.setEdgeEnabled(from,to,edgeID,enabled);
	dbg_checkGraph();
#endif
	return changed;
}
bool setEdgeEnabledUnchecked(int from, int to,int edgeID, bool connected){
	static int iter = 0;
	if(++iter==531){
		int a=1;
	}
	addEdge(from,to,edgeID);
	bool changed = false;

	dbg_checkGraph();
	assert( edges[edgeID].edgeID==edgeID);assert( edges[edgeID].to==to);assert( edges[edgeID].from==from);
	if(! edges[edgeID].enabled){
		//printf("%d: Enable Edge (%d->%d)\n",iter, edges[edgeID].from,edges[edgeID].to);
		edges[edgeID].enabled=true;
		changed = insert_unchecked(edgeID,connected);
	}
#ifndef NDEBUG
	dbg.setEdgeEnabled(from,to,edgeID,true);
	dbg_checkGraph();
#endif
	return changed;
}

};

#endif /* DYNAMICCONNECT_H_ */
