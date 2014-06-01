/*
 * DynamicGraph.h
 *
 *  Created on: 2013-07-15
 *      Author: sam
 */

#ifndef DYNAMICGRAPH_H_
#define DYNAMICGRAPH_H_
#include "mtl/Vec.h"
#include "GraphTheoryTypes.h"
using namespace Minisat;
#ifndef NDEBUG
#include <cstdio>
#endif
#ifndef NDEBUG
#define RECORD
#endif
template<class EdgeStatus=DefaultEdgeStatus >
class DynamicGraph{
public:
	EdgeStatus & edge_status;
	int nodes;
	int edges;
	int modifications;
	int additions;
	int deletions;
	int historyclears;
	int next_id;
	bool is_changed;
	bool allocated=false;
#ifdef RECORD
	FILE * outfile;
#endif
	struct Edge{
		int node;
		int id;
	};
	vec<vec<Edge> > adjacency;//adj list
	vec<vec<Edge> > inverted_adjacency;//adj list

	vec<vec<Edge> > adjacency_undirected;//adj list
	vec<int> weights;
	struct FullEdge{
		int from;
		int to;
		int id;
		int weight;
		FullEdge():from(-1),to(-1),id(-1),weight(1){}
		FullEdge(int from,int to, int id,int weight):from(from),to(to),id(id),weight(weight){}
	};

	vec<FullEdge> all_edges;

	struct EdgeChange{
		bool addition;
/*
		int u;//from
		int v;//top
*/
		int id;
		int mod;
		int prev_mod;
	};
	vec<EdgeChange> history;
	//vec<char> edge_status;
	DynamicGraph():edge_status(*new EdgeStatus()), nodes(0),edges(0),modifications(0),additions(0),deletions(0),historyclears(0),next_id(0),is_changed(true){
		allocated=true;
#ifdef RECORD
		outfile=nullptr;
#endif
	}
	DynamicGraph(EdgeStatus & _status):edge_status(_status), nodes(0),edges(0),modifications(0),additions(0),deletions(0),historyclears(0),next_id(0),is_changed(true){
#ifdef RECORD
		outfile=nullptr;
#endif
		allocated=false;
	}
	~DynamicGraph(){
		if(allocated){
			delete(&edge_status);
		}
	}
	void addNodes(int n){
		for(int i = 0;i<n;i++)
			addNode();
	}

	bool hasEdge(int from, int to)const{
		for(int i = 0;i< adjacency[from].size();i++){
			if(adjacency[from][i].node==to && edgeEnabled(adjacency[from][i].id)){
				return true;
			}
		}
		return false;
	}
	bool hasEdgeUndirected(int from, int to)const{
		for(int i = 0;i< adjacency_undirected[from].size();i++){
			if(adjacency_undirected[from][i].node==to && edgeEnabled(adjacency_undirected[from][i].id)){
				return true;
			}
		}
		return false;
	}
	int addNode(){

		adjacency.push();//adj list
		adjacency_undirected.push();
		inverted_adjacency.push();
		modifications++;
		additions=modifications;
		deletions=modifications;
		clearHistory();
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"node %d\n",nodes);
				fflush(outfile);
			}
#endif
		return nodes++;
	}

	bool edgeEnabled(int edgeID)const{
		assert(edgeID<edge_status.size());
		return edge_status[edgeID];
	}
	bool isEdge(int edgeID)const{
		return edgeID<all_edges.size() && all_edges[edgeID].id ==edgeID;
	}
	//Instead of actually adding and removing edges, tag each edge with an 'enabled/disabled' label, and just expect reading algorithms to check and respect that label.
	void addEdge(int from, int to, int id=-1, int weight=1){
		assert(from<nodes);
		assert(to<nodes);
		assert(from>=0);
		assert(to>=0);
		if(id<0){
			id = next_id++;
		}else{
			if(id>=next_id){
				next_id=id+1;
			}
		}

		edges=id+1;
		adjacency[from].push({to,id});
		adjacency_undirected[from].push({to,id});
		adjacency_undirected[to].push({from,id});
		edge_status.growTo(id+1);
		inverted_adjacency[to].push({from,id});
		all_edges.growTo(id+1);
		all_edges[id]={from,to,id,weight};
		weights.growTo(id+1,0);
		weights[id]=weight;
		//weights.push(weight);
		modifications++;
		additions=modifications;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"edge %d %d %d %d\n",from,to,weight, id+1);
				fflush(outfile);
			}
#endif
//		history.push({true,id,modifications});
		enableEdge(from,to,id);//default to enabled
	}
	int nEdgeIDs(){
		return all_edges.size();
	}

	int getWeight(int edgeID){
		return all_edges[edgeID].weight;
	}
	FullEdge getEdge(int id){
		return all_edges[id];
	}
	void enableEdge(int id){
		enableEdge(all_edges[id].from,all_edges[id].to,id);
	}
	void disableEdge(int id){
		disableEdge(all_edges[id].from,all_edges[id].to,id);
	}
	void enableEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		assert(isEdge(id));
		if(edge_status[id]!=true){
			edge_status.setStatus(id,true);

			modifications++;
			additions=modifications;
			history.push({true,id,modifications,additions});
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"%d\n", id+1);
				fflush(outfile);
			}
#endif
		}
	}

	bool undoEnableEdge( int id){
		assert(id>=0);
		assert(id<edge_status.size());
		assert(isEdge(id));
		if(!history.size())
			return false;

		if(history.last().addition  && history.last().id==id && history.last().mod==modifications){
			edge_status.setStatus(id,false);
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"-%d\n", id+1);
				fflush(outfile);
			}
#endif
			modifications--;
			additions = history.last().prev_mod;
			history.pop();
			return true;
		}
		return false;
	}

	void disableEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		assert(isEdge(id));
		if(edge_status[id]!=false){
			edge_status.setStatus(id,false);
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"-%d\n", id+1);
				fflush(outfile);
			}
#endif

			modifications++;

			history.push({false,id,modifications,deletions});
			deletions=modifications;
		}
	}

	bool undoDisableEdge( int id){
		assert(id>=0);
		assert(id<edge_status.size());
		assert(isEdge(id));
		if(!history.size())
			return false;

		if(!history.last().addition  && history.last().id==id && history.last().mod==modifications){
			edge_status.setStatus(id,true);
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"%d\n", id+1);
				fflush(outfile);
			}
#endif
			modifications--;
			deletions = history.last().prev_mod;
			history.pop();
			return true;
		}
		return false;
	}

	void drawFull(){
#ifndef NDEBUG
			printf("digraph{\n");
			for(int i = 0;i<nodes;i++){
				printf("n%d\n", i);
			}

			for(int i = 0;i<adjacency.size();i++){
				for(int j =0;j<adjacency[i].size();j++){
				int id  =adjacency[i][j].id;
				int u = adjacency[i][j].node;
				const char * s = "black";
				if( edgeEnabled(id))
					s="blue";
				else
					s="red";

				printf("n%d -> n%d [label=\"v%d\",color=\"%s\"]\n", i,u, id, s);
				}
			}
			printf("}\n");
#endif
		}
	//Removes _all_ edges (from, to)
	/*void removeEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		{
			vec<Edge>& adj= adjacency[from];
			int i,j = 0;
			for(i = 0;i<adj.size();i++){
				if(adj[i]==to){
					edges--;
				}else{
					adj[j++]=adj[i];
				}
			}
			adj.shrink(i-j);
		}

		modifications++;
		deletions=modifications;
		history.push({false,from,to,id,modifications});
	}*/

	bool rewindHistory(int steps){

		int cur_modifications = modifications;
		for(int i = 0;i<steps;i++){
			EdgeChange & e = history.last();
			if(e.addition){
				if(!undoEnableEdge(e.id)){
					return false;
				}
			}else{
				if(!undoDisableEdge(e.id)){
					return false;
				}
			}

		}
		assert(modifications==cur_modifications-steps);
		return true;
	}

	int getCurrentHistory(){
		return modifications;
	}

	void clearHistory(){
		if(history.size()>1000){
			history.clear();
			historyclears++;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"clearHistory\n");
				fflush(outfile);
			}
#endif
		}
	}
	//force a new modification
	void invalidate(){
		modifications++;
		additions=modifications;
		modifications++;
		deletions=modifications;

#ifdef RECORD
			if(outfile){
				fprintf(outfile,"invalidate\n");
				fflush(outfile);
			}
#endif
	}

	void markChanged(){
		is_changed=true;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"markChanged\n");
				fflush(outfile);
			}
#endif
	}
	bool changed(){
		return is_changed;
	}

	void clearChanged(){
		is_changed=false;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"clearChanged\n");
				fflush(outfile);
			}
#endif
	}
};


#endif /* DYNAMICGRAPH_H_ */
