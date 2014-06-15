/*
 * DynamicGraph.h
 *
 *  Created on: 2013-07-15
 *      Author: sam
 */

#ifndef DYNAMICGRAPH_H_
#define DYNAMICGRAPH_H_
#include <vector>
#include "graph/GraphTheoryTypes.h"

#ifndef NDEBUG
#define RECORD
#include <cstdio>
#endif
namespace dgl{
class DynamicGraph{
public:
	std::vector<bool> edge_status;
	int nodes;
	int edges;
	int modifications;
	int additions;
	int deletions;
	int historyclears;
	int next_id;
	bool is_changed;
	//bool allocated=false;
#ifdef RECORD
	FILE * outfile;
#endif
	struct Edge{
		int node;
		int id;
	};
	std::vector<std::vector<Edge> > adjacency;//adj list
	std::vector<std::vector<Edge> > inverted_adjacency;//adj list
	std::vector<std::vector<Edge> > adjacency_undirected;//adj list
	std::vector<int> weights;
	struct FullEdge{
		int from;
		int to;
		int id;
		int weight;
		FullEdge():from(-1),to(-1),id(-1),weight(1){}
		FullEdge(int from,int to, int id,int weight):from(from),to(to),id(id),weight(weight){}
	};

	std::vector<FullEdge> all_edges;

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
	std::vector<EdgeChange> history;

	DynamicGraph():nodes(0),edges(0),modifications(0),additions(0),deletions(0),historyclears(0),next_id(0),is_changed(true){
		//allocated=true;
#ifdef RECORD
		outfile=nullptr;
#endif
	}
/*	DynamicGraph(): nodes(0),edges(0),modifications(0),additions(0),deletions(0),historyclears(0),next_id(0),is_changed(true){
#ifdef RECORD
		outfile=nullptr;
#endif
		//allocated=false;
	}*/
	~DynamicGraph(){
	/*	if(allocated){
			delete(&edge_status);
		}*/
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

		adjacency.push_back({});//adj list
		adjacency_undirected.push_back({});
		inverted_adjacency.push_back({});
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

		edges=next_id;
		adjacency[from].push_back({to,id});
		adjacency_undirected[from].push_back({to,id});
		adjacency_undirected[to].push_back({from,id});
		if(edge_status.size()<=id)
			edge_status.resize(id+1);
		inverted_adjacency[to].push_back({from,id});
		if(all_edges.size()<=id)
			all_edges.resize(id+1);
		all_edges[id]={from,to,id,weight};
		if(weights.size()<=id)
			weights.resize(id+1,0);
		weights[id]=weight;
		//weights.push_back(weight);
		modifications++;
		additions=modifications;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"edge %d %d %d %d\n",from,to,weight, id+1);
				fflush(outfile);
			}
#endif
//		history.push_back({true,id,modifications});
		enableEdge(from,to,id);//default to enabled
	}
	int nEdgeIDs(){
		assert(edges==all_edges.size());
		return edges;//all_edges.size();
	}
	std::vector<int> & getWeights(){
		return weights;
	}
	int getWeight(int edgeID){
		return weights[edgeID];
		//return all_edges[edgeID].weight;
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
			edge_status[id]=true;
			//edge_status.setStatus(id,true);

			modifications++;
			additions=modifications;
			history.push_back({true,id,modifications,additions});
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

		if(history.back().addition  && history.back().id==id && history.back().mod==modifications){
			//edge_status.setStatus(id,false);
			edge_status[id]=false;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"-%d\n", id+1);
				fflush(outfile);
			}
#endif
			modifications--;
			additions = history.back().prev_mod;
			history.pop_back();
			return true;
		}
		return false;
	}

	void disableEdge(int from, int to, int id){
		assert(id>=0);
		assert(id<edge_status.size());
		assert(isEdge(id));
		if(edge_status[id]!=false){
			//edge_status.setStatus(id,false);
			edge_status[id]=false;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"-%d\n", id+1);
				fflush(outfile);
			}
#endif

			modifications++;

			history.push_back({false,id,modifications,deletions});
			deletions=modifications;
		}
	}

	bool undoDisableEdge( int id){
		assert(id>=0);
		assert(id<edge_status.size());
		assert(isEdge(id));
		if(!history.size())
			return false;

		if(!history.back().addition  && history.back().id==id && history.back().mod==modifications){
			//edge_status.setStatus(id,true);
			edge_status[id]=true;
#ifdef RECORD
			if(outfile){
				fprintf(outfile,"%d\n", id+1);
				fflush(outfile);
			}
#endif
			modifications--;
			deletions = history.back().prev_mod;
			history.pop_back();
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
			std::vector<Edge>& adj= adjacency[from];
			int i,j = 0;
			for(i = 0;i<adj.size();i++){
				if(adj[i]==to){
					edges--;
				}else{
					adj[j++]=adj[i];
				}
			}
			adj.resize(j);
		}

		modifications++;
		deletions=modifications;
		history.push_back({false,from,to,id,modifications});
	}*/

	bool rewindHistory(int steps){

		int cur_modifications = modifications;
		for(int i = 0;i<steps;i++){
			EdgeChange & e = history.back();
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


};
#endif /* DYNAMICGRAPH_H_ */
