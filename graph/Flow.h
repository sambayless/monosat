/*
 * Flow.h
 *
 *  Created on: 2013-04-15
 *      Author: sam
 */

#ifndef FLOW_H_
#define FLOW_H_

#include "mtl/Vec.h"
using namespace Minisat;
/*void BBFS(int v, int u,  vec<vec<int > > & backward_adj,vec<int> residual, vec<int> distance){
	static vec<char> seen;
	seen.clear();
	seen.growTo(backward_adj.size());
	seen[v]=true;
	static vec<int> w;
	w.clear();
	w.push(v);
	int n = backward_adj.size();//number of nodes

	while(w.size()){
		int x = w.last();
		w.pop();
		if(x==u)
			break;
		for (int i = 0;i<backward_adj[x].size();i++){
			int y = backward_adj[x][i];
			//this is an edge from y to x
			if (!seen[y] && residual[y]>0 && distance[y]>=n){
				seen[y]=true;
				w.push(y);
			}
		}
	}
}*/




void BFS(int v, int u,  vec<vec<int > > & adj,vec<int> & affected_out){
	static vec<char> in_w;
	static vec<char> in_affected;
	affected_out.clear();
	in_w.clear();
	in_w.growTo(adj.size());
	in_w[v]=true;

	in_affected.clear();
	in_affected.growTo(adj.size());
	in_affected[v]=true;

	static vec<int> w;
	w.clear();
	w.push(v);
	int n = adj.size();//number of nodes

	while(w.size()){
		int x = w.last();
		w.pop();
		if(x==u)
			break;
		for (int i = 0;i<adj[x].size();i++){
			int y = adj[x][i];
			//this is an edge from y to x
			if (!in_affected[y]){
				in_affected[y]=true;
				affected_out.push(y);
				if(!in_w[y]){
					in_w[y]=true;
					w.push(y);
				}
			}
		}
	}
}

vec<vec<bool> > f;//flow
vec<int> height;
vec<int> excess;
vec<bool> seen;

int n;
void push(int u, int v){
	int send =!f[u][v];
	if(send>excess[u]){
		send=excess[u];
	}
	f[u][v]=true;
	f[v][u]=false;
	excess[u]-=send;
	excess[v]+=send;
}

void relabel(int u,vec<vec<int > > & adj){
	int min_height = -1;
	//for (int v = 0;v<n;v++){
	for(int i = 0;i<adj[u].size();i++){
		int v = adj[u][i];
		if(!f[u][v]){//if capacities are all boolean, then we can change this into an iteration through the adj list for u
			if(min_height<0 || height[v]<min_height){
				min_height=height[v];
			}
			height[u]=min_height+1;
		}
	}
	//}
}

void discharge(int u,vec<vec<int > > & adj){
	while (excess[u] > 0){
		 if (seen[u] < n){
			 int v = seen[u];
			 if (!f[u][v] && (height[u] > height[v]))
				 push(u, v);
			 else
				 seen[u] += 1;
		 }else{ // we have checked all neighbours. must relabel
			 relabel(u,adj);
			 seen[u] = 0;
		 }
	}
}
void moveToFront(int i, vec<int> &A) {
              int temp = A[i];
              int n;
              for (n = i; n > 0; n--){
                      A[n] = A[n-1];
              }
              A[0] = temp;
      }
int relabel_to_front(vec<vec<int > > & adj,int s,int t,){
	//from wikipedia
	int n = adj.size();
	 seen.growTo(n);
	 excess.growTo(n);
	 height.growTo(n);
	f.growTo(n);
	//residual capacity from u to v is C[u][v] - F[u][v]
	vec<int> nodelist;

	for(int i = 0;i<n;i++){
		if(i!=s && i != t)
			nodelist.push(i);
	}

	height[s] = n;//    longest path from source to sink is less than n long
	excess[s] = 0xF0F0F0F0; // # send as much flow as possible to neighbours of source
	for(int v = 0;v<n;v++)
		 push(s, v);

    int p = 0;
    while (p < n - 2) {
            int u = nodelist[p];
            int old_height = height[u];
            discharge(u,adj);
            if (height[u] > old_height) {
                    moveToFront(p,nodelist);
                    p=0;
            }
            else
				p += 1;
    }

	 int flow = 0;
	 for(int i=0;i<f[s].size();i++)
		 flow+=f[s][i];
	 return flow;


}

/*
void push(int u, int v){
	int send = c[u][v]-f[u][v];
	if(send>excess[u]){
		send=excess[u];
	}
	f[u][v]+=send;
	f[v][u]-=send;
	excess[u]-=send;
	excess[v]+=send;
}

void relabel(int u){
	int min_height = -1;
	for (int v = 0;v<n;v++){
		if(c[u][v] - f[u][v]>0){//if capacities are all boolean, then we can change this into an iteration through the adj list for u
			if(min_height<0 || height[v]<min_height){
				min_height=height[v];
			}
			height[u]=min_height+1;
		}
	}
}

void discharge(int u){
	while (excess[u] > 0){
		 if (seen[u] < n){
			 int v = seen[u];
			 if ((c[u][v] - f[u][v] > 0) && (height[u] > height[v]))
				 push(u, v);
			 else
				 seen[u] += 1;
		 }else{ // we have checked all neighbours. must relabel
			 relabel(u);
			 seen[u] = 0;
		 }
	}
}
void moveToFront(int i, vec<int> &A) {
              int temp = A[i];
              int n;
              for (n = i; n > 0; n--){
                      A[n] = A[n-1];
              }
              A[0] = temp;
      }
int relabel_to_front(vec<vec<int > > & adj,int s,int t,){
	//from wikipedia
	int n = adj.size();
	 seen.growTo(n);
	 excess.growTo(n);
	 height.growTo(n);
	f.growTo(n);
	//residual capacity from u to v is C[u][v] - F[u][v]
	vec<int> nodelist;

	for(int i = 0;i<n;i++){
		if(i!=s && i != t)
			nodelist.push(i);
	}

	height[s] = n;//    longest path from source to sink is less than n long
	excess[s] = 0xF0F0F0F0; // # send as much flow as possible to neighbours of source
	for(int v = 0;v<n;v++)
		 push(s, v);

    int p = 0;
    while (p < n - 2) {
            int u = nodelist[p];
            int old_height = height[u];
            discharge(u);
            if (height[u] > old_height) {
                    moveToFront(p,nodelist);
                    p=0;
            }
            else
				p += 1;
    }

	 int flow = 0;
	 for(int i=0;i<f[s].size();i++)
		 flow+=f[s][i];
	 return flow;


}
*/




#endif /* FLOW_H_ */
