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

#ifndef FLOW_H_
#define FLOW_H_

#include <vector>

namespace dgl{
void BFS(int v, int u,  std::vector<std::vector<int > > & adj,std::vector<int> & affected_out){
	static std::vector<char> in_w;
	static std::vector<char> in_affected;
	affected_out.clear();
	in_w.clear();
	in_w.resize(adj.size());
	in_w[v]=true;

	in_affected.clear();
	in_affected.resize(adj.size());
	in_affected[v]=true;

	static std::vector<int> w;
	w.clear();
	w.push_back(v);
	int n = adj.size();//number of nodes

	while(w.size()){
		int x = w.back();
		w.pop_back();
		if(x==u)
			break;
		for (int i = 0;i<adj[x].size();i++){
			int y = adj[x][i];
			//this is an edge from y to x
			if (!in_affected[y]){
				in_affected[y]=true;
				affected_out.push_back(y);
				if(!in_w[y]){
					in_w[y]=true;
					w.push_back(y);
				}
			}
		}
	}
}
private:
std::vector<std::vector<bool> > f;//flow
std::vector<int> height;
std::vector<int> excess;
std::vector<bool> seen;

int n;
public:
void push_back(int u, int v){
	int send =!f[u][v];
	if(send>excess[u]){
		send=excess[u];
	}
	f[u][v]=true;
	f[v][u]=false;
	excess[u]-=send;
	excess[v]+=send;
}

void relabel(int u,std::vector<std::vector<int > > & adj){
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

void discharge(int u,std::vector<std::vector<int > > & adj){
	while (excess[u] > 0){
		 if (seen[u] < n){
			 int v = seen[u];
			 if (!f[u][v] && (height[u] > height[v]))
				 push_back(u, v);
			 else
				 seen[u] += 1;
		 }else{ // we have checked all neighbours. must relabel
			 relabel(u,adj);
			 seen[u] = 0;
		 }
	}
}
void moveToFront(int i, std::vector<int> &A) {
              int temp = A[i];
              int n;
              for (n = i; n > 0; n--){
                      A[n] = A[n-1];
              }
              A[0] = temp;
      }
int relabel_to_front(std::vector<std::vector<int > > & adj,int s,int t,){
	//from wikipedia
	int n = adj.size();
	 seen.resize(n);
	 excess.resize(n);
	 height.resize(n);
	f.resize(n);
	//residual capacity from u to v is C[u][v] - F[u][v]
	std::vector<int> nodelist;

	for(int i = 0;i<n;i++){
		if(i!=s && i != t)
			nodelist.push_back(i);
	}

	height[s] = n;//    longest path from source to sink is less than n long
	excess[s] = 0xF0F0F0F0; // # send as much flow as possible to neighbours of source
	for(int v = 0;v<n;v++)
		 push_back(s, v);

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
void push_back(int u, int v){
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
				 push_back(u, v);
			 else
				 seen[u] += 1;
		 }else{ // we have checked all neighbours. must relabel
			 relabel(u);
			 seen[u] = 0;
		 }
	}
}
void moveToFront(int i, std::vector<int> &A) {
              int temp = A[i];
              int n;
              for (n = i; n > 0; n--){
                      A[n] = A[n-1];
              }
              A[0] = temp;
      }
int relabel_to_front(std::vector<std::vector<int > > & adj,int s,int t,){
	//from wikipedia
	int n = adj.size();
	 seen.resize(n);
	 excess.resize(n);
	 height.resize(n);
	f.resize(n);
	//residual capacity from u to v is C[u][v] - F[u][v]
	std::vector<int> nodelist;

	for(int i = 0;i<n;i++){
		if(i!=s && i != t)
			nodelist.push_back(i);
	}

	height[s] = n;//    longest path from source to sink is less than n long
	excess[s] = 0xF0F0F0F0; // # send as much flow as possible to neighbours of source
	for(int v = 0;v<n;v++)
		 push_back(s, v);

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

};


#endif /* FLOW_H_ */
