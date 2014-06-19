
#ifndef ALLPAIRS_H_
#define ALLPAIRS_H_

#include <vector>
namespace dgl{
class AllPairs{
public:

	struct NullStatus{
		void setReachable(int u, int v, bool reachable){

		}
		void setMininumDistance(int u,int v,bool reachable,int distance){

		}

	};
	static NullStatus nullStatus;


	virtual ~AllPairs(){};

	virtual void addSource(int s)=0;

	virtual void update( )=0;


	virtual bool connected_unsafe(int from,int t)=0;
	virtual bool connected_unchecked(int from,int t)=0;
	virtual bool connected(int from,int t)=0;
	virtual int distance(int from,int t)=0;
	virtual int distance_unsafe(int from,int t)=0;

	virtual void getPath(int source, int to, std::vector<int> & path_store)=0;
};
};

#endif
